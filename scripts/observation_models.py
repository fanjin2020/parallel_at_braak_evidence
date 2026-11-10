"""Observation diagnostics and fixed exploratory weighted trajectory sensitivities."""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import longitudinal_observation as previous
from analysis_utils import NotEstimable, csv, dump, checkpoint, bootstrap_people
from observation_statistics import dispersion, annotate_effects


def formula(name):
    common = 'C(period) + ' + previous.COV + ' + tau01 + last_cdr + years_since_last + duration_years'
    if name == 'full_history':
        return common + ' + baseline_cdrsb + prior_visits'
    if name == 'reduced_history':
        return common
    raise ValueError(name)


def effective_n(w):
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / np.dot(w, w)) if len(w) and np.dot(w, w) > 0 else np.nan


def quantiles(values):
    a = np.asarray(values, float)
    a = a[np.isfinite(a)]
    if not len(a):
        return dict(n=0)
    return dict(n=len(a), **{k: float(v) for k, v in zip(
        ['minimum', 'q01', 'q05', 'q25', 'median', 'q75', 'q95', 'q99', 'maximum'],
        np.quantile(a, [0, .01, .05, .25, .5, .75, .95, .99, 1]))})


def coverage(periods):
    rows = []
    for grouping in ['all', 'diagnosis']:
        d = periods.copy()
        d['group'] = 'all' if grouping == 'all' else d.diagnosis
        for (period, group), g in d.groupby(['period', 'group'], observed=True):
            rows.append(dict(grouping=grouping, group=group, period=int(period),
                people=int(g.rid.nunique()), observed_people=int(g.observed.sum()),
                visits=int(g.visits.sum()), median_window_years=float(g.duration_years.median()),
                empirical_support='all_zero' if g.observed.max() == 0 else
                'all_one' if g.observed.min() == 1 else 'variable',
                interpretation='Observation counts do not establish scheduled visits, survival or reasons for absence'))
    return rows


def fit_observation(periods, model):
    if periods.duplicated(['rid', 'period']).any():
        raise ValueError('Duplicated participant-window')
    d = periods.reset_index(drop=True).copy()
    by = d.groupby('period')['observed'].agg(['min', 'max'])
    variable = by.index[by['min'] != by['max']]
    d['support'] = d.period.map({k: 'variable' if r['min'] != r['max'] else
        'empirical_all_one' if r['min'] == 1 else 'empirical_all_zero' for k, r in by.iterrows()})
    active = d.period.isin(variable)
    fit_data = d.loc[active]
    if len(fit_data) < 40 or fit_data.observed.nunique() != 2:
        raise NotEstimable('Insufficient variable observation windows')
    X = patsy.dmatrix(formula(model), fit_data, return_type='dataframe', NA_action='raise')
    rank = int(np.linalg.matrix_rank(X))
    if rank != X.shape[1]:
        raise NotEstimable(f'Observation design rank {rank}/{X.shape[1]}; no terms automatically removed')
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        f = sm.GLM(fit_data.observed, X, family=sm.families.Binomial()).fit(maxiter=200)
    pred = np.asarray(f.predict(X), float)
    separated = any('Separation' in type(w.message).__name__ for w in ws)
    if not f.converged or separated or not np.isfinite(pred).all() or not np.isfinite(f.params).all():
        raise NotEstimable('Observation model failed convergence, separation or finite-value checks')
    if np.any((pred <= 0) | (pred >= 1)):
        raise NotEstimable('Numerically degenerate observation probabilities; no clipping applied')
    d['probability'] = np.nan
    d.loc[active, 'probability'] = pred
    d.loc[d.support.eq('empirical_all_one'), 'probability'] = 1.0
    d['marginal'] = d.period.map(d.groupby('period')['observed'].mean())
    # No probability of one is assigned to windows containing no observed outcomes.
    d['raw_weight'] = d.marginal / d.probability
    info = dict(model=model, n_variable_windows=int(active.sum()), design_rank=rank,
        design_columns=X.shape[1], converged=bool(f.converged), separation_warning=separated,
        n_zero_periods=int((by['max'] == 0).sum()), n_one_periods=int((by['min'] == 1).sum()),
        **previous.probability_diagnostics(pred, fit_data.observed))
    info['original_rule_passed'] = not (info['n_probability_below'] or info['n_probability_above'])
    return d, info


def map_weights(periods, long, trim=None):
    lookup = periods.set_index(['rid', 'period']).raw_weight
    raw = np.array([1.0 if r.period == -1 else lookup.loc[(r.rid, r.period)] for r in long.itertuples()], float)
    follow = long.period.ne(-1).to_numpy()
    if not follow.any() or not np.isfinite(raw).all() or np.any(raw <= 0):
        raise NotEstimable('Nonfinite or nonpositive observed weights')
    cuts = np.quantile(raw[follow], trim) if trim is not None else [raw[follow].min(), raw[follow].max()]
    used = raw.copy()
    used[follow] = np.clip(raw[follow], *cuts)
    return raw, used, cuts


def weight_diagnostics(long, raw, used, cuts):
    follow = long.period.ne(-1).to_numpy()
    person = pd.DataFrame({'rid': long.loc[follow, 'rid'].to_numpy(), 'weight': raw[follow]}).groupby('rid').weight.sum()
    share = person / person.sum()
    return dict(raw_followup_ess=effective_n(raw[follow]), used_followup_ess=effective_n(used[follow]),
        raw_effective_people=effective_n(person), raw_max_person_share=float(share.max()),
        raw_top_five_people_share=float(share.nlargest(5).sum()), raw_weight_max=float(raw[follow].max()),
        used_weight_max=float(used[follow].max()), trim_low=float(cuts[0]), trim_high=float(cuts[1]),
        clipped_observations=int(np.sum(raw != used)), n_followup_visits=int(follow.sum()),
        n_people=int(long.rid.nunique()), n_observations=len(long))


def support_gate(diag, plan):
    reasons = []
    if diag['raw_effective_people'] < plan['minimum_effective_people']:
        reasons.append('raw effective participant count below fixed reporting safeguard')
    if diag['raw_max_person_share'] > plan['maximum_person_weight_share']:
        reasons.append('one participant exceeds fixed raw weight concentration safeguard')
    return '; '.join(reasons)


def balance_rows(periods, model, scenario, cuts):
    target = periods.loc[periods.support.ne('empirical_all_zero')].copy()
    observed = target.loc[target.observed.eq(1)].copy()
    w = observed.raw_weight.to_numpy()
    # One record per observed window, using the same cutoffs as the outcome fit.
    w = np.clip(w, *cuts)
    targets = ['age10', 'sex_male', 'education_c', 'apoe4_carrier', 'tau01',
               'centiloid10', 'hip1000', 'baseline_cdrsb', 'last_cdr', 'prior_visits', 'years_since_last']
    rows = []
    for period, tg in target.groupby('period'):
        ix = observed.period.eq(period).to_numpy()
        og = observed.loc[ix]
        if not len(og):
            continue
        for field in targets:
            scale = dispersion(tg[field])
            sd = scale['target_sd']
            mean = float(tg[field].mean())
            ow = float(np.average(og[field], weights=w[ix]))
            rows.append(dict(model=model, scenario=scenario, period=int(period), variable=field,
                target_mean=mean, observed_mean=float(og[field].mean()), weighted_mean=ow,
                unweighted_smd=(float(og[field].mean())-mean)/sd if scale['smd_applicable'] else np.nan,
                weighted_smd=(ow-mean)/sd if scale['smd_applicable'] else np.nan,
                **scale,
                target_windows=len(tg), observed_windows=len(og),
                scope='Within-period measured-covariate diagnostic; target excludes all-zero periods'))
    return rows


def fixed_estimates(periods, long, plan, diagnostics=False):
    results, fit_rows, prob_rows, weight_rows, balance, panels = [], [], [], [], [], {}
    def point(weights, name):
        fit, design = previous.fit_gee(long, weights=weights)
        ans = previous.change(fit, design, long, name)
        ans['delta_ci_low'] = ans.pop('ci_low')
        ans['delta_ci_high'] = ans.pop('ci_high')
        return ans
    results.append(dict(**point(None, 'unweighted'), status='estimated', original_rule_passed=True))
    for model in plan['observation_models']:
        try:
            pp, info = fit_observation(periods, model)
        except (NotEstimable, np.linalg.LinAlgError) as e:
            fit_rows.append(dict(model=model, status='not_estimable', reason=str(e)))
            for scenario in plan['weight_scenarios']:
                results.append(dict(model=model+'__'+scenario, status='not_estimable', reason=str(e)))
            continue
        fit_rows.append(dict(**info, status='fit_completed', reason=''))
        if diagnostics:
            for group, g in pp.loc[pp.support.eq('variable')].groupby('observed'):
                prob_rows.append(dict(model=model, observed=int(group), **quantiles(g.probability)))
            panels[model] = pp
        for scenario, trim in plan['weight_scenarios'].items():
            label = model+'__'+scenario
            raw, used, cuts = map_weights(pp, long, trim)
            wd = weight_diagnostics(long, raw, used, cuts)
            reason = support_gate(wd, plan)
            if diagnostics:
                weight_rows.append(dict(model=model, scenario=scenario, **wd,
                    **{'used_'+k:v for k,v in quantiles(used[long.period.ne(-1)]).items()}))
                balance.extend(balance_rows(pp, model, scenario, cuts))
            if reason:
                results.append(dict(model=label, status='support_limited', reason=reason,
                    original_rule_passed=info['original_rule_passed']))
                continue
            try:
                ans = point(used, label)
                results.append(dict(**ans, status='estimated', original_rule_passed=info['original_rule_passed'],
                    support_warning=bool(info['n_zero_periods'] or not info['original_rule_passed']),
                    interpretation='Exploratory weighted sensitivity, not a replacement primary result or proof of positivity'))
            except (NotEstimable, np.linalg.LinAlgError) as e:
                results.append(dict(model=label, status='not_estimable', reason=str(e), original_rule_passed=info['original_rule_passed']))
    return results, fit_rows, prob_rows, weight_rows, balance, panels


def bootstrap_sample(periods, long, seed):
    sample = bootstrap_people(long, np.random.default_rng(seed))
    groups = dict(tuple(periods.groupby('rid', sort=False)))
    blocks = []
    for row in sample[['rid', 'original_rid']].drop_duplicates().itertuples():
        p = groups[row.original_rid].copy()
        p['rid'] = row.rid
        blocks.append(p)
    return pd.concat(blocks, ignore_index=True), sample


def run(root, out, plan, diagnostic=False):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    people, common, cdr = previous.common_data(root)
    cutoff, cutoff_source = previous.ls.cutoff_date(root/'data/private/adni/raw')
    periods, long = previous.observation_periods(common, cdr, cutoff)
    if len(periods) != len(periods[['rid', 'period']].drop_duplicates()):
        raise ValueError('Duplicate windows')
    if int(periods.visits.sum()) != int(long.period.ne(-1).sum()):
        raise ValueError('Window visit counts do not match analysis rows')
    print(f'Common sample: {common.rid.nunique()} participants, {len(common)} observations; {len(periods)} windows.', flush=True)
    points, models, probs, weights, balance, panels = fixed_estimates(periods, long, plan, True)
    prior = pd.read_csv(root/'revision/longitudinal/common_sample_trajectory_contrasts.csv')
    ref = prior.loc[prior.model.eq('fractional_logit_MRI') & prior.year.eq(2)]
    if len(ref) != 1 or not np.isclose(points[0]['estimate'], float(ref.estimate.iloc[0]), rtol=1e-6, atol=1e-9):
        raise ValueError('Primary model does not reproduce retained baseline result; weighting not run')
    checks = [dict(check='primary_effect_reproduction', actual=points[0]['estimate'], expected=float(ref.estimate.iloc[0]), passed=True),
        dict(check='common_people', actual=int(common.rid.nunique()), expected=int(ref.n_standardized_people.iloc[0]), passed=int(common.rid.nunique())==int(ref.n_standardized_people.iloc[0])),
        dict(check='followup_window_mapping', actual=int(periods.visits.sum()), expected=int(long.period.ne(-1).sum()), passed=True)]
    if not all(c['passed'] for c in checks):
        raise ValueError('Primary sample mismatch')
    csv(out/'checks.csv', checks)
    eligible=pd.DataFrame(people)
    eligible['in_common']=eligible.rid.isin(common.rid)
    flow=[dict(stage='baseline_eligible',n_people=len(eligible)),
          dict(stage='has_required_followup',n_people=int(eligible.included.sum())),
          dict(stage='common_model',n_people=int(eligible.in_common.sum()))]
    csv(out/'selection_flow.csv',flow)
    contrasts=[]
    for selection in ['included','in_common']:
        for col in ['age_at_tau','sex_male','education_years','apoe4_carrier','tau_temporal_suvr','baseline_cdrsb','centiloids','hippocampus_icv']:
            if col not in eligible:continue
            values=pd.to_numeric(eligible[col],errors='coerce')
            a=values.loc[eligible[selection]].dropna();b=values.loc[~eligible[selection]].dropna()
            sd=np.sqrt((a.var()+b.var())/2) if len(a)>1 and len(b)>1 else np.nan
            contrasts.append(dict(selection=selection,variable=col,selected_n=len(a),unselected_n=len(b),
                selected_mean=a.mean(),unselected_mean=b.mean(),smd=(a.mean()-b.mean())/sd if sd>0 else np.nan,
                scope='Descriptive initial-selection comparison; no selection weighting or mortality inference'))
    csv(out/'selection_comparison.csv',contrasts)
    csv(out/'coverage.csv', coverage(periods));csv(out/'observation_models.csv', models)
    csv(out/'probability_summary.csv', probs);csv(out/'weight_summary.csv', weights)
    csv(out/'balance.csv', balance)
    csv(out/'point_estimates.csv', points)
    old = next((d for d in models if d['model']=='full_history'), {})
    dump(out/'original_stopping_rule.json', dict(status='passed' if old.get('original_rule_passed') else 'not_passed',
        diagnostics=old, inference_from_original_rule='No weighted estimate is claimed under the original rule when it fails'))
    # Participant identifiers and individual weights stay outside the submission export.
    private = out/'private';private.mkdir(exist_ok=True)
    for model, pp in panels.items():
        raw, _, _ = map_weights(pp, long)
        d = long[['rid','period','time_years']].copy();d['raw_weight']=raw
        csv(private/(model+'_weights_RESTRICTED.csv'), d)
    draws = []
    B = int(plan['bootstrap'])
    for b in range(B):
        def one():
            pp, ll = bootstrap_sample(periods, long, plan['seed']+b)
            try:
                rows = fixed_estimates(pp, ll, plan)[0]
                return [dict(replicate=b+1, **r) for r in rows]
            except (NotEstimable, np.linalg.LinAlgError) as e:
                return [dict(replicate=b+1, model=r['model'], status='not_estimable', reason=str(e)) for r in points]
        draws.extend(checkpoint(out, 'observation_boot_'+str(b+1), one))
        dump(out/'progress.json', dict(stage='bootstrap', completed=b+1, requested=B, diagnostic_only=diagnostic))
        if (b+1)%10==0 or b+1==B:
            print(f'Observation bootstrap {b+1}/{B}; all fixed specifications evaluated.', flush=True)
    draw = pd.DataFrame(draws)
    csv(out/'bootstrap.csv', draw)
    for row in points:
        values = pd.to_numeric(draw.loc[draw.model.eq(row['model']) & draw.status.eq('estimated'), 'estimate'], errors='coerce').dropna() if 'estimate' in draw else pd.Series(dtype=float)
        ok = len(values)>=plan['minimum_valid_fraction']*B and not diagnostic
        row.update(bootstrap_requested=B, bootstrap_attempted=B, bootstrap_valid=len(values),
            ci_low=float(values.quantile(.025)) if ok and row['status']=='estimated' else np.nan,
            ci_high=float(values.quantile(.975)) if ok and row['status']=='estimated' else np.nan,
            interval='participant percentile bootstrap; both models and trimming re-estimated', diagnostic_only=diagnostic)
        if row['status']=='estimated':
            row['status']='diagnostic_only' if diagnostic else 'bootstrap_stable' if ok else 'bootstrap_unstable'
    corrected, counts = annotate_effects(pd.DataFrame(points), draw)
    points = corrected.to_dict('records')
    csv(out/'effects.csv', corrected)
    csv(out/'bootstrap_summary.csv', counts)
    dump(out/'specification.json', dict(plan=plan, observation_formula={m:formula(m) for m in plan['observation_models']},
        outcome_formula=previous.RHS, administrative_cutoff=str(cutoff), cutoff_source=cutoff_source,
        n_people=int(common.rid.nunique()), n_observations=len(common), diagnostic_only=diagnostic,
        raw_weight='period-specific observed proportion / fitted probability; baseline weight one; repeated visits share a window weight',
        history='Observed CDR history at or before the left boundary; current-window outcome excluded',
        all_zero_periods='No probability or outcome is imputed. Their absence of observed support is retained.',
        selection='Conditions on the existing common sample; does not recover the 336-person baseline population'))
    return panels, points


