"""Delta-adjusted two-year endpoint sensitivity within documented survivors."""
from datetime import timedelta
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import patsy
from scipy.special import expit, logit
from scipy.stats import t
import statsmodels.api as sm
from common import NotEstimable, finite_matrix, write_csv, write_json, sha256

RHS = ('age10 + sex_male + education_c + apoe4_carrier + baseline_cdrsb'
       ' + C(diagnosis) + centiloid10 + hip1000 + tau01')


def survival_class(upper, cutoff, visit_dates, death_date=None, conflict=False):
    if conflict:
        return 'conflicting_evidence'
    if death_date and any(d > death_date for d in visit_dates if d <= cutoff):
        return 'conflicting_evidence'
    if death_date and death_date <= upper:
        return 'documented_death'
    if (death_date and upper < death_date <= cutoff) or any(upper <= d <= cutoff for d in visit_dates):
        return 'documented_alive'
    return 'unknown'


def prepare(base, out, plan, vital_file=None):
    import analyze_tau_incremental as inc
    import longitudinal as core
    import longitudinal_sensitivity as ls
    raw=base/'data/private/adni/raw'
    inputs=[base/'data/private/adni/adni_pet_aligned_master.csv',
            *[raw/name for name in ['CDR.csv','DXSUM.csv','UCBERKELEY_AMY_6MM.csv','UCSFFSX7.csv','DATA_DOWNLOADED_DATE.csv']]]
    if vital_file:inputs.append(Path(vital_file))
    write_json(out/'MNAR_input_provenance.json',dict(
        inputs=[dict(file=p.name,bytes=p.stat().st_size,sha256=sha256(p)) for p in inputs],
        download_performed=False,participant_records_exported=False))
    print('Reading retained ADNI master and clinical visits',flush=True)
    master,_=core.load_master(base/'data/private/adni/adni_pet_aligned_master.csv')
    ids={p['rid'] for p in master}
    cdr,_=core.load_visits(raw/'CDR.csv',ids,'cdr')
    dx,_=core.load_visits(raw/'DXSUM.csv',ids,'dx')
    eligible,_,_=core.make_cohort(master,cdr)
    print('Matching retained amyloid and MRI baseline measurements',flush=True)
    amy,_=inc.match_centiloids(master,raw/'UCBERKELEY_AMY_6MM.csv')
    mri,_=inc.match_mri_dates(master,raw/'UCSFFSX7.csv')
    people=inc.enrich(eligible,amy,mri,dx)
    d = pd.DataFrame(people)
    if d.empty:raise NotEstimable('No A-positive participants with eligible baseline clinical data')
    needed = inc.PREDICTOR_COLUMNS
    good = d[needed].notna().all(axis=1) & d.diagnosis.ne('Unavailable') & d.mri_match_status.eq('matched')
    flow = [dict(stage='baseline_A_positive_eligible', n=len(d)),
            dict(stage='complete_baseline_predictors', n=int(good.sum()))]
    d = inc.scaled(d.loc[good].copy()).reset_index(drop=True)
    if d.empty:
        write_csv(out/'MNAR_sample_flow.csv',flow)
        raise NotEstimable('No participants with the complete fixed baseline predictor set')
    cutoff, source = ls.cutoff_date(base/'data/private/adni/raw')
    vital = {}
    if vital_file:
        v = pd.read_csv(vital_file, dtype=str).fillna('')
        required = {'rid', 'last_known_alive_date', 'death_date', 'source_verified', 'source'}
        if not required.issubset(v) or v.rid.duplicated().any():
            raise ValueError('Vital-status input requires unique rid and the documented template columns')
        for r in v.to_dict('records'):
            if r['source_verified'].lower() not in ('true', '1') or not r['source'].strip():
                raise ValueError('Every supplied vital-status row needs verified source provenance')
            for k in ['last_known_alive_date', 'death_date']:
                value = core.iso_date(r[k]) if r[k] else None
                if r[k] and value is None:
                    raise ValueError('Vital-status dates must be complete ISO dates')
                r[k] = value
            key=core.rid(r['rid'])
            if key is None or key in vital:raise ValueError('Invalid or equivalent duplicate RID in verified vital-status data')
            vital[key] = r
    rows = []
    for p in d.to_dict('records'):
        key = p['rid']
        target = p['baseline_date'] + timedelta(days=365.25 * plan['years'])
        left, right = target-timedelta(days=plan['window_days']), target+timedelta(days=plan['window_days'])
        visits = {day: values for day, values in cdr.get(key, {}).items() if day <= cutoff}
        valid_dates = [day for day, values in visits.items() if len(values) == 1]
        record = vital.get(key, {})
        alive_date = record.get('last_known_alive_date')
        if alive_date:
            valid_dates.append(alive_date)
        alive = survival_class(right, cutoff, valid_dates, record.get('death_date'))
        selected, reason = core.nearest(visits, target, plan['window_days'])
        outcome = selected[1] if selected else np.nan
        if selected and selected[0] <= max(p['baseline_date'], p['anchor']):
            raise ValueError('Endpoint does not follow baseline/PET')
        rows.append(dict(**p, target_date=target, window_end=right,
            administrative_opportunity=right <= cutoff, survival=alive,
            outcome_cdrsb=outcome, endpoint_status=reason,
            endpoint_offset_days=(selected[0]-target).days if selected else np.nan))
    all_rows = pd.DataFrame(rows)
    private = out/'private'
    write_csv(private/'eligibility.csv', all_rows)
    opportunity = all_rows.administrative_opportunity
    flow.append(dict(stage='full_administrative_opportunity', n=int(opportunity.sum())))
    for label, g in all_rows.loc[opportunity].groupby('survival'):
        flow.append(dict(stage=label, n=len(g)))
    selected = all_rows.loc[opportunity & all_rows.survival.eq('documented_alive')].reset_index(drop=True)
    flow += [dict(stage='selected_survivor_target', n=len(selected)),
             dict(stage='observed_two_year_endpoint', n=int(selected.outcome_cdrsb.notna().sum())),
             dict(stage='missing_two_year_endpoint', n=int(selected.outcome_cdrsb.isna().sum()))]
    write_csv(out/'MNAR_sample_flow.csv', flow)
    write_json(out/'MNAR_estimand.json', dict(**plan, cutoff=str(cutoff), cutoff_source=source,
        survival_evidence='A valid CDR assessment on/after the window end, or an explicitly verified vital-status table. ADSL blank death flags and end-of-study dates are not treated as survival evidence.',
        population_boundary='Conditions on documented survival and complete baseline predictors. Does not identify outcomes after death, account for unknown survival, or replace the 185-person trajectory model.'))
    return selected


def rubin(q, u, df_complete):
    q, u = np.asarray(q, float), np.asarray(u, float)
    if len(q) < 2 or len(q) != len(u) or not np.isfinite(q).all() or not np.isfinite(u).all() or np.any(u < 0):
        raise ValueError('Invalid pooling inputs')
    m = len(q); between = float(q.var(ddof=1)); within = float(u.mean())
    total = within + (1+1/m)*between
    if total <= 0:
        raise NotEstimable('Zero pooled variance')
    fraction = (1+1/m)*between/total
    old_df = (m-1)/fraction**2 if fraction > 0 else np.inf
    observed_df = (df_complete+1)/(df_complete+3)*df_complete*(1-fraction)
    df = 1/(1/old_df+1/observed_df)
    estimate = float(q.mean()); se = float(np.sqrt(total)); crit = float(t.ppf(.975, df))
    return dict(estimate=estimate, se=se, ci_low=estimate-crit*se, ci_high=estimate+crit*se,
        nominal_p=float(2*t.sf(abs(estimate/se), df)), df=float(df),
        within_variance=within, between_variance=between, missing_variance_fraction=fraction,
        monte_carlo_se=float(np.sqrt(between/m)), imputations=m)


def pmm_draws(X, y, m, k, seed):
    observed = np.isfinite(y); miss = ~observed
    xo = finite_matrix(X[observed]); ym = logit((y[observed]+.25)/18.5)
    if len(xo) <= xo.shape[1]+3:
        raise NotEstimable('Insufficient observed outcomes for the imputation design')
    beta = np.linalg.lstsq(xo, ym, rcond=None)[0]
    residual = ym-xo@beta; rss = float(residual@residual); df = len(xo)-xo.shape[1]
    if rss <= 1e-10:
        raise NotEstimable('Degenerate outcome distribution for PMM')
    chol = np.linalg.cholesky(np.linalg.inv(xo.T@xo))
    donor_prediction = xo@beta
    rng = np.random.default_rng(seed)
    for i in range(m):
        variance = rss/rng.chisquare(df)
        sampled = beta + np.sqrt(variance)*(chol@rng.normal(size=len(beta)))
        prediction = X[miss]@sampled
        completed = y.copy()
        for j, pred in zip(np.flatnonzero(miss), prediction):
            distances = np.abs(donor_prediction-pred)
            # Shuffle before sorting so tied predictions do not always select the same donors.
            order = rng.permutation(len(distances))
            choices = order[np.argsort(distances[order], kind='stable')[:min(k, len(order))]]
            completed[j] = y[observed][rng.choice(choices)]
        yield completed


def fractional_hc3(X, y_fraction, mu):
    X=np.asarray(X,float); y_fraction=np.asarray(y_fraction,float); mu=np.asarray(mu,float)
    weights=mu*(1-mu)
    try:bread=np.linalg.inv(X.T@(weights[:,None]*X))
    except np.linalg.LinAlgError as e:raise NotEstimable('Singular fractional-logit information matrix') from e
    leverage=weights*np.einsum('ij,jk,ik->i',X,bread,X)
    if not np.isfinite(leverage).all() or np.any(leverage>=1-1e-10) or np.any(leverage < -1e-8):
        raise NotEstimable('Invalid or unit leverage for HC3 covariance')
    scores=X*((y_fraction-mu)/(1-leverage))[:,None]
    covariance=bread@(scores.T@scores)@bread
    return (covariance+covariance.T)/2


def effect(X, high, low, y, scale):
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        fit = sm.GLM(y/18, X, family=sm.families.Binomial()).fit(maxiter=200)
    if (not fit.converged or not np.isfinite(fit.params).all() or
        any('Separation' in type(w.message).__name__ for w in ws)):
        raise NotEstimable('Fractional-logit convergence/separation failure')
    ph, pl = expit(high@fit.params), expit(low@fit.params)
    q = float(18*np.mean(ph-pl)*scale)
    grad = 18*np.mean(high*(ph*(1-ph))[:, None]-low*(pl*(1-pl))[:, None], axis=0)*scale
    # GLM HC3 is computed explicitly; the generic statsmodels HC3 label can fall back to HC0.
    covariance=fractional_hc3(X,y/18,expit(X@fit.params))
    u = float(grad@covariance@grad)
    if not np.isfinite(u) or u <= 0:
        raise NotEstimable('Invalid contrast variance')
    return q, u


def analyze(d, out, plan, seed, diagnostic=False):
    y = d.outcome_cdrsb.to_numpy(float)
    if len(d) == 0:
        raise NotEstimable('No participants with documented survival through the full outcome window')
    X = patsy.dmatrix(RHS, d, return_type='dataframe', NA_action='raise')
    matrix = finite_matrix(X); observed = np.isfinite(y); missing = ~observed
    if observed.sum() < max(plan['minimum_observed'], 3*X.shape[1]):
        raise NotEstimable('Too few observed endpoints for the fixed covariate model')
    if not missing.any():
        raise NotEstimable('No missing endpoints in the documented-survivor target; MNAR not applicable')
    if np.any((y[observed]<0)|(y[observed]>18)):
        raise ValueError('CDR-SB outside [0,18]')
    q1, median, q3 = d.tau_temporal_suvr.quantile([.25, .5, .75])
    iqr = float(q3-q1)
    if iqr <= 0:
        raise NotEstimable('Zero tau IQR')
    hi, lo = d.copy(), d.copy()
    hi['tau01'], lo['tau01'] = (q3-1.34)/.1, (q1-1.34)/.1
    high = np.asarray(patsy.build_design_matrices([X.design_info], hi)[0])
    low = np.asarray(patsy.build_design_matrices([X.design_info], lo)[0])
    scenarios = [(a, b) for a in plan['delta_points'] for b in plan['tau_gradient_points']]
    estimates = {s: [] for s in scenarios}; variances = {s: [] for s in scenarios}
    clipped = {s: 0 for s in scenarios}; shifts = {s: [] for s in scenarios}
    failures = []
    m = 5 if diagnostic else plan['imputations']
    centered = (d.tau_temporal_suvr.to_numpy()-median)/iqr
    for i, mar in enumerate(pmm_draws(matrix, y, m, plan['pmm_donors'], seed), 1):
        for delta, gamma in scenarios:
            scenario = (delta, gamma)
            raw = mar[missing]+delta+gamma*centered[missing]
            completed = mar.copy(); completed[missing] = np.clip(raw, 0, 18)
            clipped[scenario] += int(np.sum(raw != completed[missing]))
            shifts[scenario].append(float(np.mean(completed[missing]-mar[missing])))
            try:
                estimate, variance = effect(matrix, high, low, completed, .1/iqr)
                estimates[scenario].append(estimate); variances[scenario].append(variance)
            except NotEstimable as e:
                failures.append(dict(imputation=i, delta=delta, tau_gradient=gamma, reason=str(e)))
        if i == 1 or i % 10 == 0 or i == m:
            print(f'MNAR imputations {i}/{m}; fixed scenarios {len(scenarios)}', flush=True)
    results = []
    for s in scenarios:
        row = dict(delta_points=s[0], tau_gradient_points=s[1], n_target=len(d), n_observed=int(observed.sum()),
            n_missing=int(missing.sum()), tau_iqr=iqr, requested_imputations=m, valid_imputations=len(estimates[s]),
            clipped_fraction=clipped[s]/(m*missing.sum()), mean_actual_shift=float(np.mean(shifts[s])),
            diagnostic_only=diagnostic)
        if len(estimates[s]) == m:
            row.update(status='estimated', **rubin(estimates[s], variances[s], len(d)-X.shape[1]))
        else:
            row.update(status='failed_imputations', reason='No incomplete-imputation pooling')
        results.append(row)
    write_csv(out/'MNAR_scenarios.csv', results)
    write_csv(out/'MNAR_fit_failures.csv', pd.DataFrame(failures, columns=['imputation','delta','tau_gradient','reason']))
    write_json(out/'MNAR_model.json', dict(formula=RHS, family='fractional logit; explicit leverage-adjusted HC3 sandwich covariance',
        covariance='Bread=(X\'WX)^-1; scores=x_i*(y_i/18-mu_i)/(1-h_i), h_i=w_i*x_i\'Bread*x_i; covariance=Bread*sum(scores*scores\')*Bread',
        imputation='Bayesian linear PMM on logit((CDR-SB+0.25)/18.5), type-1 matching, five observed donors',
        adjustment='Only missing outcomes shifted: delta + gamma*(tau - target median)/target IQR; truncated to [0,18]',
        pooling='Rubin scalar contrast pooling with Barnard-Rubin finite-sample degrees of freedom',
        reporting='All fixed scenarios retained; this is a conditional two-year endpoint sensitivity, not a fitted death/MNAR selection mechanism.'))
    return pd.DataFrame(results)
