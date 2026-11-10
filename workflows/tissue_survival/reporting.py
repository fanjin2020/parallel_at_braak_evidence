"""Validate aggregate datafill outputs and prepare one set of submission figures."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import textwrap

import numpy as np
import pandas as pd
from result_schema import CELL_NAMES, STATES, definition

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TARGETS = ['MTHFD2','OXTR','HOPX','CXCL10','IRF1','CD14','TAP1','PTPN6',
           'TNFAIP2','ICAM1','HAMP','ITPKC','JAK3','GMFG','HCK','CEBPB']
MODELS = ['age_sex_batch_RIN', 'plus_PMI']
SOURCE_CSV = {
    'ADNI_event_variable_dictionary.csv', 'ADNI_ignored_record_counts.csv',
    'ADNI_table_inventory.csv', 'Attrition_reason_counts.csv',
    'GSE174367_donors.csv', 'GSE174367_donor_cell_counts.csv',
    'GSE174367_exact_contrast_status.csv', 'GSE174367_GO_all_models.csv',
    'GSE174367_model_dimensions.csv', 'GSE174367_model_status.csv',
    'GSE174367_target_results.csv', 'MNAR_fit_failures.csv',
    'MNAR_sample_flow.csv', 'MNAR_scenarios.csv', 'Survival_classification_changes.csv',
}
SOURCE_JSON = {'Analysis_plan.json', 'Environment.json', 'ADNI_vital_provenance.json',
               'GSE174367_source_provenance.json', 'MNAR_model.json'}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')


def write_csv(path, rows, columns=None):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False, encoding='utf-8-sig')


def bh(p):
    p = np.asarray(p, dtype=float)
    if not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError('BH input must contain finite probabilities')
    order = np.argsort(p, kind='stable')
    ranked = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(1, ranked)
    return out


def truth(values):
    x = pd.Series(values).astype(str).str.casefold()
    if not x.isin(['true', 'false', '1', '0']).all():
        raise ValueError('Unrecognized boolean values')
    return x.isin(['true', '1']).to_numpy()


def validate(final, plan):
    final = Path(final); checks = []; notes = []
    def check(name, ok, detail=''):
        checks.append(dict(check=name, passed=bool(ok), detail=str(detail)))
    def close(name, x, y, tol=1e-10):
        x, y = np.asarray(x, float), np.asarray(y, float)
        ok = x.shape == y.shape and np.isfinite(x).all() and np.isfinite(y).all()
        err = float(np.max(np.abs(x-y))) if ok and x.size else 0.0
        check(name, ok and err <= tol, f'max_absolute_difference={err:.3g}')
    def qcheck(name, frame, col):
        close(name, frame[col], bh(frame.PValue.fillna(1)))
    def tests(frame, name, estimate='logFC'):
        tested = frame.status.eq('tested')
        check(name+'_untested_P_missing', frame.loc[~tested, 'PValue'].isna().all())
        p = frame.loc[tested, 'PValue'].to_numpy(float)
        check(name+'_tested_P_valid', np.isfinite(p).all() and np.all((p >= 0) & (p <= 1)))
        if estimate in frame:
            t = frame.loc[tested]
            check(name+'_CI_order', np.isfinite(t[[estimate,'CI_L','CI_R']].to_numpy(float)).all()
                  and ((t.CI_L <= t[estimate]) & (t[estimate] <= t.CI_R)).all())

    mpath = final/'MNAR_scenarios.csv'
    if mpath.is_file():
        from scipy.stats import t as student_t
        m = pd.read_csv(mpath); spec = plan['mnar']
        expected = {(float(d),float(g)) for d in spec['delta_points'] for g in spec['tau_gradient_points']}
        check('MNAR_scenario_grid', len(m) == len(expected) and set(zip(m.delta_points,m.tau_gradient_points)) == expected)
        check('MNAR_all_estimated', m.status.eq('estimated').all())
        check('MNAR_counts', (m.n_observed + m.n_missing == m.n_target).all()
              and m.n_target.nunique() == 1 and m.n_observed.nunique() == 1 and m.n_missing.nunique() == 1)
        diagnostic = truth(m.diagnostic_only)
        expected_m = np.where(diagnostic, 5, spec['imputations'])
        check('MNAR_complete_imputations', (m.requested_imputations == expected_m).all()
              and (m.valid_imputations == m.requested_imputations).all())
        check('MNAR_ranges', m.tau_iqr.gt(0).all() and m.clipped_fraction.between(0,1).all())
        good = m.loc[m.status.eq('estimated')]
        if len(good):
            check('MNAR_pooled_count', (good.imputations == good.requested_imputations).all())
            total = good.within_variance + (1+1/good.imputations)*good.between_variance
            check('MNAR_variances_and_df', good.within_variance.ge(0).all() and good.between_variance.ge(0).all()
                  and total.gt(0).all() and good.df.gt(0).all())
            close('MNAR_pooled_SE', good.se, np.sqrt(total))
            close('MNAR_variance_fraction', good.missing_variance_fraction, (1+1/good.imputations)*good.between_variance/total)
            close('MNAR_MCSE', good.monte_carlo_se, np.sqrt(good.between_variance/good.imputations))
            crit = student_t.ppf(.975, good.df)
            close('MNAR_CI_low', good.ci_low, good.estimate-crit*good.se)
            close('MNAR_CI_high', good.ci_high, good.estimate+crit*good.se)
            close('MNAR_P', good.nominal_p, 2*student_t.sf(np.abs(good.estimate/good.se),good.df))
        failures = pd.read_csv(final/'MNAR_fit_failures.csv')
        check('MNAR_no_fit_failures', failures.empty)
        f = pd.read_csv(final/'MNAR_sample_flow.csv')
        check('MNAR_unique_flow_stages', not f.stage.duplicated().any())
        flow = f.set_index('stage')['n'].to_dict()
        close('MNAR_target_count', [flow['selected_survivor_target']], [m.n_target.iloc[0]])
        close('MNAR_observed_count', [flow['observed_two_year_endpoint']], [m.n_observed.iloc[0]])
        close('MNAR_missing_count', [flow['missing_two_year_endpoint']], [m.n_missing.iloc[0]])
        survival = ['documented_alive','documented_death','documented_death_interval_before_window_end',
                    'documented_death_by_window_end','unknown','unknown_death_timing','conflicting_evidence']
        close('MNAR_survival_partition', [sum(flow.get(k,0) for k in survival)], [flow['full_administrative_opportunity']])
        close('MNAR_survivor_target', [flow.get('documented_alive',0)], [flow['selected_survivor_target']])
        if (final/'Survival_classification_changes.csv').is_file():
            changes = pd.read_csv(final/'Survival_classification_changes.csv')
            c = changes.loc[truth(changes.administrative_opportunity)].groupby('updated')['n'].sum().to_dict()
            check('MNAR_transition_partition', all(c.get(k,0) == flow.get(k,0) for k in survival))
        if diagnostic.any(): notes.append('Diagnostic imputations are present; do not use as final inference.')

    tpath = final/'GSE174367_target_results.csv'
    if tpath.is_file():
        tar = pd.read_csv(tpath); go = pd.read_csv(final/'GSE174367_GO_all_models.csv',low_memory=False)
        status = pd.read_csv(final/'GSE174367_model_status.csv')
        planned = {(c,m) for c in CELL_NAMES for m in MODELS}
        check('single_nucleus_model_grid', len(status) == len(planned) and set(zip(status.cell_type,status.model)) == planned)
        check('single_nucleus_donor_counts', (status.n_low+status.n_high == status.n_donors).all())
        estimated = status.loc[status.status.eq('estimated')]
        check('single_nucleus_group_thresholds', estimated.n_low.ge(6).all() and estimated.n_high.ge(6).all())
        check('target_unique_keys', not tar.duplicated(['cell_type','model','gene_symbol']).any())
        check('GO_unique_keys', not go.duplicated(['cell_type','model','go_id']).any())
        check('target_fixed_family', len(tar) == 224)
        check('GO_fixed_family', len(go) == 67564)
        go_sets = [set(g.go_id) for _,g in go.groupby(['cell_type','model'])]
        check('GO_same_retained_family', len(go_sets)==14 and all(x==go_sets[0] for x in go_sets))
        for c,m in sorted(planned):
            sub = tar.loc[tar.cell_type.eq(c) & tar.model.eq(m)]
            check(f'target_16_{c}_{m}', len(sub) == 16 and set(sub.gene_symbol) == set(TARGETS))
            subgo = go.loc[go.cell_type.eq(c) & go.model.eq(m)]
            check(f'GO_4826_{c}_{m}', len(subgo) == 4826 and subgo.go_id.nunique() == 4826)
            qcheck(f'GO_BH_{c}_{m}', subgo, 'q_fixed_GO_family')
            parent = status.loc[status.cell_type.eq(c) & status.model.eq(m)]
            if len(parent) == 1 and parent.status.iloc[0] != 'estimated':
                check(f'no_tests_without_model_{c}_{m}', not sub.status.eq('tested').any() and not subgo.status.eq('tested').any())
        tests(tar,'target'); tests(go,'GO')
        qcheck('target_BH_224',tar,'q_all_planned_targets')
        prim = tar.cell_type.eq('MG') & tar.model.eq('age_sex_batch_RIN')
        qcheck('target_primary_BH_16',tar.loc[prim],'q_primary_MG')
        check('primary_q_only_in_MG',tar.loc[~prim,'q_primary_MG'].isna().all())
        qcheck('GO_BH_67564',go,'q_all_cell_models')
        tested = tar.loc[tar.status.eq('tested')]
        check('target_direction_labels', np.array_equal(truth(tested.direction_match), tested.logFC.lt(0).to_numpy()))
        check('GO_tested_size',go.loc[go.status.eq('tested'),'n_genes'].between(10,300).all())
        d = pd.read_csv(final/'GSE174367_donors.csv'); cells = pd.read_csv(final/'GSE174367_donor_cell_counts.csv')
        check('donor_identifiers', not d.donor_id.duplicated().any() and not cells.donor_id.duplicated().any()
              and set(d.donor_id) == set(cells.donor_id))
        check('cell_counts_nonnegative_integer', (cells[list(CELL_NAMES)].ge(0) & cells[list(CELL_NAMES)].mod(1).eq(0)).all().all())
        dims = pd.read_csv(final/'GSE174367_model_dimensions.csv')
        for r in dims.itertuples():
            g = go.loc[go.cell_type.eq(r.cell_type) & go.model.eq(r.model)]
            check(f'GO_test_count_{r.cell_type}_{r.model}', r.GO_tested == g.status.eq('tested').sum())
        if status.status.ne('estimated').any(): notes.append('Some single-nucleus models are not estimable; retain their reasons and empty estimates.')
        notes.append('The primary target panel is retained regardless of significance; GO follow-ups do not replace it.')

    interaction_status = final/'GSE174367_interaction_status.csv'
    if interaction_status.is_file():
        s = pd.read_csv(interaction_status)
        check('interaction_two_models',len(s)==2 and set(s.model)==set(MODELS))
        check('interaction_paired_counts',(s.n_low+s.n_high==s.n_donors).all()
              and (s.n_pseudobulks==2*s.n_donors).all())
        good=s.loc[s.status.eq('estimated')]
        check('interaction_eligibility',good.n_low.ge(6).all() and good.n_high.ge(6).all()
              and good.independent_design_df.ge(5).all() and good.block_correlation.between(-1,1,inclusive='neither').all())
    for table, size, qcol, estimate in [
        ('GSE174367_interaction_targets.csv',16,'q_fixed_targets','interaction_logFC'),
        ('GSE174367_interaction_GO.csv',4826,'q_fixed_GO_family','interaction_logFC'),
    ]:
        if not (final/table).is_file(): continue
        frame = pd.read_csv(final/table)
        check(table+'_family_size',len(frame)==2*size)
        key = 'gene_symbol' if size==16 else 'go_id'
        check(table+'_unique_keys',not frame.duplicated(['model',key]).any())
        for model in MODELS:
            sub = frame.loc[frame.model.eq(model)]
            check(table+'_'+model+'_family',len(sub)==size)
            qcheck(table+'_'+model+'_BH',sub,qcol)
        qcheck(table+'_global_BH',frame,'q_all_interaction_models')
        tests(frame,table,estimate)
        notes.append('OPC–INH interaction is an exploratory follow-up selected after inspection of the preceding results.')
    if not mpath.is_file() and not tpath.is_file():
        check('at_least_one_analysis',False,'No recognized result module found.')
    return dict(passed=all(x['passed'] for x in checks),checks=checks,notes=list(dict.fromkeys(notes)))


def figures(final, destination=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder = Path(destination) if destination else final/'Figures'; folder.mkdir(parents=True,exist_ok=True)
    captions = []
    def save(fig,name):
        name={'Figure_MNAR':'FigureS4','Figure_GSE174367_MG':'FigureS19'}[name]
        for ext in ['pdf','png','svg']:
            fig.savefig(folder/(name+'.'+ext),dpi=300,bbox_inches='tight')
        plt.close(fig)
    if (final/'GSE174367_target_results.csv').is_file():
        d = pd.read_csv(final/'GSE174367_target_results.csv')
        d = d.loc[d.cell_type.eq('MG') & d.model.eq('age_sex_batch_RIN')].sort_values('gene_symbol').reset_index(drop=True)
        tested = d.status.eq('tested')
        fig,ax = plt.subplots(figsize=(9,7.5))
        for i,row in d.iterrows():
            if tested.iloc[i]:
                ax.errorbar(row.logFC,i,xerr=[[row.logFC-row.CI_L],[row.CI_R-row.logFC]],fmt='o',
                            color='#235789',capsize=3,markersize=4)
            else:
                ax.text(.03,i,'Not tested',transform=ax.get_yaxis_transform(),ha='left',va='center',color='.45',fontsize=9)
            label = f'{row.q_primary_MG:.3g}' if tested.iloc[i] else '—'
            ax.text(1.02,i,label,transform=ax.get_yaxis_transform(),va='center',fontsize=9)
        ax.text(1.02,1.025,'q (16 targets)',transform=ax.transAxes,fontsize=9)
        ax.axvline(0,color='.6',lw=.8); ax.set_yticks(range(len(d)),d.gene_symbol); ax.invert_yaxis()
        ax.set_xlabel('Braak V–VI minus I–II: log2 fold change (marginal 95% CI)')
        ax.set_title('GSE174367: donor-level microglial pseudobulk')
        ntested = int(tested.sum()); nsig = int((tested & d.q_primary_MG.lt(.05)).sum())
        note = f'{ntested}/{len(d)} targets tested; {nsig} with primary-family q < 0.05. Untested targets have no effect estimate.'
        fig.text(.02,.02,note,fontsize=9); fig.tight_layout(rect=[0,.055,.87,1]); save(fig,'Figure_GSE174367_MG')
        captions.append('Figure_GSE174367_MG. Fixed 16-gene target panel in donor-level microglial pseudobulk. '
            'The model adjusts age, sex, experimental batch and RIN. Points and bars show log2 fold changes '
            'and marginal 95% confidence intervals. The displayed q values use the fixed 16-target BH family, '
            'with untested hypotheses counted as P=1 only for adjustment. Not-tested rows are not zero effects. '
            + note + ' Targets absent from the fitted background may be low-expression, absent or ambiguously mapped; '
            'a legacy combined status alone does not distinguish these causes.')
    if (final/'MNAR_scenarios.csv').is_file():
        allrows = pd.read_csv(final/'MNAR_scenarios.csv')
        d = allrows.loc[allrows.status.eq('estimated')]
        if len(d):
            groups = list(d.groupby('tau_gradient_points',sort=True))
            fig,axes = plt.subplots(1,len(groups),figsize=(4*len(groups),4.6),sharey=True,squeeze=False)
            for ax,(gamma,g) in zip(axes.ravel(),groups):
                g = g.sort_values('delta_points'); ax.axhline(0,color='.6',lw=.8)
                ax.errorbar(g.delta_points,g.estimate,yerr=[g.estimate-g.ci_low,g.ci_high-g.estimate],fmt='o-',capsize=3)
                ax.set_title(f'Tau-gradient shift γ = {gamma:g}'); ax.set_xlabel('Overall shift δ (CDR-SB points)')
            axes[0,0].set_ylabel('Two-year CDR-SB contrast per 0.1 tau SUVR\n(marginal 95% CI)')
            fig.suptitle('Endpoint sensitivity among documented survivors',fontsize=12)
            foot = ('Only missing outcomes shifted by δ + γ × (tau − median)/IQR, then bounded to [0,18].\n'
                    'Time origin: selected baseline CDR date; outcome window: 2 years ±180 days. Not an annual slope.')
            fig.text(.02,.015,foot,fontsize=9);fig.tight_layout(rect=[0,.14,1,.94]);save(fig,'Figure_MNAR')
            n = int(d.n_target.iloc[0]); observed = int(d.n_observed.iloc[0]); missing = int(d.n_missing.iloc[0])
            captions.append('Figure_MNAR. Two-year CDR-SB endpoint sensitivity in '+str(n)+' A+ participants with '
                f'complete baseline predictors, full administrative opportunity and documented survival through the outcome-window end '
                f'({observed} observed and {missing} missing endpoints). The clock begins at the selected baseline CDR assessment. '
                'Only missing outcomes receive δ + γ × (tau − target median)/target IQR, followed by truncation to [0,18]. '
                'The scenario table records clipped_fraction and mean_actual_shift; nominal shifts need not equal actual shifts. '
                'Contrasts across the target tau IQR are normalized per 0.1 SUVR. Bars are marginal 95% intervals from Rubin pooling. '
                'The parameters are sensitivity assumptions, not estimated missingness mechanisms. Unknown survival and deaths are '
                'excluded and counted; this does not replace the repeated-measures trajectory analysis.')
    return captions

