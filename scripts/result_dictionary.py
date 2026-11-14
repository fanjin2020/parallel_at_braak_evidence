"""Definitions whose meaning depends on the exported source record."""

COMMON_PERIODS = {'SD040', 'FIG_Figure4_support', 'FIG_FigureS20_periods'}
PERIOD_SCOPES = COMMON_PERIODS | {'SD135', 'FIG_FigureS9_support', 'FIG_FigureS10_period_counts'}
CLOCK = ['Time origin used to calculate elapsed follow-up', 'category',
         'baseline_CDR: selected baseline CDR date; tau_PET: aligned tau PET date; elapsed days divided by 365.25']
GROUPING = ['Variable used to stratify the follow-up summary', 'category',
            'all; diagnosis; tau_group. Strata from different grouping variables overlap']
COMMON_GROUP = ['Stratum in the common complete-covariate longitudinal sample', 'category',
                'all; CN/MCI/Dementia; low tau <= sample Q25, high tau >= sample Q75, middle otherwise; quantiles calculated in the common sample']
MODEL = ['Named adjustment or outcome-model specification', 'label',
         'Labels identify fitted specifications; they do not imply prespecification. See the model formula and manuscript Methods']


def result_definition(scope, field):
    """Return a scoped definition, or None without inferring missing semantics."""
    if scope == 'SD037':
        meanings = {
            'grouping': GROUPING,
            'group': COMMON_GROUP,
            'clock': CLOCK,
            'scope': ['Follow-up records used to derive each participant\'s last observed time', 'category',
                      'all_available: all eligible later CDR records; analysis_3years: eligible records within three years after tau PET. The baseline-CDR span can exceed three years'],
            'n': ['Number of distinct participants with at least one eligible follow-up in this stratum and scope', 'participants',
                  'Each participant contributes one last-follow-up time, not one row per visit'],
            'median': ['Median of participant-level last-follow-up times', 'years', 'Calculated within grouping, group, scope and clock'],
            'q25': ['25th percentile of participant-level last-follow-up times', 'years', 'Calculated within grouping, group, scope and clock'],
            'q75': ['75th percentile of participant-level last-follow-up times', 'years', 'Calculated within grouping, group, scope and clock'],
            'maximum': ['Maximum participant-level last-follow-up time', 'years', 'Calculated within grouping, group, scope and clock'],
        }
        return meanings.get(field)
    if scope in PERIOD_SCOPES:
        if field == 'period':
            if scope == 'SD135':
                return ['CDR assessment interval since selected baseline CDR in the clinical-stage sensitivity models', 'years',
                        'baseline: t=0; (0,1]y: 0<t<=1; (1,2]y: 1<t<=2; >2y: t>2. Follow-up is limited to three years after tau PET, so the baseline-CDR span can exceed three years']
            if scope == 'FIG_FigureS9_support':
                return ['CDR assessment interval since selected baseline CDR in the larger longitudinal sensitivity sample', 'years',
                        'baseline: t=0; 0_to_1: 0<t<=1; 1_to_2: 1<t<=2; over_2: t>2. Includes baseline; follow-up is limited to three years after tau PET']
            population = 'common complete-covariate longitudinal sample' if scope in COMMON_PERIODS else 'larger longitudinal sensitivity sample'
            return ['Observed follow-up interval in the '+population+', using the time origin in clock', 'years',
                    '0_to_1: 0<t<=1; 1_to_2: 1<t<=2; over_2: t>2. Baseline excluded; records restricted to three years after tau PET. With baseline_CDR, over_2 is not a fixed 2-3-year interval']
        if field == 'clock': return CLOCK
        if field == 'grouping': return GROUPING
        if field == 'group' and scope in COMMON_PERIODS: return COMMON_GROUP
        if field == 'people':
            return ['Distinct participants with an assessment in the indicated interval and stratum', 'participants',
                    'A participant may occur in several intervals; interval counts must not be added to infer total unique participants']
        if field == 'visits':
            return ['Number of CDR assessments in the indicated interval and stratum', 'assessments',
                    'Multiple assessments from one participant are counted; baseline is included only when period is explicitly baseline']
    if scope == 'QC_091' and field == 'log2FC':
        return ['Adjusted GSE44770 AD-minus-control expression contrast', 'log2 expression',
                'AD minus control; adjusted for age, sex, PMI, tissue pH, RIN and batch. Not a Braak V-VI-minus-0 contrast']
    if scope == 'FIG_FigureS15_source':
        return {
            'log2FC': ['Adjusted expression contrast for the cohort and model identified in the same row', 'log2 expression',
                       'GSE44770 / AD_vs_control: AD minus control, adjusted for age, sex, PMI, pH, RIN and batch. GSE106241: Braak V-VI minus 0, adjusted for age, sex and APOE, with categorical RNA quality added in age_sex_APOE_RNAquality'],
            'cohort': ['Independent public-tissue dataset supplying the estimate', 'GEO accession', 'GSE44770 or GSE106241; contrasts differ by cohort'],
            'model': ['External-cohort contrast and covariate specification', 'label',
                      'AD_vs_control: GSE44770 AD-control; age_sex_APOE and age_sex_APOE_RNAquality: GSE106241 Braak V-VI-0, without/with categorical RNA-quality adjustment'],
        }.get(field)
    if scope == 'SD112' and field == 'model':
        return ['Longitudinal outcome family and cumulative adjustment set', 'label',
                'Gaussian or fractional-logit; demographic, amyloid, stage, MRI. The fully adjusted fractional-logit model was designated the primary exploratory model post hoc; all eight marginal models use the same participants']
    return None


def refresh_rows(rows):
    """Update verified definitions without adding variables or changing field keys."""
    changed = []
    for row in rows:
        meaning = result_definition(row['scope'], row['field']) if row['role'] == 'exported_result' else None
        if meaning is None and row['role'] == 'exported_result' and row['field'] == 'model' and 'Labels identify planned models' in row['coding_or_missing']:
            meaning = MODEL
        if meaning:
            old = {k: row[k] for k in ['definition','unit','coding_or_missing','definition_status']}
            row.update(definition=meaning[0], unit=meaning[1], coding_or_missing=meaning[2], definition_status='documented')
            if any(row[k] != old[k] for k in old): changed.append({'scope': row['scope'], 'field': row['field']})
    return changed
