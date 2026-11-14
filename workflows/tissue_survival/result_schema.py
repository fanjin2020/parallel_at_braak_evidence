"""Definitions for exported aggregate results; unregistered fields are errors."""

CELL_NAMES = {
    'MG': 'Microglia', 'ASC': 'Astrocytes', 'ODC': 'Oligodendrocytes',
    'OPC': 'Oligodendrocyte precursor cells', 'EX': 'Excitatory neurons',
    'INH': 'Inhibitory neurons', 'PER.END': 'Pericyte/endothelial class',
}

# Each entry has a definition, unit, and interpretation of missing values.
FIELDS = {
    'table': ('Exported table or authorized source-table name.', 'text', 'Required.'),
    'file': ('Repository-relative file name; no participant identifier.', 'path', 'Required.'),
    'sha256': ('SHA-256 content digest of the named file.', 'hexadecimal digest', 'Required.'),
    'scope': ('Whether a record describes executed reporting code or an available analysis-code snapshot.', 'category', 'Required.'),
    'rows': ('Number of records read from this source table; zero when unavailable.', 'records', 'Required.'),
    'package': ('R package from which the authorized source table was exported.', 'text', 'Package provenance not supplied.'),
    'version': ('Version of the source package recorded during export.', 'version string', 'Package provenance not supplied.'),
    'reason': ('Recorded exclusion, non-estimability, failure, or source disposition explanation; interpret within the named table.', 'text', 'No reason needed for successful records.'),
    'n': ('Number of records or participants in the named stage/category; see table-specific definition.', 'count', 'Required.'),
    'n_participants': ('Distinct participants with the named disposition reason; a participant may have more than one reason.', 'participants', 'Required.'),
    'stage': ('Cohort-selection step or mutually exclusive survival category within the administrative-opportunity cohort.', 'category', 'Required.'),
    'previous': ('Survival classification using the retained clinical-visit evidence.', 'category', 'Required.'),
    'updated': ('Survival classification after adding authorized event-table evidence.', 'category', 'Required.'),
    'administrative_opportunity': ('Whether the full two-year outcome window ended by the fixed data cutoff.', 'boolean', 'Required.'),
    'donor_id': ('Public GSE174367 donor identifier; the independent unit for single-nucleus inference.', 'identifier', 'Required.'),
    'age': ('Age of the GSE174367 donor as recorded in the public metadata.', 'years', 'Unavailable in public metadata; excluded if required by the model.'),
    'Sex': ('Sex label in the public donor metadata: F or M.', 'category', 'Unavailable in public metadata.'),
    'sex_male': ('Male-sex indicator derived from the public metadata: M=1, F=0.', '0/1', 'Unavailable in public metadata.'),
    'batch': ('Public GSE174367 experimental batch label, modeled as a factor.', 'category', 'Unavailable in public metadata.'),
    'rin': ('Public GSE174367 RNA integrity number; it does not supply RIN for GSE131617.', 'RIN score', 'Unavailable in public metadata.'),
    'pmi': ('Public GSE174367 postmortem interval; required only in the PMI sensitivity model.', 'hours', 'Donor omitted from the PMI-complete sensitivity model.'),
    'Tangle.Stage': ('Original public donor tangle-stage label.', 'category', 'Unavailable in public metadata.'),
    'braak': ('Numeric donor Braak stage parsed from Tangle.Stage.', 'stage 0–6', 'Unavailable in public metadata.'),
    'cell_type': ('Author-assigned GSE174367 cell class; codes are defined in the cell-count table.', 'category', 'Required.'),
    'model': ('age_sex_batch_RIN adjusts age, sex, batch and RIN; plus_PMI also adjusts PMI in complete donors.', 'category', 'Required.'),
    'status': ('Analysis state, not a significance flag; see Result_status_codes.csv.', 'category', 'Required.'),
    'n_low': ('Independent eligible donors at Braak I–II in this model.', 'donors', 'Required.'),
    'n_high': ('Independent eligible donors at Braak V–VI in this model.', 'donors', 'Required.'),
    'n_donors': ('Independent eligible donors used or evaluated for this model; not the number of nuclei.', 'donors', 'Required.'),
    'n_stage_0': ('Number of Braak-0 donors in the independent dataset.', 'donors', 'Required.'),
    'contrast': ('Explicit stage or interaction comparison being tested.', 'text', 'Required.'),
    'formula': ('Covariate design used or evaluated for the model.', 'formula', 'Design could not be constructed.'),
    'primary': ('Whether this is the predeclared primary MG model.', 'boolean', 'Required.'),
    'genes_tested': ('Expression features retained in the genome-wide model after expression filtering.', 'features', 'Model was not estimable.'),
    'genomewide_features': ('Expression features actually fitted in the full expression background.', 'features', 'Model was not estimable.'),
    'unique_GO_mapped_features': ('Fitted features with unambiguous one-to-one Entrez mapping used for GO ranking.', 'features', 'Model was not estimable.'),
    'GO_tested': ('GO sets with 10–300 measured members that were actually tested; the BH family still includes all 4,826 sets.', 'sets', 'Model was not estimable.'),
    'gene_symbol': ('Public gene symbol; ambiguous symbols are not selected for target tests.', 'symbol', 'No unambiguous symbol.'),
    'logFC': ('Adjusted log2 expression fold change: Braak V–VI minus I–II.', 'log2 fold change', 'Target was not tested or model was not estimable.'),
    'CI_L': ('Lower marginal 95% confidence limit for the fitted log2 contrast; not a simultaneous/FDR interval.', 'log2 units', 'Target was not tested.'),
    'CI_R': ('Upper marginal 95% confidence limit for the fitted log2 contrast; not a simultaneous/FDR interval.', 'log2 units', 'Target was not tested.'),
    'PValue': ('Unadjusted test P value; missing tests enter only the fixed BH calculation as P=1.', 'probability', 'Hypothesis was not tested; not an observed P=1.'),
    'FDR_genomewide': ('BH-adjusted P value across all expression features fitted in this model.', 'adjusted probability', 'Target was not tested.'),
    'q_primary_MG': ('BH adjustment over the 16 fixed MG primary targets; untested hypotheses count as P=1.', 'adjusted probability', 'Not the primary MG model.'),
    'q_all_planned_targets': ('BH adjustment over 16 targets × 7 cell classes × 2 models (224); untested hypotheses count as P=1.', 'adjusted probability', 'Required for planned target rows.'),
    'discovery_direction': ('Direction of the corresponding bulk target result; down for the fixed panel.', 'category', 'Discovery direction unavailable.'),
    'direction_match': ('For tested targets, whether the effect sign agrees with discovery; this does not establish significance.', 'boolean', 'Target was not tested.'),
    'go_id': ('GO identifier from the retained 4,826-set annotation.', 'identifier', 'Required.'),
    'n_genes': ('Uniquely mapped measured members of the GO set; zero for models not estimated.', 'genes', 'Required.'),
    'direction': ('Direction of the competitive gene-set result from cameraPR.', 'Up/Down', 'GO set or model was not tested.'),
    'q_fixed_GO_family': ('BH adjustment over all 4,826 prespecified GO sets within the model; untested P=1.', 'adjusted probability', 'Required for planned GO rows.'),
    'q_all_cell_models': ('BH adjustment over 4,826 GO sets × 7 cell classes × 2 models (67,564); untested P=1.', 'adjusted probability', 'Required for planned GO rows.'),
    'delta_points': ('Assumed overall shift added only to missing two-year CDR-SB outcomes before truncation.', 'CDR-SB points', 'Required.'),
    'tau_gradient_points': ('Assumed extra missing-outcome shift per target tau IQR: gamma × (tau−median)/IQR, before truncation.', 'CDR-SB points per tau IQR', 'Required.'),
    'n_target': ('Documented-survivor endpoint target with complete baseline predictors and full administrative opportunity.', 'participants', 'Required.'),
    'n_observed': ('Target participants with an observed CDR-SB endpoint in the two-year ±180-day window.', 'participants', 'Required.'),
    'n_missing': ('Target participants whose two-year endpoint is imputed.', 'participants', 'Required.'),
    'tau_iqr': ('75th minus 25th percentile of temporal tau SUVR in the target population.', 'SUVR', 'Required; must be positive.'),
    'requested_imputations': ('Number of imputations planned for this run; diagnostic runs are explicitly flagged.', 'imputations', 'Required.'),
    'valid_imputations': ('Imputations with a valid fitted contrast for this scenario; incomplete sets are not pooled.', 'imputations', 'Required.'),
    'clipped_fraction': ('Fraction of missing-endpoint draws truncated to [0,18], across requested imputations.', 'proportion', 'No missing endpoints or analysis not run.'),
    'mean_actual_shift': ('Mean change from the MAR imputed endpoint after scenario adjustment and [0,18] truncation.', 'CDR-SB points', 'No missing endpoints or analysis not run.'),
    'diagnostic_only': ('True for reduced diagnostic runs; such runs are not final inferential results.', 'boolean', 'Required.'),
    'estimate': ('Adjusted two-year CDR-SB level contrast across target tau IQR, normalized per 0.1 SUVR; not an annual slope.', 'CDR-SB points per 0.1 SUVR', 'Scenario was not successfully pooled.'),
    'se': ('Square root of pooled within- and between-imputation contrast variance.', 'CDR-SB points per 0.1 SUVR', 'Scenario was not successfully pooled.'),
    'ci_low': ('Lower 95% confidence limit from pooled contrast and Barnard–Rubin t degrees of freedom.', 'CDR-SB points per 0.1 SUVR', 'Scenario was not successfully pooled.'),
    'ci_high': ('Upper 95% confidence limit from pooled contrast and Barnard–Rubin t degrees of freedom.', 'CDR-SB points per 0.1 SUVR', 'Scenario was not successfully pooled.'),
    'nominal_p': ('Two-sided nominal P value for the pooled contrast; scenarios are sensitivity assumptions, not a discovery family.', 'probability', 'Scenario was not successfully pooled.'),
    'df': ('Barnard–Rubin finite-sample degrees of freedom for the pooled scalar contrast.', 'degrees of freedom', 'Scenario was not successfully pooled.'),
    'within_variance': ('Mean estimated contrast variance across successful imputations.', '(CDR-SB points per 0.1 SUVR)^2', 'Scenario was not successfully pooled.'),
    'between_variance': ('Sample variance of the estimated contrasts across imputations.', '(CDR-SB points per 0.1 SUVR)^2', 'Scenario was not successfully pooled.'),
    'missing_variance_fraction': ('(1+1/m) times between-imputation variance divided by total pooled variance.', 'proportion', 'Scenario was not successfully pooled.'),
    'monte_carlo_se': ('Monte Carlo standard error of the pooled point estimate: sqrt(between_variance/m).', 'CDR-SB points per 0.1 SUVR', 'Scenario was not successfully pooled.'),
    'imputations': ('Number of successful imputations actually pooled; must equal the requested count.', 'imputations', 'Scenario was not successfully pooled.'),
    'imputation': ('Index of the imputation in which a model failed.', 'index', 'No failure rows when all fits succeed.'),
    'delta': ('Overall shift for a failed MNAR fit; same definition as delta_points.', 'CDR-SB points', 'No failure rows when all fits succeed.'),
    'tau_gradient': ('Tau-dependent shift for a failed MNAR fit; same definition as tau_gradient_points.', 'CDR-SB points per tau IQR', 'No failure rows when all fits succeed.'),
    'PHASE': ('ADNI phase to which the source data-dictionary entry applies.', 'category', 'Not specified in the authorized dictionary.'),
    'TBLNAME': ('Source table name in the authorized ADNI data dictionary.', 'text', 'Not specified in source.'),
    'FLDNAME': ('Source variable name in the authorized ADNI data dictionary.', 'text', 'Not specified in source.'),
    'TEXT': ('Verbatim variable description from the authorized source dictionary.', 'text', 'Not specified in source.'),
    'CODE': ('Verbatim phase-specific source coding; not inferred from an observed numeric value.', 'text', 'No coding list supplied.'),
    'UNITS': ('Units stated in the authorized source dictionary.', 'text', 'No unit supplied; not inferred.'),
    'analysis_group': ('Cell/model combination or endpoint module summarized in this row.', 'text', 'Not applicable.'),
    'n_results': ('Number of result rows with this state in the named table and analysis group.', 'rows', 'Required.'),
    'meaning': ('Meaning of the recorded state, independent of statistical significance.', 'text', 'Required.'),
    'is_tested': ('Whether this state represents an estimated statistical test.', 'boolean', 'Required.'),
    'allowed_inference': ('Interpretation allowed by this result state.', 'text', 'Required.'),
    'pair': ('Fixed exploratory cell comparison: OPC minus INH.', 'category', 'Required.'),
    'interaction_logFC': ('Difference in stage effects: (V–VI minus I–II in OPC) minus (V–VI minus I–II in INH).', 'log2 difference-in-differences', 'Interaction was not tested.'),
    'q_fixed_targets': ('BH over 16 fixed targets within this interaction model; untested P=1.', 'adjusted probability', 'Required for planned targets.'),
    'q_all_interaction_models': ('BH over both interaction models: 32 target hypotheses or 9,652 GO hypotheses, depending on table; untested P=1.', 'adjusted probability', 'Required for planned rows.'),
    'block_correlation': ('Genome-wide consensus within-donor correlation used in the blocked limma fit.', 'correlation', 'Interaction model was not fitted.'),
    'n_pseudobulks': ('Number of donor-by-cell pseudobulk samples in the paired interaction model, two per donor.', 'samples', 'Required.'),
    'independent_design_df': ('Number of paired donors minus the rank of the donor-level covariate design.', 'degrees of freedom', 'Design not constructed.'),
}
for code, name in CELL_NAMES.items():
    FIELDS[code] = (f'Annotated {name.lower()} nuclei contributed by this donor to pseudobulk aggregation.', 'nuclei', 'Required; zero means no annotated nuclei.')

STATES = {
    'tested': ('A valid statistical test was fitted; significance requires the specified q value.', True, 'Report estimate, uncertainty and the declared multiple-testing family.'),
    'estimated': ('The model or scalar contrast was estimated.', True, 'Estimation alone does not imply significance.'),
    'not_estimable': ('Eligibility, sample-size, design, or fit requirements were not met; inspect reason.', False, 'No effect estimate or absence-of-effect conclusion.'),
    'model_not_estimable': ('The parent expression model could not be fitted.', False, 'No gene-set test.'),
    'unmeasured_filtered_or_ambiguous': ('Legacy combined status: target not retained in the fitted expression background, absent, or ambiguously mapped.', False, 'Not tested; this label alone does not prove the gene was absent from the platform.'),
    'absent_from_input': ('Target symbol absent from input gene annotation.', False, 'Not tested.'),
    'ambiguous_gene_symbol': ('More than one input feature has the target symbol.', False, 'Not tested; no arbitrary feature selected.'),
    'filtered_low_expression': ('Unambiguous target was removed by the full-background expression filter.', False, 'Not tested; not a zero effect.'),
    'insufficient_measured_members': ('GO set has fewer than 10 or more than 300 measured, uniquely mapped members.', False, 'Not tested; contributes P=1 only for fixed-family BH bookkeeping.'),
    'failed_imputations': ('At least one requested imputation did not yield a valid contrast; incomplete pooling is prohibited.', False, 'Do not interpret as a completed MNAR scenario.'),
    'not_analyzed': ('The contrast was not part of the fitted analysis.', False, 'No result available.'),
    'read': ('Authorized local source table was read.', False, 'Source inventory only.'),
    'not_available': ('Requested source table was not available locally.', False, 'No data were invented or downloaded.'),
}


def definition(table, field):
    if table.startswith('GSE174367_interaction') and field in ('CI_L', 'CI_R', 'direction'):
        if field == 'direction':
            return ('Direction of enrichment in the OPC-minus-INH stage-interaction ranking; not the direction within one cell type.', 'Up/Down', 'GO test not performed.')
        return (('Lower' if field == 'CI_L' else 'Upper') + ' marginal 95% confidence limit for the log2 difference-in-differences.', 'log2 difference-in-differences', 'Interaction was not tested.')
    if field == 'n':
        if table == 'ADNI_ignored_record_counts.csv':
            return ('Number of source-table records ignored for this recorded reason.', 'source records', 'Required.')
        return ('Number of participants in the named stage or classification-transition cell.', 'participants', 'Required.')
    if field not in FIELDS:
        raise ValueError(f'No field definition registered: {table}:{field}')
    return FIELDS[field]
