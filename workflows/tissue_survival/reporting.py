"""Final survivor and single-nucleus figures."""
from pathlib import Path
import pandas as pd


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

