"""The six clinical-support panels used in the final numbered figure set."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from figure_primitives import at_display_label

COLORS=['#28678e','#b55a38','#749079','#8f6da5']


def read(root, folder, name):
    return pd.read_csv(Path(root)/folder/(name+'.csv'))


def forest(ax,df,labels,estimate='estimate',low='ci_low',high='ci_high',color=COLORS[0]):
    y=np.arange(len(df));values=df[estimate].to_numpy(float)
    ax.errorbar(values,y,xerr=np.vstack([values-df[low],df[high]-values]),fmt='o',color=color,capsize=2,ms=3,lw=1)
    ax.axvline(0,color='.55',lw=.7);ax.set_yticks(y,[at_display_label(v) for v in labels]);ax.invert_yaxis();ax.margins(y=.10)


def followup_support(root,save):
    d=read(root,'longitudinal','monthly_visit_support')
    if 'clock' not in d:raise ValueError('Follow-up descriptions require explicit clocks')
    d=d.loc[d.clock.eq('baseline_CDR')].copy()
    fig,axes=plt.subplots(2,1,figsize=(7.2,4.6),sharex=True,layout='constrained')
    for ax,grouping,title in zip(axes,['diagnosis','tau_group'],['a  By baseline diagnosis','b  By baseline tau burden']):
        s=d.loc[d.grouping.eq(grouping)];n_months=int(s.month.max());bottom=np.zeros(n_months)
        for i,(group,g) in enumerate(s.groupby('group')):
            values=g.set_index('month').reindex(range(1,n_months+1)).visits.fillna(0).to_numpy()
            ax.bar((np.arange(1,n_months+1)-.5)/12,values,width=.075,bottom=bottom,color=COLORS[i%4],label=group)
            bottom+=values
        ax.axvspan(0,1,color='#f0d8b5',alpha=.35,zorder=-1);ax.set_ylabel('Observed visits')
        ax.set_ylim(0,max(1,bottom.max())*1.2)
        ax.set_title(title,loc='left',fontweight='bold');ax.legend(frameon=False,ncol=4,fontsize=7)
    axes[1].set_xlabel('Years since selected baseline CDR');axes[1].set_xlim(0,n_months/12+.04)
    periods=read(root,'longitudinal','observed_visit_periods')
    first=periods.loc[periods.grouping.eq('all')&periods.clock.eq('baseline_CDR')&periods.period.eq('0_to_1')]
    if len(first)!=1:raise ValueError('Ambiguous first-year support')
    first=first.iloc[0]
    axes[0].text(.04,.88,f'First year\n{int(first.people)} participants, {int(first.visits)} visits',transform=axes[0].transAxes,fontsize=8,va='top')
    save(fig,'FigureS10',{'source':d,'period_counts':periods},'Show the empirical visit distribution without smoothing sparse early observations.',
         'Monthly counts of actual post-baseline assessments, within the PET+3-year analysis horizon. Clock is selected baseline CDR, matching the trajectory model; assessments can extend beyond year 3 on this axis when baseline CDR precedes PET. A person can contribute in more than one month. Shaded region is the first year; bars are observed visits. Low/high tau are <=Q25/>=Q75 of the baseline-eligible sample. PET-clock descriptions, unique participant counts by period, near-year windows and durations are separate tables.')


def regional_tau(root,save):
    d=read(root,'clinical','cross_sectional_results')
    s=d.loc[d.family.eq('regional_tau')].copy()
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.4),layout='constrained')
    names={'tau_temporal_suvr_per01':'Temporal meta-ROI','tau_pvc_per01':'Temporal meta-ROI, PVC','tau_entorhinal_per01':'Entorhinal','tau_mtl_per01':'Medial temporal','tau_temporal_neo_per01':'Temporal neocortex','tau_neocortex_per01':'Neocortex excluding MTL'}
    for ax,y,letter in zip(axes,['hippocampus_icv','cdrsb'],['a','b']):
        z=s.loc[s.outcome.eq(y)]
        forest(ax,z,[names[x] for x in z.analysis],estimate='beta_outcome_sd',low='ci_low_outcome_sd',high='ci_high_outcome_sd')
        ax.set_title(letter+'  '+('Hippocampus/ICV' if y=='hippocampus_icv' else 'CDR-SB'),loc='left',fontweight='bold')
        ax.set_xlabel('Outcome SD per 0.1 SUVR')
    save(fig,'FigureS11',{'source':s},'Compare continuous tau definitions and partial-volume correction.',
         'A+ complete cases shared by every ROI within each outcome, adjusted for age, sex, education, APOE and Centiloid. HC3 95% CIs; BH q values across the named regional-tau family are in source data. ROI scales and outcome samples differ; standardized outcome SD is sample-specific. The regional measures are research definitions, without external biological-stage calibration.')


def tau_distribution(root,save):
    d=read(root,'clinical','tau_threshold_histogram')
    fig,ax=plt.subplots(figsize=(7.2,2.8),layout='constrained')
    for i,(group,g) in enumerate(d.groupby('at_stage')):
        ax.step((g.lower+g.upper)/2,g.n,where='mid',label=at_display_label(group),color=COLORS[i])
    ax.axvline(1.34,color='.35',linestyle='--',lw=1,label='Operational cutoff 1.34')
    ax.set_xlabel('Temporal meta-region flortaucipir SUVR');ax.set_ylabel('Participants');ax.legend(frameon=False,fontsize=7,ncol=3)
    save(fig,'FigureS12',{'source':d},'Describe observed tau values around the classification cutoff.',
         'Fixed 0.05-SUVR bins covering all finite values. Groups use PET-Core amyloid and temporal tau classification. Counts within 0.02, 0.05 and 0.10 of the operational 1.34 cutoff are reported separately.')


def visit_process(root,save):
    d=read(root,'visit_process','visit_process_results')
    fig,ax=plt.subplots(figsize=(7.2,3.0),layout='constrained')
    for i,(model,label,color) in enumerate([('unweighted_windows','Unweighted',COLORS[0]),('history_weighted_windows','History weighted',COLORS[1])]):
        s=d.loc[d.model.eq(model)].sort_values('year')
        ax.errorbar(s.year+(i-.5)*.08,s.estimate,yerr=np.vstack([s.estimate-s.ci_low,s.ci_high-s.estimate]),fmt='o',ms=4,capsize=3,color=color,label=label)
    support=d.loc[d.model.eq('unweighted_windows')].set_index('year')
    ax.set_xticks([1,2,3],[f'Year {y}\n{int(support.loc[y,"n_observed"])} observed' for y in [1,2,3]])
    ax.axhline(0,color='.5',lw=.7);ax.set_ylabel('Conditional CDR-SB level contrast\nper 0.1 SUVR normalized IQR');ax.legend(frameon=False)
    ax.set_title(f'Target population: {int(support.n_target.iloc[0])} participants',loc='left',fontsize=9)
    save(fig,'FigureS16',{'source':d},'Sensitivity to measured time-varying observation history.',
         'Fixed non-overlapping annual windows relative to tau PET, +/-6 months, requiring administrative opportunity for all three windows. Same target population before/after weighting; one nearest observation per window. Fractional-logit model conditional on baseline CDR, diagnosis, demographics and amyloid. Bars are percentile participant-bootstrap CIs refitting observation and outcome models; requested/valid replicate counts and observed counts are in source data. The year-3 window extends to 3.5 years.')


def regional_staging(root,save):
    d=read(root,'clinical','biological_stage_counts')
    counts=d.groupby('biological_stage')['n'].sum();order=[k for k in ['amyloid_negative','A','B','C','D','discordant_neocortical_without_MTL','Unavailable'] if k in counts]
    labels={'amyloid_negative':'Amyloid negative','discordant_neocortical_without_MTL':'Discordant pattern','Unavailable':'Unclassified'}
    fig,ax=plt.subplots(figsize=(7.2,2.6),layout='constrained');ax.bar(range(len(order)),counts.reindex(order),color=COLORS[0]);ax.set_xticks(range(len(order)),[labels.get(k,'Stage '+k) for k in order],rotation=15,ha='right');ax.set_ylabel('Participants')
    for i,k in enumerate(order):ax.text(i,counts[k]+3,str(counts[k]),ha='center',fontsize=7)
    ax.margins(y=.15)
    save(fig,'FigureS17',{'source':d},'Operational-stage assignments, including discordant and unclassified participants.',
         'Exploratory AA-inspired PET staging using outcome-independent two-component mixture cutoffs estimated in the Berkeley processing pipeline. External stage calibration is unavailable. MTL: bilateral entorhinal; neocortex: volume-weighted CTX parcels excluding entorhinal and parahippocampal cortex. Threshold bootstrap diagnostics and all stage counts are reported separately.')


def stage_uncertainty(root,save):
    d=read(root,'clinical','stage_propagated_intervals')
    if d.empty or not {'ci_low','ci_high'}.issubset(d):
        raise ValueError('Stage uncertainty results are unavailable')
    outcomes=list(d.outcome.unique())
    fig,axes=plt.subplots(1,len(outcomes),figsize=(7.2,3.8),squeeze=False)
    for i,(ax,outcome) in enumerate(zip(axes.flat,outcomes)):
        s=d.loc[d.outcome.eq(outcome)].dropna(subset=['estimate','ci_low','ci_high'])
        y=np.arange(len(s));scale=1e4 if outcome=='hippocampus_icv' else 1
        ax.hlines(y,s.ci_low*scale,s.ci_high*scale,color='#447c9a')
        ax.plot(s.estimate*scale,y,'o',color='#447c9a')
        ax.set_yticks(y,s.contrast.str.replace('_vs_',' vs ',regex=False));ax.invert_yaxis()
        ax.set_title(chr(97+i)+'  '+{'cdrsb':'CDR-SB','hippocampus_icv':'Hippocampus/ICV'}.get(outcome,outcome),loc='left')
        ax.axvline(0,color='#777777',lw=.7)
        ax.set_xlabel('Adjusted stage difference'+(r' ($\times 10^{-4}$)' if scale!=1 else ''))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    fig.tight_layout()
    save(fig,'FigureS19',{'source':d},'Uncertainty in operational regional-stage contrasts.',
         'Exploratory regional-stage differences. Each participant-bootstrap draw re-estimates PET cutoffs, assigns stages and fits the outcome model. Bars show percentile 95% intervals; intervals are withheld when fewer than 90% of requested refits succeed. External clinical-stage validation is unavailable.')
