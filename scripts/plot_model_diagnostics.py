"""Revision figures from aggregate model outputs; no participant data are exported."""
import argparse
import json
from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pet_staging_covariates import output_dir

COLORS=['#28678e','#b55a38','#749079','#8f6da5']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
                     'axes.spines.right':False,'axes.spines.top':False,
                     'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})


def forest(ax,df,labels,estimate='estimate',low='ci_low',high='ci_high',color=COLORS[0]):
    y=np.arange(len(df));values=df[estimate].to_numpy(float)
    ax.errorbar(values,y,xerr=np.vstack([values-df[low],df[high]-values]),fmt='o',color=color,capsize=2,ms=3,lw=1)
    ax.axvline(0,color='.55',lw=.7);ax.set_yticks(y,labels);ax.invert_yaxis();ax.margins(y=.10)


def run(root,out,save_hook=None):
    records=[];missing=[]
    def read(folder,name):
        path=root/folder/(name+'.csv')
        if not path.exists():missing.append(dict(component=folder+'/'+name,status='source_not_available'));return None
        return pd.read_csv(path)
    def save(fig,name,tables,claim,caption):
        if save_hook is not None:
            save_hook(fig,name,tables,claim,caption)
            records.append(dict(figure=name,claim=claim,caption=caption,sources=list(tables)))
            return
        fig.canvas.draw()
        for ext in ['pdf','svg','png']:
            fig.savefig(out/(name+'.'+ext),dpi=220,bbox_inches='tight',pad_inches=.16)
        for label,table in tables.items():table.to_csv(out/(name+'_'+label+'.csv'),index=False)
        records.append(dict(figure=name,claim=claim,caption=caption,sources=list(tables)))
        plt.close(fig)
    d=read('longitudinal','common_sample_trajectory_contrasts')
    if d is not None:
        d=d.loc[d.year.eq(2)&d.family.isin(['gaussian','fractional_logit'])&~d.model.str.contains('conditional')].copy()
        labels={x:y for x,y in zip(['demographic','amyloid','stage','MRI'],['Demographics','+ Centiloid','+ diagnosis','+ hippocampus/ICV'])}
        fig,axes=plt.subplots(1,2,figsize=(7.2,3.1),sharex=True,layout='constrained')
        for ax,family,letter in zip(axes,['gaussian','fractional_logit'],['a','b']):
            s=d.loc[d.family.eq(family)]
            forest(ax,s,[labels[x] for x in s.adjustment],color=COLORS[family=='fractional_logit'])
            ax.set_title(letter+'  '+family.replace('_',' ').capitalize(),loc='left',fontweight='bold')
            ax.set_xlabel('Two-year CDR-SB change contrast\nper 0.1 SUVR normalized IQR')
        save(fig,'Revision_01_trajectory',{'source':d},'Compare model-dependent two-year changes on one common sample.',
             'Same complete-case participants and tau IQR in both families. Time starts at the selected baseline CDR assessment, not tau PET. Points are standardized marginal change contrasts; bars are nominal 95% delta-method CIs from participant-clustered GEE covariance. Fractional-logit + MRI is the primary revision exploratory trajectory model. This is not an individual prognostic rate; analytic selection is not included in CIs.')
    d=read('longitudinal','monthly_visit_support')
    if d is not None:
        if 'clock' not in d:raise ValueError('Rebuild follow-up descriptions with explicit clocks')
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
        periods=read('longitudinal','observed_visit_periods')
        if periods is None:raise ValueError('Exact first-year participant counts required')
        first=periods.loc[periods.grouping.eq('all')&periods.clock.eq('baseline_CDR')&periods.period.eq('0_to_1')]
        if len(first)!=1:raise ValueError('Ambiguous first-year support')
        first=first.iloc[0]
        axes[0].text(.04,.88,f'First year\n{int(first.people)} participants, {int(first.visits)} visits',transform=axes[0].transAxes,fontsize=8,va='top')
        save(fig,'Revision_02_followup_support',{'source':d,'period_counts':periods},'Show the empirical visit distribution without smoothing sparse early observations.',
             'Monthly counts of actual post-baseline assessments, within the original PET+3-year analysis horizon. Clock is selected baseline CDR, matching the trajectory model; assessments can therefore extend beyond year 3 on this axis when baseline CDR precedes PET. A person can contribute in more than one month. Shaded region is the first year; bars are not predictions. Low/high tau are <=Q25/>=Q75 of the baseline-eligible sample. PET-clock descriptions, exact unique participant counts by period, near-year windows and durations are separate tables; do not mix the two time origins.')
    d=read('clinical','cross_sectional_results')
    if d is not None:
        s=d.loc[d.family.eq('regional_tau')].copy()
        fig,axes=plt.subplots(1,2,figsize=(7.2,3.4),layout='constrained')
        names={'tau_temporal_suvr_per01':'Temporal meta-ROI','tau_pvc_per01':'Temporal meta-ROI, PVC','tau_entorhinal_per01':'Entorhinal','tau_mtl_per01':'Medial temporal','tau_temporal_neo_per01':'Temporal neocortex','tau_neocortex_per01':'Neocortex excluding MTL'}
        for ax,y,letter in zip(axes,['hippocampus_icv','cdrsb'],['a','b']):
            z=s.loc[s.outcome.eq(y)]
            forest(ax,z,[names[x] for x in z.analysis],estimate='beta_outcome_sd',low='ci_low_outcome_sd',high='ci_high_outcome_sd')
            ax.set_title(letter+'  '+('Hippocampus/ICV' if y=='hippocampus_icv' else 'CDR-SB'),loc='left',fontweight='bold')
            ax.set_xlabel('Outcome SD per 0.1 SUVR')
        save(fig,'Revision_03_regional_tau',{'source':s},'Compare continuous tau definitions and partial-volume correction.',
             'A+ complete cases shared by every ROI within each outcome, adjusted for age, sex, education, APOE and Centiloid. HC3 95% CIs; BH q values across the named regional-tau family are in the source CSV. ROI scales and outcome samples differ; standardized outcome SD is sample-specific. These research ROIs do not implement validated AA biological staging.')
    d=read('clinical','tau_threshold_histogram')
    if d is not None:
        fig,ax=plt.subplots(figsize=(7.2,2.8),layout='constrained')
        for i,(group,g) in enumerate(d.groupby('at_stage')):
            ax.step((g.lower+g.upper)/2,g.n,where='mid',label=group,color=COLORS[i])
        ax.axvline(1.34,color='.35',linestyle='--',lw=1,label='Operational cutoff 1.34')
        ax.set_xlabel('Temporal meta-region flortaucipir SUVR');ax.set_ylabel('Participants');ax.legend(frameon=False,fontsize=7,ncol=3)
        save(fig,'Revision_04_tau_distribution',{'source':d},'Describe all observed tau values around the original classification cutoff.',
             'Fixed 0.05-SUVR bins covering all finite values; no smoothing or threshold optimization. Groups inherit the original PET-Core amyloid and temporal tau classification. This plot does not validate the 1.34 cutoff. Counts within 0.02, 0.05 and 0.10 of the cutoff are in population_description.csv.')
    d=read('transcriptomics','non_UTY_gene_sensitivities')
    if d is not None:
        if 'OXTR1' in set(d.gene_symbol):raise ValueError('OXTR1 label needs source verification; do not silently rename')
        genes=sorted(d.gene_symbol.unique())
        fig,axes=plt.subplots(1,3,figsize=(7.2,max(3.3,len(genes)*.24)),sharey=True,sharex=True,layout='constrained')
        for ax,model,letter in zip(axes,['original_covariates','male_only','female_only'],['a','b','c']):
            z=d.loc[d.model.eq(model)].set_index('gene_symbol').reindex(genes).reset_index()
            forest(ax,z,genes,estimate='log2FC');ax.set_xlabel('Log2 fold change')
            ax.set_title(letter+'  '+{'original_covariates':'All donors','male_only':'Male donors','female_only':'Female donors'}[model],loc='left',fontweight='bold')
            for tick in ax.get_yticklabels():tick.set_fontstyle('italic')
        save(fig,'Revision_05_gene_sex_sensitivity',{'source':d},'Display actual sex-restricted refits of the original non-UTY targets.',
             'Braak V–VI versus 0, adjusted for age, PMI, APOE and region; sex included only when variable. Donor-blocked limma; moderated 95% CIs. Each panel uses its own eligible donors. Targets were selected in the discovery cohort; this is a sensitivity analysis, not independent validation. Genome-wide and selected-family exploratory q values are in source data. UTY is reported separately.')
        regional=d.loc[d.model.eq('region_interaction')]
        fig,ax=plt.subplots(figsize=(7.2,max(3.5,len(genes)*.29)),layout='constrained')
        for j,region in enumerate(['EC','FC','TC']):
            z=regional.loc[regional.region.eq(region)].set_index('gene_symbol').reindex(genes)
            y=np.arange(len(genes))+(j-1)*.22
            ax.errorbar(z.log2FC,y,xerr=np.vstack([z.log2FC-z.ci_low,z.ci_high-z.log2FC]),fmt=['o','s','^'][j],ms=3,capsize=1,color=COLORS[j],label=region)
        ax.set_yticks(range(len(genes)),genes,fontstyle='italic');ax.invert_yaxis();ax.axvline(0,color='.5',lw=.8)
        ax.set_xlabel('Region-specific late-Braak log2 fold change');ax.legend(frameon=False,ncol=3,loc='best')
        save(fig,'Revision_06_gene_regions',{'source':regional},'Distinguish directional agreement from equal effect magnitude across regions.',
             'Braak-by-region cell-means model adjusted for age, sex, PMI and APOE, with donor correlation. Moderated 95% CIs. A joint two-contrast interaction F test with genome-wide BH correction is reported in region_heterogeneity.csv. Overlapping intervals or nonsignificant interaction do not establish equal regional effects.')
    d=read('pathways','current_GO_significant_grouped')
    if d is not None and len(d):
        d=d.loc[d.display_representative.eq(True)].sort_values(['q','p','go_id']).head(12).copy()
        labels=['\n'.join(textwrap.wrap(n,48)) for n in d['name']]
        fig,ax=plt.subplots(figsize=(7.2,max(3.5,len(d)*.45)),layout='constrained')
        x=-np.log10(d.q.clip(lower=np.finfo(float).tiny))
        ax.barh(range(len(d)),x,color=[COLORS[0] if x=='Down' else COLORS[1] for x in d.direction])
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(color=COLORS[0],label='Lower in late Braak'),Patch(color=COLORS[1],label='Higher in late Braak')],frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,-.10),fontsize=7)
        ax.set_yticks(range(len(d)),labels);ax.invert_yaxis();ax.set_xlabel('-log10(BH q), all current eligible GO BP terms')
        save(fig,'Revision_07_current_GO',{'source':d},'Show representative current GO associations after removing redundant display groups.',
             'Up to 12 lowest-q representatives of significant full-model GO sets, after greedy Jaccard >=0.5 display grouping. All current 10–300-gene sets are tested before BH correction; obsolete terms are excluded before testing. Blue: Down; orange: Up. Full results and overlap memberships remain available. This is annotation-based enrichment, not cell activation or pathway activity.')
    d=read('external','external_target_results')
    if d is not None:
        d=d.sort_values('gene_symbol').copy()
        fig,ax=plt.subplots(figsize=(7.2,max(3.5,len(d)*.27)),layout='constrained')
        for i,row in enumerate(d.itertuples()):
            if pd.isna(row.log2FC):
                ax.text(.98,i,'Not measured',ha='right',va='center',transform=ax.get_yaxis_transform(),fontsize=7,color='.4')
                continue
            color=COLORS[0] if row.direction_agrees else COLORS[1]
            ax.errorbar(row.log2FC,i,xerr=[[row.log2FC-row.ci_low],[row.ci_high-row.log2FC]],fmt='o',ms=3,capsize=2,color=color)
        ax.set_yticks(range(len(d)),d.gene_symbol,fontstyle='italic');ax.invert_yaxis();ax.axvline(0,color='.5',lw=.8)
        ax.set_xlabel('External AD versus control log2 fold change')
        from matplotlib.lines import Line2D
        ax.legend(handles=[Line2D([0],[0],marker='o',color=COLORS[0],lw=0,label='Same direction'),Line2D([0],[0],marker='o',color=COLORS[1],lw=0,label='Opposite direction')],frameon=False,ncol=2,loc='lower right')
        save(fig,'Revision_08_external_genes',{'source':d},'Report every preselected target, including opposite effects and unmeasured genes.',
             'Independent GSE44770 HBTRC prefrontal AD/control contrast, adjusted for age, sex, PMI, tissue pH, RIN and batch. Moderated 95% CIs; BH correction across all 16 original non-UTY targets, including an unmeasured target in the family size. Colors denote agreement with discovery direction, not significance. This contrast differs from Braak V–VI versus 0 and is not an exact replication or a link to ADNI mechanisms.')
    d=read('cell_composition','non_UTY_composition_sensitivity')
    if d is not None:
        genes=sorted(d.gene_symbol.unique());fig,axes=plt.subplots(1,2,figsize=(7.2,5.0),sharex=True,sharey=True,layout='constrained')
        reference=d.loc[d.model.eq('unadjusted_composition')].set_index('gene_symbol').reindex(genes)
        for ax,model,title in zip(axes,['CA_joint_quantile','NG_joint_quantile'],['a  Temporal-cortex reference','b  Prefrontal-cortex reference']):
            adjusted=d.loc[d.model.eq(model)].set_index('gene_symbol').reindex(genes)
            for source,color,offset,label in [(reference,'.55',-.12,'Base model'),(adjusted,COLORS[0],.12,'Composition adjusted')]:
                ax.errorbar(source.log2FC,np.arange(len(genes))+offset,xerr=np.vstack([source.log2FC-source.ci_low,source.ci_high-source.log2FC]),fmt='o',ms=2.5,capsize=2,color=color,label=label)
            ax.set_title(title,loc='left',fontweight='bold');ax.axvline(0,color='.5',lw=.7);ax.set_yticks(range(len(genes)),genes,fontstyle='italic')
            ax.set_xlabel('Late Braak versus 0 log2 fold change');ax.legend(frameon=False,fontsize=6.5,loc='lower left')
        axes[0].invert_yaxis()
        save(fig,'Revision_09_cell_composition',{'source':d},'Compare original target-gene effects before and after reference-derived composition adjustment.',
             'Six-class dtangle reference-based estimates from published CA and NG single-nucleus profiles; joint-quantile variant shown. Original 17 target genes excluded from marker selection. All 16 non-UTY genes shown, 71 donors and 213 tissue samples. Donor-blocked moderated 95% CIs are conditional on estimated composition and original gene selection; native-scale sensitivities remain in source data. These fractions are cross-platform estimates, not validated cell counts or functional states.')
    d=read('visit_process','visit_process_results')
    if d is not None and len(d):
        fig,ax=plt.subplots(figsize=(7.2,3.0),layout='constrained')
        for i,(model,label,color) in enumerate([('unweighted_windows','Unweighted',COLORS[0]),('history_weighted_windows','History weighted',COLORS[1])]):
            s=d.loc[d.model.eq(model)].sort_values('year')
            ax.errorbar(s.year+(i-.5)*.08,s.estimate,yerr=np.vstack([s.estimate-s.ci_low,s.ci_high-s.estimate]),fmt='o',ms=4,capsize=3,color=color,label=label)
        support=d.loc[d.model.eq('unweighted_windows')].set_index('year')
        ax.set_xticks([1,2,3],[f'Year {y}\n{int(support.loc[y,"n_observed"])} observed' for y in [1,2,3]])
        ax.axhline(0,color='.5',lw=.7);ax.set_ylabel('Conditional CDR-SB level contrast\nper 0.1 SUVR normalized IQR');ax.legend(frameon=False)
        ax.set_title(f'Target population: {int(support.n_target.iloc[0])} participants',loc='left',fontsize=9)
        save(fig,'Revision_10_visit_process',{'source':d},'Show sensitivity to measured time-varying observation history, without claiming missing-not-at-random correction.',
             'Fixed non-overlapping annual windows relative to tau PET, +/-6 months, requiring administrative opportunity for all three windows. Same target population before/after weighting; one nearest observation per window. Fractional-logit model conditional on baseline CDR, diagnosis, demographics and amyloid. Bars are percentile participant-bootstrap CIs refitting observation and outcome models; requested/valid replicate counts and actual observed counts are in the source table. Year-3 window extends to 3.5 years and differs from the original 3-year analysis.')
    d=read('clinical','biological_stage_counts')
    if d is not None:
        counts=d.groupby('biological_stage')['n'].sum();order=[k for k in ['amyloid_negative','A','B','C','D','discordant_neocortical_without_MTL','Unavailable'] if k in counts]
        labels={'amyloid_negative':'Amyloid negative','discordant_neocortical_without_MTL':'Discordant pattern','Unavailable':'Unclassified'}
        fig,ax=plt.subplots(figsize=(7.2,2.6),layout='constrained');ax.bar(range(len(order)),counts.reindex(order),color=COLORS[0]);ax.set_xticks(range(len(order)),[labels.get(k,'Stage '+k) for k in order],rotation=15,ha='right');ax.set_ylabel('Participants')
        for i,k in enumerate(order):ax.text(i,counts[k]+3,str(counts[k]),ha='center',fontsize=7)
        ax.margins(y=.15)
        save(fig,'Revision_11_exploratory_staging',{'source':d},'Show all operational-stage assignments, including discordant and unclassified participants.',
             'Exploratory AA-inspired PET staging using outcome-independent two-component mixture cutoffs estimated in the current Berkeley processing pipeline. This is not an externally validated implementation of AA staging. MTL: bilateral entorhinal; neocortex: volume-weighted CTX parcels excluding entorhinal and parahippocampal cortex. Threshold bootstrap diagnostics and all stage counts are reported separately; no outcome-based threshold selection.')
    if save_hook is None:
        (out/'figure_captions.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
        pd.DataFrame(missing or [dict(component='requested_figure_sources',status='ok')]).to_csv(out/'figure_status.csv',index=False)
        print('Diagnostic figures exported:',len(records),'Missing source tables:',len(missing),flush=True)


