"""Shared clinical and trajectory panels for the single publication figure set."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib import colors
from matplotlib.ticker import MaxNLocator

COLORS=['#28678e','#b55a38','#749079','#8f6da5']
LABELS={'hippocampus_icv':'Hippocampus / ICV','cdrsb':'CDR-SB','fdg_metaroi':'FDG MetaROI'}

class Figures:
    def __init__(self,root,out):
        self.root=root;self.out=out;self.out.mkdir(parents=True,exist_ok=True)
        if any(self.out.iterdir()):raise ValueError('Use an empty final figure directory')
        self.source_out=root/'reporting/figure_sources';self.source_out.mkdir(parents=True,exist_ok=True)
        self.records=[]
    def read(self,path):
        p=self.root/path
        if not p.is_file():raise FileNotFoundError(f'Current-run figure input missing: {path}')
        return pd.read_csv(p)
    def save(self,fig,name,tables,claim,interval='See current-run source model specification'):
        raise NotImplementedError('Use PublicationFigures through the unified export entry point')
    def forest(self,ax,d,label,beta='beta',low='ci95_lower',high='ci95_upper',scale=1):
        x=d[beta].to_numpy(float)*scale;lo=d[low].to_numpy(float)*scale;hi=d[high].to_numpy(float)*scale
        if not np.isfinite(np.c_[x,lo,hi]).all():raise ValueError('Nonfinite forest estimates')
        y=np.arange(len(d))
        ax.hlines(y,lo,hi,color=COLORS[0],lw=1.1);ax.plot(x,y,'o',color=COLORS[0],ms=4)
        ax.axvline(0,color='.5',lw=.7);ax.set_yticks(y,label);ax.invert_yaxis();ax.margins(y=.2)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
    def original(self):
        a=self.read('results/adni/adni_results.csv')
        master=self.read('data/private/adni/adni_pet_aligned_master.csv')
        counts=master.at_stage.value_counts().reindex(['A-T-','A-T+','A+T-','A+T+'],fill_value=0).rename_axis('at_stage').reset_index(name='n')
        fig,ax=plt.subplots(figsize=(5.4,3.5),layout='constrained')
        bars=ax.bar(counts.at_stage,counts.n,color=[COLORS[0],'#b9b9b9',COLORS[2],COLORS[1]])
        ax.bar_label(bars,padding=3);ax.set_ylabel('Participants');ax.margins(y=.15)
        self.save(fig,'FigureS1',{'source':counts},'Show the complete PET-aligned four-state cohort.')
        for name,analysis in [('Figure2','confirmatory_180day'),('FigureS2','exploratory_FDG'),('FigureS4','APOE_stratified_display')]:
            d=a.loc[a.analysis.eq(analysis)].copy();outcomes=list(d.outcome.unique())
            fig,axs=plt.subplots(1,len(outcomes),figsize=(7.2,3.5),layout='constrained',squeeze=False)
            for i,(ax,y) in enumerate(zip(axs.flat,outcomes)):
                z=d.loc[d.outcome.eq(y)]
                labels=z.contrast.astype(str)
                if name=='FigureS4':labels=labels+'; APOE '+z.apoe4_status.astype(str)
                scale=1e4 if y=='hippocampus_icv' else 1
                self.forest(ax,z,labels,scale=scale)
                ax.set_title((chr(97+i)+'  ' if len(outcomes)>1 else '')+LABELS.get(y,y),loc='left')
                ax.set_xlabel('Adjusted difference'+(r' ($\times 10^{-4}$)' if scale!=1 else ''))
            self.save(fig,name,{'source':d},analysis,'Points and HC3 95% confidence intervals; BH families remain in source CSV.')
        d=a.loc[a.analysis.eq('APOE_joint_interaction')]
        fig,ax=plt.subplots(figsize=(5.4,3.5),layout='constrained')
        ax.bar([LABELS[y] for y in d.outcome],-np.log10(d.p_value),color=COLORS[0])
        ax.axhline(-np.log10(.05),color='.4',ls='--');ax.set_ylabel('−log10 joint-test P')
        self.save(fig,'FigureS3',{'source':d},'Joint APOE-by-A/T interaction tests.')
        self.flow(master,counts)
    def flow(self,master,counts):
        m=self.read('data/public/gse131617/sample_manifest.csv')
        original=self.read('results/extensions/longitudinal/longitudinal_results.csv')
        o=original.loc[original.analysis.eq('primary_Aplus_180d')].iloc[0]
        common=self.read('revision/longitudinal/common_sample_trajectory_contrasts.csv')
        n=common.n_standardized_people.unique()
        if len(n)!=1:raise ValueError('Inconsistent common trajectory population')
        target=self.read('data/contracts/gene_symbols.csv')
        go=self.read('targeted/qc_joint/GO_all_models.csv');full=go.loc[go.model.eq('baseline')]
        retained=self.read('cel_qc/retained_sample_manifest.csv')
        support=self.read('enhancement/clinical/common_observed_visit_periods.csv')
        post=support.loc[support.grouping.eq('all')&support.group.eq('all')&support.clock.eq('baseline_CDR'),'visits'].sum()
        observations=int(n[0]+post)
        rows=[dict(component='PET-aligned participants',n=len(master)),
              dict(component='Three-state participants',n=int(master.at_stage.ne('A-T+').sum())),
              dict(component='Larger-sample sensitivity participants',n=int(o.n_participants)),
              dict(component='Larger-sample sensitivity observations',n=int(o.n_observations)),
              dict(component='Common longitudinal participants',n=int(n[0])),
              dict(component='Common longitudinal observations',n=observations),
              dict(component='Tissue donors',n=m.donor_key.nunique()),
              dict(component='Provided tissue arrays',n=len(m)),dict(component='QC-retained tissue arrays',n=len(retained)),
              dict(component='Technical exclusions',n=len(m)-len(retained)),
              dict(component='Non-Y selected targets',n=int(target.gene_symbol.ne('UTY').sum())),
              dict(component='GO tested',n=len(full)),dict(component='GO q<0.05',n=int(full.q.lt(.05).sum()))]
        fig,ax=plt.subplots(figsize=(7.2,6.6));ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)
        boxes=[
            (.02,.74,.52,.19,f'ADNI: {len(master)} PET-aligned participants\nThree-state cohort: {rows[1]["n"]}'),
            (.02,.45,.52,.19,f'Longitudinal cohort\n{int(o.n_participants)} people / {int(o.n_observations)} observations\nLarger-sample sensitivity models'),
            (.02,.14,.52,.20,f'Common-model subset\n{int(n[0])} people / {observations} observations\nPrimary exploratory trajectory model'),
            (.60,.68,.38,.25,f'GSE131617: {m.donor_key.nunique()} donors\n{len(m)} provided arrays\n{len(m)-len(retained)} technical exclusion; {len(retained)} retained\nBraak V–VI versus 0'),
            (.60,.14,.38,.42,f'{int(target.gene_symbol.ne("UTY").sum())} non-Y selected targets\nUTY reported separately\n{len(full)} GO sets tested\n{int(full.q.lt(.05).sum())} with q < 0.05\n(QC-passed baseline)\n\nTechnical and external-cohort checks')]
        for x,y,w,h,label in boxes:
            ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.008',fc='#edf3f5',ec='#647783',lw=.8))
            ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=8.7,linespacing=1.5)
        for x,y,yy in [(.28,.74,.64),(.28,.45,.34),(.79,.68,.56)]:
            ax.annotate('',xy=(x,yy),xytext=(x,y),arrowprops=dict(arrowstyle='->',color='#647783'))
        ax.text(.5,.04,'Independent datasets; no participant-level PET–tissue linkage.',ha='center',fontsize=8)
        self.save(fig,'Figure1',{'counts':pd.DataFrame(rows)},'Trace distinct clinical and tissue analyses without implying linked individuals.')
    def trajectories(self):
        d=self.read('revision/longitudinal/common_sample_trajectory_contrasts.csv')
        d=d.loc[d.year.eq(2)&d.adjustment.isin(['demographic','amyloid','stage','MRI'])].copy()
        s=self.read('enhancement/clinical/common_observed_visit_periods.csv')
        s=s.loc[s.grouping.eq('all')&s.group.eq('all')]
        fig=plt.figure(figsize=(7.2,6.5),layout='constrained');g=fig.add_gridspec(2,2,height_ratios=[1.2,1])
        ax=fig.add_subplot(g[0,:]);adjustments=['demographic','amyloid','stage','MRI']
        for j,(family,shift) in enumerate([('fractional_logit',-.12),('gaussian',.12)]):
            z=d.loc[d.family.eq(family)].set_index('adjustment').reindex(adjustments)
            if z.estimate.isna().any():raise ValueError('Missing common-model contrast')
            y=np.arange(4)+shift
            ax.hlines(y,z.ci_low,z.ci_high,color=COLORS[j]);ax.plot(z.estimate,y,'o',color=COLORS[j],label=family.replace('_',' '))
        ax.set_yticks(range(4),['Demographics','+ Amyloid','+ Clinical stage','+ MRI']);ax.invert_yaxis()
        ax.axvline(0,color='.5',lw=.7);ax.set_xlabel('Two-year CDR-SB contrast per 0.1 SUVR')
        ax.legend(frameon=False,ncol=2,loc='lower left',bbox_to_anchor=(0,1))
        ax.set_title(f'Common sample: {int(d.n_standardized_people.iloc[0])} participants',loc='right',fontsize=8)
        ax.text(-.12,1.07,'a',transform=ax.transAxes,fontweight='bold',fontsize=10)
        for col,clock in enumerate(['baseline_CDR','tau_PET']):
            ax=fig.add_subplot(g[1,col])
            z=s.loc[s.clock.eq(clock)].set_index('period').reindex(['0_to_1','1_to_2','over_2'])
            bars=ax.bar(range(3),z.visits,color=COLORS[0])
            ax.bar_label(bars,labels=[f'{int(r.visits)} visits\n{int(r.people)} people' for r in z.itertuples()],padding=4,fontsize=7)
            ax.margins(y=.3);ax.set_xticks(range(3),['0–1','1–2','>2']);ax.set_xlabel('Years from '+clock.replace('_',' '))
            ax.set_ylabel('Observed follow-up assessments')
            ax.text(-.12,1.04,chr(98+col),transform=ax.transAxes,fontweight='bold',fontsize=10)
        self.save(fig,'Figure4',{'model':d,'support':s},'Common-sample model dependence alongside actual follow-up support.',
                  'Nominal delta-method 95% CIs using participant-clustered covariance; descriptive support counts.')
        old=self.read('results/extensions/stage_sensitivity/standardized_stage_change.csv')
        periods=self.read('revision/longitudinal/observed_visit_periods.csv')
        p=periods.loc[periods.grouping.eq('all')&periods.group.eq('all')&periods.clock.eq('baseline_CDR')]
        fig,axs=plt.subplots(1,2,figsize=(7.2,4),layout='constrained')
        label_map={'gaussian_original':'Gaussian','gaussian_stage':'Gaussian + stage',
                   'fractional_logit_original':'Fractional logit','fractional_logit_stage':'Fractional logit + stage'}
        labels=old.model.map(label_map)
        if labels.isna().any():raise ValueError('Unlabeled longitudinal sensitivity model')
        self.forest(axs[0],old,labels+' / year '+old.year.astype(str),beta='estimate',low='ci_low',high='ci_high')
        axs[0].set_title('a  Larger-sample model contrasts',loc='left');axs[1].set_title('b  Observed follow-up',loc='left')
        axs[0].set_xlabel('CDR-SB contrast per 0.1 SUVR')
        base=int(old.n_standardized_people.iloc[0])
        shown=pd.concat([pd.DataFrame([dict(period='baseline',people=base,visits=base)]),p[['period','people','visits']]],ignore_index=True)
        bars=axs[1].bar(range(len(shown)),shown.visits,color=COLORS[0])
        axs[1].bar_label(bars,labels=[str(int(n))+' people' for n in shown.people],padding=3,fontsize=7)
        axs[1].set_xticks(range(len(shown)),['Baseline','0–1','1–2','>2']);axs[1].margins(y=.2)
        axs[1].set_ylabel('Actual assessments');axs[1].set_xlabel('Years since baseline CDR')
        self.save(fig,'FigureS9',{'models':old,'support':shown},'Larger-sample longitudinal sensitivities and observed support.')
        monthly=self.read('enhancement/clinical/common_monthly_visit_support.csv');monthly=monthly.loc[monthly.grouping.eq('all')]
        fig,axs=plt.subplots(2,1,figsize=(7.2,5),layout='constrained')
        for i,(ax,clock) in enumerate(zip(axs,['baseline_CDR','tau_PET'])):
            z=monthly.loc[monthly.clock.eq(clock)].sort_values('month')
            ax.bar((z.month-.5)/12,z.visits,width=.075,color=COLORS[0]);ax.axvspan(0,1,color='#edf1f1',zorder=-1)
            first=s.loc[s.clock.eq(clock)&s.period.eq('0_to_1')].iloc[0]
            clock_label={'baseline_CDR':'Baseline CDR','tau_PET':'Tau PET'}[clock]
            ax.set_title(f'{chr(97+i)}  {clock_label}; first year: {int(first.visits)} visits / {int(first.people)} people',loc='left',fontsize=9)
            ax.set_ylabel('Observed assessments')
        axs[1].set_xlabel('Years from the indicated time origin')
        self.save(fig,'FigureS20',{'monthly':monthly,'periods':s},'Monthly actual follow-up using two explicit clocks.')


