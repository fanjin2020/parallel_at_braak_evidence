"""One QC-aligned figure set and compact source archive; input runs remain read-only."""
import csv as csvlib
import io
import json
from pathlib import Path
import shutil
import textwrap
import zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from figure_primitives import Figures, COLORS
from figure_support import followup_support, regional_tau, tau_distribution
from figure_support import visit_process, regional_staging, stage_uncertainty
from export_io import check_public,read_csv_exact,csv_bytes
from analysis_utils import csv,dump,sha


MODELS=['baseline','scan_date','joint_CA_joint_quantile','joint_NG_joint_quantile']
NAMES=['Baseline','Scan date','Scan date + CA','Scan date + NG']
CODE=Path(__file__).resolve().parents[1]

def composition_label(model):
    if model=='baseline':return 'Baseline'
    if model=='scan_date':return 'Scan date'
    prefix='Scan date + ' if model.startswith('joint_') else ''
    ref='CA' if '_CA_' in model else 'NG'
    scale='joint quantile' if model.endswith('_joint_quantile') else 'native log scale'
    return prefix+ref+' composition ('+scale+')'

def prepare(root,review,out):
    """Construct display/source tables, without changing statistical estimates."""
    root=Path(root);review=Path(review);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    go=pd.read_csv(root/'targeted/qc_joint/GO_all_models.csv')
    terms=pd.read_csv(root/'revision/go_annotation/go_terms.csv')
    members=pd.read_csv(root/'revision/go_annotation/go_gene_membership.csv',dtype=str)
    reps=pd.read_csv(root/'cel_qc/fixed_gene_representatives.csv',dtype=str)
    sets={key:set(g.entrez_id)&set(reps.entrez_id) for key,g in members.groupby('go_id',sort=False)}
    d=go.loc[go.model.eq('baseline')].merge(terms[['go_id','name']],on='go_id',validate='one_to_one')
    if len(d)!=4826 or d.go_id.duplicated().any():raise ValueError('Incomplete QC GO family')
    d=d.loc[d.q.lt(.05)].sort_values(['q','p','go_id'],kind='stable').reset_index(drop=True)
    representatives=[];groups=[];flags=[]
    for row in d.itertuples():
        group=next((i for i,rep in enumerate(representatives,1) if len(sets[row.go_id]&sets[rep])/len(sets[row.go_id]|sets[rep])>=.5),None)
        is_rep=group is None
        if is_rep:representatives.append(row.go_id);group=len(representatives)
        groups.append(group);flags.append(is_rep)
    d['overlap_group']=groups;d['display_representative']=flags
    csv(out/'QC_GO_grouped.csv',d);csv(out/'QC_GO_display_12.csv',d.loc[d.display_representative].head(12))
    targets=pd.read_csv(root/'targeted/qc_joint/target_effects_full_background.csv')
    target=targets.loc[targets.model.isin(MODELS)&targets.gene_symbol.ne('UTY')].copy()
    if len(target)!=64:raise ValueError('Expected 16 non-Y targets across four QC models')
    target['population']='QC-retained 212 arrays / 71 donors';target['target_scope']='16 non-Y selected targets'
    csv(out/'QC_nonY_targets.csv',target)
    effects=pd.read_csv(review/'tissue/target_effects.csv')
    male=effects.loc[effects.kind.eq('male')&effects.gene_symbol.eq('UTY')].copy()
    if len(male)!=4:raise ValueError('Missing male-only UTY models')
    diagnostics=pd.read_csv(review/'tissue/model_diagnostics.csv')
    male=male.merge(diagnostics.loc[diagnostics.kind.eq('male'),['model','n_samples','n_donors']],on='model',validate='one_to_one')
    male['population']='QC-retained male-only arrays';male['target_scope']='Y-linked target, reported separately'
    csv(out/'UTY_male_only.csv',male)
    for name,kinds in [('QC_sex_subgroups',['male','female']),('QC_age_sensitivity',['reference','age_overlap','age_spline']),('QC_region_targets',['feature_region'])]:
        csv(out/(name+'.csv'),effects.loc[effects.kind.isin(kinds)&effects.gene_symbol.ne('UTY')])
    rows=[]
    for model in MODELS:
        z=targets.loc[targets.model.eq(model)]
        for scope,part in [('all_17_including_UTY',z),('non_Y_16',z.loc[z.gene_symbol.ne('UTY')]),('UTY_pooled_diagnostic_only',z.loc[z.gene_symbol.eq('UTY')])]:
            rows.append(dict(model=model,scope=scope,tested=len(part),q05=int(part.q.lt(.05).sum()),samples=212,donors=71))
    csv(out/'QC_target_counts.csv',rows)
    baseline=target.loc[target.model.eq('baseline'),['gene_symbol','log2FC','q']].rename(columns={'log2FC':'QC_baseline_log2FC','q':'QC_baseline_q'})
    for cohort,path in [('GSE44770','revision/external/external_target_results.csv'),('GSE106241','enhancement/external/external_target_results.csv')]:
        external=pd.read_csv(root/path)
        external=external.drop(columns=['discovery_direction','direction_agrees'],errors='ignore')
        external=external.merge(baseline,on='gene_symbol',how='left',validate='many_to_one')
        external['QC_direction_agrees']=np.where(external.log2FC.notna(),np.sign(external.log2FC)==np.sign(external.QC_baseline_log2FC),pd.NA)
        external['comparison_scope']='Related AD/control context, not a matched Braak replication' if cohort=='GSE44770' else 'Independent temporal-cortex Braak comparison; region and models differ'
        csv(out/(cohort+'_QC_gene_comparison.csv'),external)
    raw=root/'data/private/adni/raw'
    date=pd.read_csv(raw/'DATA_DOWNLOADED_DATE.csv',dtype=str)
    provenance=dict(snapshot_dates=date.to_dict('records'),PET_tables=[dict(table=name,sha256=sha(raw/name)) for name in ['UCBERKELEY_AMY_6MM.csv','UCBERKELEY_TAU_6MM.csv']],
        scope='Authorized local CSV snapshot; package and PET processing-release metadata are not inferred from filenames.')
    dump(out/'ADNI_snapshot_provenance.json',provenance)

class PublicationFigures(Figures):
    def __init__(self,root,review,new,out):
        self.root=Path(root);self.review=Path(review);self.new=Path(new);self.out=Path(out)
        self.out.mkdir(parents=True,exist_ok=True)
        if any(self.out.iterdir()):raise ValueError('Final figures must be written to an empty directory')
        self.source_out=self.new/'figure_sources';self.source_out.mkdir(parents=True,exist_ok=True)
        self.records=[]
        self.numbering={f'FigureS{old}':f'FigureS{i}' for i,old in enumerate([12,8,22,23,1,11,17,19,2,3,4,5,21,6,13,7,18,14,24,20,16,9,10,15],1)}
    def table(self,name):return pd.read_csv(self.new/'publication_tables'/name)
    def reviewed(self,name):return pd.read_csv(self.review/name)
    def save(self,fig,name,tables,claim,interval='See source model specification'):
        if any(r['figure']==name for r in self.records):raise ValueError('Duplicate final figure '+name)
        if name=='Figure1':
            for item in fig.findobj(matplotlib.text.Text):
                item.set_text(item.get_text().replace('Technical and external-cohort checks','Technical and external-\ncohort checks'))
        final_name=self.numbering.get(name,name)
        fig.canvas.draw()
        for ext in ['pdf','svg','png']:
            fig.savefig(self.out/(final_name+'.'+ext),dpi=360,bbox_inches='tight',pad_inches=.16,facecolor='white')
        plt.close(fig)
        files=[]
        for label,data in tables.items():
            file=name+'_'+label+'.csv'
            data.to_csv(self.source_out/file,index=False);files.append(file)
        self.records.append(dict(figure=name,final_figure=final_name,claim=claim,source_data=files,interval=interval))
        print(final_name+' generated',flush=True)
    def genes(self):
        targets=self.read('targeted/qc_joint/target_effects_full_background.csv')
        non=targets.loc[targets.model.eq('baseline')&targets.gene_symbol.ne('UTY')].sort_values('log2FC')
        y=np.arange(len(non))
        fig,axs=plt.subplots(1,2,figsize=(8.4,max(4.8,.32*len(non))),layout='constrained',gridspec_kw={'width_ratios':[1,1.15]})
        axs[0].hlines(y,non.ci_low,non.ci_high,color=COLORS[0],lw=1)
        axs[0].scatter(non.log2FC,y,facecolors=[COLORS[0] if q<.05 else 'white' for q in non.q],edgecolors=COLORS[0],s=26,zorder=3)
        axs[0].axvline(0,color='.5',lw=.6)
        axs[0].set_yticks(y,non.gene_symbol,fontstyle='italic');axs[0].invert_yaxis()
        axs[0].set_xlabel('Log2 fold change (95% CI)');axs[0].set_title('a  QC-passed target effects',loc='left')
        axs[0].text(0,-.12,'Filled: full-feature q < 0.05',transform=axs[0].transAxes,fontsize=7)
        go=self.read('targeted/qc_joint/GO_all_models.csv')
        models=['baseline','scan_date','joint_CA_joint_quantile','joint_NG_joint_quantile']
        labels=['Baseline','Scan date','Scan date + CA','Scan date + NG']
        counts=[]
        for model in models:
            z=go.loc[go.model.eq(model)]
            if z.empty:raise ValueError('Missing QC GO model: '+model)
            for direction in ['Down','Up']:
                counts.append(dict(model=model,direction=direction,n=int((z.q.lt(.05)&z.direction.eq(direction)).sum()),tested=len(z)))
        counts=pd.DataFrame(counts);left=np.zeros(4)
        for i,direction in enumerate(['Down','Up']):
            values=counts.loc[counts.direction.eq(direction)].set_index('model').loc[models,'n'].to_numpy()
            axs[1].barh(range(4),values,left=left,color=COLORS[i],label=direction+'-ranked')
            left+=values
        for i,n in enumerate(left):axs[1].text(n+max(left)*.02,i,str(int(n)),va='center',fontsize=8)
        axs[1].set_yticks(range(4),labels);axs[1].invert_yaxis();axs[1].set_xlim(0,max(left)*1.15)
        axs[1].set_xlabel('GO sets with within-model q < 0.05')
        axs[1].set_title('b  QC-passed GO associations',loc='left');axs[1].legend(frameon=False,loc='lower right')
        axs[1].xaxis.set_major_locator(MaxNLocator(integer=True,nbins=4))
        self.save(fig,'Figure3',{'targets':non,'GO_counts':counts},'QC-passed target effects and full-family GO support.',
                  'Nominal moderated 95% CIs; target q uses the full feature family. GO counts use all eligible sets per model.')
        f=self.read('targeted/qc_joint/features_baseline.csv');t=self.read('data/contracts/gene_symbols.csv')
        f['id']=f.id.astype(str);t['feature_id']=t.feature_id.astype(str)
        f=f.merge(t[['feature_id','gene_symbol']],left_on='id',right_on='feature_id',how='left',validate='one_to_one')
        fig,ax=plt.subplots(figsize=(7.2,4.6),layout='constrained')
        y=-np.log10(f.q.clip(lower=1e-300));ax.scatter(f.log2FC,y,s=5,color='#bdc6ca',alpha=.65,rasterized=True)
        for symbol,marker,color in [('non_Y','o',COLORS[0]),('UTY','D',COLORS[1])]:
            keep=f.gene_symbol.notna()& (f.gene_symbol.eq('UTY') if symbol=='UTY' else f.gene_symbol.ne('UTY'))
            ax.scatter(f.loc[keep,'log2FC'],y[keep],s=24,marker=marker,color=color,label='UTY (descriptive)' if symbol=='UTY' else '16 non-Y selected targets')
        ax.axhline(-np.log10(.05),ls=':',color='.4',lw=.8);ax.axvline(0,color='.5',lw=.7)
        ax.set(xlabel='QC baseline: Braak V–VI versus 0 log2 fold change',ylabel='−log10 full-feature q',title='212 retained arrays; 22,011 tested features')
        ax.legend(frameon=False,fontsize=7)
        self.save(fig,'FigureS5',{'source':f},'QC-retained full-feature landscape with fixed target highlighting.','Full-feature BH q; targets were selected before these sensitivity fits.')
        d=self.reviewed('tissue/influence_summary.csv');d=d.loc[d.gene_symbol.ne('UTY')]
        fig,axes=plt.subplots(2,2,figsize=(8.4,8.2),layout='constrained',sharex=True)
        genes=sorted(d.gene_symbol.unique());offsets={'omit_donor':-.13,'omit_scan_date':.13}
        for ax,model,label in zip(axes.flat,MODELS,NAMES):
            for j,(omission,shift) in enumerate(offsets.items()):
                z=d.loc[d.model.eq(model)&d.omission.eq(omission)].set_index('gene_symbol').loc[genes]
                y=np.arange(16)+shift;ax.hlines(y,z.min_effect,z.max_effect,color=COLORS[j],lw=1.5)
                ax.plot(z.reference_effect,y,'o',color=COLORS[j],ms=2.6,label='Omit one donor' if j==0 else 'Omit one scan date')
            ax.set_yticks(range(16),genes,fontstyle='italic');ax.invert_yaxis();ax.axvline(0,lw=.6,color='.5');ax.set_title(label,loc='left');ax.set_xlabel('Log2 fold change');ax.xaxis.set_major_locator(MaxNLocator(4))
        axes[0,0].legend(frameon=False,fontsize=7,loc='lower left')
        self.save(fig,'FigureS6',{'source':d},'QC target influence ranges across all donor and scan-date omissions.',
                  'Lines are min–max estimates over leave-one-out fits, NOT confidence intervals; points are full-sample references.')
    def imported(self):
        for panel in [followup_support,regional_tau,tau_distribution,visit_process,regional_staging]:
            panel(self.root/'revision',self.save)
        stage_uncertainty(self.root/'enhancement',self.save)
    def qc_supplements(self):
        z=self.table('QC_GO_display_12.csv').sort_values(['q','p','go_id'],kind='stable')
        fig,ax=plt.subplots(figsize=(7.2,5.7),layout='constrained')
        labels=[textwrap.fill(f'{r.go_id}  {r.name}',48) for r in z.itertuples()]
        ax.barh(range(len(z)),-np.log10(z.q),color=[COLORS[0] if x=='Down' else COLORS[1] for x in z.direction])
        ax.set_yticks(range(len(z)),labels,fontsize=7);ax.invert_yaxis();ax.set_xlabel('−log10 full-family q (QC baseline)')
        self.save(fig,'FigureS7',{'source':z},'Nonredundant QC-baseline GO associations; deterministic q, P, GO-ID ordering.',
                  'Twelve representatives after Jaccard grouping; not twelve newly tested hypotheses.')
        d=self.table('QC_age_sensitivity.csv');genes=sorted(d.gene_symbol.unique())
        fig,axes=plt.subplots(2,2,figsize=(8.4,8.4),layout='constrained',sharex=True)
        for ax,model,label in zip(axes.flat,MODELS,NAMES):
            for j,(kind,title) in enumerate([('reference','Linear age'),('age_overlap','Age overlap'),('age_spline','Age spline')]):
                z=d.loc[d.model.eq(model)&d.kind.eq(kind)].set_index('gene_symbol').loc[genes]
                y=np.arange(16)+(j-1)*.2;ax.hlines(y,z.ci_low,z.ci_high,color=COLORS[j],lw=.6)
                ax.plot(z.effect,y,'o',color=COLORS[j],ms=2.6,label=title)
            ax.set_yticks(range(16),genes,fontstyle='italic');ax.invert_yaxis();ax.axvline(0,color='.5',lw=.6);ax.set_title(label,loc='left');ax.set_xlabel('Late Braak versus 0 log2 fold change');ax.xaxis.set_major_locator(MaxNLocator(4))
        axes[0,0].legend(frameon=False,fontsize=6.5,loc='lower left')
        self.save(fig,'FigureS8',{'source':d},'QC age sensitivity in full, age-overlap, and spline-adjusted samples.','Moderated nominal 95% CIs; restricted age range changes the population.')
        d=pd.read_csv(self.review/'sex_interaction/target_sex_interactions.csv')
        fig,axes=plt.subplots(2,2,figsize=(8.4,8),layout='constrained',sharex=True)
        genes=sorted(d.target_symbol.unique())
        for ax,model,label in zip(axes.flat,MODELS,NAMES):
            z=d.loc[d.model.eq(model)].set_index('target_symbol').reindex(genes)
            if z.interaction_late_effect.isna().any():raise ValueError('Incomplete formal sex-interaction target panel')
            y=np.arange(16);ax.hlines(y,z.interaction_late_ci_low,z.interaction_late_ci_high,color=COLORS[0],lw=.8)
            ax.scatter(z.interaction_late_effect,y,s=22,edgecolors=COLORS[0],facecolors=[COLORS[0] if q<.05 else 'white' for q in z.interaction_late_q],zorder=3)
            ax.set_yticks(y,genes,fontstyle='italic');ax.invert_yaxis();ax.axvline(0,color='.5',lw=.6);ax.set_title(label,loc='left');ax.set_xlabel('Male − female late-Braak effect\n(log2 expression scale)');ax.xaxis.set_major_locator(MaxNLocator(4))
        self.save(fig,'FigureS13',{'source':d},'Formal Braak-by-sex interaction, not a comparison of subgroup significance.',
                  'Nominal moderated 95% CIs; filled symbols: one-df interaction q<0.05 across all 22,011 features. Omnibus three-df tests remain in source data. UTY excluded from this panel.')
        d=self.table('QC_region_targets.csv');d=d.loc[d.model.eq('baseline')]
        fig,axes=plt.subplots(1,3,figsize=(8.4,5.7),layout='constrained',sharex=True)
        for ax,region in zip(axes,['EC','FC','TC']):
            z=d.loc[d.region.eq(region)].set_index('gene_symbol').reindex(genes)
            self.forest(ax,z,genes,beta='effect',low='ci_low',high='ci_high');ax.set_title(region);ax.set_xlabel('Log2 fold change')
        self.save(fig,'FigureS14',{'source':d},'Region-specific effects from the QC baseline region-interaction model.','Nominal moderated 95% CIs; formal region tests use the full background.')
        external=[]
        for cohort in ['GSE44770','GSE106241']:
            z=self.table(cohort+'_QC_gene_comparison.csv');z['cohort']=cohort
            if 'model' not in z:z['model']='AD_vs_control'
            external.append(z)
        d=pd.concat(external,ignore_index=True);panels=list(d.groupby(['cohort','model'],sort=False))
        values=np.full((16,len(panels)),np.nan)
        for j,(_,part) in enumerate(panels):
            z=part.set_index('gene_symbol').reindex(genes);values[:,j]=np.where(z.log2FC.notna(),np.sign(z.log2FC*z.QC_baseline_log2FC),np.nan)
        fig,ax=plt.subplots(figsize=(7.2,5.6),layout='constrained')
        cmap=matplotlib.colors.ListedColormap([COLORS[1],'#d8d8d8',COLORS[0]]);cmap.set_bad('#f1f1f1')
        ax.imshow(values,cmap=cmap,vmin=-1,vmax=1,aspect='auto')
        labels={'AD_vs_control':'AD versus control','age_sex_APOE':'Age + sex + APOE','age_sex_APOE_RNAquality':'Age + sex + APOE\n+ RNA quality'}
        ax.set_yticks(range(16),genes,fontstyle='italic');ax.set_xticks(range(len(panels)),[key[0]+'\n'+labels.get(key[1],key[1].replace('_',' ')) for key,_ in panels],fontsize=7)
        for i in range(16):
            for j in range(len(panels)):ax.text(j,i,'NA' if np.isnan(values[i,j]) else ('+' if values[i,j]>0 else '−'),ha='center',va='center',color='white' if np.isfinite(values[i,j]) else '.5',fontsize=8)
        ax.set_title('Direction relative to QC baseline: + same; − opposite; NA unmeasured',fontsize=9)
        self.save(fig,'FigureS15',{'source':d},'External gene directions relative to the QC baseline, retaining unmeasured and discordant results.',
                  'Direction symbols are not significance labels. Full external effect estimates and distinct q families remain in source data.')
        d=self.read('targeted/qc_joint/GO_all_models.csv');rows=[]
        for model,g in d.groupby('model',sort=False):
            for direction in ['Down','Up']:rows.append(dict(model=model,direction=direction,n=int((g.q.lt(.05)&g.direction.eq(direction)).sum()),tested=len(g)))
        counts=pd.DataFrame(rows);models=list(counts.model.unique());fig,ax=plt.subplots(figsize=(7.2,5.4),layout='constrained');left=np.zeros(len(models))
        for j,direction in enumerate(['Down','Up']):
            n=counts.loc[counts.direction.eq(direction)].set_index('model').loc[models,'n'].to_numpy();ax.barh(range(len(models)),n,left=left,color=COLORS[j],label=direction+'-ranked');left+=n
        labels=[composition_label(x) for x in models]
        ax.set_yticks(range(len(models)),labels,fontsize=7);ax.invert_yaxis();ax.set_xlabel('GO sets with within-model q < 0.05 (all 4,826 tested)');ax.legend(frameon=False,fontsize=7)
        self.save(fig,'FigureS18',{'counts':counts},'All ten QC composition/preprocessing model specifications; no significance-selected model display.')
    def qc(self):
        d=self.table('QC_nonY_targets.csv');male=self.table('UTY_male_only.csv');genes=sorted(d.gene_symbol.unique())
        fig=plt.figure(figsize=(7.5,8),layout='constrained');grid=fig.add_gridspec(2,1,height_ratios=[4,1]);ax=fig.add_subplot(grid[0])
        for j,model in enumerate(MODELS):
            z=d.loc[d.model.eq(model)].set_index('gene_symbol').loc[genes];y=np.arange(16)+(j-1.5)*.17
            ax.hlines(y,z.ci_low,z.ci_high,color=COLORS[j],lw=.8);ax.plot(z.log2FC,y,'o',color=COLORS[j],ms=3,label=NAMES[j])
        ax.set_yticks(range(16),genes,fontstyle='italic');ax.invert_yaxis();ax.axvline(0,color='.5',lw=.7)
        ax.set(title='a  Sixteen non-Y targets; 212 retained arrays',xlabel='Late Braak versus 0 log2 fold change');ax.legend(frameon=False,ncol=2,fontsize=7)
        ax=fig.add_subplot(grid[1]);z=male.set_index('model').loc[MODELS]
        self.forest(ax,z,NAMES,beta='effect',low='ci_low',high='ci_high');ax.set(title='b  UTY, male-only analysis',xlabel='Late Braak versus 0 log2 fold change')
        self.save(fig,'FigureS21',{'nonY':d,'UTY_male_only':male},'QC pooled non-Y estimates and a separate male-only UTY panel.','Moderated nominal 95% CIs; full-feature BH within each model and population.')
    def calibration(self):
        d=self.reviewed('prediction/calibration_bins.csv');u=self.reviewed('prediction/calibration_uncertainty.csv')
        panels=list(d.component.unique());fig,axes=plt.subplots(2,2,figsize=(7.5,6.2),layout='constrained')
        for ax,component in zip(axes.flat,panels):
            g=d.loc[d.component.eq(component)&d.family.eq('fractional_logit')];lim=max(1,float(max(g.predicted.max(),g.observed.max()))*1.12)
            ax.plot([0,lim],[0,lim],color='.5',ls=':',lw=.8)
            for model,color,marker,label in [('clinical_amyloid_mri',COLORS[2],'o','Without tau'),('clinical_amyloid_mri_tau',COLORS[0],'s','With tau')]:
                z=g.loc[g.model.eq(model)];ax.plot(z.predicted,z.observed,color=color,marker=marker,lw=1,label=label)
            title=component.replace('landmark_','Day ').replace('_participant',' | Participant CV').replace('_site',' | Site CV')
            n=int(g.loc[g.model.eq('clinical_amyloid_mri')].n.sum())
            ax.set(title=title+f'\nn = {n}',xlabel='Mean predicted CDR-SB',ylabel='Mean observed CDR-SB',xlim=(0,lim),ylim=(0,lim));ax.set_aspect('equal');ax.legend(frameon=False,fontsize=7)
        self.save(fig,'FigureS22',{'bins':d,'uncertainty':u},'Held-out calibration using participant- or site-grouped internal validation.',
                  'Five equal-count descriptive bins; no bin-level CI. Refitted bootstrap intervals and valid counts are in source data. Day 0/180 specify permitted baseline information timing.')
    def run(self):
        plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.spines.right':False,'axes.spines.top':False,'pdf.fonttype':42,'svg.fonttype':'none'})
        self.original();self.genes();self.trajectories();self.imported();self.qc_supplements();self.qc();self.calibration()
        expected={f'Figure{i}' for i in range(1,5)}|{f'FigureS{i}' for i in range(1,23)}
        if {r['figure'] for r in self.records}!=expected:raise ValueError('Final figure inventory mismatch')
        dump(self.new/'figure_manifest.json',self.records)
        return self.records
