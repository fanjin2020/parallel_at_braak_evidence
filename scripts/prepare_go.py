"""Build current human GO BP sets from versioned ontology and GAF annotations."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
from collections import defaultdict
from datetime import datetime,timezone
from functools import lru_cache
from pathlib import Path

GAF='https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz'
GENE='https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz'


def ontology(path):
    terms={}; current=None; version=''
    with path.open(encoding='utf-8') as f:
        for line in f:
            line=line.strip()
            if line.startswith('data-version:'):version=line.split(': ',1)[1]
            if line=='[Term]':
                current={'parents':[],'alt':[],'obsolete':False};continue
            if line.startswith('['):current=None;continue
            if current is None:continue
            if line.startswith('id: GO:'):
                current['id']=line[4:];terms[current['id']]=current
            elif line.startswith('name: '):current['name']=line[6:]
            elif line.startswith('namespace: '):current['namespace']=line[11:]
            elif line.startswith('alt_id: '):current['alt'].append(line[8:])
            elif line.startswith('is_a: '):current['parents'].append(line.split()[1])
            elif line.startswith('relationship: part_of '):current['parents'].append(line.split()[2])
            elif line=='is_obsolete: true':current['obsolete']=True
    return terms,version


def run(out,rank_file,annotation_url=GAF,ontology_url=None,ontology_file=None,annotation_file=None,gene_info_file=None):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'go_gene_membership.csv').exists():raise ValueError('Output already built; choose a new directory')
    gaf=out/'human.gaf.gz';gene=out/'human_gene_info.gz';obo=out/'go-basic.obo'
    print('Reading matched local GO release and human gene dictionary',flush=True)
    for local,target,url in [(annotation_file,gaf,annotation_url),(gene_info_file,gene,GENE)]:
        if local:
            if Path(local).resolve()!=target.resolve():shutil.copy2(local,target)
        else:raise FileNotFoundError('The local GO and gene annotation snapshots are required.')
    headers=[]
    with gzip.open(gaf,'rt',encoding='utf-8') as f:
        for line in f:
            if not line.startswith('!'):break
            headers.append(line.strip())
    version_line=next((x for x in headers if x.lower().startswith('!go-version:')),'')
    match=re.search(r'(\d{4}-\d{2}-\d{2})',version_line)
    if ontology_url is None:
        if not match:raise ValueError('No ontology version in GAF header; provide --ontology-url for a documented matching release')
        ontology_url='https://release.geneontology.org/'+match[1]+'/ontology/go-basic.obo'
    if ontology_file:
        if Path(ontology_file).resolve()!=obo.resolve(): shutil.copy2(ontology_file,obo)
    else:
        raise FileNotFoundError('Supply the cached OBO file matching the GAF release.')
    terms,version=ontology(obo)
    if match and match[1] not in version:
        raise ValueError('GAF ontology release and supplied OBO differ; do not silently combine')
    alias={alt:k for k,v in terms.items() for alt in v['alt']}
    @lru_cache(maxsize=None)
    def ancestors(go):
        if go not in terms or terms[go]['obsolete']:return frozenset()
        result={go}
        for parent in terms[go]['parents']:result.update(ancestors(alias.get(parent,parent)))
        return frozenset(result)
    symbols=defaultdict(set);gene_info=[]
    with gzip.open(gene,'rt',encoding='utf-8') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['#tax_id']!='9606':continue
            symbols[r['Symbol']].add(r['GeneID'])
            gene_info.append(dict(entrez_id=r['GeneID'],gene_symbol=r['Symbol'],chromosome=r['chromosome']))
    with rank_file.open(encoding='utf-8-sig') as f:
        universe={r['entrez_id'] for r in csv.DictReader(f) if r['model']=='full'}
    sets=defaultdict(set);counts=defaultdict(int)
    with gzip.open(gaf,'rt',encoding='utf-8') as f:
        for line in f:
            if line.startswith('!'):continue
            a=line.rstrip('\n').split('\t');counts['annotation_rows']+=1
            if len(a)<15:raise ValueError('Malformed GAF row')
            if a[8]!='P' or 'NOT' in a[3].split('|') or a[12].split('|')[0]!='taxon:9606':continue
            ids=symbols.get(a[2],set())
            if len(ids)!=1:counts['unmapped_or_ambiguous_symbol']+=1;continue
            gene_id=next(iter(ids))
            if gene_id not in universe:continue
            go=alias.get(a[4],a[4])
            for ancestor in ancestors(go):
                if terms[ancestor].get('namespace')=='biological_process':sets[ancestor].add(gene_id)
    sets={k:v for k,v in sets.items() if 10<=len(v)<=300}
    if not sets:raise ValueError('No eligible GO sets mapped')
    for name,rows,fields in [
        ('go_gene_membership.csv',[dict(go_id=k,entrez_id=g) for k in sorted(sets) for g in sorted(sets[k])],['go_id','entrez_id']),
        ('go_terms.csv',[dict(go_id=k,name=terms[k]['name'],n_genes=len(sets[k])) for k in sorted(sets)],['go_id','name','n_genes']),
        ('human_gene_dictionary.csv',gene_info,['entrez_id','gene_symbol','chromosome'])]:
        with (out/name).open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    hashes={}
    for path in [gaf,gene,obo]:
        with path.open('rb') as source: hashes[path.name]=hashlib.file_digest(source,'sha256').hexdigest()
    provenance=dict(retrieved_utc=datetime.now(timezone.utc).isoformat(),annotation_url=annotation_url,ontology_url=ontology_url,
                    ontology_version=version,gaf_header=headers,gene_dictionary_url=GENE,
                    supplied_local_files=bool(ontology_file or annotation_file or gene_info_file),
                    sha256=hashes,
                    universe_genes=len(universe),tested_GO_sets=len(sets),annotation_counts=dict(counts),
                    mapping='Unique current official gene symbol to Entrez; ambiguous/unmapped excluded, no synonym guessing',
                    rules='Human BP; exclude NOT and obsolete terms before testing; propagate is_a and part_of; include all evidence codes; 10-300 genes in measured universe')
    snapshot_file=out/'quickgo_snapshot/quickgo_snapshot.json'
    if snapshot_file.exists():provenance['quickgo_graph_snapshot']=json.loads(snapshot_file.read_text(encoding='utf-8'))
    (out/'annotation_versions.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print('GO BP sets:',len(sets),'; ontology:',version,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--gene-ranks',type=Path,required=True);p.add_argument('--annotation-url',default=GAF);p.add_argument('--ontology-url');p.add_argument('--ontology-file',type=Path)
    p.add_argument('--annotation-file',type=Path);p.add_argument('--gene-info-file',type=Path)
    a=p.parse_args();run(a.output_dir,a.gene_ranks,a.annotation_url,a.ontology_url,a.ontology_file,a.annotation_file,a.gene_info_file)
