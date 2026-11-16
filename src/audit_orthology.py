"""Cross-check each ST6 dog gene with Ensembl dog-to-human orthology."""
import gzip
import json
import re
import time
from io import StringIO
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'notebook/results'
CACHE=ROOT/'data/downloaded/orthology'

def fetch(symbol):
    CACHE.mkdir(parents=True,exist_ok=True)
    path=CACHE/f'{symbol}.json'
    if path.exists():
        return json.loads(path.read_text())
    url=f'https://rest.ensembl.org/homology/symbol/canis_lupus_familiaris/{symbol}'
    for attempt in range(5):
        try:
            r=requests.get(url,params={'target_species':'homo_sapiens','type':'orthologues','sequence':'none'},headers={'Accept':'application/json'},timeout=45)
            if r.ok:
                value=r.json()
                path.write_text(json.dumps(value))
                return value
        except (requests.RequestException,ValueError):
            pass
        time.sleep(2*(attempt+1))
    return None

def biomart(names):
    path=CACHE/'dog_human_biomart.tsv'
    if not path.exists():
        query=('<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">'
               '<Dataset name="clfamiliaris_gene_ensembl" interface="default">'
               f'<Filter name="external_gene_name" value="{",".join(names)}"/>'
               '<Attribute name="external_gene_name"/><Attribute name="ensembl_gene_id"/>'
               '<Attribute name="hsapiens_homolog_ensembl_gene"/>'
               '<Attribute name="hsapiens_homolog_orthology_type"/></Dataset></Query>')
        r=requests.get('https://www.ensembl.org/biomart/martservice',params={'query':query},timeout=120)
        r.raise_for_status()
        if not r.text.startswith('Gene name\t'):
            raise ValueError(r.text[:200])
        path.write_text(r.text)
    return pd.read_csv(path,sep='\t',dtype=str).fillna('')

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    dog=pd.read_excel(ROOT/'data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx',sheet_name='ST6',header=1)
    dog=dog[dog.Species.eq('Dog')].copy()
    names=dog['Gene*'].tolist()
    st5=pd.read_excel(ROOT/'data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx',sheet_name='ST5',header=1)
    st5_names={part.strip() for value in st5['Genes (NetColoc)'].dropna() for part in str(value).split(',')}
    st5_names.discard('no genes')
    pd.DataFrame({'Gene*':sorted(set(names)|st5_names)}).assign(in_st5=lambda x:x['Gene*'].isin(st5_names),in_st6=lambda x:x['Gene*'].isin(names)).to_csv(OUT/'st5_st6_gene_audit.csv',index=False)
    symbol_ids={}
    with gzip.open(ROOT/'data/downloaded/annotation/gencode.v50.annotation.gtf.gz','rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            a=line.split('\t')
            if len(a)<9 or a[2]!='gene' or 'gene_type "protein_coding"' not in a[8]:continue
            attrs=dict(re.findall(r'(\w+) "([^"]+)"',a[8]))
            symbol_ids.setdefault(attrs['gene_name'],set()).add(attrs['gene_id'].split('.')[0])
    table=biomart(names+['NALF1'])
    byname={name:table[table['Gene name'].eq(name)] for name in names+['NALF1']}
    missing=[s for s in names if byname[s].empty]
    records={}
    with ThreadPoolExecutor(max_workers=4) as ex:
        futures={ex.submit(fetch,s):s for s in missing}
        for f in as_completed(futures):
            records[futures[f]]=f.result()
    out=[]
    for row in dog.itertuples(index=False):
        symbol=row[0]
        expected='NALF1' if symbol=='FAM155A' else symbol
        ids=symbol_ids.get(expected,set())
        if symbol=='FAM155A': ids={'ENSG00000204442'} if 'ENSG00000204442' in ids else set()
        response=records.get(symbol)
        hom=[] if response is None else [h for datum in response.get('data',[]) for h in datum.get('homologies',[]) if h.get('target',{}).get('species')=='homo_sapiens']
        biomart_rows=byname['NALF1' if symbol=='FAM155A' and byname[symbol].empty else symbol]
        orth_ids=set(biomart_rows.loc[biomart_rows['Human homology type'].eq('ortholog_one2one'),'Human gene stable ID']) | {h['target']['id'] for h in hom if h.get('type')=='ortholog_one2one'}
        all_ids=set(biomart_rows['Human gene stable ID']) | {h['target']['id'] for h in hom}
        all_ids.discard('')
        toga=str(row[2]).strip()
        if toga!='1:1': status='excluded_non_one_to_one'
        elif biomart_rows.empty and response is None: status='orthology_unavailable'
        elif len(ids)!=1: status='symbol_unresolved'
        elif ids==orth_ids: status='confirmed'
        else: status='orthology_conflict'
        out.append({'Gene*':symbol,'TOGA':toga,'human_symbol':expected,'human_gene_id':next(iter(ids)) if len(ids)==1 else '', 'ensembl_one_to_one_ids':';'.join(sorted(orth_ids)),'ensembl_all_human_ids':';'.join(sorted(all_ids)),'mapping_status':status,'ensembl_url':f'https://rest.ensembl.org/homology/symbol/canis_lupus_familiaris/{symbol}?target_species=homo_sapiens;type=orthologues'})
    df=pd.DataFrame(out)
    df.to_csv(OUT/'orthology_audit.csv',index=False)
    print(df.mapping_status.value_counts().to_string())
    print(df.loc[df.mapping_status.ne('confirmed'),['Gene*','human_symbol','human_gene_id','ensembl_one_to_one_ids','mapping_status']].to_string(index=False))

if __name__=='__main__':main()
