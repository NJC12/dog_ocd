"""Compare human OCD GWAS signal around dog GWAS ortholog TSSs.

Run from the project root: python3 src/compare_ocd.py
"""
import gzip
import hashlib
import heapq
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "notebook/results"
CHAIN = Path.home() / "data/liftover/hg19ToHg38.over.chain"
GWAS = DATA / "downloaded/gwas/daner_OCD_full_wo23andMe_190522.gz"
GTF = DATA / "downloaded/annotation/gencode.v50.annotation.gtf.gz"
SEED = 20261013
PLOT_MAX_LOGP = 15


def chain_segments(path):
    """Nonoverlapping hg19 0-based segments; overlapping chains are ambiguous."""
    events = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("chain "):
                a = line.split()
                tchr, ts, qchr, qs, qsize, strand = a[2], int(a[5]), a[7], int(a[10]), int(a[8]), a[9]
                continue
            a = list(map(int, line.split()))
            size = a[0]
            if tchr.startswith('chr') and tchr[3:] in {str(i) for i in range(1,23)}:
                # q coordinates in the chain are relative to its strand.
                if strand == '+':
                    offset = qs - ts
                    sign = 1
                else:
                    offset = qsize - qs - 1 + ts
                    sign = -1
                block = (tchr[3:], qchr.removeprefix('chr'), sign, offset)
                events[tchr[3:]].append((ts, 1, block))
                events[tchr[3:]].append((ts + size, -1, block))
            if len(a) == 3:
                ts += size + a[1]
                qs += size + a[2]
    out = {}
    for chrom, ev in events.items():
        ev.sort(key=lambda x: x[0])
        active = Counter()
        seg = []
        prev = None
        for pos, group in __import__('itertools').groupby(ev, key=lambda x:x[0]):
            if prev is not None and pos > prev and active:
                mapping = next(iter(active)) if sum(active.values()) == 1 else None
                seg.append((prev, pos, mapping))
            for _, typ, block in group:
                active[block] += typ
                if active[block] == 0:
                    del active[block]
            prev = pos
        out[chrom] = seg
    return out


def lift_positions(chrom, bp, segments):
    """Map 1-based GRCh37 SNP BP to 1-based GRCh38; status 0/1/2 = unmapped/mapped/ambiguous."""
    seg = segments.get(str(chrom), [])
    starts = np.array([x[0] for x in seg], dtype=np.int64)
    ends = np.array([x[1] for x in seg], dtype=np.int64)
    idx = np.searchsorted(starts, np.asarray(bp, dtype=np.int64)-1, side='right')-1
    valid = idx >= 0
    safe = np.maximum(idx, 0)
    valid &= np.asarray(bp)-1 < ends[safe] if len(ends) else False
    status = np.zeros(len(bp), dtype=np.int8)
    mapped = np.zeros(len(bp), dtype=np.int64)
    dest = np.full(len(bp), '', dtype='<U3')
    if not len(seg):
        return dest, mapped, status
    signs=np.array([s[2][2] if s[2] else 0 for s in seg],dtype=np.int8)
    offsets=np.array([s[2][3] if s[2] else 0 for s in seg],dtype=np.int64)
    targets=np.array([s[2][1] if s[2] else '' for s in seg],dtype='<U3')
    status[valid]=np.where(signs[safe[valid]]==0,2,1)
    good=status==1
    dest[good]=targets[safe[good]]
    mapped[good]=signs[safe[good]]*(np.asarray(bp)[good]-1)+offsets[safe[good]]+1
    return dest, mapped, status


def read_genes(path):
    rows = []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            a = line.rstrip().split('\t')
            if len(a) != 9 or a[2] != 'gene' or a[0].removeprefix('chr') not in {str(i) for i in range(1,23)} or 'gene_type "protein_coding"' not in a[8]:
                continue
            attrs = dict(__import__('re').findall(r'(\w+) "([^"]+)"', a[8]))
            rows.append(dict(chrom=a[0].removeprefix('chr'), start=int(a[3]), end=int(a[4]), strand=a[6], gene_id=attrs['gene_id'].split('.')[0], symbol=attrs['gene_name'], tss=int(a[3]) if a[6]=='+' else int(a[4])))
    return pd.DataFrame(rows)


def read_dog(path):
    d = pd.read_excel(path, sheet_name='ST6', header=1)
    d = d[d['Species'].eq('Dog')].copy()
    d['TOGA'] = d['TOGA (Human:Dog)**'].astype(str).str.strip()
    return d


def read_gwas(segments):
    expected = '9ab61cb7412c42a0ba207856af6da0f8'
    h = hashlib.md5()
    with open(GWAS, 'rb') as f:
        for b in iter(lambda:f.read(8*1024*1024), b''):
            h.update(b)
    if h.hexdigest() != expected:
        raise ValueError('GWAS MD5 mismatch')
    batches = defaultdict(list)
    audit = Counter()
    samples = []
    for chunk in pd.read_csv(GWAS, sep='\t', compression='gzip', usecols=['CHR','SNP','BP','P'], chunksize=300000, dtype={'CHR':str,'SNP':str,'BP':str,'P':str}):
        audit['total'] += len(chunk)
        for chrom, sub in chunk.groupby('CHR'):
            if chrom not in {str(i) for i in range(1,23)}:
                audit['nonautosomal'] += len(sub)
                continue
            bp = pd.to_numeric(sub.BP, errors='coerce').to_numpy(dtype=float)
            p = pd.to_numeric(sub.P, errors='coerce').to_numpy(dtype=float)
            valid = np.isfinite(bp) & (bp >= 1) & (bp == np.floor(bp)) & np.isfinite(p) & (p > 0) & (p <= 1)
            audit['invalid_fields'] += int((~valid).sum())
            if not valid.any():
                continue
            sub = sub.iloc[np.flatnonzero(valid)]
            bp = bp[valid].astype(np.int64)
            p = p[valid]
            dest, pos, status = lift_positions(chrom, bp, segments)
            audit['unmapped'] += int((status==0).sum())
            audit['ambiguous'] += int((status==2).sum())
            good = (status==1) & np.isin(dest,list(map(str,range(1,23)))) & (pos>=1)
            audit['cross_chromosome_mapped'] += int((good & (dest!=chrom)).sum())
            audit['mapped_nonautosomal_or_invalid'] += int(((status==1) & ~good).sum())
            audit['mapped'] += int(good.sum())
            if len(samples)<20:
                ix=np.flatnonzero(good)[:20-len(samples)]
                samples.extend([dict(chrom=chrom, grch38_chrom=str(dest[i]), grch37=int(bp[i]), grch38=int(pos[i]), snp=sub.SNP.iloc[i]) for i in ix])
            for target in np.unique(dest[good]):
                take=good & (dest==target)
                batches[target].append(pd.DataFrame({'pos':pos[take], 'p':p[take], 'snp':sub.SNP.to_numpy()[take]}))
    bychr = {}
    for chrom, parts in batches.items():
        d=pd.concat(parts,ignore_index=True).sort_values(['pos','p'],kind='mergesort')
        bychr[chrom]=(d.pos.to_numpy(),d.p.to_numpy(),d.snp.to_numpy())
    return bychr, dict(audit), samples


def window_min(tss, positions, pvalues, snps, width=1000000):
    left=max(1,tss-width); right=tss+width
    lo=np.searchsorted(positions,left,'left'); hi=np.searchsorted(positions,right,'right')
    if lo==hi:
        return (np.nan,'',np.nan,0)
    j=lo+np.argmin(pvalues[lo:hi])
    return (pvalues[j],snps[j],positions[j],hi-lo)


def rank_and_permute(d, n=10000, p_col='min_p', count_col='n_snps'):
    covered=d[d[count_col].gt(0)].copy()
    covered['snp_stratum']=covered.groupby('chrom')[count_col].transform(lambda x:pd.qcut(x.rank(method='first'),min(5,len(x)),labels=False))
    covered['signal_rank']=covered.groupby(['chrom','snp_stratum'])[p_col].transform(lambda x:x.rank(ascending=True,method='average',pct=True))
    # Circular shifts act on the ordered gene labels within each chromosome.
    chroms={c:x.sort_values('tss') for c,x in covered.groupby('chrom')}
    rng=np.random.default_rng(SEED)
    out={}
    for name,maskcol in [('all','dog'),('significant','dog_significant')]:
        observed=np.mean(covered.loc[covered[maskcol],'signal_rank'])
        null=np.empty(n)
        for k in range(n):
            vals=[]
            for c,x in chroms.items():
                labels=np.roll(x[maskcol].to_numpy(),rng.integers(len(x)))
                vals.extend(x.signal_rank.to_numpy()[labels])
            null[k]=np.mean(vals)
        out[name]={'n':int(covered[maskcol].sum()),'observed_mean_rank':float(observed),'null_mean':float(null.mean()),'empirical_one_sided_p':float((1+(null<=observed).sum())/(n+1)),'permutations':n}
    return covered, out


def plot(d, path, cad=None):
    windows=(('±1 Mb','min_p','n_snps'),('±500 kb','min_p_500kb','n_snps_500kb'))
    groups=(('All covered protein-coding genes',None,'#888888'),
            ('Dog GWAS significant','dog_significant','#d55e00'),
            ('Dog GWAS suggestive','dog_suggestive','#0072b2'))
    datasets=[('OCD',d)] + ([('CAD',cad)] if cad is not None else [])
    fig,axes=plt.subplots(2*len(datasets),3,figsize=(14.5,4*len(datasets)),sharex=True,sharey=False)
    bins=np.linspace(0,PLOT_MAX_LOGP,31)
    for trait,frame in datasets:
        trait_row=0 if trait=='OCD' else 2
        trait_max_density=0.
        for offset,(distance,p_col,count_col) in enumerate(windows):
            for col,(label,mask_col,color) in enumerate(groups):
                mask=frame[count_col].gt(0)
                if mask_col is not None:
                    mask &= frame[mask_col]
                x=-np.log10(np.maximum(frame.loc[mask,p_col].to_numpy(dtype=float),np.nextafter(0,1)))
                # Pool stronger signals in the final bin so every panel has
                # the same 0–15 x-axis without silently dropping those genes.
                x=np.minimum(x,PLOT_MAX_LOGP)
                ax=axes[trait_row+offset,col]
                heights,_,_=ax.hist(x,bins=bins,density=True,histtype='step',linewidth=2,color=color)
                trait_max_density=max(trait_max_density,float(heights.max()))
                ax.set_title(f'{trait} {distance} · {label}\n(n={len(x):,})')
                ax.set_xlim(0,PLOT_MAX_LOGP)
                ax.set_xticks([0,5,10,PLOT_MAX_LOGP])
                ax.set_xticklabels(['0','5','10',f'≥{PLOT_MAX_LOGP}'])
                ax.tick_params(axis='x',labelbottom=True)
        for ax in axes[trait_row:trait_row+2].flat:
            ax.set_ylim(0,trait_max_density*1.08)
    fig.supxlabel('−log₁₀(minimum GWAS p-value near TSS)')
    fig.supylabel('Gene density')
    fig.tight_layout(rect=(0.01,0.02,0.96,1));fig.savefig(path,dpi=200);plt.close(fig)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    seg=chain_segments(CHAIN)
    genes=read_genes(GTF)
    dog=read_dog(ROOT/'data/downloaded/supplementary/NATURE_SupplementaryTables.xlsx')
    orth=pd.read_csv(OUT/'orthology_audit.csv')
    dog=dog.merge(orth.drop(columns='TOGA'),on='Gene*',how='left',validate='one_to_one')
    cached=DATA/'generated/mapped_grch38'
    audit_path=OUT/'gwas_liftover_audit.json'
    if audit_path.exists() and len(list(cached.glob('chr*.parquet')))==22:
        bychr={}
        for chrom in map(str,range(1,23)):
            x=pd.read_parquet(cached/f'chr{chrom}.parquet')
            bychr[chrom]=(x.pos.to_numpy(),x.p.to_numpy(),x.snp.to_numpy())
        info=json.loads(audit_path.read_text())
        audit,samples=info['counts'],info['samples']
    else:
        bychr,audit,samples=read_gwas(seg)
    for item in samples:
        back=[]
        for source_chrom, pieces in seg.items():
            for start,end,block in pieces:
                if block is None or block[1]!=item['grch38_chrom']:
                    continue
                _,_,sign,offset=block
                zero=sign*(item['grch38']-1-offset)
                if start<=zero<end:
                    back.append((source_chrom,int(zero+1)))
        item['reverse_matches']=back
        item['reverse_unique']=back==[(item['chrom'],item['grch37'])]
    if not all(x['reverse_unique'] for x in samples):
        raise ValueError('Liftover reverse check failed')
    rows=[]
    rows_500kb=[]
    for row in genes.itertuples(index=False):
        pos,p,snp=bychr.get(row.chrom,(np.array([],dtype=int),np.array([]),np.array([])))
        rows.append(window_min(row.tss,pos,p,snp,width=1000000))
        rows_500kb.append(window_min(row.tss,pos,p,snp,width=500000))
    genes[['min_p','lead_snp','lead_pos_grch38','n_snps']]=pd.DataFrame(rows,index=genes.index)
    genes[['min_p_500kb','lead_snp_500kb','lead_pos_grch38_500kb','n_snps_500kb']]=pd.DataFrame(rows_500kb,index=genes.index)
    if not genes.n_snps_500kb.le(genes.n_snps).all():
        raise ValueError('A 500 kb window contains more SNPs than its 1 Mb window')
    if not genes.loc[genes.n_snps_500kb.gt(0),'min_p_500kb'].ge(genes.loc[genes.n_snps_500kb.gt(0),'min_p']).all():
        raise ValueError('A 500 kb minimum p-value is below its 1 Mb minimum')
    selected=dog[(dog.TOGA=='1:1') & dog.mapping_status.eq('confirmed')]
    ids=set(selected.human_gene_id)
    sigids=set(selected.loc[selected['dogCD significance'].eq('significant'),'human_gene_id'])
    sugids=set(selected.loc[selected['dogCD significance'].eq('suggestive'),'human_gene_id'])
    genes['dog']=genes.gene_id.isin(ids)
    genes['dog_significant']=genes.gene_id.isin(sigids)
    genes['dog_suggestive']=genes.gene_id.isin(sugids)
    ranked,stats=rank_and_permute(genes)
    ranked_500kb,stats_500kb=rank_and_permute(genes,p_col='min_p_500kb',count_col='n_snps_500kb')
    genes=genes.merge(ranked[['gene_id','snp_stratum','signal_rank']],on='gene_id',how='left',validate='one_to_one')
    genes=genes.merge(ranked_500kb[['gene_id','snp_stratum','signal_rank']].rename(columns={'snp_stratum':'snp_stratum_500kb','signal_rank':'signal_rank_500kb'}),on='gene_id',how='left',validate='one_to_one')
    genes.to_csv(OUT/'per_gene_ocd.csv',index=False)
    dog.to_csv(OUT/'dog_ortholog_decisions.csv',index=False)
    cad_path=OUT/'per_gene_cad.csv'
    cad=pd.read_csv(cad_path) if cad_path.exists() else None
    plot(genes,OUT/'ocd_p_distributions.png',cad=cad)
    summary={'gwas':audit,'genes_total':len(genes),'dog_eligible':len(selected),
             'windows':{
                 '1mb':{'width_each_side_bp':1000000,'genes_covered':int(genes.n_snps.gt(0).sum()),'dog_covered':int((genes.dog & genes.n_snps.gt(0)).sum()),'dog_significant_covered':int((genes.dog_significant & genes.n_snps.gt(0)).sum()),'tests':stats},
                 '500kb':{'width_each_side_bp':500000,'genes_covered':int(genes.n_snps_500kb.gt(0).sum()),'dog_covered':int((genes.dog & genes.n_snps_500kb.gt(0)).sum()),'dog_significant_covered':int((genes.dog_significant & genes.n_snps_500kb.gt(0)).sum()),'tests':stats_500kb}},
             'liftover_sample':samples}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    main()
