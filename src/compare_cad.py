"""Compare harmonised CAD GWAS signal near the established dog OCD genes.

Run after compare_ocd.py: python3 src/compare_cad.py
"""
import hashlib
import json
from collections import Counter

import numpy as np
import pandas as pd

try:
    from .compare_ocd import OUT, ROOT, plot, rank_and_permute, window_min
except ImportError:  # direct execution from src/
    from compare_ocd import OUT, ROOT, plot, rank_and_permute, window_min

CAD = ROOT / 'data/downloaded/gwas/GCST90132314.h.tsv.gz'
CAD_MD5 = '82e19c1817c52d5a94c8f42d9f48ccd0'
AUTOSOMES = {str(i) for i in range(1, 23)}
WINDOWS = (('1mb', 1000000, ''), ('500kb', 500000, '_500kb'))


def verify_checksum(path):
    digest = hashlib.md5()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    actual = digest.hexdigest()
    if actual != CAD_MD5:
        raise ValueError(f'CAD checksum mismatch: {actual}')
    return actual


def read_cad_gene_windows(path, genes, chunksize=300000):
    """Scan the chromosome-sorted GRCh38 file, retaining one chromosome at a time."""
    by_chrom = {chrom: group for chrom, group in genes.groupby('chrom')}
    results = {key: {} for key, _, _ in WINDOWS}
    audit = Counter()
    current_chrom = None
    chunks = []
    last_pos = 0
    seen_chroms = set()

    def finish_chrom():
        if current_chrom is None or not chunks:
            return
        all_rows = pd.concat(chunks, ignore_index=True)
        positions = all_rows.pos.to_numpy(dtype=np.int64)
        pvalues = all_rows.p.to_numpy(dtype=float)
        snps = all_rows.snp.to_numpy(dtype=str)
        for row in by_chrom.get(current_chrom, pd.DataFrame()).itertuples(index=False):
            for key, width, _ in WINDOWS:
                results[key][row.gene_id] = window_min(row.tss, positions, pvalues, snps, width)

    columns = ['chromosome', 'base_pair_location', 'p_value', 'rsid',
               'effect_allele', 'other_allele']
    for chunk in pd.read_csv(path, sep='\t', compression='gzip', usecols=columns,
                             dtype=str, chunksize=chunksize, keep_default_na=False):
        audit['total'] += len(chunk)
        for chrom, sub in chunk.groupby('chromosome', sort=False):
            chrom = str(chrom)
            if chrom not in AUTOSOMES:
                audit['nonautosomal'] += len(sub)
                continue
            if chrom != current_chrom:
                finish_chrom()
                chunks = []
                last_pos = 0
                if chrom in seen_chroms:
                    raise ValueError(f'CAD chromosome {chrom} is not grouped contiguously')
                seen_chroms.add(chrom)
                current_chrom = chrom
            bp = pd.to_numeric(sub.base_pair_location, errors='coerce').to_numpy(dtype=float)
            p = pd.to_numeric(sub.p_value, errors='coerce').to_numpy(dtype=float)
            valid = (np.isfinite(bp) & (bp >= 1) & (bp == np.floor(bp)) &
                     np.isfinite(p) & (p >= 0) & (p <= 1))
            audit['invalid_fields'] += int((~valid).sum())
            if not valid.any():
                continue
            bp = bp[valid].astype(np.int64)
            p = p[valid]
            if (bp[0] < last_pos or np.any(np.diff(bp) < 0)):
                raise ValueError(f'CAD chromosome {chrom} is not position-sorted')
            last_pos = int(bp[-1])
            sub = sub.iloc[np.flatnonzero(valid)]
            rsid = sub.rsid.to_numpy(dtype=str)
            fallback = (np.char.add(np.char.add(np.char.add(chrom + ':', bp.astype(str)), ':'),
                                   np.char.add(sub.effect_allele.to_numpy(dtype=str),
                                               np.char.add(':', sub.other_allele.to_numpy(dtype=str)))))
            snp = np.where((rsid != '') & (rsid != 'NA'), rsid, fallback)
            chunks.append(pd.DataFrame({'pos': bp, 'p': p, 'snp': snp}))
            audit['valid_autosomal'] += len(bp)
            audit['p_zero'] += int((p == 0).sum())
    finish_chrom()
    if audit['valid_autosomal'] == 0:
        raise ValueError('No valid autosomal CAD variants')
    out = genes[['gene_id']].copy()
    for key, _, suffix in WINDOWS:
        vals = [results[key].get(gid, (np.nan, '', np.nan, 0)) for gid in out.gene_id]
        out[[f'min_p{suffix}', f'lead_snp{suffix}',
             f'lead_pos_grch38{suffix}', f'n_snps{suffix}']] = pd.DataFrame(vals, index=out.index)
    return out, dict(audit)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    checksum = verify_checksum(CAD)
    ocd = pd.read_csv(OUT / 'per_gene_ocd.csv', dtype={'chrom': str})
    score, audit = read_cad_gene_windows(CAD, ocd[['gene_id', 'chrom', 'tss']])
    cad = ocd[['gene_id', 'symbol', 'chrom', 'start', 'end', 'strand', 'tss',
               'dog', 'dog_significant', 'dog_suggestive']].merge(
                   score, on='gene_id', validate='one_to_one')
    if not cad.n_snps_500kb.le(cad.n_snps).all():
        raise ValueError('A CAD 500 kb window contains more variants than its 1 Mb window')
    covered = cad.n_snps_500kb.gt(0)
    if not cad.loc[covered, 'min_p_500kb'].ge(cad.loc[covered, 'min_p']).all():
        raise ValueError('A CAD 500 kb minimum p-value is below its 1 Mb minimum')
    summary = {'source': 'GCST90132314.h.tsv.gz', 'md5': checksum,
               'coordinate_system': '1-based GRCh38', 'gwas': audit,
               'genes_total': len(cad), 'dog_eligible': int(cad.dog.sum()), 'windows': {}}
    for key, width, suffix in WINDOWS:
        ranked, tests = rank_and_permute(cad, p_col=f'min_p{suffix}', count_col=f'n_snps{suffix}')
        cad = cad.merge(ranked[['gene_id', 'snp_stratum', 'signal_rank']].rename(
            columns={'snp_stratum': f'snp_stratum{suffix}',
                     'signal_rank': f'signal_rank{suffix}'}), on='gene_id',
            how='left', validate='one_to_one')
        covered = cad[f'n_snps{suffix}'].gt(0)
        summary['windows'][key] = {
            'width_each_side_bp': width,
            'genes_covered': int(covered.sum()),
            'dog_covered': int((cad.dog & covered).sum()),
            'dog_significant_covered': int((cad.dog_significant & covered).sum()),
            'tests': tests}
    cad.to_csv(OUT / 'per_gene_cad.csv', index=False)
    (OUT / 'cad_summary.json').write_text(json.dumps(summary, indent=2))
    plot(ocd, OUT / 'ocd_p_distributions.png', cad=cad)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
