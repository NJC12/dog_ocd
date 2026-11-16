"""Score the supplied GRCh37 psychiatric GWAS in GRCh38 TSS windows.

The reader accepts whitespace-delimited files in any chromosome/position order.
This reader uses the same liftover and inclusive windows as the OCD analysis.
"""

import gzip
import hashlib
import json
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from .compare_ocd import chain_segments, lift_positions, window_min


PSYCHIATRIC_GWAS = {
    'BIP': {
        'name': 'Bipolar disorder',
        'filename': 'bip2024_multianc_no23andMe.gz',
        'md5': '6615ee8a7c6443df124848d17ba1336e',
        'source': 'https://figshare.com/articles/dataset/bip2024/27216117?file=49760775',
    },
    'ADHD': {
        'name': 'ADHD',
        'filename': 'ADHD2022_iPSYCH_deCODE_PGC.meta.gz',
        'md5': 'ea57120c720e880e34f67f57483dcbb4',
        'source': 'https://figshare.com/articles/dataset/adhd2022/22564390',
    },
    'MDD': {
        'name': 'Major depressive disorder',
        'filename': 'pgc-mdd2025_no23andMe_div_v3-49-46-01.tsv.gz',
        'md5': '919829e44688fe7b0d745199a7893746',
        'source': 'https://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255?file=51487010',
        'columns': {'#CHROM': 'CHR', 'POS': 'BP', 'ID': 'SNP', 'PVAL': 'P'},
    },
    'SCZ': {
        'name': 'Schizophrenia',
        'filename': 'PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv.gz',
        'md5': 'e186621c9c290605ad7e19968bfc3dfc',
        'source': 'https://figshare.com/articles/dataset/scz2022/19426775?file=34517861',
        'columns': {'CHROM': 'CHR', 'POS': 'BP', 'ID': 'SNP', 'PVAL': 'P'},
    },
}
AUTOSOMES = {str(i) for i in range(1, 23)}
WINDOWS = (('1mb', 1_000_000, ''), ('500kb', 500_000, '_500kb'))


def file_md5(path):
    digest = hashlib.md5()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_psychiatric_chunks(path, columns=None, chunksize=300_000, nrows=None):
    """Read data after any leading ## metadata, preserving a #CHROM header."""
    columns = columns or {name: name for name in ['CHR', 'SNP', 'BP', 'P']}
    skiprows = 0
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if not line.startswith('##'):
                break
            skiprows += 1
    for chunk in pd.read_csv(path, sep=r'\s+', usecols=list(columns),
                             skiprows=skiprows, dtype=str, keep_default_na=False,
                             chunksize=chunksize, nrows=nrows):
        yield chunk.rename(columns=columns)


def read_psychiatric_gene_windows(path, genes, segments, chunksize=300_000,
                                 nrows=None, columns=None):
    """Validate, lift, sort and score a whitespace-delimited CHR/SNP/BP/P file.

    Reject nonpositive/nonfinite p-values as in the original OCD reader. Keep
    every valid mapped row (no MAF, INFO, LD, or ancestry filtering is added).
    ``nrows`` is only for a pilot and must be None for the full analysis.
    """
    batches = defaultdict(list)
    counts = Counter({key: 0 for key in (
        'total', 'nonautosomal', 'invalid_fields', 'p_zero', 'valid_autosomal',
        'unmapped', 'ambiguous', 'mapped_nonautosomal_or_invalid',
        'cross_chromosome_mapped', 'mapped')})
    samples = []
    for chunk in read_psychiatric_chunks(path, columns, chunksize, nrows):
        counts['total'] += len(chunk)
        for chrom, sub in chunk.groupby('CHR', sort=False):
            if chrom not in AUTOSOMES:
                counts['nonautosomal'] += len(sub)
                continue
            bp = pd.to_numeric(sub.BP, errors='coerce').to_numpy(dtype=float)
            p = pd.to_numeric(sub.P, errors='coerce').to_numpy(dtype=float)
            valid = (np.isfinite(bp) & (bp >= 1) & (bp == np.floor(bp)) &
                     np.isfinite(p) & (p > 0) & (p <= 1))
            counts['p_zero'] += int((p == 0).sum())
            counts['invalid_fields'] += int((~valid).sum())
            counts['valid_autosomal'] += int(valid.sum())
            if not valid.any():
                continue
            sub = sub.iloc[np.flatnonzero(valid)]
            bp = bp[valid].astype(np.int64)
            p = p[valid]
            dest, pos, status = lift_positions(chrom, bp, segments)
            counts['unmapped'] += int((status == 0).sum())
            counts['ambiguous'] += int((status == 2).sum())
            good = (status == 1) & np.isin(dest, list(AUTOSOMES)) & (pos >= 1)
            counts['mapped_nonautosomal_or_invalid'] += int(((status == 1) & ~good).sum())
            counts['cross_chromosome_mapped'] += int((good & (dest != chrom)).sum())
            counts['mapped'] += int(good.sum())
            for i in np.flatnonzero(good)[:max(0, 20 - len(samples))]:
                samples.append({'chrom': chrom, 'grch37': int(bp[i]),
                                'grch38_chrom': str(dest[i]), 'grch38': int(pos[i]),
                                'snp': sub.SNP.iloc[i]})
            snps = sub.SNP.to_numpy()
            for target in np.unique(dest[good]):
                take = good & (dest == target)
                batches[target].append(pd.DataFrame({
                    'pos': pos[take], 'p': p[take], 'snp': snps[take]}))
    if counts['mapped'] == 0:
        raise ValueError(f'No valid mapped autosomal variants in {path}')
    assert counts['total'] == (counts['nonautosomal'] + counts['invalid_fields'] +
                               counts['valid_autosomal'])
    assert counts['valid_autosomal'] == (counts['mapped'] + counts['unmapped'] +
        counts['ambiguous'] + counts['mapped_nonautosomal_or_invalid'])

    results = {key: {} for key, _, _ in WINDOWS}
    for chrom, chrom_genes in genes.groupby('chrom'):
        parts = batches.pop(str(chrom), [])
        if not parts:
            continue
        variants = pd.concat(parts, ignore_index=True).sort_values(
            ['pos', 'p'], kind='mergesort')
        positions = variants.pos.to_numpy(dtype=np.int64)
        pvalues = variants.p.to_numpy(dtype=float)
        snps = variants.snp.to_numpy()
        for gene in chrom_genes.itertuples(index=False):
            for key, width, _ in WINDOWS:
                results[key][gene.gene_id] = window_min(
                    gene.tss, positions, pvalues, snps, width)
    scores = genes[['gene_id']].copy()
    for key, _, suffix in WINDOWS:
        values = [results[key].get(gid, (np.nan, '', np.nan, 0)) for gid in scores.gene_id]
        scores[[f'min_p{suffix}', f'lead_snp{suffix}',
                f'lead_pos_grch38{suffix}', f'n_snps{suffix}']] = pd.DataFrame(
                    values, index=scores.index)
    return scores, {'counts': dict(counts), 'samples': samples}


def load_psychiatric_gene_windows(trait, root, genes, chain, output, rebuild=False):
    """Reuse verified gene scores or build them; never reuse stale permutations."""
    spec = PSYCHIATRIC_GWAS[trait]
    path = root / 'data/downloaded/gwas' / spec['filename']
    checksum = file_md5(path)
    if checksum != spec['md5']:
        raise ValueError(f'{trait} checksum mismatch: {checksum}')
    chain_checksum = file_md5(chain)
    score_path = output / f'per_gene_{trait.lower()}.csv'
    audit_path = output / f'{trait.lower()}_summary.json'
    if not rebuild and score_path.exists() and audit_path.exists():
        summary = json.loads(audit_path.read_text())
        if (summary.get('md5') == checksum and
                summary.get('chain_md5') == chain_checksum and
                summary.get('scores_md5') == file_md5(score_path)):
            frame = pd.read_csv(score_path, dtype={'chrom': str})
            return frame, summary, False
    print(f'Building {trait} gene windows from {path.name}', flush=True)
    scores, audit = read_psychiatric_gene_windows(
        path, genes, chain_segments(chain), columns=spec.get('columns'))
    frame = genes.merge(scores, on='gene_id', validate='one_to_one')
    output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(score_path, index=False)
    summary = {'source': spec['source'], 'filename': spec['filename'],
               'md5': checksum, 'chain_md5': chain_checksum,
               'scores_md5': file_md5(score_path),
               'input_coordinates': '1-based GRCh37 (hg19)',
               'coordinate_system': '1-based GRCh38', 'gwas': audit,
               'genes_total': len(frame), 'dog_eligible': int(frame.dog.sum()),
               'windows': {}}
    audit_path.write_text(json.dumps(summary, indent=2))
    return frame, summary, True
