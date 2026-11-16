"""Validate and lift full OCD GWAS; cache mapped SNPs for repeated analysis."""
import json

import pandas as pd

from compare_ocd import DATA, OUT, CHAIN, chain_segments, read_gwas

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    mapped,audit,samples=read_gwas(chain_segments(CHAIN))
    folder=DATA/'generated/mapped_grch38'
    folder.mkdir(exist_ok=True)
    for chrom,(pos,p,snp) in mapped.items():
        pd.DataFrame({'pos':pos,'p':p,'snp':snp}).to_parquet(folder/f'chr{chrom}.parquet',index=False)
    (OUT/'gwas_liftover_audit.json').write_text(json.dumps({'counts':audit,'samples':samples},indent=2))
    print(json.dumps(audit,indent=2))

if __name__=='__main__':main()
