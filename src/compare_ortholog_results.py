"""Compare saved 78-gene Ensembl and 83-gene TOGA notebook results."""
from pathlib import Path
import pandas as pd
from toga_gene_set import adjust_ocd_tests

root = Path(__file__).resolve().parents[1] / 'notebook/results/ipynb'
new = root / 'toga'
mw = adjust_ocd_tests(pd.read_csv(root / 'mann_whitney_results.csv')).merge(
    pd.read_csv(new / 'mann_whitney_results.csv'),
    on=['trait', 'window', 'dog_set'], suffixes=('_78', '_83'), validate='one_to_one')
mw.to_csv(new / 'ortholog_comparison_mann_whitney.csv', index=False)
print('ALL MANN WHITNEY')
columns = ['trait', 'window', 'dog_set', 'probability_dog_stronger_78',
           'probability_dog_stronger_83', 'one_sided_p_78', 'one_sided_p_83',
           'holm_adjusted_p_78', 'holm_adjusted_p_83']
print(mw[columns].to_string(index=False))
print('\nHOLM SIGNIFICANCE CHANGES')
print(mw.loc[mw.holm_adjusted_p_78.lt(.05) != mw.holm_adjusted_p_83.lt(.05), columns].to_string(index=False))
print('\nNOMINAL SIGNIFICANCE CHANGES')
print(mw.loc[mw.one_sided_p_78.lt(.05) != mw.one_sided_p_83.lt(.05), columns].to_string(index=False))
shifts = pd.read_csv(root / 'circular_shift_results.csv').merge(
    pd.read_csv(new / 'circular_shift_results.csv'),
    on=['Trait', 'Window', 'Dog set'], suffixes=('_78', '_83'), validate='one_to_one')
shifts.to_csv(new / 'ortholog_comparison_circular_shifts.csv', index=False)
print('\nALL CIRCULAR SHIFTS')
print(shifts[['Trait', 'Window', 'Dog set', 'Empirical p_78', 'Empirical p_83']].to_string(index=False))
print('\nCIRCULAR SIGNIFICANCE CHANGES')
print(shifts.loc[shifts['Empirical p_78'].lt(.05) != shifts['Empirical p_83'].lt(.05),
                  ['Trait', 'Window', 'Dog set', 'Empirical p_78', 'Empirical p_83']].to_string(index=False))
print('\nOCD THRESHOLDS')
print(pd.read_csv(new / 'ocd_threshold_counts.csv').to_string(index=False))
restored = pd.read_csv(new / 'per_gene_ocd.csv')
print('\nRESTORED OCD GENES')
print(restored.loc[restored.symbol.isin(['EMC2', 'PPA2', 'PLXNA2', 'PROSER3', 'SH3GL3']),
                   ['symbol', 'min_p', 'min_p_500kb']].to_string(index=False))
