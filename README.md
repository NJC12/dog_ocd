# Dog OCD genes and human GWAS signal

The main analysis compares human GWAS associations near dog compulsive-disorder
gene orthologs across OCD and five control traits. The current dog set contains
83 ST6 TOGA one-to-one human orthologs; the Human set contains all 48 ST6 Human
genes. A separate constraint analysis uses the original ST6 groups. The OCD
comparisons are exploratory because ST6 selection used human OCD results.

## Run the notebook

Open [dog_ocd_analysis.ipynb](./dog_ocd_analysis.ipynb) at the project root.
From this directory:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r config/requirements-notebook.txt
.venv/bin/python src/execute_analysis_notebook.py
```

The runner tests the code, executes all default cells, and saves the notebook.
Default execution reuses established gene-window scores and verified permutation
caches. PDF plots and TSV tables, each with a matching TXT explanation, are
exported to `plots and tables/`. Scientific tables and other analysis outputs
remain under `notebook/results/`.

## Project layout

| Location | Contents |
| --- | --- |
| `data/downloaded/` | Original external inputs, source metadata and API responses; each file has a provenance/contents TXT companion. |
| `data/generated/` | Lifted OCD SNP caches, extracted metadata, checksum registry and file-move manifest. |
| `notebook/results/` | Generated analysis results, established scoring inputs and historical comparisons. |
| `plots and tables/` | Current notebook PDF/TSV exports and matching explanations. |
| `src/` | Analysis scripts, reusable modules, notebook runner and validation tools. |
| `tests/` | Focused scientific fixture tests. |
| `config/` | Project-local environment requirements. |
| `docs/` | Detailed analysis instructions and cleanup candidates for review. |

The UCSC chain and Agarwal posterior inputs are still read from the existing
shared `~/data/` locations; their required paths and contents are documented in
[data/README.md](./data/README.md). Moving this project to another machine requires
providing these shared inputs as well.

## Guides

1. [Data organization and provenance](./data/README.md)
2. [Every source script and its notebook inclusion](./src/README.md)
3. [Detailed analysis and reproduction guide](./docs/analysis-guide.md)
4. [Files proposed for removal, pending your decision](./docs/cleanup-candidates.md)

Run the fixture tests independently with:

```bash
.venv/bin/python -m unittest discover -s tests -v
```

The public repository is [NJC12/dog_ocd](https://github.com/NJC12/dog_ocd).
The project environment, Python/Jupyter caches, and bulk GWAS/GENCODE downloads
are excluded from Git. Download those inputs from the sources documented in
their tracked provenance TXT companions under `data/downloaded/`; the local
copies remain available in this workspace. Generated results and notebook
exports are included.

The All by All burden CSV and derived burden outputs are retained locally and
excluded from the public repository. The underlying associations can be viewed
at https://allbyall.researchallofus.org/ by searching for phenotype code
`obsessive_compulsive_disorder` (ancestry `meta`, pLoF, MAF < 0.001). See the
[input publication notes](./docs/publication-review.md) for source restrictions
and the supplementary workbook's verified bioRxiv source.
