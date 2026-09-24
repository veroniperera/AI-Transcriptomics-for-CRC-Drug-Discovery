# AI-guided identification of novel colorectal cancer inhibitors

Transcriptomics-based biomarker identification (GSE156451, TCGA-COAD, GDSC2, DepMap) followed by ML-based virtual screening for SHMT2 inhibitors.

## Setup

```bash
conda env create -f environment.yml
conda activate crc-pipeline
Rscript R/install_packages.R
```

GNINA is not on conda-forge; install the release binary from https://github.com/gnina/gnina/releases.

R 4.5, Python 3.12. Docking needs `smina`, `gnina`, `vina` (1.2.7) and `obabel` on PATH. Seed 42 throughout. Run every script from the repository root.

## Pipeline

| Step | Script | Output |
|---|---|---|
| 1 | `R/01_preprocessing_dea.R` | QC figures, limma DEGs (`results/tables/DEG_results_*.csv`) |
| 2 | `R/02_wgcna.R` | modules, module-trait table, hub genes (`hub_genes_ranked.csv`) |
| 3 | `R/03_tf_inference.R` | DoRothEA/VIPER TF activity, TF-hub overlap, Fisher test |
| 4 | `R/04_tcga_validation.R` | external validation (`TableS1_external_validation.csv`) |
| 5 | `R/05_gdsc2_pharmacogx.R` | gene-drug associations, Fig. 4B-F |
| 6 | `R/06_depmap_dependency.R` | DepMap dependency counts |
| 7 | `python/01_curate_chembl.py` | `data/processed/shmt2_curated.csv` |
| 8 | `python/02_features_split.py` | 54 descriptors, fingerprints, scaffold-disjoint 80/20 split |
| 9 | `python/03_benchmark.py` | 10 models, Optuna (6 trials), Tables 2-3, Fig. 9 EF table |
| 10 | `python/04_yrandomization_ad.py` | Y-randomization (Table 4), applicability domain, final model |
| 11 | `python/05_screening_library.py` | Sorafenib similarity library, ranking, top-50 docking set |
| 12 | `python/06_docking.py` | Smina + GNINA docking, Table 5 |
| 13 | `python/07_dude_enrichment.py dock` / `analyze` | DUD-E ROC-AUC and EF with bootstrap CIs |

```bash
./run_all.sh
python python/07_dude_enrichment.py dock --actives actives_final.ism --decoys decoys_final.ism
```

`run_all.sh` ends with the DUD-E analysis on the bundled `results/dude/dude_docking_SMALLBOX.csv`. Re-run the `dock` command first to regenerate that file. The DepMap step needs `CRISPRGeneEffect.csv` in `data/raw/`.

## Values to check after a run

| Quantity | Manuscript |
|---|---|
| Genes after filtering | 16,047 |
| DEGs (up / down) | 449 (252 / 197) |
| Soft-threshold R2 at beta = 15 / 16 | 0.807 / 0.817 |
| Modules | 10 |
| Hub candidates (turquoise / blue) | 1,005 / 592 |
| Curated ChEMBL compounds (active / inactive) | 235 (128 / 107) |
| Scaffold split (train / test) | 188 / 47 |
| DUD-E compounds docked (actives) | 507 (9) |
| DUD-E ROC-AUC (LE) | 0.707 |

## Configuration

Analysis settings are defined in `R/00_setup.R` (`QUANTILE_NORM`, `EBAYES_TREND`, `WGCNA_TOP_MAD`, `SOFT_POWER`, thresholds) and `python/common.py` (box, pIC50 cut-off, seed).

## Notes

- DUD-E ranking uses ligand efficiency, LE = -Vina affinity / heavy atoms. Compounds without a Vina score are excluded.
- Descriptor list (`python/features.py`) follows the 44 2D + 10 3D scheme in the manuscript.
