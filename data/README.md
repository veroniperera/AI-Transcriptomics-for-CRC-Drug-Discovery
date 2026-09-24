# Data

Nothing here is versioned. Each source is downloaded by the scripts unless noted.

| Source | Accession / file | Used by |
|---|---|---|
| GEO | GSE156451 (supplementary CPM/TMM count matrix) | `R/01_preprocessing_dea.R` |
| TCGA-COAD | STAR counts via TCGAbiolinks | `R/04_tcga_validation.R` |
| GDSC2 | PharmacoGx `GDSC_2020(v2-8.2)` | `R/05_gdsc2_pharmacogx.R` |
| DepMap | `CRISPRGeneEffect.csv` (Public 29Q1, Chronos), place in `data/raw/` manually | `R/06_depmap_dependency.R` |
| ChEMBL | target `CHEMBL4295747`, IC50, binding assays | `python/01_curate_chembl.py` |
| PDB | 7BYI | `python/06_docking.py` |
| DUD-E | actives/decoys generated at https://dude.docking.org/ (`actives_final.ism`, `decoys_final.ism`) | `python/07_dude_enrichment.py` |

`results/dude/dude_docking_SMALLBOX.csv` (510 submitted, 507 docked) is the docking output used for Fig. 13.
