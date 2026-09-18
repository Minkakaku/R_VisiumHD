# Repository instructions

- This is the R/Seurat Visium HD main workflow, defaulting to Space Ranger 8 µm bins.
- Put reusable functions in `src/R/`; shared defaults in root `config/`; project samples, thresholds, paths, and annotations in `projects/<project>/config/`.
- Keep numbered project workflow scripts short and append-only. Do not add notebooks or a `run_all` across QC and annotation checkpoints.
- Save RDS intermediates under project `data/processed/`, references under `data/reference/`, and final outputs under numbered project `results/` directories.
- Save each formal figure as matching PNG/PDF and keep its supporting table and `README_results_CN.md` in the same analysis directory.
- Keep Xenium and historical scripts in `source_catalog/` or `legacy/original/`; never source them from the formal workflow and preserve original bytes.
- Marker analysis uses projected full bins and the full `Spatial.008um` assay. ProjectData must use PCA with feature loadings; Harmony may guide sketch clustering but is not the projection reduction.
- Do not treat Visium HD bins as segmented cells, BANKSY domains as cell types, Moran's I as group DE, or bins as biological replicates.
- Keep concise Chinese comments above formal functions and important algorithm blocks.
