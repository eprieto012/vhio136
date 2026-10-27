# vhio136 — single-cell RNA + TCR analysis

Analysis code for patient vhio136 (RIO136): 10x Genomics scRNA-seq + scTCR-seq from two sequencing batches (`rootpath`, `vhio`), ambient-RNA correction, QC, MAGIC imputation, ProjecTILs/STACAS annotation and integration, and signature scoring.

The computational environment is fully pinned:

| Layer | Pinned by |
|---|---|
| R 4.4.0 + system libraries | `renv_vhio136/env_lock/conda_renv_vhio136.explicit.txt` (conda, linux-64) |
| 347 R packages (CRAN, Bioconductor 3.20, GitHub commits) | `renv_vhio136/renv.lock` (renv 1.3.1) |
| Python 3.11 for MAGIC (`magic-impute` 3.0.0) | `renv_vhio136/env_lock/conda_magic_env.explicit.txt` + `magic_env.pip.txt` |

## Repository layout

```
scripts/        interactive-R.sh, run_rmd.sh (SLURM launcher)
  single-cell/  R Markdown analyses (m01 → m05)
  tcr-library/  TCR library analyses
  bulk/         bulk analyses
setup/          setup.R (paths, seeds, threads, Python), libraries.R, colors.R, functions.R
setup/refs/     ProjecTILs reference maps (downloaded automatically, not tracked)
renv_vhio136/   renv project: renv.lock, activate.R, settings.json, env_lock/
data/           input data (not tracked, see below)
results/        outputs (not tracked)
```

## 1. Recreate the environment

Requirements: Linux x86_64, [miniforge/conda](https://github.com/conda-forge/miniforge), git.

```bash
git clone git@github.com:eprieto012/vhio136.git
cd vhio136

# R 4.4.0 and system libraries
conda create -p ./envs/renv_vhio136 --file renv_vhio136/env_lock/conda_renv_vhio136.explicit.txt

# Python environment for MAGIC
conda create -p ./envs/magic_env --file renv_vhio136/env_lock/conda_magic_env.explicit.txt
./envs/magic_env/bin/pip install -r renv_vhio136/env_lock/magic_env.pip.txt

# Tell the launcher and setup.R where the environments are
export ENV_REPRO="$PWD/envs/renv_vhio136"
export MAGIC_PYTHON="$PWD/envs/magic_env/bin/python"

# Start R inside the renv project and install the exact package versions
./scripts/interactive-R.sh
```

```r
renv::restore()   # installs the 347 packages recorded in renv.lock
renv::status()    # should report: "No issues found"
```

`interactive-R.sh` puts the conda R first on `PATH`, unsets `R_LIBS*` and starts R from `renv_vhio136/`, whose `.Rprofile` activates renv, so only the renv library is used. Always start R through it.

### Starting R

**Option A — VS Code (interactive, with plot viewer).** Requires the [R extension](https://marketplace.visualstudio.com/items?itemName=REditorSupport.r) (REditorSupport). Point the extension to the launcher in `.vscode/settings.json`:

```json
{
  "r.rterm.linux": "/absolute/path/to/vhio136/scripts/interactive-R.sh",
  "r.sessionWatcher": true,
  "r.plot.backend": "standard"
}
```

Then `Ctrl+Shift+P` → **R: Create R Terminal**. The first time, accept the notification asking to install the helper package `sess` (installed in the extension's own folder, not in the renv library). Run `.Rmd` chunks with `Ctrl+Enter`; plots appear in the R Plot panel.

**Option B — terminal (manual, or long batch runs).**

```bash
./scripts/interactive-R.sh            # interactive R session
tmux new -s vhio136              # recommended for long renders: survives disconnections
```

A session started this way is not connected to the VS Code plot viewer. To connect it, run `Ctrl+Shift+P` → **R: Attach External R Session (Copy command)** and paste the copied command into the R console.

Check that the right environment is active in either case:

```r
renv::project()   # .../vhio136/renv_vhio136
.libPaths()[1]    # .../renv_vhio136/renv/library/...
```

The `.yml` files in `env_lock/` are human-readable equivalents of the explicit specs, useful on non-Linux systems (versions approximate).

## 2. Input data

Raw FASTQ files will be deposited at EGA (accession: **TBD**). The pipeline starts from Cell Ranger outputs, organised as:

```
data/
├── rootpath/<sample>/          # 8 libraries: pre_*/post2_d42_* (bulk1-3, pd1)
│   ├── sample_filtered_feature_bc_matrix/   (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
│   ├── raw_feature_bc_matrix/
│   ├── filtered_contig_annotations.csv
│   └── clonotypes.csv
├── vhio/<sample>/              # 4 libraries: post1_d14_*, post6_d140_* (bulk, pd1); same structure
├── outs/                       # auxiliary tables
│   ├── pre_post2_library_clonotypes.rds
│   ├── tcr_lib.txt              # TCR library (TCRb, TCR.ID)
│   ├── tcr_lib_annotated.txt
│   ├── p136-Tfrozen_TCRB in frame.tsv
│   ├── pbl_raw_RIO136.tsv       # bulk TCRβ repertoire of PBLs (Fig. 3f)
│   ├── color_TCR.IDs.txt
│   ├── signatures.tsv
│   ├── yossef_reactivity.tsv    # reactive TCRs reported by Yossef et al.
│   ├── heatmap_genes.csv        # gene categories for the iTAG heatmap (m04)
│   ├── vdjdb/vdjdb_rio136_{tra,trb,paired}.tsv
│   └── iedb/receptor_table_export_1713292478.xlsx
└── yossef-cohort/               # Cell Ranger outputs of Yossef et al. (dbGaP), one folder per sample (only these folders)
```

ProjecTILs human reference maps (`CD4T_human_ref_v1.rds`, `CD8T_human_ref_v1.rds`) are downloaded from figshare into `setup/refs/` the first time `m02` runs.

## 3. Run the analyses

Run the R Markdown files in `scripts/single-cell/` in order, either on SLURM from the project root (`sbatch -J <name> scripts/run_rmd.sh scripts/single-cell/<file>.Rmd`) or from R started with `interactive-R.sh`:

| Step | Script | Main outputs |
|---|---|---|
| 1a | `m01_data_load_rootpath_latest.Rmd` | `rootpath/seurat/` (`soupX_rootpath.rds`, `SCT_rootpath.rds`), `rootpath/plots/` |
| 1b | `m01_data_load_vhio_latest.Rmd` | `vhio/seurat/` (`soupX_vhio.rds`, `SCT_vhio.rds`), `vhio/plots/` |
| 1c | `m01_data_load_yossef_latest.Rmd` | `results/yossef-cohort/seurat/` (`soupX_yossef.rds`, `int_filt_annot.rds`), `results/yossef-cohort/plots/` (QC) |
| 2 | `m02_data_integration_latest.Rmd` | `integrated/seurat/` (`cd8_int*.rds`, `cd4_int.rds`) |
| 3 | `m03_plots_analyses_latest.Rmd` | `integrated/plots/manuscript/` (Fig. 3, S7–S9), `integrated/traj/` |
| 4 | `m04_signatures_latest.Rmd` | `integrated/markers/` (iTAG signatures and plots), updated `integrated/seurat/cd8_int_nn.rds` |
| 5 | `m05_signatures_yossef_latest.Rmd` | `results/yossef-cohort/seurat/int_filt_annot_signatures.rds`, `results/yossef-cohort/plots/` (iTAG performance in Yossef et al.) |

Steps 1a–1c are independent; 5 needs 1c and 4. Outputs of the Yossef et al. cohort go to `results/yossef-cohort/`, all others under `results/vhio136/`:

```
results/vhio136/
├── rootpath/          # batch 1 (m01_data_load_rootpath)
│   ├── seurat/        # SoupX/SCT objects and QC tables
│   └── plots/         # QC plots
├── vhio/              # batch 2 (m01_data_load_vhio)
│   ├── seurat/
│   └── plots/
└── integrated/        # both batches (m02-m04)
    ├── seurat/        # integrated CD4/CD8 objects, Yossef et al. object with scores
    ├── plots/         # manuscript figures (m03)
    ├── traj/          # clonotype phenotype trajectories (tables + png)
    └── markers/       # iTAG signature gene sets and plots (m04)
```

```r
rmarkdown::render("../scripts/single-cell/m01_data_load_rootpath_latest.Rmd", knit_root_dir = getwd())
```

## Reproducibility notes

- Seeds are fixed in `setup/setup.R` (`SEED = 123`, `RNGkind("L'Ecuyer-CMRG")`) and in Python (NumPy, `random`).
- BLAS/OpenMP threads are limited to 1. Results should be identical on the same platform; small numerical differences (e.g. UMAP coordinates) are possible on different CPUs.
- `renv.lock` records `future` 1.70.0 with `future.apply` 1.21.0, which declares `future >= 1.80.0`. This is the combination the analyses were run with and it loads correctly; it is kept on purpose.

