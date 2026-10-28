# vhio136 — single-cell RNA + TCR analysis

Analysis code for patient vhio136 (RIO136): 10x Genomics scRNA-seq + scTCR-seq from two sequencing batches (`rootpath`, `vhio`), ambient-RNA estimation (QC only), QC, MAGIC imputation, ProjecTILs/STACAS annotation and integration, and signature scoring, with validation in the Yossef et al. cohort.

The computational environment is fully pinned:

| Layer | Pinned by |
|---|---|
| R 4.4.0 + system libraries | `renv-vhio136/env-lock/conda_renv_vhio136.explicit.txt` (conda, linux-64) |
| 347 R packages (CRAN, Bioconductor 3.20, GitHub commits) | `renv-vhio136/renv.lock` (renv 1.3.1) |
| Python 3.11 for MAGIC (`magic-impute` 3.0.0) | `renv-vhio136/env-lock/conda_magic_env.explicit.txt` + `magic_env.pip.txt` |

## Repository layout

```
scripts/        interactive-R.sh, run-rmd.sh (SLURM launcher)
  single-cell/  R Markdown analyses (m01 → m05)
  tcr-library/  TCR library analyses
setup/          setup.R (paths, seeds, threads, Python), libraries.R, colors.R, functions.R, functions-yossef.R
setup/refs/     ProjecTILs reference maps (not tracked, see section 2)
renv-vhio136/   renv project: renv.lock, activate.R, settings.json, env-lock/
data/           input data (only data/aux/ is tracked, see section 2)
results/        outputs (not tracked; created automatically by run-rmd.sh)
logs/           SLURM logs (not tracked; created automatically by run-rmd.sh)
```

## 1. Recreate the environment

Requirements: Linux x86_64, [miniforge/conda](https://github.com/conda-forge/miniforge), git.

```bash
git clone git@github.com:eprieto012/vhio136.git
cd vhio136

# R 4.4.0 and system libraries
conda create -p ./envs/renv_vhio136 --file renv-vhio136/env-lock/conda_renv_vhio136.explicit.txt

# Python environment for MAGIC
conda create -p ./envs/magic_env --file renv-vhio136/env-lock/conda_magic_env.explicit.txt
./envs/magic_env/bin/pip install -r renv-vhio136/env-lock/magic_env.pip.txt

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

`interactive-R.sh` puts the conda R first on `PATH`, unsets `R_LIBS*` and starts R from `renv-vhio136/`, whose `.Rprofile` activates renv, so only the renv library is used. Always start R through it.

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
renv::project()   # .../vhio136/renv-vhio136
.libPaths()[1]    # .../renv-vhio136/renv/library/...
```

The `.yml` files in `env-lock/` are human-readable equivalents of the explicit specs, useful on non-Linux systems (versions approximate).

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
├── bulk/                       # bulk TCRβ repertoires (not tracked)
│   ├── blood/d<day>.tsv         # PBLs, one file per timepoint
│   └── tumor/p136-Tfrozen_TCRB in frame.tsv   # tumor biopsy (m02)
├── aux/                        # auxiliary input tables (path set as aux_dir in setup/setup.R; tracked in git, except bulk-repertoire/)
│   ├── tcr-library/
│   │   ├── tcr_library.txt                       # TCR library (TCRb, TCR.ID) (m03)
│   │   ├── tcr_library_annotated.txt             # library with reactivity and selection criteria (m02)
│   │   ├── tcr_library_clonotypes_pre_post2.rds  # library clonotypes, used to keep library cells at doublet removal (m01)
│   │   └── tcr_library_colors.txt                # colour and emergence timepoint per reactive TCR.ID (m03)
│   ├── bulk-repertoire/
│   │   └── pbl_trb_timecourse_vhio136.tsv        # processed bulk TCRβ repertoire of PBLs, all timepoints (m03, Fig. 3f; not tracked)
│   ├── viral-databases/
│   │   ├── vdjdb/vdjdb_vhio136_{tra,trb,paired}.tsv
│   │   └── iedb/iedb_receptor_export_2024-04-16.xlsx  # (m02)
│   ├── signatures/
│   │   └── published_signatures.tsv              # Oliveira, Lowery, Yossef, Hanada, TRTPBL gene sets (m02, m05)
│   └── yossef-cohort/
│       └── yossef_reactive_tcrs.tsv              # reactive TCRs reported by Yossef et al. (only reactive TCRs; patients with data in yossef-cohort/) (m01_yossef)
└── yossef-cohort/               # Cell Ranger outputs of Yossef et al. (dbGaP), one folder per sample (only these folders)
```

ProjecTILs human reference maps (`CD4T_human_ref_v1.rds`, `CD8T_human_ref_v1.rds`) must be placed in `setup/refs/` before running `m01_data_load_yossef` (CD8T) and `m02` (CD4T, CD8T). They are not downloaded automatically: `m02` stops if they are missing. Get them from the [ProjecTILs reference atlases](https://github.com/carmonalab/ProjecTILs#reference-atlases).

## 3. Run the analyses

Run the R Markdown files in `scripts/single-cell/` from the project root with the SLURM launcher. Called directly (not through `sbatch`), `run-rmd.sh` creates `logs/` and `results/` if they do not exist and submits itself; any option before the `.Rmd` is passed to `sbatch`:

```bash
cd /path/to/vhio136
scripts/run-rmd.sh -J m01_rootpath scripts/single-cell/m01_data_load_rootpath_latest.Rmd
scripts/run-rmd.sh -J m02_integration --mem=200gb scripts/single-cell/m02_data_integration_latest.Rmd
```

- `logs/` must exist before `sbatch` is called (SLURM opens the log file before the job starts and does not create its directory). Calling `sbatch -J <name> scripts/run-rmd.sh <file>.Rmd` directly also works, but only if `logs/` already exists.
- Folders are only created when missing: `logs/` accumulates one `<name>_<jobid>.log` per job, and outputs in `results/` are overwritten when a script is re-run (files a script no longer produces are not removed).
- Resources default to the `#SBATCH` header of `run-rmd.sh` (25 CPUs, 100 GB). `m02` peaks at ~160 GB, so run it with `--mem=200gb`.
- Each script ends by printing `"Plots saved successfully."`: if this line is not at the end of the log, the run did not finish.

| Step | Script | Main outputs |
|---|---|---|
| 1a | `m01_data_load_rootpath_latest.Rmd` | `results/vhio136/rootpath/seurat/` (`soupX_rootpath.rds`, `SCT_rootpath.rds`), `rootpath/plots/` (QC) |
| 1b | `m01_data_load_vhio_latest.Rmd` | `results/vhio136/vhio/seurat/` (`soupX_vhio.rds`, `SCT_vhio.rds`), `vhio/plots/` (QC) |
| 1c | `m01_data_load_yossef_latest.Rmd` | `results/yossef-cohort/seurat/` (`soupX_yossef.rds`, `int_filt_annot.rds`, with scGate/ProjecTILs CD8 annotation), `results/yossef-cohort/plots/` (QC) |
| 2 | `m02_data_integration_latest.Rmd` | `results/vhio136/integrated/seurat/` (`cd8_int.rds`, `cd8_int_nn.rds`, `cd8_int_bulk*.rds`, `cd4_int.rds`), `integrated/plots/` |
| 3 | `m03_plots_analyses_latest.Rmd` | `results/manuscript/` (Fig. 3, S7–S9), `results/vhio136/integrated/phenotype-trajectories/` |
| 4 | `m04_signatures_latest.Rmd` | `results/vhio136/integrated/DEG-signatures/cd8_int_nn/` (iTAG gene sets and plots); overwrites `integrated/seurat/cd8_int_nn.rds` with the iTAG scores |
| 5 | `m05_signatures_yossef_latest.Rmd` | `results/yossef-cohort/seurat/int_filt_annot_signatures.rds`, `results/yossef-cohort/signatures/` (all cells + `*_cd8_nn.pdf`) |

Order: 1a–1c are independent and can run in parallel; 2 needs 1a and 1b; 3 and 4 need 2; 5 needs 1c and 4.

**Manual checkpoint after step 2.** The CD8+ non-naïve clusters are named by hand in `m02` (`recode` of `STACAS_snn_res.0.5` into `manual.annotation`). Whenever the integration changes, check the cluster markers printed by `m02` and update that mapping before running steps 3–5; otherwise phenotypes are mislabelled downstream without any error.

If `m02` is re-run, re-run `m04` afterwards: `m02` rewrites `cd8_int_nn.rds` without the iTAG scores that `m04` adds.

```
results/
├── manuscript/                      # manuscript figures (m03)
├── vhio136/
│   ├── rootpath/{seurat,plots}/     # batch 1 (m01_data_load_rootpath)
│   ├── vhio/{seurat,plots}/         # batch 2 (m01_data_load_vhio)
│   └── integrated/                  # both batches (m02-m04)
│       ├── seurat/                  # integrated CD4/CD8 objects
│       ├── plots/                   # integration QC (m02)
│       ├── phenotype-trajectories/  # clonotype phenotype trajectories: plots/, tables/ (m03)
│       └── DEG-signatures/cd8_int_nn/
│           ├── markers/             # iTAG gene sets (.txt, itag_signatures.rds)
│           └── plots/{markers,yield}/
└── yossef-cohort/
    ├── seurat/                      # Yossef et al. objects (m01, m05)
    ├── plots/                       # QC (m01)
    └── signatures/                  # iTAG performance in Yossef et al. (m05)
```

To run a script interactively instead, start R with `interactive-R.sh` (it starts in `renv-vhio136/`, so project paths are `../setup`, `../results`, ...) and render it:

```r
rmarkdown::render("../scripts/single-cell/m01_data_load_rootpath_latest.Rmd", knit_root_dir = getwd())
```

## Reproducibility notes

- Seeds are fixed in `setup/setup.R` (`SEED = 123`, `RNGkind("L'Ecuyer-CMRG")`) and in Python (NumPy, `random`).
- BLAS/OpenMP threads are limited to 1. Results should be identical on the same platform; small numerical differences (e.g. UMAP coordinates) are possible on different CPUs.
- `renv.lock` records `future` 1.70.0 with `future.apply` 1.21.0, which declares `future >= 1.80.0`. This is the combination the analyses were run with and it loads correctly; it is kept on purpose.

