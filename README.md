# Longitudinal multiomic profiling and functional TCR screening in patient VHIO-136

This repository contains the code used to find the tumor-reactive T cells of patient VHIO-136 and to follow them over the course of treatment. It combines four types of data: bulk TCRβ sequencing of blood samples taken over time and of a tumor biopsy taken before treatment, flow cytometry of blood T cells, a functional screen of 1,001 candidate T-cell receptors (TCRs), and single-cell RNA and TCR sequencing.

## Workflow

The analysis has four steps:

1. **Tracking the TCR repertoire over time** (section 4). We follow each T-cell clonotype (the T cells that share one TCR) in blood at 12 timepoints and in the baseline tumor biopsy, and measure how T-cell phenotypes change by flow cytometry. From these data and the single-cell TCRs of two early timepoints (pre-treatment and post2, day 42), we chose 1,001 TCRs to test: clonotypes that were abundant, expanding over time, shared with the tumor, or found in the PD-1<sup>hi</sup> T cells. The reason each TCR was chosen, and the result of its test, are listed in `resources/tcr-library/TCR-lib-VHIO136-annotated.txt`.
2. **Screening the TCR library** (section 5). The 1,001 TCRs are introduced into T cells from healthy donors, which are then co-cultured with the patient's own tumor cell line. T cells that recognize the tumor turn on the activation marker 4-1BB and are sorted. A TCR is called reactive when it is enriched in the 4-1BB+ cells compared with the 4-1BB− cells. The same screen is run against an irrelevant tumor cell line (specificity control) and against B cells that present the targets of the positive-control TCRs included in the library (positive control). Pipeline: sequencing reads → counts per TCR → differential enrichment with edgeR → volcano plots.
3. **Single-cell analysis** (section 6). We process and integrate the single-cell RNA and TCR data of four timepoints: pre-treatment and post2, sequenced by RootPath (`rootpath`), and post1 (day 14) and post6 (day 140), sequenced in-house (`vhio`). Matching each cell's TCR to the screening results tells us which cells are tumor-reactive, so we can follow them over time and describe their gene-expression states. We then derive gene signatures that change with time and test whether they also identify tumor-reactive T cells in an independent cohort (Yossef et al.).
4. **Combining bulk and single-cell data** (section 7). We draw the blood frequency of each tumor-reactive clonotype over time, colored by its single-cell phenotype, to show which clonotypes expand, persist or disappear during treatment and how their state changes.

```
1. long-tracking   script-bulk-tcrb ──────────────────┐
                   script-flow-phenotype              ├─► 1,001 candidate TCRs
                   Cell Ranger rootpath (pre, post2) ─┘

2. tcr-library     run-hits.sh ─► script-hit-identification-TCL ─┬─► script-volcano-plots-TCL
                                                                 └─► reactive TCRs ──────┐
                                                                                         │
3. single-cell     m01-rootpath ─┐                                                       │
                   m01-vhio ─────┴─► m02 ◄───────────────────────────────────────────────┘
                                     └─► m03 ─┬─► m04 ─► m05 ◄── m01-yossef
                                              │
4. integration                                └─► script-bulk-sc-phenotype
```

The computational environment is fully pinned:

| Layer                                                                                                              | Pinned by                                                                                                                     |
| ------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------- |
| R 4.4.0 + system libraries                                                                                         | `setup/env-lock/conda_renv_vhio136.explicit.txt` (conda, linux-64)                                                          |
| 348 R packages (CRAN, Bioconductor 3.20, GitHub commits)                                                           | `renv-vhio136/renv.lock` (renv 1.3.1)                                                                                       |
| Python 3.11 for MAGIC (`magic-impute` 3.0.0)                                                                     | `setup/env-lock/conda_magic_env.explicit.txt` + `magic_env.pip.txt`                                                       |
| TCR library screen alignment: fastp 0.23.2, cutadapt 4.0, seqkit 2.14.0, bowtie2 2.5.5, samtools 1.25 (Python 3.9) | `setup/env-lock/conda_tcrlib_hits.explicit.txt` (conda, linux-64)                                                           |
| Cell Ranger 7.1.0 (vhio136) and 9.0.0 (Yossef et al.), 10x references                                              | `scripts/single-cell/cellranger/` (versions and references in each script; not installed by the steps below, see section 2) |

## Repository layout

```
vhio136/
├── scripts/                               # interactive-R.sh, run-rmd.sh (SLURM launcher, section 3)
│   ├── long-tracking/                     # longitudinal PBL follow-up (sections 4 and 7)
│   │   ├── script-bulk-tcrb.Rmd           # bulk TCRβ repertoire
│   │   ├── script-flow-phenotype.Rmd      # flow-cytometry phenotyping
│   │   └── script-bulk-sc-phenotype.Rmd   # bulk trajectories colored by single-cell phenotype
│   ├── tcr-library/                       # TCR library screen (section 5)
│   │   ├── hits/                          # FASTQ → counts per library construct (align_sample.py, run-hits.sh)
│   │   └── rmd/                           # counts → edgeR hit calling and volcano plots
│   └── single-cell/                       # single-cell analysis (section 6)
│       ├── cellranger/                    # Cell Ranger multi: vhio136/ (rootpath + vhio), yossef-cohort/ (FASTQ → data/)
│       └── rmd/                           # R Markdown analyses (m01 → m05)
├── setup/                                 # setup.R (paths, seeds, threads, Python), libraries.R, colors.R, functions*.R
│   ├── env-lock/                          # conda/pip locks: R base, MAGIC (Python) and TCR library alignment tools
│   └── refs/                              # ProjecTILs reference maps (not tracked, see section 2)
├── renv-vhio136/                          # renv project: renv.lock, activate.R, settings.json
├── resources/                             # reference tables and flow-cytometry exports (tracked, see section 2)
├── data/                                  # sequencing-derived input data (not tracked, see section 2)
├── results/                               # outputs (not tracked; created by run-rmd.sh)
└── logs/                                  # SLURM logs (not tracked; created by run-rmd.sh)
```

## 1. Recreate the environment

Requirements: Linux x86_64 with SLURM, [miniforge/conda](https://github.com/conda-forge/miniforge) and git. The single-cell steps need up to ~130 GB of RAM (`m02`); the other steps run with less than 100 GB.

```bash
git clone git@github.com:eprieto012/vhio136.git
cd vhio136

# R 4.4.0 and system libraries
conda create -p ./envs/renv_vhio136 --file setup/env-lock/conda_renv_vhio136.explicit.txt

# Python environment for MAGIC
conda create -p ./envs/magic_env --file setup/env-lock/conda_magic_env.explicit.txt
./envs/magic_env/bin/pip install -r setup/env-lock/magic_env.pip.txt

# Alignment tools for the TCR library screen (fastp, cutadapt, seqkit, bowtie2, samtools)
conda create -p ./envs/tcrlib_hits --file setup/env-lock/conda_tcrlib_hits.explicit.txt

# Tell the launchers and setup.R where the environments are
export ENV_REPRO="$PWD/envs/renv_vhio136"
export MAGIC_PYTHON="$PWD/envs/magic_env/bin/python"
export ENV_HITS="$PWD/envs/tcrlib_hits"

# Start R inside the renv project and install the exact package versions
./scripts/interactive-R.sh
```

```r
renv::restore()   # installs the 348 packages recorded in renv.lock
renv::status()    # should report: "No issues found"
```

`interactive-R.sh` puts the conda R first on `PATH`, unsets `R_LIBS*` and starts R from `renv-vhio136/`, whose `.Rprofile` activates renv, so only the renv library is used. Always start R through it. The `.yml` files in `setup/env-lock/` are human-readable equivalents of the explicit specs, useful on non-Linux systems (versions approximate).

### Starting R

Use VS Code with the [R extension](https://marketplace.visualstudio.com/items?itemName=REditorSupport.r) (REditorSupport). Point the extension to the project launcher in `.vscode/settings.json`:

```json
{
  "r.rterm.linux": "/absolute/path/to/vhio136/scripts/interactive-R.sh",
  "r.sessionWatcher": true,
  "r.plot.backend": "standard"
}
```

Then `Ctrl+Shift+P` → **R: Create R Terminal**. The first time, accept the notification asking to install the helper package `sess` (installed in the extension's own folder, not in the renv library). Run `.Rmd` chunks with `Ctrl+Enter`; plots appear in the R Plot panel. Full runs of a script go through SLURM instead (section 3).

To check that the project environment is active:

```r
renv::project()   # .../vhio136/renv-vhio136
.libPaths()[1]    # .../renv-vhio136/renv/library/...
```


## 2. Input data

| Data                                                    | Where                         | Access     |
| ------------------------------------------------------- | ----------------------------- | ---------- |
| Single-cell RNA/TCR sequencing (FASTQ), VHIO-136        | EGA, accession **TBD**        | controlled |
| TCR library screen (FASTQ)                              | EGA, accession **TBD**        | controlled |
| Bulk TCRβ repertoires (Adaptive), blood and tumor       | EGA, accession **TBD**        | controlled |
| Flow-cytometry exports                                  | `resources/flow-phenotype/`   | this repository |
| Single-cell data of Yossef et al.                       | dbGaP, accession **TBD**      | controlled |

The single-cell analysis starts from Cell Ranger outputs (below); the TCR library screen starts from its FASTQ files (section 5).

### From FASTQ to Cell Ranger outputs

`data/single-cell/` was generated with `cellranger multi`. The configuration templates and the run script are in `scripts/single-cell/cellranger/`:

| Output folder                       | Script folder      | Cell Ranger | GEX / VDJ reference                       |
| ----------------------------------- | ------------------ | ----------- | ----------------------------------------- |
| `data/single-cell/rootpath/`      | `vhio136/`       | 7.1.0       | GRCh38-2020-A / GRCh38-alts-ensembl-7.0.0 |
| `data/single-cell/vhio/`          | `vhio136/`       | 7.1.0       | GRCh38-2020-A / GRCh38-alts-ensembl-7.0.0 |
| `data/single-cell/yossef-cohort/` | `yossef-cohort/` | 9.0.0       | GRCh38-2024-A / GRCh38-alts-ensembl-7.1.0 |

Cell Ranger is run once per library with these settings. The outputs used in the analysis (the filtered and raw feature-barcode matrices, `filtered_contig_annotations.csv` and `clonotypes.csv`) are placed in `data/single-cell/<folder>/<library>/`, where `m01` reads them. How to run it is described in the header of `run-cellranger-multi.sh`.

Cell Ranger and its references are not included. The references are available from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger/downloads); Cell Ranger was run from the public Docker image `litd/docker-cellranger`, at these exact versions:

```bash
singularity build cellranger-7.1.0.sif docker://litd/docker-cellranger@sha256:5b27ba5674db54cb9820d912cef3ba841dc582f137fba7c8353a91881ae1f469  # vhio136
singularity build cellranger-9.0.0.sif docker://litd/docker-cellranger@sha256:91fe8317edca1ae48f73738e521b5baf872d2b3be0fd1b4db2707f50ebac22b2  # Yossef et al.
```

### Layout of data/

Sequencing-derived input data, not tracked in git:

```
data/
├── bulk-tcrb/                  # bulk TCRβ repertoires (Adaptive)
│   ├── blood/d<day>.tsv         # PBLs, one file per timepoint (script-bulk-tcrb, m03 Fig. 3f; every .tsv here is read as a timepoint)
│   └── tumor/p136-Tfrozen_TCRB_in_frame.tsv   # tumor biopsy (script-bulk-tcrb, m02)
└── single-cell/                # Cell Ranger outputs
    ├── rootpath/<sample>/      # 8 libraries: pre_*/post2_d42_* (bulk1-3, pd1)
    │   ├── sample_filtered_feature_bc_matrix/   (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
    │   ├── raw_feature_bc_matrix/
    │   ├── filtered_contig_annotations.csv
    │   └── clonotypes.csv
    ├── vhio/<sample>/          # 4 libraries: post1_d14_*, post6_d140_* (bulk, pd1); same structure
    └── yossef-cohort/<sample>/ # Yossef et al. (dbGaP), 12 libraries; same structure
```

### Layout of resources/

Small tables read by the scripts, tracked in git (path set as `resources_dir` in `setup/setup.R`):

```
resources/
├── flow-phenotype/                                # FlowJo exports of the PBL phenotyping (script-flow-phenotype, section 4.2)
│   ├── VHIO136-PD1-cat.csv                        # PD-1 hi/int/neg counts per timepoint
│   ├── VHIO136-CD4-CD8-panel{1,2}.csv             # marker frequencies per PD-1 level
│   └── VHIO136-reanalysis-{CD4,CD8}pan{1,2}.csv   # FlowJo gate exports per marker-positive population
├── tcr-library/
│   ├── TCR-lib-VHIO136-library-TCRb-clono.txt     # TCR library (TCRb, TCR.ID) (script-bulk-tcrb, m03)
│   ├── TCR-lib-VHIO136-annotated.txt              # library with selection criteria and reactivity: written by script-bulk-tcrb, read by m02
│   ├── TCR-lib-VHIO136-clonotypes-pre-post2.rds   # library clonotypes, used to keep library cells at doublet removal (m01)
│   ├── TCR-lib-VHIO136-AAV.fasta                  # library sequences (3011 constructs); bowtie2 index built next to it in bowtie2_indices/ (not tracked)
│   └── TCR-lib-VHIO136-screen-samples.tsv         # the 26 screen samples: FASTQ ↔ analysis ID (section 5)
├── viral-databases/
│   ├── vdjdb/vdjdb-vhio136-{tra,trb,paired}.tsv   # (script-bulk-tcrb, m02)
│   └── iedb/iedb-receptor-export-2024-04-16.xlsx  # (m02)
├── signatures/
│   └── published-signatures.tsv                   # Oliveira, Lowery, Yossef, Hanada, TRTPBL gene sets (m02, m05)
└── yossef-cohort/
    └── yossef-reactive-tcrs.tsv                   # reactive TCRs reported by Yossef et al. (only reactive TCRs; patients with data in single-cell/yossef-cohort/) (m01-yossef)
```

ProjecTILs human reference maps (`CD4T_human_ref_v1.rds`, `CD8T_human_ref_v1.rds`) must be placed in `setup/refs/` before running `m01-data-load-yossef` (CD8T) and `m02` (CD4T, CD8T). They are not downloaded automatically: `m02` stops if they are missing. Get them from the [ProjecTILs reference atlases](https://github.com/carmonalab/ProjecTILs#reference-atlases).

## 3. Running the R Markdown scripts

Every `.Rmd` (sections 4–7) is submitted to SLURM with `run-rmd.sh`, from the project root:

```bash
cd /path/to/vhio136
scripts/run-rmd.sh -J m01-rootpath scripts/single-cell/rmd/m01-data-load-rootpath.Rmd
scripts/run-rmd.sh -J m02-integration --mem=200gb scripts/single-cell/rmd/m02-data-integration.Rmd
```

Options before the `.Rmd` go to `sbatch`; the default is 25 CPUs and 100 GB, except `m02`, which needs `--mem=200gb`. Outputs go to `results/` and the log of each job to `logs/<name>_<jobid>.log`.

### Full run, in order

Each job is submitted on its own, so wait for a step to finish before starting the next one. Times are those of our last full run (25 CPUs).

```bash
# 1. Longitudinal follow-up (section 4)                                    ~3 min each
scripts/run-rmd.sh -J long_tracking scripts/long-tracking/script-bulk-tcrb.Rmd
scripts/run-rmd.sh -J flow_pheno    scripts/long-tracking/script-flow-phenotype.Rmd

# 2. TCR library screen (section 5)                                        ~10 min per sample (30 CPUs)
export RAW_DIR=/path/to/fastq
scripts/tcr-library/hits/run-hits.sh screen-1
scripts/tcr-library/hits/run-hits.sh screen-2
scripts/run-rmd.sh -J tcl_hits    scripts/tcr-library/rmd/script-hit-identification-TCL.Rmd   # < 1 min
scripts/run-rmd.sh -J tcl_volcano scripts/tcr-library/rmd/script-volcano-plots-TCL.Rmd        # < 1 min

# 3. Single-cell analysis (section 6); the three m01 can run in parallel
scripts/run-rmd.sh -J m01-rootpath scripts/single-cell/rmd/m01-data-load-rootpath.Rmd          # ~1.7 h
scripts/run-rmd.sh -J m01-vhio     scripts/single-cell/rmd/m01-data-load-vhio.Rmd              # ~1.4 h
scripts/run-rmd.sh -J m01-yossef   scripts/single-cell/rmd/m01-data-load-yossef.Rmd            # ~3.2 h
scripts/run-rmd.sh -J m02 --mem=200gb scripts/single-cell/rmd/m02-data-integration.Rmd         # ~1.6 h
scripts/run-rmd.sh -J m03 scripts/single-cell/rmd/m03-plots-analyses.Rmd                       # ~3 min
scripts/run-rmd.sh -J m04 scripts/single-cell/rmd/m04-signatures.Rmd                           # ~6 min
scripts/run-rmd.sh -J m05 scripts/single-cell/rmd/m05-signatures-yossef.Rmd                    # ~3 min

# 4. Bulk × single-cell integration (section 7)                            ~1 min
scripts/run-rmd.sh -J bulk_sc_pheno scripts/long-tracking/script-bulk-sc-phenotype.Rmd
```

## 4. Longitudinal follow-up of PBLs

### 4.1 Bulk TCRβ repertoire

`script-bulk-tcrb.Rmd` follows the T-cell clonotypes in the patient's blood over the 12 sampled timepoints: which ones are shared with the tumor, which expand or appear, which have a known viral specificity (VDJdb), and how similar the repertoire is between timepoints. It also annotates each library TCR with why it was selected and the result of its test.

```bash
scripts/run-rmd.sh -J long_tracking scripts/long-tracking/script-bulk-tcrb.Rmd
```

- **Inputs:** `data/bulk-tcrb/`, `data/single-cell/rootpath/` and `resources/`.
- **Outputs:** figures in `results/long-tracking/repertoire/` and `resources/tcr-library/TCR-lib-VHIO136-annotated.txt`, read by `m02`. The reactive TCR IDs are written by hand in this script; if they change, re-run it and then `m02` → `m05`.

### 4.2 Flow-cytometry phenotyping

`script-flow-phenotype.Rmd` follows PD-1 levels and activation and differentiation markers of blood CD4+ and CD8+ T cells over 7 timepoints.

```bash
scripts/run-rmd.sh -J flow_pheno scripts/long-tracking/script-flow-phenotype.Rmd
```

- **Inputs:** the FlowJo exports in `resources/flow-phenotype/`.
- **Outputs:** 10 PDFs in `results/long-tracking/phenotype/`.

## 5. TCR library screen

The 1,001 TCRs were expressed in T cells of two healthy donors, co-cultured with target cells and sorted by the activation marker 4-1BB. A TCR is reactive when it is enriched in the 4-1BB+ cells. Each donor was screened twice against three targets:

- **Relevant**: the patient's own tumor cell line (VHIO-136 TCL).
- **Irrelevant**: an unrelated tumor cell line, as specificity control (NCI3998 in screen 1, VHIO-008 TCL in screen 2).
- **B cells** of patient VHIO-008, as positive control for the control TCRs included in the library.

The 26 samples are listed in `resources/tcr-library/TCR-lib-VHIO136-screen-samples.tsv`. **FASTQ names and analysis IDs differ** (e.g. `Screen1_TCR_1` is `TCR_8` in the analysis): always match samples through this table.

### 5.1 Alignment: FASTQ → counts

`run-hits.sh` counts the reads of each library TCR in every sample (one SLURM job per sample), using the `tcrlib_hits` environment (section 1).

```bash
export RAW_DIR=/path/to/fastq ENV_HITS=/path/to/envs/tcrlib_hits
scripts/tcr-library/hits/run-hits.sh screen-1
scripts/tcr-library/hits/run-hits.sh screen-2
```

- **Inputs:** the FASTQs in `<RAW_DIR>/<screen>/<sample>_{1,2}.fastq.gz` (names from the sample table) and the library sequences in `resources/tcr-library/`.
- **Outputs:** one `<sample>.counts.txt` per sample in `results/tcr-library/<screen>/<sample>/`.

### 5.2 Hit calling and volcano plots

`script-hit-identification-TCL.Rmd` compares 4-1BB+ and 4-1BB− cells for each target with edgeR. A TCR is a hit when at least 2 of its 3 synonymous variants have logFC > 1.5 and FDR ≤ 0.05. `script-volcano-plots-TCL.Rmd` plots the results and highlights the tumor-specific hits.

```bash
scripts/run-rmd.sh -J tcl_hits scripts/tcr-library/rmd/script-hit-identification-TCL.Rmd
scripts/run-rmd.sh -J tcl_volcano scripts/tcr-library/rmd/script-volcano-plots-TCL.Rmd
```

- **Inputs:** the counts (5.1) and the sample table.
- **Outputs:** `results/tcr-library/`, including `volcano.pdf` and the list of tumor-specific hits (`relevant_hits_specific.txt`).

The reactive TCRs (screen and extended validation) are then entered by hand in `script-bulk-tcrb.Rmd` (section 4.1).

## 6. Single-cell analysis

Run the scripts in `scripts/single-cell/rmd/` in this order, as in section 3. Outputs go to `results/single-cell/`.

| Step | Script                         | What it does                                                            | Needs                                |
| ---- | ------------------------------ | ----------------------------------------------------------------------- | ------------------------------------ |
| 1a   | `m01-data-load-rootpath.Rmd` | QC and normalization of the RootPath batch (pre, post2)                 | `data/single-cell/rootpath/`       |
| 1b   | `m01-data-load-vhio.Rmd`     | QC and normalization of the in-house batch (post1, post6)               | `data/single-cell/vhio/`           |
| 1c   | `m01-data-load-yossef.Rmd`   | QC and CD8+ annotation of the Yossef et al. cohort                      | `data/single-cell/yossef-cohort/`  |
| 2    | `m02-data-integration.Rmd`   | integrates both batches; annotates cells with TCR reactivity            | 1a, 1b, the library annotation (4.1) |
| 3    | `m03-plots-analyses.Rmd`     | manuscript figures (Fig. 3, S7–S9) and phenotype tables per clonotype  | 2                                    |
| 4    | `m04-signatures.Rmd`         | derives the iTAG gene signatures                                        | 2                                    |
| 5    | `m05-signatures-yossef.Rmd`  | tests the iTAG signatures in the Yossef et al. cohort (AUC per patient) | 1c, 4                                |

If you re-run `m02`, re-run `m04` afterwards (it adds the signature scores to the object `m02` writes). The CD8+ cluster names are set by hand in `m02`: if the integration changes, check them before running steps 3–5.

## 7. Bulk × single-cell integration

`script-bulk-sc-phenotype.Rmd` draws the blood trajectories of the clonotypes colored by their single-cell phenotype. Run it last: it needs the tables written by `m03`.

```bash
scripts/run-rmd.sh -J bulk_sc_pheno scripts/long-tracking/script-bulk-sc-phenotype.Rmd
```

- **Inputs:** `data/bulk-tcrb/blood/` and the `m03` tables in `results/single-cell/integrated/phenotype-trajectories/tables/`.
- **Outputs:** `results/long-tracking/phenotype-trajectories/` (`plots/`, `tables/`).

## Results

```
results/
├── long-tracking/
│   ├── repertoire/                  # bulk TCRβ figures (4.1)
│   ├── phenotype/                   # flow-cytometry figures (4.2)
│   └── phenotype-trajectories/      # bulk trajectories colored by single-cell phenotype (7)
├── tcr-library/                     # screen counts, edgeR results and volcano plots (5)
└── single-cell/                     # (6)
    ├── rootpath/, vhio/, yossef-cohort/   # Seurat objects and QC per batch (m01)
    ├── integrated/                  # integrated objects, QC, manuscript figures, phenotype tables (m02, m03)
    └── DEG-signatures/              # iTAG signatures (m04) and their performance in Yossef et al. (m05)
```

## Reproducibility notes

- TCR library counts: re-running the alignment reproduces the counts except for 1–3 reads in a few constructs per sample, because fastp's UMI deduplication is not deterministic (tested on `Screen1_TCR_1`: 30–40 of 2,790 constructs differ by 1–3 reads).
- Seeds are fixed in `setup/setup.R` (`SEED = 123`, `RNGkind("L'Ecuyer-CMRG")`) and in Python (NumPy, `random`).
- BLAS/OpenMP threads are limited to 1. Results should be identical on the same platform; small numerical differences (e.g. UMAP coordinates) are possible on different CPUs.
- `renv.lock` records `future` 1.70.0 with `future.apply` 1.21.0, which declares `future >= 1.80.0`. This is the combination the analyses were run with and it loads correctly; it is kept on purpose.

## Citation

If you use this code, please cite: **TBD** (manuscript in preparation).

## Contact

Questions about the code: Endika Prieto (endikaprieto@vhio.net). Data access is requested through EGA and dbGaP (see section 2).

## License

The code is released under the [MIT License](LICENSE). Data are not covered by it: see section 2 for their sources.
