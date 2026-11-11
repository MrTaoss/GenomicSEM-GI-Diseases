# Genomic SEM analysis of chronic gastrointestinal diseases

Code accompanying the manuscript *Resolving the latent genetic architecture and spatial enrichment patterns of common gastrointestinal diseases*.

## Scripts

- `GenomicSEM_GI_core_analysis.R` — LDSC, parallel analysis, EFA, CFA, chromosome-split validation, and factor GWAS for the gastrointestinal disease factors.
- `gsMap_GI_spatial_analysis.py` — preparation of the Harnik et al. small-intestine Visium data, gsMap spatial genetic enrichment, Cauchy combination tests, and anatomical-region analysis.

## Software

- R 4.5.2; GenomicSEM 0.0.5; lavaan 0.6-21; data.table; Matrix.
- Python 3.11; gsMap 1.73.8; anndata; numpy; pandas; scipy; Pillow.
- PLINK reference files, HapMap3 SNP lists, European-ancestry LD reference data, and GENCODE v46lift37 annotations are required for the relevant steps.

Package versions listed above are those recorded for the main analysis where available. Other dependency versions are not pinned here.

## Inputs

**Genomic SEM:** harmonized GWAS summary statistics for CD, UC, GERD, PUD, GASTD, IBS, CONST, DIVD, and HEM; CeD summary statistics for the preliminary ten-trait analysis; munged HapMap3 files; GenomicSEM SNP reference; and European LD-score files.

**gsMap:** eight adult proximal small-intestine Visium sections (P1, P5, P8, P9, P10, P11, P12, P13) from Harnik et al. (2024), including the authors' count matrices, spatial metadata, pixel classifications, and manual anatomical annotations. The gsMap stage requires the three final Q_SNP/local-LD-filtered factor GWAS files named `FUPPER.sumstats.gz`, `FBOWEL.sumstats.gz`, and `FIBD.sumstats.gz`, together with gsMap reference resources.

The input datasets and large reference files are not included in this repository. Harnik et al. source files and anatomical annotation code are available from https://github.com/yotamharnik/ST_Human_SI (reference commit: `cbccfecca66b7896e0595905ad3c01f4e6484ecf`).

## Running the scripts

Set `project_dir` in the R script, and `PROJECT_ROOT` and `GWAS_DIR` in the Python script, to local data directories. Prepare the input files in the directory structure expected by the scripts.

```bash
Rscript GenomicSEM_GI_core_analysis.R

python gsMap_GI_spatial_analysis.py --stage prepare
python gsMap_GI_spatial_analysis.py --stage gsmap
python gsMap_GI_spatial_analysis.py --stage anatomy
```

Run the three Python stages in the order shown. The anatomy stage uses the spatial-LDSC results from the preceding gsMap stage.

## Scope

The R script produces factor GWAS output before the downstream Q_SNP/local-LD filtering. The Python script takes the *filtered* factor GWAS summary statistics as input. GWAS preprocessing, downstream variant filtering, and the other analyses reported in the manuscript are not implemented in these two scripts.

The scripts require the original input datasets and specified software environments; results are not bundled in this repository.
