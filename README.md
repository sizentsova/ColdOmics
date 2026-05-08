# Multi-omics Snakemake Pipeline for Robust TF Identification in *Arabidopsis thaliana*

## Overview

This Snakemake pipeline integrates **transcriptomics, cistromics, regulomics, genomics, and climatomics** data to identify **robust transcription factors (TFs)** associated with cold stress response in *Arabidopsis thaliana*.

The workflow combines differential expression analysis, TF enrichment, regulatory network analysis, variant annotation, and climate association testing into a unified multi-layer framework.

---

## Workflow Summary

The pipeline consists of the following steps:

### 1. Transcriptomics

* Processes microarray and RNA-seq datasets
* Identifies differentially expressed genes (DEGs)
* Clusters time-series data into:

  * Early Response (ER)
  * Late Response (LR)
  * Very Late Response (VR)
* Applies binomial filtering to identify **robust DEGs**

### 2. Cistromics (TF enrichment)

* Uses external CisCross results (**manual step required**)
* Identifies TFs enriched in DEG sets

### 3. Genomics

* Extracts TF genomic regions (BED)
* Subsets variants from 1001 Genomes VCF
* Annotates variants using VEP
* Extracts genotype matrices (ALT / REF / NA)

### 4. Regulomics

* Computes TF connectivity from regulatory networks

### 5. Climatomics

* Performs Monte Carlo permutation tests
* Associates TF variants with seasonal temperature data

### 6. Integration

* Combines all omics layers
* Produces a final table of **robust multi-layer TFs**

---

## Final Outputs

```
results/
├── transcriptomics/
│   ├── robustDEG.csv
│   ├── binomial_threshold.txt
│   ├── upER.csv
│   ├── downER.csv
│   ├── upLR.csv
│   ├── downLR.csv
│   ├── upVR.csv
│   └── downVR.csv
├── cistromics/
│   └── tf_cis_enrichment.csv
├── regulomics/
│   ├── Connectivity_per_TF.csv
│   └── Connectivity_network_summary.csv
├── genomics/
│   └── tf_variants_genotype_alt.csv
├── climatomics/
│   └── Monte-Carlo_permutation_robustTF.csv
└── integrated/
    └── robust_TF_multilayer.csv
```

---

## Installation

### Requirements

* Snakemake (≥7)
* Conda

### Setup

```bash
git clone git@github.com:sizentsova/ColdOmics.git
cd ColdOmics
```

Run with automatic environment creation:

```bash
snakemake --use-conda --cores N
```

Example:

```bash
snakemake --use-conda --cores 4
```

---

## Input Data Structure

```
data/
├── transcriptomics/   # DEGs extracted from RNA-seq and microarray transcriptoms
├── cistromics/        # CisCross TF enrichment and network results
├── genomics/
│   └── 1001genomes_snp-short-indel_only_ACGTN.vcf.gz
└── climatomics/
    └── CHELSA_min_temp_ALL_1058_accessions.xlsx
```

---

## Genomic Data Requirement

The pipeline requires the Arabidopsis 1001 Genomes VCF file:

`1001genomes_snp-short-indel_only_ACGTN.vcf.gz`

### Download

Download from:

https://1001genomes.org/data/GMI-MPI/releases/current/

### Placement

Place the file in:

```
data/genomics/
```

Final path:

```
data/genomics/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
```

> The pipeline will automatically index the VCF using `bcftools`.

---

## Manual Step: CisCross TF Network Analysis

This pipeline requires manual interaction with the CisCross web tool:

https://plamorph.sysbio.ru/ciscross/FindTFnet_index.html

---

### Step 1: Prepare DEG Lists

After transcriptomics, the following files are generated:

```
results/transcriptomics/
├── upER.csv
├── downER.csv
├── upLR.csv
├── downLR.csv
├── upVR.csv
├── downVR.csv
```

Each file contains **TAIR IDs of robust DEGs**.

---

### Step 2: Upload DEG Lists to CisCross

For each response phase:

| Phase | Upload        |
| ----- | ------------- |
| ER    | upER & downER |
| LR    | upLR & downLR |
| VR    | upVR & downVR |

#### Important:

* Upload **UP and DOWN lists separately**
* Provide **only TAIR IDs (one per line)**
* Remove headers if required by the tool

---

### Step 3: Download CisCross Results

For each phase (ER, LR, VR) save separate Ciscross outputs and place them in the following directory:

```
data/cistromics/
```

Example structure:

```
data/cistromics/
├── ciscross_ER_1500p_0.05.txt
├── ciscross_LR_1500p_0.05.txt
├── ciscross_VR_1500p_0.05.txt
```
Example files are located in the `resource/ciscross_examples/` directory. These can be copied into `data/cistromics/` for testing the ColdOmics pipeline.

---

## Step 4. Signal completion

After placing all required files, create a completion marker file:

```bash
touch data/cistromics/ciscross.complete
```

This file signals that all inputs are ready and the pipeline can proceed. The snakemake pipeline will wait up to **30 minutes** for these files.

---

### Notes

* File naming must match expected patterns (ER/LR/VR, up/down)
* CisCross column names should remain unchanged
* The pipeline standardizes formatting internally

---

## Scripts

| Step              | Script                                     |
| ----------------- | ------------------------------------------ |
| Transcriptomics   | `scripts/1.Transcriptomics.R`              |
| Cistromics        | `scripts/2.Cistromics.R`                   |
| Regulomics        | `scripts/3.Regulomics.R`                   |
| TF regions        | `scripts/4.Extract_TF_regions.R`           |
| Variant genotypes | `scripts/4.Extract_variant_genotypes.R`    |
| Climatomics       | `scripts/5.Climatomics.R`                  |
| Integration       | `scripts/6.Integrate_multi-omics_layers.R` |

---


