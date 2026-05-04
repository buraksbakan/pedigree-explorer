# IBIS IBD Detection Pipeline

IBD (Identity-by-Descent) segment detection using **IBIS v1.20.9** (Seidman et al., 2020), a fast, scalable IBD detector for un-phased genotype data. IBIS uses centiMorgan (cM)-based filtering with a population genetic map.

This pipeline was developed as part of the **Pedigree Explorer** project for detection and visualisation of shared chromosomal regions in pedigree sequencing. It runs IBIS across three parameter configurations (sensitive, literature-recommended, and strict) on PLINK binary input data, then converts the output `.seg` files to BED format for visualisation in the [Pedigree Explorer GUI](../GUI/README.txt).

---

## Overview

IBIS detects IBD regions by scanning the genome for contiguous regions of high identity-by-state (IBS) between pairs of individuals. A window-based approach identifies candidate segments while tolerating a limited number of mismatches, reducing fragmentation caused by genotyping errors, isolated recombination events, or de novo mutations.

By default, IBIS detects only **IBD1 sharing** (regions where one haplotype is shared). Optional `-2` flag enables detection of IBD2 segments.

**Key features:**
- Phase-free detection (no haplotype phasing required)
- cM-based filtering using a population genetic map
- All-vs-all pairwise mode for cohort analysis
- Fast: ~7 minutes on 3,000 individuals

In the Pedigree Explorer project, IBIS is applied alongside [TRUFFLE](../Truffle) (un-phased, IBS marker-based) and [RaPID](../RaPID) (phased) to enable cross-validation of IBD calls across detection algorithms.

---

## Repository Structure

```
IBIS/
├── README.md                # This file
└── Scripts/
    └── ibis.sh              # Run IBIS across three parameter configurations
```

---

## Prerequisites

- A `qsub`-compatible job scheduler (PBS, Torque, PBSPro, or equivalent) — tested on Crescent2 at Cranfield University
- IBIS v1.20.9 (loaded via module system or installed locally)
- `bedtools` (for sorting; loaded via module system)
- `awk` (standard Linux utility)
- PLINK binary input files (`.bed`, `.bim`, `.fam`) with genetic distances annotated in the `.bim` file — typically produced by the [preprocessing pipeline](../../Preprocessing_Pipeline)

### Input File Requirements

IBIS requires **PLINK binary format** with the genetic distance column populated in the `.bim` file. The preprocessing pipeline produces the required files:

```
input_prefix.bed    # Binary genotype data
input_prefix.bim    # Variant info with genetic distances (cM)
input_prefix.fam    # Sample/family info
```

To prepare these files, use:
- `s02_Illumina_PLINK_Conversion.sh` — Converts VCF to PLINK
- `s03_Illumina_Add_Genetic_Map.sh` — Adds genetic distances using the GRCh38 Beagle genetic map

See the [preprocessing pipeline README](../Preprocessing_Pipeline) for details.

---

## Configuration

The script contains a **USER-DEFINED INPUTS** section near the top — this is the only part you need to edit:

| Variable | Description |
|----------|-------------|
| `BASE_DIR` | Your project root directory |
| `PLINK_INPUT_1` | PLINK file prefix for first dataset (no extension) |
| `PLINK_INPUT_2` | PLINK file prefix for second dataset (no extension) |
| `OUT_DIR_1` | Output directory for first dataset |
| `OUT_DIR_2` | Output directory for second dataset |

Example:

```bash
BASE_DIR="/path/to/your/project"
PLINK_INPUT_1="${BASE_DIR}/data/plink_dataset_1/illumina_filtered_mapped"
PLINK_INPUT_2="${BASE_DIR}/data/plink_dataset_2/illumina_filtered_mapped"
OUT_DIR_1="${BASE_DIR}/results/dataset_1/ibis"
OUT_DIR_2="${BASE_DIR}/results/dataset_2/ibis"
```

The script runs IBIS sequentially on both datasets; comment out unused dataset blocks if you only have one input.

---

## Parameter Configurations

The script tests **three parameter configurations** to assess the effect of detection stringency on IBD calls. These configurations vary by minimum segment length, minimum marker count, and error rate tolerance:

| Configuration | `-min_l` (cM) | `-mt` (markers) | `-errorRate` | Use case |
|---------------|---------------|-----------------|--------------|----------|
| **Sensitive** (lenient) | 3 | 50 | 0.06 | Exploratory analysis; captures shorter segments at risk of background noise |
| **Seidman** (literature-recommended) | 7 | 400 | 0.004 | Robust default for general relatedness estimation; from Seidman et al. (2020) |
| **Strict** | 10 | 200 | 0.20 | High-specificity analyses (validation, fine-mapping); may underestimate total sharing |

**IBIS flags explained:**
- `-min_l` — minimum segment length in cM
- `-mt` — minimum number of markers per segment
- `-errorRate` — tolerated genotyping error rate per segment

**Configuration choice guidance:**
- For initial discovery: use **sensitive** to identify candidate regions
- For relationship inference: use **Seidman** as the standard reference
- For high-confidence calls: use **strict** to focus on the strongest shared haplotypes

The literature-recommended configuration (Seidman et al., 2020) provides a balanced default and is generally suitable for most analyses.

---

## Usage

**Submit:**
```bash
qsub ibis.sh
```

**The script performs the following steps:**

1. Loads required modules (`ibis`, `bedtools`)
2. Runs IBIS on dataset 1 with all three configurations
3. Runs IBIS on dataset 2 with all three configurations
4. Converts all `.seg` outputs to sorted BED format

**Expected runtime:** ~10–30 minutes for typical pedigree-sized inputs (~5 samples). Larger cohorts will take longer; adjust `walltime` accordingly.

---

## Output Files

For each dataset and configuration, IBIS produces:

| File | Description |
|------|-------------|
| `{config}.seg` | Tab-delimited segment output (sample1, sample2, chrom, phys_start, phys_end, IBD_type, gen_start, gen_end, gen_length_cM, marker_count, error_count, error_density) |
| `{config}.coef` | Per-pair kinship coefficients and degree of relatedness *(if `-printCoef` enabled)* |
| `{config}.bed` | Converted BED format (produced by post-processing block in the script) |

**Example output structure:**

```
results/dataset_1/ibis/
├── lenient.seg
├── lenient.bed
├── seidman.seg
├── seidman.bed
├── strict.seg
└── strict.bed
```

---

## .seg Output Format

IBIS produces tab-delimited `.seg` files with the following columns:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | sample1 | First sample ID |
| 2 | sample2 | Second sample ID |
| 3 | chrom | Chromosome |
| 4 | phys_start_pos | Start position (bp) |
| 5 | phys_end_pos | End position (bp) |
| 6 | IBD_type | `IBD1` or `IBD2` |
| 7 | genetic_start_pos | Start position (cM) |
| 8 | genetic_end_pos | End position (cM) |
| 9 | genetic_seg_length | Segment length (cM) |
| 10 | marker_count | Number of markers in segment |
| 11 | error_count | Number of errors tolerated |
| 12 | error_density | Errors per marker |

---

## BED Conversion

The script automatically converts each `.seg` file to a four-column BED format compatible with the [Pedigree Explorer GUI](../../GUI):

```
chromosome  start_bp  end_bp  sample_pair
```

**Conversion logic:**
- BED is 0-based, half-open; PLINK/IBIS positions are 1-based
- Therefore `start = phys_start - 1`, `end = phys_end`
- Sample pair name is formatted as `sample1-sample2`
- Output is sorted by chromosome and start position for downstream tools

Example output:
```
chr1   75230111   75245678   IHCAPX8-1-IHCAPX8-2
chr1  169631899  171007200   IHCAPX8-1-IHCAPX8-2
chr2   44105199   64252800   IHCAPX8-1-IHCAPX8-2
```

---

## Important Notes

### cM vs Mb Filtering

IBIS uses centiMorgan (cM)-based thresholds via the population genetic map embedded in the `.bim` file. This is more biologically meaningful than physical-distance (Mb) filtering because segment length is determined by recombination events, which vary across the genome. As a result, IBIS automatically accounts for regional variation in recombination rate (e.g. centromeric vs subtelomeric regions).

### Detection Sensitivity

IBIS reliably detects segments **≥ 7 cM** with quality comparable to Refined IBD and GERMLINE (Seidman et al., 2020). Below this threshold, sensitivity decreases and false positive rates increase. The Seidman configuration is calibrated to this threshold.

For detecting more distant relationships (5th–6th degree relatives), the sensitive configuration may capture additional shorter segments, but with reduced specificity.

### Genotyping Error Tolerance

The `-errorRate` parameter controls how many mismatches are tolerated within a segment. Higher values are more permissive but may inflate segment counts:

- 0.004 (Seidman) — strict; suitable for high-quality data
- 0.06 (sensitive) — permissive; tolerates short stretches of mismatch
- 0.20 (strict) — counterintuitively permissive on errors but combined with longer length and high marker count thresholds

---

## Notes

- The script uses `set -euo pipefail` for safer execution — any error stops the script.
- IBIS source: https://github.com/williamslab/ibis
- Output `.seg` files for high-IBD pairs (e.g. siblings) can be large — ensure adequate storage.
- All scripts print status messages with timestamps. Check the scheduler output logs (`.o<jobid>` for PBS-style schedulers) to confirm completion.

---

## Citation

If you use this pipeline, please cite the IBIS paper:

> Seidman, D.N., Shenoy, S.A., Kim, M., Babu, R., Woods, I.G., Dyer, T.D., Lehman, D.M., Curran, J.E., Duggirala, R., Blangero, J., & Williams, A.L. (2020). Rapid, phase-free detection of long identity-by-descent segments enables effective relationship classification. *The American Journal of Human Genetics*, 106(4), 453–466.

---

## Institution Details

Cranfield University  
Supervisor: Dr Alexey Larionov  
Course Lead: Dr Maria Anastasiadi  
Support staff: Sajad Falsafi Zadeh  
Course: MSc Applied Bioinformatics 2025-26
