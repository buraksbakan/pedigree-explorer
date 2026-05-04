# RaPID IBD DETECTION PIPELINE

IBD (Identity-by-Descent) segment detection using RaPID v1.7 (Naseri et al.,2019), a fast haplotype based IBD detector designed for phased genotype data.

This pipeline was developed as part of the **Pedigree Explorer** project for detection and visualisation of shared chromosomal regions in pedigree sequencing. It runs RaPID
across three window size configurations (strict w = 75, default w = 250, and lenient w = 500) on phased VCF input data then merges per-chromosome outputs for downstream 
analysis.

---

## Overview

RaPID uses a randomised path cover algorithm over positional Burrows-Wheeler transform (PBWT) based window approach. A segment is called when a sufficient number of independent
random projections covers the agreed region.

**Key features:**
- Phased


## Repository Structure

```
RaPID/
├── README.md                # This file
├── rapid_env.yml            # Conda Environment used for running the script
└── Scripts/
    └── s01_map_interpolating.sh             # Step 1: Interpolate genetic map to VCF sites
    └── s02_rapid_3window_size.sh            # Step 2: Run RaPID across three different window sizes
```

---

## Prerequisites

- A `qsub`-compatible job scheduler (PBS, Torque, PBSPro, or equivalent) — tested on Crescent2 at Cranfield University
- RaPID v1.7 (loaded via module system)
- `python3` (for running s01_map_interpolating.sh)
- `bedtools` (for sorting; loaded via module system)
- `awk` (standard Linux utility)
- Multi-sample, chromosome split, Block-phased GT-only PacBio VCFs — produced by the [preprocessing pipeline](../Preprocessing_Pipeline)
- GRCh38 reference genetic maps - available from the RaPID GitHub repository

---

## Configuration

The script contains a **USER-DEFINED INPUTS** section near the top — this is the only part you need to edit:

| Variable | Description |
|----------|-------------|
| `BASE_DIR` | Your project root directory |
| `SCRIPT_DIR` | Script directory for interpolate_loci.py |
| `MAP_DIR` | Genetic Map directory for interpolating VCFs |
| `VCF_DIR` | VCF directory |
| `INTERPOLATED_MAP_DIR` | Chromosome specific interpolated map directory  |
| `OUTPUT_BASE` | Output Directory |

Example:

```bash
BASE_DIR="/path/to/your/project"
SCRIPT_DIR="${BASE_DIR}/script/directory/interpolate_loci.py"
MAP_DIR="${BASE_DIR}/genetic/maps"
VCF_DIR="${BASE_DIR}/path/to/split/vcfs"
INTERPOLATED_MAP_DIR="$BASE_DIR/platinum_pedigree_dataset/interpolated_maps"
OUTPUT_BASE="$BASE_DIR/platinum_pedigree_dataset/results"
```
---

## Parameter Configurations

The script tests **three parameter configurations** to assess the effect of detection stringency on IBD calls. These configurations vary by minimum segment length, minimum marker count, and error rate tolerance:

| Configuration | `-w` (window size) | `-r` (runs) | `-s`  (successes)| `-d` (min cM)| Use case |
|---------------|---------------|-----------------|--------------|----------|----------|
| **Strict** | 75 | 10 | 2 | 5 |High-specificity detection; retains only the strongest, longest shared haplotypes |
| **Default** (literature-recommended) | 250 | 10 | 2 | 5 |Balanced detection for general pedigree IBD analysis; based on the RaPID repository (ZhiGroup, 2019) |
| **Lenient** | 500 | 10 | 2 | 5 |Exploratory analysis; maximises segment recovery at the cost of increased false positives; produces biology sensible results |

---

## RaPID flags explained:

- **-w** — window size in markers; smaller windows increase sensitivity but may reduce specificity  
- **-r** — number of independent random path cover iterations  
- **-s** — minimum iterations that must agree for a segment call  
- **-d** — minimum IBD segment length in cM  
- **-i** — input VCF  
- **-g** — interpolated genetic map file  
- **-o** — output prefix  

**Note:** The lenient setting (**-w 500**) exceeds the 300-SNP maximum defined in the author-provided `parameter_estimation.py` script. Results obtained under this configuration should therefore be interpreted with caution, as they fall outside the developer-recommended range. However, in our analysis, this setting provided improved sensitivity for detecting longer IBD segments.

---


