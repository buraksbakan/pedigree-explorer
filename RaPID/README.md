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
- Haplotype-based detection using a randomised path cover algorithm over a PBWT structure.
- Phased detection of IBD segments (phase-set PS tags must be removed for PacBio block-phased VCFs).
- cM-based filtering using a per-chromosome interpolated population genetic map.
- All-vs-all pairwise IBD detection across all samples in the input VCF.
- Per-chromosome processing with merged genome-wide output


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

### Input File Requirements

RaPID requires a **phased, GT-only VCF** and a per-chromosome interpolated genetic map

> **Important:** Phase-set (PS) tags must be removed from PacBio VCFs before running RaPID. Retaining PS information alongside block-phased genotypes causes segmentation faults. GT-only VCFs are produced during the [preprocessing pipeline](../Preprocessing_Pipeline) stage.

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

**Note**:The lenient setting (**-w 500**) exceeds the 300-SNP maximum defined in the author-provided `parameter_estimation.py` script. Results obtained under this configuration should therefore be interpreted with caution, as they fall outside the developer-recommended range. However, in our analysis, this setting provided improved sensitivity for detecting longer IBD segments.

---
## Usage

**Submit Step 1:**
```bash
qsub s01_map_interpolating.sh
```
**The script performs the following steps:**

1. Activates conda enviroment (rapid_env) and loads required modules
2. Interpolates GRCh38 genetic map to VCF site positions for each chromosome (chr1–22)

**Submit Step 2:**
```bash
qsub s02_rapid_3window_size.sh
```
**The script performs the following steps:**

3. Runs RaPID per chromosome for each window size configuration.
4. Merges per-chromosome results into a single file per window size

**Expected runtime:** 
- ~1 minute for genetic map interpolation  
- ~1 minute for running RaPID on two samples  

For large-scale datasets, see:  
Naseri, A., Liu, X., Tang, K. et al. (2019). *RaPID: ultra-fast, powerful, and accurate detection of segments identical by descent (IBD) in biobank-scale cohorts*. Genome Biology, 20, 143.  
https://doi.org/10.1186/s13059-019-1754-8

---

## Output Files

For each window size, per-chromosome results are merged:

**Example output structure:**

```
results/
├── rapid_w75
├── rapid_w250
├── rapid_w500
├── rapid_w75_merged
├── rapid_w250_merged
├── rapid_w500_merged

```
**Note:**RaPID with -r 10 produces numbered intermediate files (chr1/your_chromosome_seperated_vcf/0_results_sorted.max.gz...chr1/your_chromosome_seperated_vcf/9_results_sorted.max.gz) alongside the final **results.max.gz**. The merge steps
targets **results.max.gz** only

---

## Output Format

RaPID produces tab-delimited `results.max.gz` for each chromomosome within the given VCF input with the following columns. These seperate outputs are merged using RaPID/Scripts/s02_rapid_3window_size.sh :

| Column | Field | Description |
|--------|-------|-------------|
| 1 | chr_name | Chromosome |
| 2 | sample_id1 | First sample ID |
| 3 | sample_id2 | Second sample ID |
| 4 | hap_id1 | Haplotype ID of the first sample |
| 5 | hap_id2 | Haplotype ID of the second sample |
| 6 | starting_pos_genomic | Start position of IBD segment (bp) |
| 7 | ending_pos_genomic | End position of IBD segment (bp) |
| 8 | genetic_length | Genetic Length of the IBD segment (cM) |
| 9 | starting_site | Starting site of SNPs in the VCF file |
| 10 | ending_site | Ending site of SNPs in the VCF file |

**Output:** Files are converted to sorted BED format for visualisation in the Pedigree Explorer (GUI) 

```
chromosome  start_bp  end_bp  sample_pair
```

**Conversion logic:**
- chr_name column stays the same
- 2nd and 3rd column of BED file `start_bp = starting_pos_genomic`, `end_bp = ending_pos_genomic`
- Sample pair name is formatted as `sample_id1-sample_id2`
- Output is sorted by chromosome and start position for downstream tools
 
---

## Notes

- Process chromosomes individually; concatenating maps before interpolation causes RaPID to process only the last chromosome
- The merge step targets `results.max.gz` only — not the numbered intermediate files
- The author-provided `parameter_estimation.py` script produced division-by-zero errors and negative window-size values in this analysis; window sizes were therefore set manually
- - RaPID source and genetic maps: https://github.com/ZhiGroup/RaPID


---
## Citation

If you use this pipeline, please cite the RaPID paper:

> Naseri, A., Liu, X., Tang, K. et al. RaPID: ultra-fast, powerful, and accurate detection of segments identical by descent (IBD) in biobank-scale cohorts. Genome Biol 20, 143 (2019). https://doi.org/10.1186/s13059-019-1754-8.

---

## Institution Details

Cranfield University  
Supervisor: Dr Alexey Larionov  
Course Lead: Dr Maria Anastasiadi  
Support staff: Sajad Falsafi Zadeh  
Course: MSc Applied Bioinformatics 2025-26

