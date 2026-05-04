# TRUFFLE IBD Detection Pipeline

IBD (Identity-by-Descent) segment detection using **TRUFFLE v1.38** (Dimitromanolakis et al., 2019), a fast, phase-free IBD detector designed for un-phased VCF data.

This pipeline includes scripts for default analysis, parameter sensitivity testing, IBS marker threshold optimization, and conversion of TRUFFLE output to BED format for visualization.

---

## Overview

TRUFFLE detects long stretches of consecutive Identity-By-State (IBS) markers between pairs of individuals to infer IBD segments. It distinguishes:
- **IBD1** — regions where at least one allele is shared
- **IBD2** — regions where both alleles are shared

Unlike phase-based methods (BEAGLE Refined IBD, HapIBD), TRUFFLE does not require haplotype phasing or a genetic map, making it suitable for rapid analysis of moderately sized cohorts.

---

## Repository Structure

```
TRUFFLE/
├── README.md                            # This file
└── Scripts/
    ├── truffle_default.sh               # Run with default parameters
    ├── truffle_L_sensitivity.sh         # Test multiple -L values (length sensitivity)
    ├── truffle_ibs_marker_filter.sh     # Optimized IBS marker thresholds
    └── truffle_to_bed.sh                # Convert .segments output → BED format
```

---

## Prerequisites

- A `qsub`-compatible job scheduler (PBS, Torque, PBSPro, or equivalent) — tested on Crescent2
- Singularity or Apptainer
- A built container (`.sif` file) containing TRUFFLE — see the [preprocessing pipeline](../Preprocessing_Pipeline/README.md)
- A bgzipped, tabix-indexed multi-sample VCF (`.vcf.gz` + `.tbi`) — typically the output of the [preprocessing pipeline](../Preprocessing_Pipeline/README.md)
- `awk` (for BED conversion — standard Linux utility)

### Input File Requirements

TRUFFLE requires a **single multi-sample VCF** containing all individuals to be compared. The VCF should be:
- bgzipped (`.vcf.gz`)
- tabix-indexed (`.tbi`)
- Filtered for autosomal biallelic SNPs (use the [preprocessing pipeline](.//Preprocessing_Pipeline/README.md))
- GRCh38 (or any consistent reference — TRUFFLE is reference-agnostic)

---

## Configuration

All scripts contain a **USER INPUTS** section near the top — this is the only part you need to edit. The variables present vary by script, but the three common ones are:

| Variable | Description |
|----------|-------------|
| `base_folder` | Your project root directory |
| `container` | Full path to your built `.sif` container file |
| `input_vcf` | Path to your input multi-sample VCF |

Example:
```bash
base_folder="/path/to/your/project"
container="/path/to/your/container.sif"
input_vcf="/path/to/your/data/input.vcf.gz"
```

Do **not** modify anything below the `# DO NOT EDIT BELOW` line in each script.

---

## Recommended Workflow

For most projects, run scripts in this order:

1. **Baseline analysis** (`truffle_default.sh`) — get a baseline IBD detection with default parameters
2. **Sensitivity check** (`truffle_L_sensitivity.sh`) — assess robustness across `-L` values *(optional but recommended)*
3. **Threshold optimization** (`truffle_ibs_marker_filter.sh`) — *only if* default detection is insufficient (e.g., for distant relatives)
4. **Format conversion** (`truffle_to_bed.sh`) — convert final output to BED for GUI visualization

---

## Step 1 — Default Run (`truffle_default.sh`)

Runs TRUFFLE with default parameters (`-L 1.0`, equivalent to a minimum 5 Mb for IBD1 and 2 Mb for IBD2). Use this as a baseline before any optimization.

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `container` | Path to the built `.sif` container |
| `input_vcf` | Path to bgzipped input VCF |
| `output_dir` | Output directory |
| `output_prefix` | Prefix for output files |

**Output:**
- `{output_dir}/truffle_default.segments` — IBD segment list

**Submit:**
```bash
qsub truffle_default.sh
```

**Interpretation:**
- For close relatives (siblings, parent-offspring): default parameters typically recover ~50% sharing as expected — no further optimization needed.
- For distant relatives (cousins or beyond): default parameters may be too stringent and detect very few or no segments — proceed to Step 3.

---

## Step 2 — L Parameter Sensitivity (`truffle_L_sensitivity.sh`)

Tests TRUFFLE across multiple `-L` values (1.0, 1.5, 2.0, 2.5, 3.0). The `-L` parameter is a **length multiplier** applied to TRUFFLE's default minimum segment lengths:

**Higher L = more stringent detection**, fewer segments reported.

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `container` | Path to the built `.sif` container |
| `input_vcf` | Path to bgzipped input VCF |
| `output_dir` | Output directory |
| `L_values` | Array of L values to test |

**Output:**
- `{output_dir}/truffle_L1.0.segments`
- `{output_dir}/truffle_L1.5.segments`
- `{output_dir}/truffle_L2.0.segments`
- ... etc.

**Submit:**
```bash
qsub truffle_L_sensitivity.sh
```

---

## Step 3 — IBS Marker Threshold Optimization (`truffle_ibs_marker_filter.sh`)

For datasets where default parameters fail to detect IBD (typically distant relatives or whole-genome sequencing data), this script applies **direct IBS marker count thresholds** for finer control over detection sensitivity.

The `--ibs1markers` and `--ibs2markers` parameters control the minimum number of consecutive IBS markers required to call an IBD1 or IBD2 segment. Lower values = more sensitive detection.

**Why this is needed:**

TRUFFLE's defaults were optimized for ~400,000-marker SNP array data (Dimitromanolakis et al., 2019). On WGS data with millions of markers, the implicit physical-length thresholds can become disproportionately large, requiring manual recalibration. This is particularly important for distant relatives whose IBD segments are inherently shorter due to extensive recombination.

**Suggested thresholds:**

| Relationship | Suggested IBS1 / IBS2 markers | Notes |
|--------------|-------------------------------|-------|
| Siblings, parent-offspring | Use defaults (Step 1) | Short stretches are common; defaults work |
| First cousins | `7000 / 1500` | ~5 Mb / 2 Mb minimum at WGS density |
| First cousins once removed | `5000 / 1000` | Detects shorter segments |
| Second cousins or more distant | `4000 / 800` | Very permissive — risk of false positives |

**Empirical calibration is recommended:** if a known relationship exists in your dataset, tune thresholds so the detected IBD percentage matches the expected value for that relationship.

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `container` | Path to the built `.sif` container |
| `input_vcf` | Path to bgzipped input VCF |
| `output_dir` | Output directory |
| `ibs1_markers` | Minimum consecutive IBS1 markers |
| `ibs2_markers` | Minimum consecutive IBS2 markers |

**Output:**
- `{output_dir}/truffle_ibs{ibs1}_{ibs2}.segments`

**Submit:**
```bash
qsub truffle_ibs_marker_filter.sh
```

---

## Step 4 — Convert to BED (`truffle_to_bed.sh`)

Converts TRUFFLE `.segments` output files to BED format suitable for GUI(GUI/README.txt). Optionally filters segments by minimum length to remove low-confidence calls.

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `input_dir` | Folder containing `.segments` files |
| `output_dir` | Output folder for BED files |
| `input_pattern` | Glob pattern for files to convert (e.g., `truffle_L*.segments`) |
| `min_length_mb` | Minimum segment length to keep (default: 2.0 Mb) |

**Output BED format:**

```
chromosome  start_bp  end_bp  sample_pair  IBD_type
```

Example:
```
chr1   169631900   171007200   IHCAPX8-1_IHCAPX8-2   IBD2
chr1    39997400    46929600   IHCAPX8-1_IHCAPX8-2   IBD1
chr2    44105200    64252800   IHCAPX8-1_IHCAPX8-2   IBD1
```

**Submit:**
```bash
qsub truffle_to_bed.sh
```

---

## Output File Format

TRUFFLE produces tab-delimited `.segments` files with the following columns:

| Column | Field | Description |
|--------|-------|-------------|
| 1 | TYPE | `IBD1` or `IBD2` |
| 2 | ID1 | First sample ID |
| 3 | ID2 | Second sample ID |
| 4 | CHROM | Chromosome |
| 5 | VARSTART | Start variant index |
| 6 | VAREND | End variant index |
| 7 | POS | Start position (Mbp) |
| 8 | (label) | `Mbp` |
| 9 | LENGTH | Segment length (Mbp) |
| 10 | (label) | `Mbp` |
| 11 | NMARKERS | Number of markers in segment |

---

## Calculating Percentage of Genome Shared

To convert TRUFFLE segment lengths to a percentage of the genome:

```
% Genome = (Σ segment lengths in Mb) / 2,881 Mb × 100
```

Where 2,881 Mb is the total autosomal genome size (GRCh38, chromosomes 1–22), excluding sex chromosomes and mitochondrial DNA — consistent with TRUFFLE's autosomal-only analysis.

This metric enables direct comparison with theoretical kinship values:

| Relationship | Expected % shared |
|--------------|-------------------|
| Full siblings | ~50% |
| Parent-offspring | ~50% |
| First cousins | ~12.5% |
| First cousins once removed | ~6.25% |
| Second cousins | ~3.13% |
| Second cousins once removed | ~1.56% |

*Reference: ISOGG (2024) Autosomal DNA statistics — https://isogg.org/wiki/Autosomal_DNA_statistics*

---

## Important Notes

### Detection Power Limitations

From the TRUFFLE paper:
- Detection power is **>80%** for segments **>5 Mb**
- Detection power drops to **<5%** for segments **<4 Mb**

This is a fundamental sensitivity limit for distant relationships where most IBD segments are shorter than 4 Mb.

### No Genetic Map Filtering

Unlike IBIS and RaPID, TRUFFLE operates on physical distance (Mbp) only, not genetic distance (cM). This contrasts with cM-based methods that automatically adapt to local recombination rate variation. Consider cross-validating distant-relative findings with IBIS or RaPID.

### Singularity Bind Path

Scripts use `--bind /` for Singularity. Adjust this if your HPC requires a more specific bind (e.g., `--bind /scratch`, `--bind /data`, `--bind /home`).

---

## Citation

If you use this pipeline, please cite:

> Dimitromanolakis, A., Paterson, A.D., & Sun, L. (2019). Fast and accurate shared segment detection and relatedness estimation in un-phased genetic data via TRUFFLE. *The American Journal of Human Genetics*, 105(1), 78–88.

TRUFFLE source: https://github.com/adimitromanolakis/truffle

---

## Institution Details

Cranfield University  
Supervisor: Dr Alexey Larionov  
Course Lead: Dr Maria Anastasiadi  
Support staff: Sajad Falsafi Zadeh  
Course: MSc Applied Bioinformatics 2025-26
