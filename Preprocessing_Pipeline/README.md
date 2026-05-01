## Data Preparation

A preprocessing pipeline for preparing short-read (Illumina) and long-read (PacBio) genotype data for identity-by-descent (IBD) analysis with [IBIS](https://github.com/williamslab/ibis), [RaPID](https://github.com/ZhiGroup/RaPID) and [Truffle](https://github.com/adimitromanolakis/truffle-website)

The pipeline is designed for HPC clusters with a `qsub`-compatible (PBS-style) job scheduler (tested on Crescent2) and uses [Singularity/Apptainer](https://apptainer.org/) containers to manage software dependencies.

### Pipeline Overview

- **Illumina (unphased):** Filter → PLINK conversion → Genetic map annotation  
- **PacBio (phased):** Per-sample filter → Merge → Chromosome split

### Repository Structure

```
Preprocessing_Pipeline/
├── Scripts/
│   ├── Unphased/               # Illumina (short-read) pipeline
│   │   ├── s01_Illumina_Filter.sh          # VCF filtering with bcftools
│   │   ├── s02_Illumina_PLINK_Conversion.sh # VCF → PLINK bed/bim/fam
│   │   └── s03_Illumina_Add_Genetic_Map.sh  # Add genetic distances to BIM
│   ├── Phased/                 # PacBio (long-read) pipeline
│   │   ├── 01_PacBio_Filter.sh             # Per-sample filter + merge
│   │   └── 02_PacBio_Split.sh              # Split merged VCF by chromosome
│   └── Add_Map_Plink/
│       └── add-map-plink.pl                # Perl script for genetic map annotation
├── config/
│   ├── Simple_Container.def    # Singularity definition (latest tool versions 01/05/26)
│   └── Complete_Container.def  # Singularity definition (pinned tool versions)
└── examples/
    ├── Illumina/               # Example PBS job logs for Illumina pipeline
    └── PacBio/                 # Example PBS job logs for PacBio pipeline
```

### Prerequisites

- A `qsub`-compatible job scheduler (PBS, Torque, PBSPro, or equivalent)
- Singularity or Apptainer
- A built container (`.sif` file) — see [Building the Container](#building-the-container)
- bgzipped and tabix-indexed input VCF files (`.vcf.gz` + `.tbi`) — see [Input File Requirements](#input-file-requirements) below
- GRCh38 genetic recombination map files (Illumina pipeline only) — available from the [Beagle genetic maps page](https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/)

---

## Input File Requirements

The two pipelines expect different input VCF formats.

### Illumina pipeline

Requires a **single multi-sample VCF** (`.vcf.gz`) containing all individuals. This is the standard output from array genotyping or joint-calling pipelines (e.g. GATK GenotypeGVCFs). The file must be bgzipped and tabix-indexed:

```bash
bcftools view -Oz -o output.vcf.gz input.vcf
bcftools index -t output.vcf.gz
```

### PacBio pipeline

Requires **one VCF per sample** (`.vcf.gz` + `.tbi`), all placed in the same folder.

If your data is in a single multi-sample VCF, split it into per-sample files using `bcftools`:

```bash
# Extract sample names and write one VCF per sample
bcftools query -l cohort.vcf.gz | while read sample; do
    bcftools view -s "$sample" -Oz -o "${sample}.vcf.gz" cohort.vcf.gz
    bcftools index -t "${sample}.vcf.gz"
done
```

This will generate:

```
sample1.vcf.gz
sample1.vcf.gz.tbi
sample2.vcf.gz
sample2.vcf.gz.tbi
...
```

Point `input_vcf_folder` in `01_PacBio_Filter.sh` at the directory containing these files.

---

## Building the Container

Two container definition files are provided in `Preprocessing_Pipeline/config/`:

| File | Description |
|------|-------------|
| `Simple_Container.def` | Builds with the latest available tool versions |
| `Complete_Container.def` | Builds with pinned versions (recommended for reproducibility) |

Both containers include: **bcftools**, **HTSlib**, **PLINK2**, **IBIS**, and **RaPID v1.7**.

Build with:

```bash
singularity build my_container.sif Complete_Container.def
```

---

## Configuration

Each script contains a **USER INPUTS** section near the top — this is the only part you need to edit. The variables present vary by script, but the three common ones are:

| Variable | Description |
|----------|-------------|
| `base_folder` | Your project root directory |
| `container` | Full path to your built `.sif` container file |
| `input_vcf` / `input_vcf_folder` | Path to your input VCF file or folder |

Example:

```bash
base_folder="/path/to/your/project"
container="/path/to/your/container.sif"
input_vcf="/path/to/your/data/input.vcf.gz"
```

Do **not** modify anything below the `# DO NOT EDIT BELOW` line in each script.

---

## Illumina (Unphased) Pipeline

Run these three scripts in order.

### Step 1 — Filter VCF (`s01_Illumina_Filter.sh`)

Filters the input VCF to biallelic autosomal SNPs passing quality thresholds.

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `input_vcf` | Path to bgzipped input VCF (`.vcf.gz`) |
| `container` | Path to the built `.sif` container |
| `min_qual` | Minimum QUAL score (default: 20) |

**Output:**
- `data/processed/Illumina_Preprocessing/clean_vcf/illumina_filtered.vcf.gz`
- bcftools stats files (raw and filtered)

**Submit:**
```bash
qsub s01_Illumina_Filter.sh
```

---

### Step 2 — Convert to PLINK (`s02_Illumina_PLINK_Conversion.sh`)

Converts the filtered VCF to PLINK binary format, applying genotype (`--geno 0.1`) and missingness (`--mind 0.1`) filters.

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `container` | Path to the built `.sif` container |

**Output:**
- `data/processed/Illumina_Preprocessing/plink/final_plink_files/illumina_filtered.{bed,bim,fam}`

**Submit:**
```bash
qsub s02_Illumina_PLINK_Conversion.sh
```

---

### Step 3 — Add Genetic Map (`s03_Illumina_Add_Genetic_Map.sh`)

Annotates the PLINK BIM file with genetic distances (cM) using `add-map-plink.pl`, then exports a sorted, indexed VCF.

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `container` | Path to the built `.sif` container |
| `perl_script` | Path to `add-map-plink.pl` |
| `map_folder` | Path to GRCh38 genetic map files (`.map`) |

> **Downloading genetic map files:** GRCh38 PLINK-format genetic maps (University of Washington / Beagle) can be downloaded from:
> **https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/**
> Choose `plink.GRCh38.map.zip` (chromosomes named `chr1`, `chr2`, … → use `chr_in_chrom_field`) or the equivalent without the `chr` prefix (→ use `no_chr_in_chrom_field`).

> **Chromosome naming:** use `chr_in_chrom_field` if your BIM uses `chr1`, `chr2`, … or `no_chr_in_chrom_field` if it uses `1`, `2`, …

**Output:**
- `data/processed/Illumina_Preprocessing/plink/mapped_plink_files/illumina_filtered_mapped.{bed,bim,fam}`
- `data/processed/Illumina_Preprocessing/plink/mapped_plink_vcf/illumina_filtered_mapped.vcf.gz` (sorted + CSI indexed)

**Submit:**
```bash
qsub s03_Illumina_Add_Genetic_Map.sh
```

---

## PacBio (Phased) Pipeline

Run these two scripts in order.

### Step 1 — Filter and Merge (`01_PacBio_Filter.sh`)

Filters each per-sample PacBio VCF to biallelic autosomal PASS SNPs, then merges them into a single VCF.

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `input_vcf_folder` | Folder containing per-sample `.vcf.gz` files |
| `container` | Path to the built `.sif` container |
| `min_qual` | Minimum QUAL score (default: 20) |

**Output:**
- `data/processed/PacBio_Preprocessing/merged_vcf/pacbio.merged.filtered.vcf.gz`
- Per-sample and merged bcftools stats files

**Submit:**
```bash
qsub 01_PacBio_Filter.sh
```

---

### Step 2 — Split by Chromosome (`02_PacBio_Split.sh`)

Strips all FORMAT fields except GT, then splits the merged VCF into one file per autosome (chr1–chr22).

**Edit the USER INPUTS section:**

| Variable | Description |
|----------|-------------|
| `base_folder` | Root project directory |
| `container` | Path to the built `.sif` container |

**Output:**
- `data/processed/PacBio_Preprocessing/split/pacbio.merged.filtered.gt_only.chr{1..22}.vcf.gz`
- Per-chromosome bcftools count files

**Submit:**
```bash
qsub 02_PacBio_Split.sh
```

---
## Expected Output Directory Structure

After running all scripts, your `data/processed/` directory will contain:

```
data/processed/
├── Illumina_Preprocessing/
│   ├── clean_vcf/
│   ├── plink/
│   │   ├── final_plink_files/
│   │   ├── mapped_plink_files/
│   │   └── mapped_plink_vcf/
│   ├── genetic_map/
│   └── stats/
└── PacBio_Preprocessing/
    ├── clean_vcf/
    ├── merged_vcf/
    ├── filtered_for_split/
    ├── split/
    ├── counts/
    └── stats/
```
## Notes

- All scripts print a **VALIDATION** section at the end of each job. Check the scheduler output logs (`.o<jobid>` for PBS-style schedulers) to confirm steps completed successfully. Example logs are provided in `Preprocessing_Pipeline/examples/`.
- `add-map-plink.pl` is sourced from the [IBIS repository](https://github.com/williamslab/ibis) and is included here for convenience.
- Scripts use `--bind /mnt/beegfs` for Singularity/Apptainer. Replace this with your HPC storage path (e.g. `/scratch`, `/data`, `/home`) if different.

---

## Institution Details

Cranfield University  
Supervisor: Dr Alexey Larionov  
Course Lead: Dr Maria Anastasiadi  
Support staff: Sajad Falsafi Zadeh  
Course: MSc Applied Bioinformatics 2025-26

