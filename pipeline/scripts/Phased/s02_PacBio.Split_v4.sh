#!/bin/bash
# s02_split_v4
#
# Cameron Brown 10Apr2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s06_split_v4

# PBS directives
#---------------

#PBS -N s02_split_v4
#PBS -l nodes=1:ncpus=16
#PBS -l walltime=00:30:00
#PBS -q half_hour
#PBS -m abe
#PBS -M cameron.brown.944@cranfield.ac.uk
#PBS -j oe
#PBS -v "CUDA_VISIBLE_DEVICES="
#PBS -W sandbox=PRIVATE
#PBS -k n

ln -s "$PWD" "$PBS_O_WORKDIR/$PBS_JOBID"

# Change to working directory
cd "$PBS_O_WORKDIR"

# Stop at runtime errors
set -e

# Threads
threads="${PBS_NCPUS:-${NCPUS:-1}}"

# Load module
module use /apps/modules/all
module load BCFtools/1.18-GCC-12.3.0

# Folders
base_folder="/mnt/beegfs/project/Alexey_Larionov/IBD-2026"
pipeline_folder="${base_folder}/data/processed/PacBio_Preprocessing"

split_folder="${pipeline_folder}/split"
counts_folder="${pipeline_folder}/counts"
filtered_folder="${pipeline_folder}/filtered_for_split"

input_vcf="${pipeline_folder}/merged_vcf/pacbio.merged.clean.vcf.gz"
filtered_vcf="${filtered_folder}/pacbio.merged.clean.gt_ps_only.vcf.gz"

# Prepare output folders
mkdir -p "${split_folder}"
mkdir -p "${counts_folder}"
mkdir -p "${filtered_folder}"

echo "----------------------------------------"
echo "STEP 1: Remove unneeded FORMAT fields"
date
echo "Input VCF: $(basename "${input_vcf}")"
echo "Filtered VCF: $(basename "${filtered_vcf}")"
echo "Threads used: ${threads}"
echo "----------------------------------------"
echo ""

# Keep only GT and PS in FORMAT
bcftools annotate \
  --threads ${threads} \
  -x ^FORMAT/GT,^FORMAT/PS \
  "${input_vcf}" \
  -Oz -o "${filtered_vcf}"

# Index filtered VCF
bcftools index -t "${filtered_vcf}"

echo "Written: $(basename "${filtered_vcf}")"
echo ""

echo "----------------------------------------"
echo "VALIDATION: FORMAT cleanup"
date
echo "----------------------------------------"

if [ ! -f "${filtered_vcf}" ]; then
  echo "FAIL: Filtered VCF not created"
  exit 1
fi
echo "PASS: Filtered VCF created"

if [ ! -f "${filtered_vcf}.tbi" ]; then
  echo "FAIL: Filtered VCF index not created"
  exit 1
fi
echo "PASS: Filtered VCF indexed"

format_count=$(bcftools view -h "${filtered_vcf}" | grep '^##FORMAT=' | wc -l)
gt_count=$(bcftools view -h "${filtered_vcf}" | grep '^##FORMAT=<ID=GT,' | wc -l)
ps_count=$(bcftools view -h "${filtered_vcf}" | grep '^##FORMAT=<ID=PS,' | wc -l)

if [ "${gt_count}" -ne 1 ]; then
  echo "FAIL: GT FORMAT field missing"
  exit 1
fi
echo "PASS: GT FORMAT field present"

if [ "${ps_count}" -ne 1 ]; then
  echo "FAIL: PS FORMAT field missing"
  exit 1
fi
echo "PASS: PS FORMAT field present"

if [ "${format_count}" -ne 2 ]; then
  echo "FAIL: Unexpected FORMAT fields remain"
  exit 1
fi
echo "PASS: Only GT and PS FORMAT fields remain"

echo ""
echo "----------------------------------------"
echo "STEP 2: Split filtered VCF by chromosome"
date
echo "----------------------------------------"
echo ""

# Split by chromosome
for i in {1..22}
do
  chr="chr${i}"
  out_vcf="${split_folder}/pacbio.merged.clean.gt_ps_only.${chr}.vcf.gz"

  echo "Processing ${chr}..."

  bcftools view \
    --threads ${threads} \
    -r "${chr}" \
    "${filtered_vcf}" \
    -Oz -o "${out_vcf}"

  bcftools index -t "${out_vcf}"

  bcftools +counts "${out_vcf}" > "${counts_folder}/counts_${chr}.txt"

  echo "Written: $(basename "${out_vcf}")"
done

echo ""
echo "----------------------------------------"
echo "VALIDATION: Split chromosome VCFs"
date
echo "----------------------------------------"

total_split_variants=0
filtered_variant_count=$(bcftools view -H "${filtered_vcf}" | wc -l)

for i in {1..22}
do
  chr="chr${i}"
  split_vcf="${split_folder}/pacbio.merged.clean.gt_ps_only.${chr}.vcf.gz"

  echo "Validating ${chr}..."

  if [ ! -f "${split_vcf}" ]; then
    echo "FAIL: Missing split VCF for ${chr}"
    exit 1
  fi
  echo "PASS: ${chr} VCF created"

  if [ ! -f "${split_vcf}.tbi" ]; then
    echo "FAIL: Missing index for ${chr}"
    exit 1
  fi
  echo "PASS: ${chr} VCF indexed"

  chrom_check=$(bcftools query -f '%CHROM\n' "${split_vcf}" | sort -u)

  if [ "${chrom_check}" != "${chr}" ]; then
    echo "FAIL: ${chr} VCF contains unexpected chromosome(s): ${chrom_check}"
    exit 1
  fi
  echo "PASS: ${chr} VCF contains only ${chr}"

  split_format_count=$(bcftools view -h "${split_vcf}" | grep '^##FORMAT=' | wc -l)
  split_gt_count=$(bcftools view -h "${split_vcf}" | grep '^##FORMAT=<ID=GT,' | wc -l)
  split_ps_count=$(bcftools view -h "${split_vcf}" | grep '^##FORMAT=<ID=PS,' | wc -l)

  if [ "${split_gt_count}" -ne 1 ] || [ "${split_ps_count}" -ne 1 ] || [ "${split_format_count}" -ne 2 ]; then
    echo "FAIL: ${chr} VCF does not contain only GT and PS FORMAT fields"
    exit 1
  fi
  echo "PASS: ${chr} VCF retains only GT and PS"

  split_variant_count=$(bcftools view -H "${split_vcf}" | wc -l)

  if [ "${split_variant_count}" -eq 0 ]; then
    echo "WARNING: ${chr} contains 0 variants"
  else
    echo "PASS: ${chr} contains ${split_variant_count} variants"
  fi

  total_split_variants=$((total_split_variants + split_variant_count))
done

echo "----------------------------------------"
echo "Filtered VCF total variants: ${filtered_variant_count}"
echo "Sum of split VCF variants: ${total_split_variants}"

if [ "${filtered_variant_count}" -ne "${total_split_variants}" ]; then
  echo "FAIL: Variant totals do not match after splitting"
  exit 1
fi
echo "PASS: Variant totals match after splitting"

echo "----------------------------------------"
echo "STEP COMPLETE: FORMAT cleanup and chromosome split finished"
date
echo "----------------------------------------"

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"