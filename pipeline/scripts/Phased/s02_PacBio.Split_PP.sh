#!/bin/bash
# s02_split_PP
#
# Cameron Brown 10Apr2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s02_split_PP

# PBS directives
#---------------

#PBS -N s02_split_PP
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

cd "$PBS_O_WORKDIR"
set -e

# Threads
threads="${PBS_NCPUS:-${NCPUS:-1}}"

# Folders
base_folder="/mnt/beegfs/project/Alexey_Larionov/IBD-2026"
pipeline_folder="${base_folder}/data/processed/PacBio_Preprocessing_PP"

split_folder="${pipeline_folder}/split"
counts_folder="${pipeline_folder}/counts"
filtered_folder="${pipeline_folder}/filtered_for_split"

input_vcf="${pipeline_folder}/merged_vcf/pacbio.merged.clean.vcf.gz"
filtered_vcf="${filtered_folder}/pacbio.merged.clean.gt_only.vcf.gz"

# Container
container="${base_folder}/github_repo/pedigree-explorer/pipeline/config/Simple_Container.sif"

# Prepare output folders
mkdir -p "${split_folder}" "${counts_folder}" "${filtered_folder}"

# Check inputs
if [ ! -f "${input_vcf}" ]; then
  echo "ERROR: Input VCF not found: ${input_vcf}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

echo "----------------------------------------"
echo "STEP 1: Keep only GT FORMAT field"
date
echo "----------------------------------------"

singularity exec --bind /mnt/beegfs "${container}" \
  bcftools annotate \
    --threads "${threads}" \
    -x ^FORMAT/GT \
    "${input_vcf}" \
    -Oz -o "${filtered_vcf}"

singularity exec --bind /mnt/beegfs "${container}" \
  bcftools index -t "${filtered_vcf}"

echo "Written: $(basename "${filtered_vcf}")"
echo ""

echo "----------------------------------------"
echo "VALIDATION: FORMAT cleanup"
date
echo "----------------------------------------"

[ -f "${filtered_vcf}" ] || { echo "FAIL: No filtered VCF"; exit 1; }
[ -f "${filtered_vcf}.tbi" ] || { echo "FAIL: No index"; exit 1; }

echo "PASS: File + index exist"

format_count=$(singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -h "${filtered_vcf}" | grep '^##FORMAT=' | wc -l)

gt_count=$(singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -h "${filtered_vcf}" | grep '^##FORMAT=<ID=GT,' | wc -l)

[ "${gt_count}" -eq 1 ] || { echo "FAIL: GT missing"; exit 1; }
[ "${format_count}" -eq 1 ] || { echo "FAIL: Extra FORMAT fields present"; exit 1; }

echo "PASS: Only GT present"

echo ""
echo "----------------------------------------"
echo "STEP 2: Split by chromosome"
date
echo "----------------------------------------"

for i in {1..22}
do
  chr="chr${i}"
  out_vcf="${split_folder}/pacbio.merged.clean.gt_only.${chr}.vcf.gz"

  echo "Processing ${chr}..."

  singularity exec --bind /mnt/beegfs "${container}" \
    bcftools view \
      --threads "${threads}" \
      -r "${chr}" \
      "${filtered_vcf}" \
      -Oz -o "${out_vcf}"

  singularity exec --bind /mnt/beegfs "${container}" \
    bcftools index -t "${out_vcf}"

  singularity exec --bind /mnt/beegfs "${container}" \
    bcftools +counts "${out_vcf}" > "${counts_folder}/counts_${chr}.txt"
done

echo ""
echo "----------------------------------------"
echo "VALIDATION: Split VCFs"
date
echo "----------------------------------------"

total_split=0
filtered_total=$(singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H "${filtered_vcf}" | wc -l)

for i in {1..22}
do
  chr="chr${i}"
  split_vcf="${split_folder}/pacbio.merged.clean.gt_only.${chr}.vcf.gz"

  echo "Checking ${chr}..."

  [ -f "${split_vcf}" ] || { echo "FAIL: Missing ${chr}"; exit 1; }
  [ -f "${split_vcf}.tbi" ] || { echo "FAIL: No index ${chr}"; exit 1; }

  chrom_check=$(singularity exec --bind /mnt/beegfs "${container}" \
    bcftools query -f '%CHROM\n' "${split_vcf}" | sort -u)

  [ "${chrom_check}" = "${chr}" ] || { echo "FAIL: Wrong chrom in ${chr}"; exit 1; }

  split_format_count=$(singularity exec --bind /mnt/beegfs "${container}" \
    bcftools view -h "${split_vcf}" | grep '^##FORMAT=' | wc -l)

  split_gt_count=$(singularity exec --bind /mnt/beegfs "${container}" \
    bcftools view -h "${split_vcf}" | grep '^##FORMAT=<ID=GT,' | wc -l)

  [ "${split_gt_count}" -eq 1 ] || { echo "FAIL: GT missing ${chr}"; exit 1; }
  [ "${split_format_count}" -eq 1 ] || { echo "FAIL: Extra FORMAT fields ${chr}"; exit 1; }

  count=$(singularity exec --bind /mnt/beegfs "${container}" \
    bcftools view -H "${split_vcf}" | wc -l)

  echo "${chr}: ${count} variants"

  total_split=$((total_split + count))
done

echo "----------------------------------------"
echo "Filtered total: ${filtered_total}"
echo "Split total:    ${total_split}"

[ "${filtered_total}" -eq "${total_split}" ] || { echo "FAIL: Count mismatch"; exit 1; }

echo "PASS: Split matches original"
echo "----------------------------------------"

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"