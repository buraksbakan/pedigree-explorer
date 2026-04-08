#!/bin/bash
# s01_PacBio_CleanMerge_v3
#
# Cameron Brown 30Mar2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s01_PacBio_CleanMerge_v3

# PBS directives
#---------------

#PBS -N s01_pacbio_cleanmerge_v3
#PBS -l nodes=1:ncpus=16
#PBS -l walltime=01:00:00
#PBS -q one_hour
#PBS -m abe
#PBS -M cameron.brown.944@cranfield.ac.uk
#PBS -j oe
#PBS -v "CUDA_VISIBLE_DEVICES="
#PBS -W sandbox=PRIVATE
#PBS -k n

ln -sf "$PWD" "$PBS_O_WORKDIR/$PBS_JOBID"

# Change to working directory
cd "$PBS_O_WORKDIR"

# Calculate number of threads
threads="${PBS_NCPUS:-${NCPUS:-1}}"

# Stop at runtime errors
set -e

# Folders and files
base_folder="/mnt/beegfs/project/Alexey_Larionov/IBD-2026"

# Main PacBio preprocessing folder
pipeline_folder="${base_folder}/data/processed/PacBio_Preprocessing"

# Subfolders
input_vcf_folder="${base_folder}/data/raw/family_8/pb"
clean_vcf_folder="${pipeline_folder}/clean_vcf"
merged_vcf_folder="${pipeline_folder}/merged_vcf"
log_folder="${pipeline_folder}/logs"

# Inputs
container="${base_folder}/github_repo/pedigree-explorer/pipeline/config/Simple_Container.sif"

# Outputs
final_vcf="${merged_vcf_folder}/pacbio.merged.clean.biallelic.nomiss.vcf.gz"
log_file="${log_folder}/s01_pacbio_cleanmerge_$(date +%Y%m%d_%H%M%S).log"

# Autosomes only
autosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

mkdir -p \
  "${clean_vcf_folder}" \
  "${merged_vcf_folder}" \
  "${log_folder}"

# Log to file and screen
exec > >(tee -i "${log_file}")
exec 2>&1

# Check inputs exist
if [ ! -d "${input_vcf_folder}" ]; then
  echo "ERROR: Input VCF folder not found: ${input_vcf_folder}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

# Start message
echo "----------------------------------------"
echo "STEP: PacBio VCF clean and merge started"
date
echo "Input folder: ${input_vcf_folder}"
echo "Final merged VCF: ${final_vcf}"
echo "Container: ${container}"
echo "PBS_NODEFILE: ${PBS_NODEFILE}"
echo "PBS_NCPUS: ${PBS_NCPUS}"
echo "NCPUS: ${NCPUS}"
echo "Threads used: ${threads}"
echo "----------------------------------------"
echo ""

# Clean each PacBio VCF
cleaned_vcfs=""

for input_vcf in "${input_vcf_folder}"/*.vcf.gz
do
  if [ ! -f "${input_vcf}" ]; then
    echo "ERROR: No input VCF files found in ${input_vcf_folder}"
    exit 1
  fi

  base_name=$(basename "${input_vcf}" .vcf.gz)
  output_vcf="${clean_vcf_folder}/${base_name}.clean.vcf.gz"

  echo "----------------------------------------"
  echo "Cleaning VCF: ${input_vcf}"
  echo "Output VCF: ${output_vcf}"
  echo "----------------------------------------"

  singularity exec --bind /mnt/beegfs "${container}" bash -c "
    bcftools view \
      --threads ${threads} \
      -m2 -M2 \
      -f PASS \
      -v snps \
      -r ${autosomes} \
      -Ou ${input_vcf} | \
    bcftools filter \
      --threads ${threads} \
      -i 'QUAL>=10' \
      -Ou | \
    bcftools sort \
      -Oz -o ${output_vcf} && \
    bcftools index -t ${output_vcf}
  "

  cleaned_vcfs="${cleaned_vcfs} ${output_vcf}"

  echo ""
  echo "Finished cleaning and indexing: ${output_vcf}"
  date
  echo ""
done

echo "----------------------------------------"
echo "STEP: Merge cleaned PacBio VCFs and apply post-merge filtering"
date
echo "----------------------------------------"
echo ""

singularity exec --bind /mnt/beegfs "${container}" bash -c "
  bcftools merge \
    --threads ${threads} \
    -Ou \
    ${cleaned_vcfs} | \
  bcftools view \
    --threads ${threads} \
    -g ^miss \
    -m2 -M2 \
    -v snps \
    -Oz -o ${final_vcf} && \
  bcftools index -t ${final_vcf}
"

echo ""
echo "Merge and post-merge filtering complete"
echo "Final merged VCF:"
echo "${final_vcf}"
date
echo ""

# Final confirmation
echo "----------------------------------------"
echo "STEP COMPLETE: PacBio VCFs cleaned, merged, and filtered successfully"
echo "Final merged biallelic no-missing file:"
echo "${final_vcf}"
echo "----------------------------------------"
date
echo ""

echo "========================================"
echo "VALIDATION: Checking final merged VCF"
date
echo "========================================"
echo ""

echo "Checking sample names..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -l "${final_vcf}"

echo ""

echo "Checking chromosomes present (should be chr1 to chr22 only)..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%CHROM\n' "${final_vcf}" | sort -u

echo ""

echo "Checking for indels (should be 0)..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -v indels "${final_vcf}" | wc -l

echo ""

echo "Checking for multiallelic sites (should be 0)..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -m3 "${final_vcf}" | wc -l

echo ""

echo "Checking FILTER distribution..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%FILTER\n' "${final_vcf}" | sort | uniq -c

echo ""

echo "Checking number of QUAL < 10..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -i 'QUAL<20' "${final_vcf}" | wc -l

echo ""

echo "Checking number of variants in final merged VCF..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H "${final_vcf}" | wc -l

echo ""

echo "Checking missing genotypes (should be 0)..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '[%GT\n]' "${final_vcf}" | grep '\./\.' | wc -l

echo ""

echo "========================================"
echo "VALIDATION COMPLETE"
date
echo "========================================"
echo ""

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"