#!/bin/bash
# s01_Illumina_Clean_v4
#
# Cameron Brown 2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s01_Illumina_Clean_v4

# PBS directives
#---------------

#PBS -N s01_illumina_clean_v4
#PBS -l nodes=1:ncpus=16
#PBS -l walltime=01:00:00
#PBS -q one_hour
#PBS -m abe
#PBS -M cameron.brown.944@cranfield.ac.uk
#PBS -j oe
#PBS -v "CUDA_VISIBLE_DEVICES="
#PBS -W sandbox=PRIVATE
#PBS -k n

ln -s "$PWD" "$PBS_O_WORKDIR/$PBS_JOBID"

# Change to working directory
cd "$PBS_O_WORKDIR"

# Calculate number of threads
threads="${PBS_NCPUS:-${NCPUS:-1}}"

# Stop at runtime errors
set -e

# Folders and files
base_folder="/mnt/beegfs/project/Alexey_Larionov/IBD-2026"

# Main Illumina preprocessing folder
pipeline_folder="${base_folder}/data/processed/Illumina_Preprocessing"

# Subfolders
clean_vcf_folder="${pipeline_folder}/clean_vcf"
plink_folder="${pipeline_folder}/plink"
plink_files_folder="${plink_folder}/final_plink_files"
plink_vcf_folder="${plink_folder}/plink_converted_vcf"
stats_folder="${pipeline_folder}/stats"
raw_stats_folder="${stats_folder}/raw_stats"
clean_stats_folder="${stats_folder}/clean_stats"

# Inputs
container="${base_folder}/github_repo/pedigree-explorer/pipeline/config/Simple_Container.sif"
input_vcf="${base_folder}/data/raw/illumina/CEPH1463.GRCh38.illumina-dragen.oa.vcf.gz"

# Outputs
output_vcf="${clean_vcf_folder}/illumina.clean.vcf.gz"
raw_stats_txt="${raw_stats_folder}/illumina.raw.stats.txt"
clean_stats_txt="${clean_stats_folder}/illumina.clean.stats.txt"

# Autosomes only
autosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

mkdir -p \
  "${clean_vcf_folder}" \
  "${plink_files_folder}" \
  "${plink_vcf_folder}" \
  "${raw_stats_folder}" \
  "${clean_stats_folder}"

# Check inputs exist
if [ ! -f "${input_vcf}" ]; then
  echo "ERROR: Input VCF not found: ${input_vcf}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

echo "----------------------------------------"
echo "STEP: Illumina VCF cleaning"
date
echo "Input VCF: $(basename "${input_vcf}")"
echo "Output VCF: $(basename "${output_vcf}")"
echo "Threads used: ${threads}"
echo "----------------------------------------"
echo ""

echo "Generating raw stats..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools stats "${input_vcf}" > "${raw_stats_txt}"
echo "Raw stats written: $(basename "${raw_stats_txt}")"
echo ""

# Filter to:
# - PASS variants only
# - biallelic sites only
# - SNPs only
# - autosomes only
# - QUAL >= 20
# - sorted compressed output
# - indexed final VCF
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
    -i 'QUAL>=20' \
    -Ou | \
  bcftools sort \
    -Oz -o ${output_vcf} && \
  bcftools index -t ${output_vcf} && \
  bcftools stats ${output_vcf} > ${clean_stats_txt}
"

echo "Filtering, sorting and indexing complete"
echo "Clean stats written: $(basename "${clean_stats_txt}")"
date
echo ""

echo "----------------------------------------"
echo "STEP COMPLETE: Illumina VCF cleaned"
echo "Output file: $(basename "${output_vcf}")"
date
echo "----------------------------------------"
echo ""

echo "----------------------------------------"
echo "VALIDATION"
date
echo "----------------------------------------"

echo "Indels:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -v indels "${output_vcf}" | wc -l

echo "----------------------------------------"

echo "Multiallelic sites:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -m3 "${output_vcf}" | wc -l

echo "----------------------------------------"

echo "Chromosomes present:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%CHROM\n' "${output_vcf}" | sort -u

echo "----------------------------------------"

echo "FILTER distribution before:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%FILTER\n' "${input_vcf}" | sort | uniq -c

echo "----------------------------------------"

echo "FILTER distribution after:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%FILTER\n' "${output_vcf}" | sort | uniq -c

echo "----------------------------------------"

echo "Variants with QUAL < 20:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -i 'QUAL<20' "${output_vcf}" | wc -l

echo "----------------------------------------"
echo "Validation complete"
date
echo "----------------------------------------"
echo ""

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"