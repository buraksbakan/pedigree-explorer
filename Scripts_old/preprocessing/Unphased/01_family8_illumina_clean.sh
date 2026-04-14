#!/bin/bash
# 01_family8_illumina_clean
#
# Cameron Brown 30Mar2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub 01_family8_illumina_clean

# PBS directives
#---------------

#PBS -N 01_family8_illumina_clean
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
container="${base_folder}/containers/plink.sif"
input_vcf="${base_folder}/data/raw/family_8/pb/2129547_Y430.GCA_000001405.deepvariant.phased.bgzip.vcf.gz"
output_folder="${base_folder}/data/raw/family_8/pb/clean"

output_vcf="${output_folder}/2129544_Y430.GCA_000001405.deepvariant.phased.clean.vcf.gz" 

# Autosomes only
autosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

# Make output folder if needed
mkdir -p "${output_folder}"

# Check inputs exist
if [ ! -f "${input_vcf}" ]; then
  echo "ERROR: Input VCF not found: ${input_vcf}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

# Start message
echo "----------------------------------------"
echo "STEP: Illumina VCF cleaning started"
date
echo "Input VCF: ${input_vcf}"
echo "Output VCF: ${output_vcf}"
echo "Container: ${container}"
echo "PBS_NODEFILE: ${PBS_NODEFILE}"
echo "PBS_NCPUS: ${PBS_NCPUS}"
echo "NCPUS: ${NCPUS}"
echo "Threads used: ${threads}"
echo "----------------------------------------"
echo ""

# Filter to:
# - PASS variants only (-f PASS)
# - biallelic sites only (-m2 -M2)
# - SNPs only (-v snps)
# - autosomes only (chr1–chr22 via -r)
# - ensure coordinate-sorted output (bcftools sort)
# - output as compressed VCF (.vcf.gz, -Oz)
# - index the final file for downstream tools
singularity exec --bind /mnt/beegfs "${container}" bash -c "
  bcftools view \
    --threads ${threads} \
    -m2 -M2 \
    -f PASS \
    -v snps \
    -r ${autosomes} \
    -Ou ${input_vcf} | \
  bcftools sort \
    -Oz -o ${output_vcf} && \
  bcftools index -t ${output_vcf}
"

echo ""
echo "Filtering, sorting and indexing complete"
date
echo ""

# Final confirmation
echo "----------------------------------------"
echo "STEP COMPLETE: Illumina VCF cleaned successfully"
echo "Output file:"
echo "${output_vcf}"
echo "----------------------------------------"
date
echo ""

echo "========================================"
echo "VALIDATION: Checking filtered VCF"
date
echo "========================================"
echo ""

echo "Checking for indels (should be 0)..."
singularity exec --bind /mnt/beegfs "${container}" \
bcftools view -H -v indels "${output_vcf}" | wc -l

echo ""

echo "Checking for multiallelic sites (should be 0)..."
singularity exec --bind /mnt/beegfs "${container}" \
bcftools view -H -m3 "${output_vcf}" | wc -l

echo ""

echo "Checking chromosomes present (should be chr1 to chr22 only)..."
singularity exec --bind /mnt/beegfs "${container}" \
bcftools query -f '%CHROM\n' "${output_vcf}" | sort -u

echo ""

echo "Checking FILTER distribution Before ..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%FILTER\n' "${input_vcf}" | sort | uniq -c
  
echo "Checking FILTER distribution After ..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%FILTER\n' "${output_vcf}" | sort | uniq -c

echo "========================================"
echo "VALIDATION COMPLETE"
date
echo "========================================"
echo ""

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"