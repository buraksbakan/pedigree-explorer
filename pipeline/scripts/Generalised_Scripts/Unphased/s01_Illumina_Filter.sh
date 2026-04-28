#!/bin/bash
# s01_Filter_Illumina
#
# Cameron Brown 2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s01_Filter_Illumina

# PBS directives
#---------------

#PBS -N s01_Filter_Illumina
#PBS -l nodes=1:ncpus=16
#PBS -l walltime=01:00:00
#PBS -q one_hour
#PBS -m abe
#PBS -M your.email@yourinstitution.ac.uk
#PBS -j oe
#PBS -v "CUDA_VISIBLE_DEVICES="
#PBS -W sandbox=PRIVATE
#PBS -k n

# ----------------------------
# INSTRUCTIONS
# ----------------------------
# 1. Edit ONLY the variables in the "USER INPUTS" section below.
# 2. Set:
#    - base_folder: your main project directory
#    - input_vcf: full path to your input VCF (.vcf.gz)
#    - container: full path to the Singularity/Apptainer container with bcftools
#    - min_qual: minimum QUAL threshold (e.g. 20 recomended)
# 3. The container definition (.def file) used to build the container is available on the project GitHub.
# 4. Do NOT modify anything below the "DO NOT EDIT" line unless you understand the pipeline.
# 5. Ensure the input VCF is bgzipped and indexed (.tbi file present).
# 6. Submit the script using:
#       qsub s01_Filter_Illumina
#
# Output:
# - Filtered VCF
# - Raw and filtered bcftools stats
# - Basic validation printed to stdout
# ----------------------------

# ----------------------------
# USER INPUTS (EDIT THESE ONLY)
# ----------------------------

base_folder="/path/to/project"
input_vcf="/path/to/input/Illumina.vcf.gz"
container="/path/to/container.sif"

# Programmable QUAL threshold
min_qual=20

# ----------------------------
# DO NOT EDIT BELOW
# ----------------------------

ln -s "$PWD" "$PBS_O_WORKDIR/$PBS_JOBID"

cd "$PBS_O_WORKDIR"

threads="${PBS_NCPUS:-${NCPUS:-1}}"

set -e

pipeline_folder="${base_folder}/data/processed/Illumina_Preprocessing"

clean_vcf_folder="${pipeline_folder}/clean_vcf"
plink_folder="${pipeline_folder}/plink"
plink_files_folder="${plink_folder}/final_plink_files"
stats_folder="${pipeline_folder}/stats"
raw_stats_folder="${stats_folder}/raw_stats"
clean_stats_folder="${stats_folder}/clean_stats"

output_vcf="${clean_vcf_folder}/illumina_filtered.vcf.gz"
raw_stats_txt="${raw_stats_folder}/illumina.raw.stats.txt"
clean_stats_txt="${clean_stats_folder}/illumina.filtered.stats.txt"

autosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

mkdir -p \
  "${clean_vcf_folder}" \
  "${plink_files_folder}" \
  "${raw_stats_folder}" \
  "${clean_stats_folder}"

if [ ! -f "${input_vcf}" ]; then
  echo "ERROR: Input VCF not found: ${input_vcf}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

echo "----------------------------------------"
echo "STEP: Illumina VCF filtering"
date
echo "Input VCF: $(basename "${input_vcf}")"
echo "Output VCF: $(basename "${output_vcf}")"
echo "Threads used: ${threads}"
echo "QUAL threshold: ${min_qual}"
echo "----------------------------------------"
echo ""

echo "Generating raw stats..."
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools stats "${input_vcf}" > "${raw_stats_txt}"
echo "Raw stats written: $(basename "${raw_stats_txt}")"
echo ""

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
    -i \"QUAL>=${min_qual}\" \
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
echo "STEP COMPLETE: Illumina VCF filtered"
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

echo "Variants with QUAL < ${min_qual}:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -i "QUAL<${min_qual}" "${output_vcf}" | wc -l

echo "----------------------------------------"
echo "Validation complete"
date
echo "----------------------------------------"
echo ""

rm -f "$PBS_O_WORKDIR/$PBS_JOBID"