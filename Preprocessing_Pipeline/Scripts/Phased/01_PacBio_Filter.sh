#!/bin/bash
# s01_PacBio_Filter
#
# Cameron Brown 10Apr2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s01_PacBio_Filter

# PBS directives
#---------------

#PBS -N s01_PacBio_Filter
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
# 
# 1. Edit ONLY the variables in the "USER INPUTS" section below.
#
# 2. Set:
#    - base_folder: your main project directory
#    - input_vcf_folder: folder containing PacBio input VCF files (.vcf.gz)
#    - container: full path to the Singularity/Apptainer container with bcftools
#    - min_qual: minimum QUAL threshold (e.g. 20 recommended)
#
# 3. The container definition (.def file) used to build the container is available on the project GitHub.
#
# 4. Do NOT modify anything below the "DO NOT EDIT" line unless you understand the pipeline.
#
# 5. Ensure input VCFs are bgzipped and indexed (.tbi files present).
#
# 6. Submit the script using:
#       qsub s01_PacBio_Filter
#
# Output:
# - Individually filtered PacBio VCFs
# - Merged filtered PacBio VCF
# - Raw and filtered bcftools stats
# - Basic validation printed to stdout
# ----------------------------

# ----------------------------
# USER INPUTS (EDIT THESE ONLY)
# ----------------------------

base_folder="/path/to/project"
input_vcf_folder="/path/to/pacbio/vcf_folder"
container="/path/to/container.sif"

# Programmable QUAL threshold
min_qual=20

# ----------------------------
# DO NOT EDIT BELOW
# ----------------------------

ln -s "$PWD" "$PBS_O_WORKDIR/$PBS_JOBID"

cd "$PBS_O_WORKDIR"

set -e

threads="${PBS_NCPUS:-${NCPUS:-1}}"

pipeline_folder="${base_folder}/data/processed/PacBio_Preprocessing"

clean_vcf_folder="${pipeline_folder}/clean_vcf"
merged_vcf_folder="${pipeline_folder}/merged_vcf"
stats_folder="${pipeline_folder}/stats"
raw_stats_folder="${stats_folder}/raw_stats"
clean_stats_folder="${stats_folder}/clean_stats"

final_vcf="${merged_vcf_folder}/pacbio.merged.filtered.vcf.gz"

autosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

mkdir -p \
  "${clean_vcf_folder}" \
  "${merged_vcf_folder}" \
  "${stats_folder}" \
  "${raw_stats_folder}" \
  "${clean_stats_folder}"

if [ ! -d "${input_vcf_folder}" ]; then
  echo "ERROR: Input VCF folder not found: ${input_vcf_folder}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

echo "----------------------------------------"
echo "STEP 1: PacBio VCF filter and merge"
date
echo "Input folder: ${input_vcf_folder}"
echo "Merged output: $(basename "${final_vcf}")"
echo "Threads used: ${threads}"
echo "QUAL threshold: ${min_qual}"
echo "----------------------------------------"
echo ""

cleaned_vcfs=""

for input_vcf in "${input_vcf_folder}"/*.vcf.gz
do
  if [ ! -f "${input_vcf}" ]; then
    echo "ERROR: No input VCF files found in ${input_vcf_folder}"
    exit 1
  fi

  base_name=$(basename "${input_vcf}")
  base_name=${base_name%.vcf.gz}
  base_name=${base_name%.bgzip}

  output_vcf="${clean_vcf_folder}/${base_name}.filtered.vcf.gz"
  raw_stats_txt="${raw_stats_folder}/${base_name}.raw.stats.txt"
  clean_stats_txt="${clean_stats_folder}/${base_name}.filtered.stats.txt"

  echo "----------------------------------------"
  echo "Processing: ${base_name}"
  echo "----------------------------------------"

  singularity exec --bind /mnt/beegfs "${container}" \
    bcftools stats "${input_vcf}" > "${raw_stats_txt}"

  echo "Raw stats written: $(basename "${raw_stats_txt}")"

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

  cleaned_vcfs="${cleaned_vcfs} ${output_vcf}"

  echo "Filtered VCF written: $(basename "${output_vcf}")"
  echo "Filtered stats written: $(basename "${clean_stats_txt}")"
  date
  echo ""
done

echo "----------------------------------------"
echo "STEP 2: Merge filtered PacBio VCFs"
date
echo "----------------------------------------"
echo ""

singularity exec --bind /mnt/beegfs "${container}" bash -c "
  bcftools merge \
    --threads ${threads} \
    --missing-to-ref \
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

echo "Merged VCF written: $(basename "${final_vcf}")"
date
echo ""

echo "----------------------------------------"
echo "STEP COMPLETE"
echo "Final file: $(basename "${final_vcf}")"
date
echo "----------------------------------------"
echo ""

echo "----------------------------------------"
echo "VALIDATION"
date
echo "----------------------------------------"

echo "Samples:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -l "${final_vcf}"

echo "----------------------------------------"

echo "Chromosomes present:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%CHROM\n' "${final_vcf}" | sort -u

echo "----------------------------------------"

echo "Indels:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -v indels "${final_vcf}" | wc -l

echo "----------------------------------------"

echo "Multiallelic sites:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -m3 "${final_vcf}" | wc -l

echo "----------------------------------------"

echo "FILTER distribution:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '%FILTER\n' "${final_vcf}" | sort | uniq -c

echo "----------------------------------------"

echo "Variants with QUAL < ${min_qual}:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H -i "QUAL<${min_qual}" "${final_vcf}" | wc -l

echo "----------------------------------------"

echo "Total variants:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools view -H "${final_vcf}" | wc -l

echo "----------------------------------------"

echo "Missing genotypes:"
singularity exec --bind /mnt/beegfs "${container}" \
  bcftools query -f '[%GT\n]' "${final_vcf}" | grep '\./\.' | wc -l

echo "----------------------------------------"
echo "Validation complete"
date
echo "----------------------------------------"
echo ""

rm -f "$PBS_O_WORKDIR/$PBS_JOBID"