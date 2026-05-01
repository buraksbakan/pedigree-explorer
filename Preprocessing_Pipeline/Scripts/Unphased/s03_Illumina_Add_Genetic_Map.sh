#!/bin/bash
# s03_Add_Genetic_Map_Illumina
#
# Cameron Brown 24Mar2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s03_Add_Genetic_Map_Illumina

# PBS directives
#---------------

#PBS -N s03_Add_Genetic_Map_Illumina
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
# 1. Run s01_Filter_Illumina and s02_Illumina_PLINK before running this script.
#
# 2. Edit ONLY the variables in the "USER INPUTS" section below.
#
# 3. Set:
#    - base_folder: your main project directory
#    - container: full path to the Singularity/Apptainer container
#    - perl_script: full path to add-map-plink.pl
#    - map_folder: full path to the genetic map files
#
# 4. A precomputed GRCh38 genetic recombination map is required.
#    Publicly available Beagle genetic maps (University of Washington) can be downloaded from:
#    https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/
#    The selected map_folder must match chromosome naming in the PLINK BIM file:
#    - use chr_in_chrom_field if chromosomes are named chr1, chr2, ...
#    - use no_chr_in_chrom_field if chromosomes are named 1, 2, ...
#    These maps are reformatted within the script and used to annotate the PLINK BIM file
#    with genetic distances (cM).
#
# 5. This script expects the PLINK files produced by s02_Illumina_PLINK:
#       ${base_folder}/data/processed/Illumina_Preprocessing/plink/final_plink_files/illumina_filtered.bed
#       ${base_folder}/data/processed/Illumina_Preprocessing/plink/final_plink_files/illumina_filtered.bim
#       ${base_folder}/data/processed/Illumina_Preprocessing/plink/final_plink_files/illumina_filtered.fam
# 
# 6. The container definition (.def file) and add-map-plink.pl script are available on the project GitHub.
# 
# 7. Do NOT modify anything below the "DO NOT EDIT" line unless you understand the pipeline.
# 
# 8. Submit the script using:
#       qsub s03_Add_Genetic_Map_Illumina
#
# Output:
# - Mapped PLINK .bed, .bim and .fam files (with genetic map added to BIM)
# - Mapped, sorted and indexed VCF (.vcf.gz + .csi)
# - Raw intermediate VCF (unsorted)
# - Basic validation printed to stdout
# ----------------------------
# ----------------------------
# USER INPUTS (EDIT THESE ONLY)
# ----------------------------

base_folder="/path/to/project"
container="/path/to/container.sif"
perl_script="/path/to/add-map-plink.pl"
map_folder="/path/to/hg38_Map/no_chr_in_chrom_field"

# ----------------------------
# DO NOT EDIT BELOW
# ----------------------------

ln -s "$PWD" "$PBS_O_WORKDIR/$PBS_JOBID"

# Change to working directory
cd "$PBS_O_WORKDIR"

# Calculate number of threads
threads="${PBS_NCPUS:-${NCPUS:-1}}"

# Stop at runtime errors
set -e

# Main Illumina preprocessing folder
pipeline_folder="${base_folder}/data/processed/Illumina_Preprocessing"

# Subfolders
plink_folder="${pipeline_folder}/plink"
plink_files_folder="${plink_folder}/final_plink_files"
mapped_plink_folder="${plink_folder}/mapped_plink_files"
mapped_plink_vcf_folder="${plink_folder}/mapped_plink_vcf"
genetic_map_folder="${pipeline_folder}/genetic_map"
converted_map_folder="${genetic_map_folder}/converted_for_perl"

# Prefixes
input_prefix="${plink_files_folder}/illumina_filtered"
output_prefix="${mapped_plink_folder}/illumina_filtered_mapped"

raw_mapped_vcf_prefix="${mapped_plink_vcf_folder}/illumina_filtered_mapped_raw"
raw_mapped_vcf="${raw_mapped_vcf_prefix}.vcf.gz"

sorted_mapped_vcf_prefix="${mapped_plink_vcf_folder}/illumina_filtered_mapped"
sorted_mapped_vcf="${sorted_mapped_vcf_prefix}.vcf.gz"

# Input files
input_bim="${input_prefix}.bim"
input_bed="${input_prefix}.bed"
input_fam="${input_prefix}.fam"

# Output files
output_bim="${output_prefix}.bim"
output_bed="${output_prefix}.bed"
output_fam="${output_prefix}.fam"

# Make output folders
mkdir -p \
  "${converted_map_folder}" \
  "${mapped_plink_folder}" \
  "${mapped_plink_vcf_folder}"

echo "----------------------------------------"
echo "STEP 1: Add genetic map to PLINK BIM"
date
echo "Input BIM: $(basename "${input_bim}")"
echo "Output BIM: $(basename "${output_bim}")"
echo "Perl script: $(basename "${perl_script}")"
echo "Threads used: ${threads}"
echo "----------------------------------------"
echo ""

# Check inputs exist
if [ ! -f "${input_bim}" ]; then
  echo "ERROR: Input BIM not found: ${input_bim}"
  echo "Run s02_Illumina_PLINK first, or check base_folder is correct."
  exit 1
fi

if [ ! -f "${input_bed}" ]; then
  echo "ERROR: Input BED not found: ${input_bed}"
  echo "Run s02_Illumina_PLINK first, or check base_folder is correct."
  exit 1
fi

if [ ! -f "${input_fam}" ]; then
  echo "ERROR: Input FAM not found: ${input_fam}"
  echo "Run s02_Illumina_PLINK first, or check base_folder is correct."
  exit 1
fi

if [ ! -f "${perl_script}" ]; then
  echo "ERROR: Perl script not found: ${perl_script}"
  exit 1
fi

if [ ! -f "${container}" ]; then
  echo "ERROR: Container not found: ${container}"
  exit 1
fi

if [ ! -d "${map_folder}" ]; then
  echo "ERROR: Genetic map folder not found: ${map_folder}"
  exit 1
fi

# Convert GRCh38 PLINK-format map files to the format expected by add-map-plink.pl
echo "Converting map files..."
for f in "${map_folder}"/*.map
do
  base_name=$(basename "${f}" .map)
  awk 'BEGIN{OFS="\t"} {print $1, $4, 0, $3}' "${f}" > "${converted_map_folder}/${base_name}.txt"
done
echo "Map conversion complete"
echo ""

# Add genetic map to BIM
echo "Running add-map-plink.pl..."
singularity exec --bind /mnt/beegfs "${container}" \
perl "${perl_script}" \
  -noheader \
  -no_set_zero \
  "${input_bim}" \
  "${converted_map_folder}"/*.txt \
  > "${output_bim}"

echo "Mapped BIM created"

# Copy BED and FAM
cp "${input_bed}" "${output_bed}"
cp "${input_fam}" "${output_fam}"

echo "BED and FAM copied"
echo ""

echo "----------------------------------------"
echo "STEP 2: Convert mapped PLINK back to VCF"
date
echo "----------------------------------------"
echo ""

# Export mapped PLINK to raw VCF
singularity exec --bind /mnt/beegfs "${container}" plink2 \
  --bfile "${output_prefix}" \
  --recode vcf bgz \
  --threads "${threads}" \
  --out "${raw_mapped_vcf_prefix}"

# Sort and index final VCF
singularity exec --bind /mnt/beegfs "${container}" bash -c "
  bcftools sort \
    -Oz -o ${sorted_mapped_vcf} ${raw_mapped_vcf} && \
  bcftools index -c ${sorted_mapped_vcf}
"

echo "Mapped VCF written: $(basename "${sorted_mapped_vcf}")"
echo ""

echo "----------------------------------------"
echo "VALIDATION"
date
echo "----------------------------------------"

echo "Checking mapped BED and FAM..."
if [ -f "${output_bed}" ] && [ -f "${output_fam}" ]; then
  echo "PASS: BED and FAM copied"
else
  echo "FAIL: BED and/or FAM missing"
  exit 1
fi
echo ""

echo "BIM line counts:"
echo "Original:"
wc -l "${input_bim}"
echo "Mapped:"
wc -l "${output_bim}"
echo ""

echo "First 5 lines of mapped BIM:"
head -5 "${output_bim}"
echo ""

echo "First 5 genetic distance values:"
awk '{print $3}' "${output_bim}" | head -5
echo ""

echo "Genetic map coverage:"
awk '{
  total++
  if ($3 == 0) zero++
}
END {
  if (total > 0)
    printf "Mapped: %.2f%% (%d/%d non-zero cM)\n", ((total-zero)/total)*100, (total-zero), total
  else
    print "No data"
}' "${output_bim}"
echo ""

echo "Checking sorted mapped VCF CSI index..."
if [ -f "${sorted_mapped_vcf}.csi" ]; then
  echo "PASS: Sorted mapped VCF CSI index found"
else
  echo "FAIL: Sorted mapped VCF CSI index missing"
  exit 1
fi
echo ""

echo "----------------------------------------"
echo "STEP COMPLETE: Genetic map added"
date
echo "----------------------------------------"
echo ""

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"