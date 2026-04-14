#!/bin/bash
# s03_Add_Genetic_Map_family8
#
# Cameron Brown 30Mar2026

# Crescent2 script
# Note: this script should be run on a compute node
# qsub s03_Add_Genetic_Map_family8

# PBS directives
#---------------

#PBS -N s03_add_genetic_map_family8
#PBS -l nodes=1:ncpus=16
#PBS -l walltime=00:30:00
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

# Load modules
module purge
module use /apps2/modules/all

# Load Perl
module load Perl

# Stop at runtime errors
set -e

# Folders and files
base_folder="/mnt/beegfs/project/Alexey_Larionov/IBD-2026"
perl_script="${base_folder}/github_repo/pedigree-explorer/Scripts/preprocessing/Perl/add-map-plink.pl"

input_prefix="${base_folder}/data/processed/illumina/plink/family_8/illumina_clean_family8"
output_prefix="${base_folder}/data/processed/illumina/plink/family_8/illumina_mapped_clean_family8"

map_folder="${base_folder}/hg38_Map/no_chr_in_chrom_field"
converted_map_folder="${base_folder}/hg38_Map/converted_for_perl"

# Input files
input_bim="${input_prefix}.bim"
input_bed="${input_prefix}.bed"
input_fam="${input_prefix}.fam"

# Output files
output_bim="${output_prefix}.bim"
output_bed="${output_prefix}.bed"
output_fam="${output_prefix}.fam"

# Make output folders
mkdir -p "${converted_map_folder}"

echo "----------------------------------------"
echo "STEP: Add genetic map to PLINK BIM"
date
echo "Input BIM: ${input_bim}"
echo "Output BIM: ${output_bim}"
echo "Perl script: ${perl_script}"
echo "Map folder: ${map_folder}"
echo "----------------------------------------"
echo ""

# Check inputs exist
if [ ! -f "${input_bim}" ]; then
  echo "ERROR: Input BIM not found: ${input_bim}"
  exit 1
fi

if [ ! -f "${input_bed}" ]; then
  echo "ERROR: Input BED not found: ${input_bed}"
  exit 1
fi

if [ ! -f "${input_fam}" ]; then
  echo "ERROR: Input FAM not found: ${input_fam}"
  exit 1
fi

if [ ! -f "${perl_script}" ]; then
  echo "ERROR: Perl script not found: ${perl_script}"
  exit 1
fi

# Convert GRCh38 PLINK-format map files to the format expected by add-map-plink.pl
echo "Converting map files..."
for f in "${map_folder}"/*.map; do
  base_name=$(basename "${f}" .map)
  awk 'BEGIN{OFS="\t"} {print $1, $4, 0, $3}' "${f}" > "${converted_map_folder}/${base_name}.txt"
done
echo "Map conversion complete"
echo ""

# Add genetic map to BIM
echo "Running add-map-plink.pl..."
perl "${perl_script}" \
  "${input_bim}" \
  "${converted_map_folder}"/*.txt \
  > "${output_bim}"
echo "Mapped BIM created"
echo ""

# Copy BED and FAM
cp "${input_bed}" "${output_bed}"
cp "${input_fam}" "${output_fam}"

echo "BED and FAM copied"
echo ""

echo "----------------------------------------"
echo "VALIDATION"
echo "----------------------------------------"
echo ""

echo "Original BIM line count:"
wc -l "${input_bim}"

echo "Mapped BIM line count:"
wc -l "${output_bim}"
echo ""

echo "First 5 lines of mapped BIM:"
head -5 "${output_bim}"
echo ""

echo "First 5 genetic distance values:"
awk '{print $3}' "${output_bim}" | head -5
echo ""

echo "Number of non-zero cM entries:"
awk '{
  total++
  if ($3 == 0) zero++
}
END {
  if (total>0)
    printf "Genetic map coverage: %.2f%% mapped (%d/%d non-zero cM)\n", ((total-zero)/total)*100, (total-zero), total
  else
    print "No data"
}' "${output_bim}"
echo ""

echo "----------------------------------------"
echo "STEP COMPLETE: Genetic map added"
date
echo "----------------------------------------"
echo ""

# Clean-up
rm -f "$PBS_O_WORKDIR/$PBS_JOBID"