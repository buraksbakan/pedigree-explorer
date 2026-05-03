#!/bin/bash
#
# s01_map_interpolating.sh
#
# Interpolates hg38 genetic map positions for each SNP in the VCF
# using interpolate_loci.py. Produces per-chromosome interpolated
# map files required as input for RaPID.
# (all paths and sample identifiers are placeholders)
#
#PBS -N map_interpolating
#PBS -l nodes=1:ncpus=12
#PBS -l walltime=06:00:00
#PBS -q six_hour
#PBS -j oe

cd "$PBS_O_WORKDIR"
set -e

source ~/.bashrc
conda activate rapid_env

echo "Starting interpolation"
date
echo ""

# =====================
# USER-DEFINED INPUTS
# =====================
BASE_DIR="/path/to/project_root"
SCRIPT_DIR="${BASE_DIR}/mapping_scripts"  # directory containing interpolate_loci.py script provided by RaPID authors
MAP_DIR="${BASE_DIR}/genetic_map"     # directory with original genetic map files provided by RaPID authors
VCF_DIR="/path/to/vcf/split"       # directory with per-chromosome multi-sample VCFs
VCF_PREFIX="sample.chr"            # filename prefix before chrN
VCF_SUFFIX=".vcf.gz"
INTERPOLATED_MAP_DIR="${BASE_DIR}/interpolated_maps"

# =====================
# INTERPOLATE PER CHR
# =====================
for CHR in {1..22}; do
    echo "  Interpolating chr${CHR}..."
    python3 "${SCRIPT_DIR}/interpolate_loci.py" \
        "${MAP_DIR}/genetic_map_GRCh38_chr${CHR}.txt" \
        "${VCF_DIR}/${VCF_PREFIX}${CHR}${VCF_SUFFIX}" \
        "${INTERPOLATED_MAP_DIR}/hg38_interpolated_map_chr${CHR}.txt"
done

echo ""
echo "Interpolation completed"
date
