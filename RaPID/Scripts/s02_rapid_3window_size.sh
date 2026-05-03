#!/bin/bash
#
# s02_rapid_3window_size.sh
#
# Runs RaPID IBD detection across multiple window sizes on
# per-chromosome multi-sample VCFs. Merges per-chromosome results into a
# single file per window size.
# (all paths and sample identifiers are placeholders)
#
#PBS -N rapid_multi_window
#PBS -l nodes=1:ncpus=16
#PBS -l walltime=06:00:00
#PBS -q six_hour
#PBS -j oe

cd "$PBS_O_WORKDIR"
set -e

module use /apps/modules/all
module load RaPID/RaPID-1.7

echo "Starting RaPID"
date
echo ""

# =====================
# USER-DEFINED INPUTS
# =====================
BASE_DIR="/path/to/project_root"
VCF_DIR="/path/to/vcf/split"       # directory with per-chromosome multi-sample VCFs
VCF_PREFIX="sample.chr"            # filename prefix before chrN
VCF_SUFFIX=".vcf.gz"
MAP_DIR="${BASE_DIR}/interpolated_maps"    # directory with interpolated map files from s01_map_interpolating.sh
OUTPUT_BASE="${BASE_DIR}/results"

# =====================
# RAPID PARAMETERS
# =====================
NUM_RUNS=10
NUM_SUCCESS=2
MIN_CM=5

# =====================
# RUN RAPID
# =====================
for WINDOW in 75 250 500; do
    MERGED_DIR="${OUTPUT_BASE}/rapid_w${WINDOW}_merged"
    mkdir -p "${MERGED_DIR}"

    for CHR in {1..22}; do
        CHR_DIR="${OUTPUT_BASE}/rapid_w${WINDOW}/chr${CHR}"
        mkdir -p "${CHR_DIR}"

        echo "  chr${CHR} w=${WINDOW}..."
        RaPID_v.1.7 \
            -i "${VCF_DIR}/${VCF_PREFIX}${CHR}${VCF_SUFFIX}" \
            -g "${MAP_DIR}/hg38_interpolated_map_chr${CHR}.txt" \
            -d "${MIN_CM}" \
            -w "${WINDOW}" \
            -r "${NUM_RUNS}" \
            -s "${NUM_SUCCESS}" \
            -o "${CHR_DIR}/rapid_chr${CHR}"
    done

    # =====================
    # MERGE RESULTS
    # =====================
    find "${OUTPUT_BASE}/rapid_w${WINDOW}" -name "results.max.gz" | sort -V | \
        xargs zcat > "${MERGED_DIR}/rapid_w${WINDOW}_merged"

    echo "  w=${WINDOW} done"
done

echo ""
echo "RaPID completed"
date
