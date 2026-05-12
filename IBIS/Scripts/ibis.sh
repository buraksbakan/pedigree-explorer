#!/bin/bash
# ibis.sh
#
# script for running IBIS on unphased Illumina PLINK data
# (all paths and sample identifiers are placeholders)

#PBS -N run_ibis
#PBS -l nodes=1:ncpus=8
#PBS -l walltime=01:00:00
#PBS -q one_hour
#PBS -j oe

cd "$PBS_O_WORKDIR"
set -euo pipefail

module use /apps/modules/all
module load ibis
module load bedtools

echo "Starting IBIS analysis"
date
echo ""

# =====================
# USER-DEFINED INPUTS
# =====================

BASE_DIR="/path/to/project_root"

PLINK_INPUT_1="${BASE_DIR}/data/plink_dataset_1/input_prefix"
PLINK_INPUT_2="${BASE_DIR}/data/plink_dataset_2/input_prefix"

OUT_DIR_1="${BASE_DIR}/results/dataset_1/ibis"
OUT_DIR_2="${BASE_DIR}/results/dataset_2/ibis"

mkdir -p "${OUT_DIR_1}" "${OUT_DIR_2}"

# =====================
# IBIS PARAMETERS
# =====================

SENSITIVE_OPTS="-min_l 3 -mt 50 -errorRate 0.06"
SEIDMAN_OPTS="-min_l 7 -mt 400 -errorRate 0.004"
STRICT_OPTS="-min_l 10 -mt 200 -errorRate 0.20"

# =====================
# RUN IBIS
# =====================

ibis -bfile "${PLINK_INPUT_1}" ${SENSITIVE_OPTS} -o "${OUT_DIR_1}/lenient"
ibis -bfile "${PLINK_INPUT_1}" ${SEIDMAN_OPTS} -o "${OUT_DIR_1}/seidman"
ibis -bfile "${PLINK_INPUT_1}" ${STRICT_OPTS} -o "${OUT_DIR_1}/strict"

ibis -bfile "${PLINK_INPUT_2}" ${SENSITIVE_OPTS} -o "${OUT_DIR_2}/lenient"
ibis -bfile "${PLINK_INPUT_2}" ${SEIDMAN_OPTS} -o "${OUT_DIR_2}/seidman"
ibis -bfile "${PLINK_INPUT_2}" ${STRICT_OPTS} -o "${OUT_DIR_2}/strict"

# =====================
# CONVERT .seg to BED
# =====================

for dir in "${OUT_DIR_1}" "${OUT_DIR_2}"; do
    for seg in "${dir}"/*.seg; do
        out="${seg%.seg}.bed"
        awk 'BEGIN{OFS="\t"} {
            a=$1; b=$2
            print "chr"$3, $4-1, $5, a"-"b
        }' "${seg}" | sort -k1,1 -k2,2n > "${out}"
    done
done

echo "IBIS analysis completed"
date
