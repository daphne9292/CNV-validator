#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="$ROOT_DIR/demo/data"
PYTHON="${PYTHON:-python}"
BENCH_BED="${BENCH_BED:-$DATA_DIR/chr22.bed}"
FASTA="${FASTA:-$DATA_DIR/chr22.fa}"
TUMOR_BAM="${TUMOR_BAM:-$DATA_DIR/WGS_IL_T_1.chr22.bwa.dedup.bam}"
NORMAL_BAM="${NORMAL_BAM:-$DATA_DIR/WGS_IL_N_1.chr22.bwa.dedup.bam}"
SNP_VCF_PREFIX="${SNP_VCF_PREFIX:-$DATA_DIR/varscan/HCC1395_WGS_IL_T_1.varscan.}"
RP_DIR="${RP_DIR:-$ROOT_DIR/demo/RP_and_SR_proof/}"
RD_DIR="${RD_DIR:-$ROOT_DIR/demo/depth_and_BAF_proof/}"
FEATURE_CSV="${FEATURE_CSV:-$ROOT_DIR/demo/RP_SR_depth_BAF_proofStatistic.csv}"
PREDICTION_CSV="${PREDICTION_CSV:-$ROOT_DIR/demo/model.prediction.csv}"
MODEL="${MODEL:-$ROOT_DIR/model/catboost_model.pkl}"
SAMPLE_NAME="${SAMPLE_NAME:-HCC1395_WGS_IL_T_1}"
MEAN_INSERT_SIZE="${MEAN_INSERT_SIZE:-413}"
STD_INSERT_SIZE="${STD_INSERT_SIZE:-106}"
THREAD_NUM="${THREAD_NUM:-1}"
EXPAND_SIZE="${EXPAND_SIZE:-500}"

for input in "$BENCH_BED" "$FASTA" "$TUMOR_BAM" "$NORMAL_BAM" "$MODEL"; do
    if [[ ! -f "$input" ]]; then
        printf 'Missing input: %s\n' "$input" >&2
        exit 1
    fi
done
mkdir -p -- "$RP_DIR" "$RD_DIR" "$(dirname -- "$FEATURE_CSV")" "$(dirname -- "$PREDICTION_CSV")"

common_args=(
    --tumor_bam "$TUMOR_BAM"
    --normal_bam "$NORMAL_BAM"
    --fasta "$FASTA"
    --bed "$BENCH_BED"
    --sample_name "$SAMPLE_NAME"
    --RP_dir "$RP_DIR"
    --RD_dir "$RD_DIR"
    --mean_insert_size "$MEAN_INSERT_SIZE"
    --std_insert_size "$STD_INSERT_SIZE"
    --expand_size "$EXPAND_SIZE"
    --snp_vcf_prefix "$SNP_VCF_PREFIX"
    --thread_num "$THREAD_NUM"
)

"$PYTHON" "$ROOT_DIR/extract.py" "${common_args[@]}" --output "$FEATURE_CSV"
# validator.py reuses the proof files created by extract.py.
"$PYTHON" "$ROOT_DIR/validator.py" "${common_args[@]}" \
    --output_feature "$FEATURE_CSV" \
    --model "$MODEL" \
    --output_model_predict "$PREDICTION_CSV"

printf 'Features: %s\nPredictions: %s\n' "$FEATURE_CSV" "$PREDICTION_CSV"
