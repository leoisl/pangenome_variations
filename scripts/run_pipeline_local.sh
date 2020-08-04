#!/usr/bin/env bash
set -eux

CORES=112
LOG_DIR=logs/

mkdir -p "$LOG_DIR"
snakemake -j "$CORES" --use-singularity >"$LOG_DIR"/pipeline.out 2>"$LOG_DIR"/pipeline.err
exit 0
