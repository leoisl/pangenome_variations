#!/usr/bin/env bash
set -eux

LOG_DIR=logs/
mkdir -p "$LOG_DIR"

snakemake --use-singularity >"$LOG_DIR"/pipeline.out 2>"$LOG_DIR"/pipeline.err
exit 0
