#!/usr/bin/env bash
set -eux

CORES=112
LOG_DIR=logs/

mkdir -p "$LOG_DIR"
snakemake --configfile config.pandora_paper_tag1.yaml -j "$CORES" --use-singularity >"$LOG_DIR"/pipeline.out 2>"$LOG_DIR"/pipeline.err
exit 0
