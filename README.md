# Pangenome Variations

A snakemake pipeline that receives as input several genomes and build pangenome SNPs between them, serving as a set of
truth variants to be found by variant callers.

It does so by comparing all the input genomes pairwisely, and thus obtaining SNPs, then deduplicating these SNPs based
on the coordinates of their alleles. In the end, we have, for each genome, truth probes from deduplicated SNPs, with
several data.

A thorough explanation of this pipeline will be done soon (either in the paper or here).

The version used in the pandora paper has tag `pandora_paper_tag1`.

# Running

## Requirements

### Dependencies

* python 3.6+;
* singularity;


### Setting up virtualenv
```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Running on the sample example:
```
unzip sample_data.zip
snakemake -j8 --use-singularity
```

## Running on the paper data:

1. `git checkout pandora_paper_tag1`

### If you want to run locally:
2. `bash scripts/run_pipeline_local.sh -j <NB_OF_THREADS> --configfile config.pandora_paper_tag1.yaml`

### If you want to run on an LSF cluster:
2. `bash scripts/submit_lsf.sh --configfile config.pandora_paper_tag1.yaml`

# TODO

Improve this README
