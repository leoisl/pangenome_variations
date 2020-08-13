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
You need `singularity` and `python 3.6+`.

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

Be sure you are on tag `pandora_paper_tag1`

If you are on a *LSF* cluster, run:
`bash scripts/submit_lsf.sh`

Otherwise, you can run locally as:
`bash scripts/run_pipeline_local.sh`

# TODO
Improve this README