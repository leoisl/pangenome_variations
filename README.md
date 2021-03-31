# Pangenome Variations

A snakemake pipeline that receives as input several genomes and builds pangenome variants between them, serving as a set of
truth variants to be found by variant callers.

It does so by comparing all the input genomes pairwisely, and thus obtaining variants, then deduplicating these variants based
on the coordinates of their alleles. In the end, we have, for each genome, truth probes from deduplicated variants, with
several data.

A thorough explanation of this pipeline can be found in [our paper](https://www.biorxiv.org/content/10.1101/2020.11.12.380378v3).

The version used in the pandora paper review has tag `pandora_paper_update_31_03_2021`.

The version used in the pandora paper submission has tag `pandora_paper_tag1`.

# Running

## Requirements

### Dependencies

* python 3.6+;
* singularity 2.4.1+;


### Setting up virtualenv
`./setup.sh`

## Running on the sample example:
```
unzip sample_data.zip
source venv/bin/activate
bash scripts/run_pipeline_local.sh -j8
```

## Running on the paper data:

1. `git checkout pandora_paper_tag1`
2. `source venv/bin/activate`

### If you want to run locally:
3. `bash scripts/run_pipeline_local.sh -j <NB_OF_THREADS> --configfile config.pandora_paper_tag1.yaml`

### If you want to run on an LSF cluster:
3. `bash scripts/submit_lsf.sh --configfile config.pandora_paper_tag1.yaml`

# Troubleshooting

If you get an error similar to this (this is an example):
```
Building DAG of jobs...
Pulling singularity image docker://leandroishilima/subsampler:pandora_paper_tag1.
WorkflowError:
Failed to pull singularity image from docker://leandroishilima/subsampler:pandora_paper_tag1:
WARNING: pull for Docker Hub is not guaranteed to produce the
WARNING: same image on repeated pull. Use Singularity Registry
WARNING: (shub://) to pull exactly equivalent images.
ERROR: Image file exists, not overwriting.

  File "/hps/nobackup2/iqbal/leandro/pandora_paper_tag1/subsampler/venv/lib/python3.7/site-packages/snakemake/deployment/singularity.py", line 88, in pull
```

pass to the running script the default location where singularity images are store.
For example, in the EBI cluster, it would be `--singularity-prefix /nfs/leia/singularity/leandro/`.
