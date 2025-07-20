# Batch Submission Using Snakemake

This workflow demonstrates how to use fpocketR with Snakemake to efficiently run pocket-finding jobs in parallel.

See the [README in the demo folder](../../demo/batch_submission_snakemake/README.md) for full details.

## Overview
- Configure batch runs in a YAML file
- Run the Snakemake workflow to process jobs in parallel
- Output and log files are organized by parameter set and sample

## Example Snakemake Command
```bash
snakemake --use-envmodules --use-conda --conda-frontend conda --cores all --rerun-incomplete --latency-wait 15 --keep-going &
```
