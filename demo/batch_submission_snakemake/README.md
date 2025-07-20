# Example Workflow: fpocketR Batch Analysis with Snakemake

This directory demonstrates how to use fpocketR with Snakemake to efficiently run pocket-finding jobs on dozens or hundreds of RNA structures in parallel. Users can copy and adapt this workflow for their own large-scale analyses.

## Purpose
Showcase a scalable, reproducible pipeline for batch RNA pocket analysis using fpocketR and Snakemake. This example is ideal for users who want to automate and parallelize pocket detection across many structures.

## Prerequisites
- [fpocketR](https://github.com/Weeks-UNC/fpocketR) installed (see main package instructions)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) installed
- Activate an environment with snakemake before running
- Python (recommended: use the version provided in the fpocketR environment)

## File Structure
- `workflow/snakefile`: Main Snakemake workflow
- `workflow/envs/fpocketR.yaml`: Conda environment used by Snakemake to run fpocketR
- `config/config.yaml`: Configuration file with parameters and sample info
- `logs/`: Directory for log files
- Output directories: `fpocketR_{parameter}/...` for results

## Workflow Steps

### 0. (Optional) Create and Activate a Conda Environment for Snakemake

If you do not already have Snakemake installed in a conda environment, you can create one as follows:

```bash
conda create -n snakemake -c conda-forge snakemake
conda activate snakemake
```

This ensures you have a clean environment with Snakemake available. You can then proceed with the workflow below.

### 1. Configure Your Batch Run
Edit `config/config.yaml` to specify:
- fpocketR parameters (e.g., m, M, D, i)
- List of samples and their ligands

```yaml
parameter:
  optimized:
    m: 3.0
    M: 5.7
    D: 1.65
    i: 42
  ....

sample:
  1ajf:
    ligand: NCO
  1am0:
    ligand: AMP
  ....
```

### 2. Run the Snakemake Workflow
Use the following command to launch the batch analysis (customize cores/jobs as needed):

```bash
snakemake --use-envmodules --use-conda --conda-frontend conda --cores all --rerun-incomplete --latency-wait 15 --keep-going &
```

This will:
- Run fpocketR jobs in parallel for each sample and parameter set
- Collect results and logs in organized output directories
- Continue running other jobs even if some fail
- Generate a summary .csv file for each parameter set

### 3. Review Results
- Output files for each sample are found in `fpocketR_{parameter}/{sample}_clean_out/`
- Summary files are generated for each parameter set (`pocketR_{parameter}/all_pocket_characteristics_{parameter}.csv`)
- Logs for each job are in `logs/`

## Notes
- This demo is designed for high-throughput, reproducible RNA pocket analysis
- You can easily scale up to hundreds of samples by editing the config file
- The workflow is robust to individual job failures (see `--keep-going`)

For more details, see the [fpocketR documentation](https://github.com/Weeks-UNC/fpocketR) and [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).
