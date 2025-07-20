# Example Workflow: Submitting fpocketR Batch Jobs with Bash

This directory demonstrates how to use a simple bash script to queue and run multiple fpocketR jobs back-to-back. This approach is useful for users who want to submit complex fpocketR commands in batches, without the need for workflow managers like Snakemake.

## Purpose
Show how to automate the submission of many fpocketR jobs using a bash script. This method is flexible for custom arguments, but runs jobs sequentially (slower than parallelized workflows).

## Prerequisites
- [fpocketR](https://github.com/Weeks-UNC/fpocketR) installed (see main package instructions)
- The fpocketR conda environment needs to be installed on your system
- Python (recommended: use the version provided in the fpocketR environment)

## File Structure
- `fpocketR_batch_file.txt`: Text file listing fpocketR command-line arguments (one job per line)
- `fpocketR_batch_submitter.sh`: Bash script to submit each job in sequence
- Output directories: Created by fpocketR for each job

## Workflow Steps

### 1. Prepare Your Batch File
List each fpocketR job as a line in `fpocketR_batch_file.txt`. Example:

```plaintext
-pdb 3e5c.pdb -o SAM-III_RS -y
-pdb 2l1v.pdb -ss 2l1v.nsd -o preQ1_RS -y
-pdb 2l1v.pdb -s 0 -ss 2l1v.nsd -o preQ1_RS_multistate -y
-pdb 8f4o_apo.pdb -al 2gdi_holo.pdb -nt 19,20,42,43 -o TPP_RS_apo_holo -y
-pdb 2gdi.pdb --chain Y --ligand TPP --dpi 10 --out TPP_RS -y
```

### 2. Run the Batch Submitter Script
Use the provided bash script to run all jobs in sequence:

```bash
bash fpocketR_batch_submitter.sh
```

This will:
- Read each line from `fpocketR_batch_file.txt`
- Run `python -m fpocketR` with the specified arguments
- Print each command before running it

### 3. Review Results
- Output files for each job are created in the specified output directories
- Check the terminal output for job progress and errors

## Notes
- This method is ideal for custom or complex fpocketR arguments
- Jobs run one after another (not in parallel)
- For parallel performance, consider using [Snakemake](../batch_submission_snakemake/README.md)

For more details, see the [fpocketR documentation](https://github.com/Weeks-UNC/fpocketR).
