# Batch Submission Using Bash Scripts

This workflow demonstrates how to use a simple bash script to queue and run multiple fpocketR jobs back-to-back.

See the [README in the demo folder](../../demo/batch_submission_bash/README.md) for full details.

## Overview
- List fpocketR jobs in a batch file
- Run the batch submitter script to process jobs sequentially
- Output files are created for each job

## Example Batch File
```plaintext
-pdb 3e5c.pdb -o SAM-III_RS -y
-pdb 2l1v.pdb -ss 2l1v.nsd -o preQ1_RS -y
-pdb 2l1v.pdb -s 0 -ss 2l1v.nsd -o preQ1_RS_multistate -y
-pdb 8f4o_apo.pdb -al 2gdi_holo.pdb -nt 19,20,42,43 -o TPP_RS_apo_holo -y
-pdb 2gdi.pdb --chain Y --ligand TPP --dpi 10 --out TPP_RS -y
```

## Run the Batch Submitter Script
```bash
bash fpocketR_batch_submitter.sh
```
