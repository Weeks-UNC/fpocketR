
# Example Workflow: fpocketR Analysis on Modeled RNA

This directory provides a step-by-step example for running fpocketR on modeled RNA structures. Users can copy and adapt this workflow for their own data. All required input files and scripts are included here.

## Purpose
Demonstrate how to merge multiple RNA model PDB files, generate secondary structure drawings, and run fpocketR multistate analysis. The example uses RNA models (R1261v1) obtained from [CASP 16](https://predictioncenter.org/casp16/index.cgi).

## Prerequisites
- [fpocketR](https://github.com/Weeks-UNC/fpocketR) installed (see main package instructions)
- Activate the fpocketR conda environment before running any commands
- Python (recommended: use the version provided in the fpocketR environment)

## File Structure
- `./scripts/merge_pdb.py`: Script to merge PDB files
- `./data/`: Directory containing CASP RNA model PDB files
- `ZTP_multistate.pdb`: Example merged multistate PDB file (output from step 1)
- `9bzc.nsd`: Example secondary structure file (generated in step 2)

## Workflow Steps

### 1. Activate the fpocketR Conda Environment
Before running any commands, activate the fpocketR environment:

```bash
source $(conda info --base)/etc/profile.d/conda.sh
conda activate fpocketR
```

### 2. Merge PDB Files
Use the provided script to merge individual PDB files into a single multistate file.

```bash
python ./scripts/merge_pdb.py ZTP ./data/ ZTP_multistate.pdb R1261v1LG408_1.pdb
```

**Arguments:**
- `sample_name` (str): Name for the sample (used in PyMOL)
- `input_dir` (path): Directory containing PDB files to merge
- `output_file` (path): Output file name (should end with .pdb or .pse)
- `intra_fit_file` (path): PDB file to use as the target state for intra_fit alignment

### 3. (Optional) Generate Secondary Structure Drawing
You can create a secondary structure drawing from a PDB file using [RNApdbee 2.0](http://rnapdbee.cs.put.poznan.pl/) and [Structure Editor](https://rna.urmc.rochester.edu/GUI/html/StructureEditor.html). 

### 4. Run fpocketR Multistate Analysis
Analyze the merged RNA structure using fpocketR:

```bash
python -m fpocketR -pdb ZTP_multistate.pdb -ss 9bzc.nsd --state 0 --knownnt 78,85,86,14 -o ZTP_multistate_out
```

**Arguments:**
- `-pdb`: Path to the merged multistate PDB file
- `-ss`: Path to the secondary structure file
- `--state`: State index to analyze
- `--knownnt`: Comma-separated list of known nucleotide indices
- `-o`: Output directory

## Notes
- All example input files are provided in this directory.
- Users are encouraged to adapt the workflow for their own RNA models and secondary structure files.

For more details, see the [fpocketR documentation](https://github.com/Weeks-UNC/fpocketR).