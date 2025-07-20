# Multistate Analysis on Modeled RNA Structures

This workflow demonstrates how to merge multiple RNA model PDB files, generate secondary structure drawings, and run fpocketR multistate analysis.

See the [README in the demo folder](../../demo/modeled_RNA/README.md) for full details.

## Overview
- Merge individual PDB files into a multistate file
- Generate secondary structure drawings
- Run fpocketR multistate analysis

## Example Command
```bash
python -m fpocketR -pdb ZTP_multistate.pdb -ss 9bzc.nsd --state 0 --knownnt 78,85,86,14 -o ZTP_multistate_out
```
