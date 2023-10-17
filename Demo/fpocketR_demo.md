
# fpocket Analysis Pipeline Demo

## Preparation

    conda activate fpocketR
    cd ~/Weeks_Lab/fpocket4/fpocketR/demo

Navigate to a working directory that contains a .pdb file(s). Secondary structure drawings such as .nsd file(s) are optional.

## Run from terminal

    python -m fpocketR -pdb 7ELR.pdb

or

    python -m fpocketR -pdb 7ELR

This argument will run fpocket, analyze pockets, and make 3D figures.

The user will be prompted to input ligand name since multiple heteroatoms are detected (input: XAN).

## Run batches of files through fpocket-R using a bash script

    bash fpocketR_batch_submitter.sh

**NOTE: Analyses run from bash scripts must specify ligand name in command line (use the -l or --ligand options).**

### Contents of shell file

> * Specify an nsd file to create 2D figures (-nsd, --nsd).
> * Specify a three character ligand residue name for holo structure analysis (-l, --ligand). 
>
>       -pdb 3E5C.pdb -nsd 3E5C.nsd -l SAM
>
> * Specify the chain id of the RNA, if not chain A (-c, --chain).
>
>       -pdb 2GDI.pdb -nsd 2GDI.nsd -l TPP -c X 
>
> * Specify upto 2 chain ids (eg. A,B) for discontiguous RNAs (-c, --chain).
> * Specify the chain containing the ligand of interest (-l, --ligand).
>
>       -pdb 1YKV.pdb -nsd 1YKV.nsd -l DAI -c A,B -lc A
>
> * Specify *fpocket* parameters to use default parameters (optimized for proteins) or your own parameters (-m, --m; -M, --M, -D, --D; -i, --i; -A, --A; -p, --p).
>
>       -pdb 3E5C.pdb -nsd 3E5C.nsd -l SAM -m 3.4 -M 6.2 -D 2.4 -i 15 -A 3 -p 0
>
> * Specify no ligand for apo structure analysis (-l, --ligand).
> * Analyze all NMR states (-s, --state).
>
>       -pdb 6MCI.pdb -nsd 6MCI.nsd -l no -s 0 -al
>
> * Specify the resolution for 3D figures to decrease render time or increase quality (-dpi, --dpi).
> * Specify a custom name for the output directory (-o, --out).
>
>       -pdb 7EZ0.pdb -nsd 7EZ0.nsd -l no -c N -dpi 50 -o group_I_intron
>

## fpocket-R options


| Input options                 | Description                                                                                        |
| :---------------------------- | :------------------------------------------------------------------------------------------------- |
| -pdb, --pdb STRING (Required) | Specify a path to a .pdb file, .cif file, or 4 charater PDB indentification code to run fpocketR.  |
| -nsd, --nsd STRING            | Specify an .nsd file or other secondary structure file for generating secondary structure figures. |



| fpocket parameter options | Description                                                                                                           |
| :------------------------ | :-------------------------------------------------------------------------------------------------------------------- |
| -m, --m FLOAT             | Sets fpocket -m flag. Specifies the minimum radius for an a-sphere. (Default: 3.0)                                    |
| -M, --M FLOAT             | Sets fpocket -M flag. Specifies the maximium radius for an a-sphere. (Default: 5.7)                                   |
| -i, --i INT               | Sets fpocket -i flag. Specifies the minimum number of a-spheres per pocket. (Default: 42)                             |
| -D, --D FLOAT             | Sets fpocket -D flag. Specifies the a-sphere clustering distance for forming pockets. (Default: 1.65)                 |
| -A, --A INT               | Sets fpocket -A flag. Specifies the number of electronegative atoms required to define a polar a-sphere (Deafult: 3). |
| -p, --p FLOAT             | Sets fpocket -p flag. Speciefies the maximum ratio of apolar a-spheres. (Default: 0)                                  |

| Output options    | Description                                                                                         |
| :---------------- | :-------------------------------------------------------------------------------------------------- |
| -o, --out STRING  | Specify name of fpocket output parent directory name. (Default: fpocket-R_out_{fpocket parameters}) |
| -n, --name STRING | Specify name prefix for fpocket_out and analysis_out subdirectories.                                |
| -y, --yes BOOLEAN | Overwrites output files and directories with same name.                                             |

| Analysis settings          | Description                                                                                                                                                                |
| :------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| -s, --state INT            | Specify a particular NMR states/model for analysis. Set to 0 for all. (Default: NONE)                                                                                      |
| -c, --chain STING          | Specify the chain(s) IDs conatining RNA (case sensitive). List upto 2 chains seperated by a comma (eg. A,B). (Default: A)                                                  |
| -l, --ligand STRING        | Specify the three character residue name of desired ligand.                                                                                                                |
| -lc, --ligandchain STRING  | Specify the chain containing the ligand of interest. (Default: same as first RNA chain.)                                                                                   |
| -off, -offset INT          | Specify the offset between the structures in the input pdb and nsd files. Manual input is required for use with .cif files. (Default: will gather offset from pdb header.) |
| -qf, --qualityfilter FLOAT | Specify the minimum fpocket score for a pocket to pass the quality filter. (Default: 0.0)                                                                                  |

| Figure settings              | Description                                                                                 |
| :--------------------------- | :------------------------------------------------------------------------------------------ |
| -dpi, --dpi INT              | Specify 3D figure resolution (dots per linear inch). (Default: 300)                         |
| -zoom, --zoom INT            | Specify zoom buffer distance to set the feild of view for 3D figures. (Default= 10)         |
| -cp, --connectpocket BOOLEAN | Visually connects pockets in 2D figures. (Default: False)                                   |
| -al, --alignligand BOOLEAN   | Align output structures to input structure. Useful for multistate analysis. (Default: True) |
