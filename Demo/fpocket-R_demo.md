
# fpocket Analysis Pipeline Demo

## Preparation

    conda activate fpocket-R
    cd ~/Weeks_Lab/fpocket4/fpocket-R_7.0/demo

Navigate to a working directory that contains a .pdb file(s). Secondary structure drawings such as .nsd file(s) are optional.

## Run from terminal

    ../fpocket-R_7.0.py -pdb 7ELR.pdb

This argument will run fpocket, analyze pockets, and make 3D figures.

The user will be prompted to input ligand name since multiple heteroatoms are detected.

## Run batches of files through fpocket-R using a bash script

    bash fpocket-R_7.0_job_submitter.sh

**Analyses run from bash scripts must specify ligand name in command (-l)**

### Contents of shell file

> * Specify an nsd file to create 2D figures (-nsd).
> * Specify a three character ligand residue name for holo structure analysis (-l). 
>
>       -pdb 3E5C.pdb -nsd 3E5C.nsd -l SAM
>
> * Specify the chain id of the RNA, if not chain A (-c).
>
>       -pdb 2GDI.pdb -nsd 2GDI.nsd -l TPP -c X 
>
> * Specify upto 2 chain ids (eg. A,B) for discontiguous RNAs (-c).
>
>       -pdb 1YKV.pdb -nsd 1YKV.nsd -l DAI -c A,B
>
> * Specify no ligand for apo structure analysis (-l).
> * Analyze all NMR states (-s).
> * Align all outputs to the input .pdb file (-al).
>   
>       -pdb 6MCI.pdb -nsd 6MCI.nsd -l no -s 0 -al
>
> * Specify the resolution for 3D figures to decrease render time or increase quality (-dpi).
> * Specify a custom name for the output directory (-o).
>
>       -pdb 7EZ0.pdb -nsd 7EZ0.nsd -l no -c N -dpi 50 -o group_I_intron
> 

## fpocket-R options


| Input options | Description                             |
| :------------- | :----------                             |
| -pdb STRING (Required)    | Specify a .pdb file to run fpocket.     |
| -nsd STRING     | Specify an .nsd file for generating secondary structure figures. |
| -a STRING       | Specify a directory contianing fpocket outputs for analysis (without running fpocket). |
| -s INT          | Specify a particular NMR states/model for analysis. Set to 0 for all. Default: NONE |

| fpocket options | Description                          |
| :-------------               | :----------                             |
| -m FLOAT                      | Sets fpocket -m flag. Specifies the minimum radius for an a-sphere. Default: 3.0 |
| -M FLOAT                      | Sets fpocket -M flag. Specifies the maximium radius for an a-sphere. Default: 5.7 |
| -i INT                        | Sets fpocket -i flag. Specifies the minimum number of a-spheres per pocket. Default: 42 |
| -D FLOAT                      | Sets fpocket -D flag. Specifies the a-sphere clustering distance for forming pockets. Default: 1.65 |
| -p FLOAT                      | Sets fpocket -p flag. Speciefies the maximum ratio of apolar a-spheres. Default: 0 |

| Output options | Description                           |
| :-------------   | :----------                           |
| -o STRING         | Specify name of fpocket output parent directory name. Default: fpocket-R_out_{fpocket parameters} |
| -n STRING         | Specify name prefix for fpocket_out and analysis_out subdirectories. |
| -y                | Overwrites output files and directories with same name. |

| Analysis settings | Description                           |
| :-------------   | :----------                           |
| -l STRING         | Specify the three character residue name of desired ligand. |
| -c STING          | Specify the chain(s) IDs conatining RNA. List upto 2 chains seperated by a comma (eg. A,B). Default: A |
| -lc STRING        | Specify the chain containing a ligand. Default: same as specified RNA chain. |
| -off INT          | Specify the offset between the structures in the input pdb and nsd files. Default: will gather offset from pdb header. |
| -qf FLOAT         | Specify the minimum fpocket score for a pocket to pass the quality filter. Default: 0.1 |
| -df FLOAT         | Specify distance filter (Angstroms) for identifying nucleotides near pockets. Default: 4.5 |

| Figure settings   | Description                           |
| :-------------   | :----------                           |
| -dpi INT          | Specify 3D figure resolution (dots per linear inch). Deafult: 300 |
| -zoom INT         | Specify zoom buffer distance to set the feild of view for 3D figures. Deafult= 10 |
| -cp               | Visually connects pockets in 2D figures. Deafult: False |
| -al               | Align output structures to input structure. Useful for multistate analysis. Default: True |
