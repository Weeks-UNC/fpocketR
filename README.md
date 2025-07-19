fpocketR
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/weeks-UNC/fpocketR/workflows/CI/badge.svg)](https://github.com/weeks-UNC/fpocketR/actions?query=workflow%3ACI)
[![PyPI version](https://img.shields.io/pypi/v/fpocketR.svg)](https://pypi.org/project/fpocketR/)

<img src="fpocketR_logo.png" alt="fpocketR logo" width="250" height="250" />

fpocketR is a modified version of [fpocket 4.0](https://github.com/Discngine/fpocket) and is optimized for finding, characterizing, and visualizing drug-like RNA-ligand binding pockets.

## Installation with Conda

### Windows Users: guide to install WSL and Ubuntu

1. fpocketR requires a Unix/Linux to run properly, this means that Windows users need to activate the Windows Subsystem for Linux (WSL). 

   * [Guide to installing WSL and Ubuntu](https://www.freecodecamp.org/news/how-to-install-wsl2-windows-subsystem-for-linux-2-on-windows-10/)

### Guide to install conda

2. Follow [guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to install conda.

   * **Windows WSL users:** Use [Linux installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).

### Install fpocketR

**NOTE for MacOS users**: fpocketR is not compatible with arm-based M1 or M2 processors (only intel (x86) processors).

3. Open your terminal and navigate to the directory where you would like to install the fpocketR package (optional).

        cd <path-to-GitHub-reposities>

    * **Tip for new Windows WSL users:** You can navigate to your Windows file system in the WSL command line by using the following command:

          cd /mnt/c/Users/<your-user-name>


4. Clone the RNAvigate and fpocketR GitHub repositories.

        git clone https://github.com/Weeks-UNC/fpocketR.git

5. Create fpocketR conda environment and install fpocketR and RNAvigate.

        cd fpocketR
        conda env create --file enviroment.yml
        conda activate fpocketR
        conda develop .

## Demo / Tutorial

[Demonstration of fpocketR usage](https://github.com/Weeks-UNC/fpocketR/blob/main/Demo/fpocketR_demo.md)

## Usage

| Input options                 | Description                                                                                      |
| :---------------------------- | :----------------------------------------------------------------------------------------------- |
| -pdb, --pdb STRING (Required) | Specify a path to a .pdb/.cif file, or a 4 charater PDB indentification code.                    |
| -nsd, --nsd STRING            | Specify an .nsd file or other secondary structure file to generate a secondary structure figure. |



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
| -n, --name STRING | Specify name prefix for fpocket_out and analysis_out subdirectories. (Default: None)                |
| -y, --yes BOOLEAN | Overwrites output files and directories with same name. (Default: False)                            |

| Analysis settings          | Description                                                                                                                                                                                                                        |
| :------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| -s, --state INT            | Specify a particular NMR states/model for analysis. Set to 0 for all. (Default: None)                                                                                                                                              |
| -c, --chain STING          | Specify the chain(s) IDs conatining RNA (case sensitive). List upto 2 chains seperated by a comma (eg. A,B). (Default: A)                                                                                                          |
| -l, --ligand STRING        | Specify the PDB ligand identification code (≤ 3 characters).                                                                                                                                                                       |
| -lc, --ligandchain STRING  | Specify the chain containing the ligand of interest. (Default: <first RNA chain>)                                                                                                                                                  |
| -off, -offset INT          | Specify offset between the starting nucleotide of the rna sequence and starting nucleotide of the PDB structure (usually = 0).<br>Manual input is required for use with .cif files. (Default: will gather offset from pdb header.) |
| -qf, --qualityfilter FLOAT | Specify the minimum fpocket score for a pocket to pass the quality filter. (Default: 0.0)                                                                                                                                          |

| Figure settings              | Description                                                                                 |
| :--------------------------- | :------------------------------------------------------------------------------------------ |
| -dpi, --dpi INT              | Specify 3D figure resolution (dots per linear inch). (Default: 300)                         |
| -zoom, --zoom INT            | Specify zoom buffer distance (Å) to set the feild of view for 3D figures. (Default: 10)     |
| -cp, --connectpocket BOOLEAN | Visually connects pockets in 2D figures. (Default: False)                                   |
| -al, --alignligand BOOLEAN   | Align output structures to input structure. Useful for multistate analysis. (Default: True) |

### Copyright

Copyright (c) 2025, Seth Veenbaas


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.11.
