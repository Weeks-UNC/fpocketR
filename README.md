# fpocket-R
Program to find drug-like RNA-ligand binding pockets.

## Installation with Conda

1. Install conda on a Unix or Linux-based OS

    **NOTE: Windows users should install Windows Subsystem for Linux (WSL) and run conda from WSL**

2. Navigate to or create a directory where you would like to install the fpocket-R and RNAvigate packages.

3. Download the RNAvigate and fpocket-R packages

        git clone https://github.com/Weeks-UNC/RNAvigate.git
        git clone https://github.com/Weeks-UNC/fpocket-R.git

4. Create fpocket-R conda eviroment

        cd fpocket-R
        conda env create -f env.yaml
        conda activate fpocket-R
        conda develop .
        cd ../RNAvigate
        conda develop .
