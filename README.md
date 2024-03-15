# fpocketR
fpocketR is a modified version of [fpocket 4.0](https://github.com/Discngine/fpocket) and is optimized for finding, characterizing, and visulizing drug-like RNA-ligand binding pockets.

## Installation with Conda

### Windows Users: guide to install WSL and Ubuntu

1. fpocketR requires a Unix/Linux to run properly, this means that Windows users need to activate the Windows Subsystem for Linux (WSL). 

   * [Guide to installing WSL and Ubuntu](https://www.freecodecamp.org/news/how-to-install-wsl2-windows-subsystem-for-linux-2-on-windows-10/)

### Guide to install conda

2. Follow [guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to install conda.

   * **Windows WSL users:** Use [Linux installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html).

### Install fpocketR and RNAviagte:

**NOTE for MacOS users**: fpocketR is only compatible with intel (x86) processors (not compatible with (arm) M1 or M2 processors) 

3. Open your terminal and navigate to the directory where you would like to install the fpocketR and RNAvigate packages (optional).

    * **Tip for new Windows (WSL) users:** You can navigate to your Windows file system in the WSL command line by using the following command:

          cd /mnt/c/Users/<your-user-name>


4. Clone the RNAvigate and fpocketR GitHub repositories. (RNAvigate is a dependency for fpocketR)

        git clone https://github.com/Weeks-UNC/RNAvigate.git
        git clone https://github.com/Weeks-UNC/fpocketR.git

5. Create fpocketR conda environment and install fpocketR and RNAvigate.

        cd fpocketR
        conda env create --file env.yml
        conda activate fpocketR
        conda develop .
        cd ../RNAvigate
        conda develop .

## fpocketR demonstration and instructions

[Demonstration of fpocketR usage](https://github.com/Weeks-UNC/fpocketR/blob/main/Demo/fpocketR_demo.md)

