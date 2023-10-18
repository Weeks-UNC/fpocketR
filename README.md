# fpocket-R
fpocket-R is a modified version of [fpocket 4.0](https://github.com/Discngine/fpocket) and is optimized for finding, characterizing, and visulizing drug-like RNA-ligand binding pockets.

## Installation with Conda

### Windows Users:

1. fpocket-R requires a Unix/Linux to run properly, this mean that Windows users need to activate the Windows Subsystem for Linux (WSL). 

   * WSL is build-in to Windows 10 and 11 and allows users to easily run conda and fpocket-R using a Linux virtual machine. 

   * [Guide to installing WSL and Ubuntu](https://www.freecodecamp.org/news/how-to-install-wsl2-windows-subsystem-for-linux-2-on-windows-10/)

   * Once WSL is activated on your computer, use WSL/Ubuntu to install conda and fpocket-R using the installation instructions for Linux.

### Linux, Unix, MacOS, and WSL users:

2. Install conda on a Unix/Linux-based system

   * [Guide to install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
        * **macOS users:** Regular installation > macOS

        * **Windows users:** User WSL/Ubuntu to install conda: Regular installation > Linux

        * **Linux users:** You can figure it out

### Linux, Unix, MacOS, and WSL users with conda installed:

**NOTE for MacOS users**: fpocketR is only compatible with intel (x86) processors (not compatible with (arm) M1 or M2 processors) 

3. Open your terminal and navigate to the directory where you would like to install the fpocket-R and RNAvigate packages (optional).

    * **Tip for new Windows (WSL) users:** You can navigate to your Windows file system in the WSL command line by using the following command:

          cd /mnt/c/Users/<your-user-name>


4. Clone the RNAvigate and fpocket-R GitHub repositories. (RNAvigate is a dependancy for fpocket-R)

        git clone https://github.com/Weeks-UNC/RNAvigate.git
        git clone https://github.com/Weeks-UNC/fpocketR.git

5. Create fpocket-R conda eviroment and install fpocket-R and RNAvigate.

        cd fpocketR
        conda env create -f env.yml
        conda activate fpocketR
        conda develop .
        cd ../RNAvigate
        conda develop .

## Demo

[Demonstration of fpocketR usage](https://github.com/Weeks-UNC/fpocketR/blob/main/Demo/fpocketR_demo.md)

