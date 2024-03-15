#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Module for controlling fpocket
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2022
#
# Version 1.0.1
#
# -----------------------------------------------------------------------------

import os
import subprocess
import shutil
from pymol import cmd


def find_pockets(
    pdb : str, chain : str, state : int, m : float, M : float, i :int,
    D : float, A : int, p : float, out : str, yes : bool
    ) -> str:
    """Pocket finding pipeline:
        - cleans pdb file to generate rna-only file
        - runs pocket prediction using fpocket
        - manages fpocket output files

    Args:
        pdb (str): Path to input .pdb file.
        chain (str): Chain identifier for desired RNA chain (default='A').
        state (int): Structural state to analyze.
        m (float): Min. a-sphere radius in angstroms (default=3.0).
        M (float): Max. a-sphere radius in angstroms (default=5.7).
        i (int): Min. number of a-spheres per pocket (default=42).
        D (float): a-sphere clustering distance in angstroms (default=1.65).
        out (str): name of fpocket output parent directory name.
        yes (boolean): Overwrite output files and directories with same name.

    Returns:
        str: path to output directory contianing fpocket outputs for analysis
    """

    # Path to clean (ligand-free) pdb file
    pdb_clean = f'{os.path.basename(pdb).rsplit(".")[0]}_clean.pdb'

    # Makes pdb_clean if it is not already a file
    if not os.path.isfile(pdb_clean):
        clean_pdb(pdb, pdb_clean)

    # Runs fpocket on cleaned pdb file.
    run_fpocket(pdb_clean, chain, state, m, M, i, D, A, p)

    # Files fpocket outputs into directories and manages overwriting.
    analysis = file_fpocket(pdb_clean, state, out, yes)
    return analysis

# -----------------------------------------------------------------------------


def clean_pdb(pdb : str, pdb_clean : str) -> None:
    """Cleans a .pdb file input and saves output as a .pdb file.
       Removes not polymer molecules (ligands) and proteins.
       Preserves modified/heteroatom RNA residues.

    Args:
        pdb (str): path to input .pdb file.
        pdb_clean (str): path to output (cleaned) .pdb file.
    """
    cmd.load(pdb)
    cmd.remove('not polymer')
    cmd.remove('byres polymer & name CA')
    cmd.alter('all', 'type="ATOM"')
    cmd.save(pdb_clean, state='0')
    cmd.reinitialize()


class MissingEnvironmentVariable(Exception):
    pass

def run_fpocket(
    pdb : str, chain :str, state : int, m : float, M : float, i :int,
    D : float, A : int, p : float
    ) -> None:
    """Detects potential binding pockets in RNA structures using fpocket.

    Args:
        pdb (str): Path to input .pdb file.
        chain (str): Chain identifier for desired RNA chain (default='A').
        state (int): Structural state to analyze.
        m (float): Min. a-sphere radius in angstroms (default=3.0).
        M (float): Max. a-sphere radius in angstroms (default=5.7).
        i (int): Min. number of a-spheres per pocket (default=42).
        D (float): a-sphere clustering distance in angstroms (default=1.65).
    """
    # Prints announcement that fpocket is searching for pockets.
    pdb_code = os.path.basename(pdb)[0:4]
    print(f'***** POCKET HUNTING {pdb_code} *****')
    # Runs fpocket bash commands
    
        
    bash_command = f'conda run -n fpocketR fpocket -f {pdb} -k {chain} -l {state} -m {m} -M {M} -i {i} -D {D} -A {A} -p {p} -w p'
            
    process = subprocess.Popen(bash_command.split(
    ), stderr=subprocess.PIPE, stdout=subprocess.PIPE, text=True)

    # Prints fpocket communications to console.
    result = process.communicate()
  
    for message in result:
        print(message)
        if 'EnvironmentLocationNotFound' in message:
            raise OSError('Unable to run fpocket because the fpocketR conda enviroment does not exist.\n'
                                   'Install the fpocketR conda enviroment.')


def file_fpocket(pdb : str, state : int, out : str, yes : bool) -> str:
    """Moves fpocket outputs into designated output directory.
       Default directory name specifies the fpocket parameters used.
       Manages overwriting files/directories if thet already exist.

    Args:
        pdb (str): Path to input .pdb file.
        state (int): Structural state to analyze.
        out (str): name of fpocket output parent directory name.
        y (boolean): Overwrites output files and directories with same name.

    Returns:
        str: path to output directory contianing fpocket outputs for analysis
    """
    # Moves fpocket output directories into a shared directory.
    source_dir = f'{pdb.rsplit(".")[0]}_out'
        
    if not os.path.isdir(source_dir):
        raise FileNotFoundError(f'fpocket output directory does not exist: {source_dir}.')

    if state is None:
        dest_dir = os.path.join(out, f'{pdb.rsplit(".")[0]}_out')
    else:
        dest_dir = os.path.join(out, f'{pdb.rsplit(".")[0]}_state{state}_out')

    if not os.path.isdir(dest_dir):
        analysis = shutil.move(source_dir, dest_dir)

    # Prompts user to overwrite an existing file with same name.
    elif os.path.isdir(dest_dir):
        if not yes:
            remove = input('A directory already exists with this name.\n'
                           f'{dest_dir}\n\n'
                           'Overwrite directory? [y/n]: ')
            print()
        if yes or remove in ('y', 'Y', 'yes', 'Yes'):
            shutil.rmtree(dest_dir)
            analysis = shutil.move(source_dir, dest_dir)
        else:
            print('Exiting program. \n'
                  'The name of the output directory can '
                  'be changed with the --name flag.')
            exit()

    # Adds state identifier number to file names.
    if state is not None:
        for file in os.listdir(dest_dir):
            file = os.path.join(dest_dir, file)
            if not os.path.isfile(file):
                continue

            head_tail = os.path.split(file)
            filename_ext = head_tail[1].split(pdb[0:-4])
            src = os.path.join(file)
            dst = os.path.join(
                head_tail[0], f'{pdb[0:-4]}_state{state}{filename_ext[1]}')

            # check if the file doesn't exist
            if not os.path.exists(dst):
                os.rename(src, dst)

    return analysis
