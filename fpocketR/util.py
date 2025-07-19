#!/usr/bin/env python3
# -----------------------------------------------------
# Utilities functions for fpocket-R
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2024
#
# Version 1.3.0
#
# -----------------------------------------------------
import os
import re
from glob import glob
from prody import *
import numpy as np


def fetch_pdb(pdb_id : str) -> str:
    ''' Downloads structure from the PDB.
    Fetches .pdb format or .cif format if the .pdb is unavailable.
    Args:
        pdb_id: 4-character alphanumeric PDB identifier
    Returns:
        str: path to .pdb (or .cif) file
    '''
    pdb_id_lower = pdb_id.lower()
    pdb_filename = fetchPDB(f'{pdb_id_lower}', compressed=False, quiet=False)
    return pdb_filename

def natsorted(filenames : list) -> list:
    """ Sorts filenames by digits in the filenames.
    Supports multiple numbers in a filename.
    Natsorted example: [
        'pocket1_atm.pdb',
        'pocket2_atm.pdb',
        '2pocket11_atm.pdb',
        'pocket10_atm.pdb',
        'pocket11_atm.pdb',
        'pocket11_atm5.pdb',
        'pocket20_atm.pdb'
    ]
    Args:
        filenames (list): List of strings to be sorted by digits.
    Returns:
        list: Sorted list, ordered by digits.
    """
    # return sorted(filenames, key=lambda x: int(re.search(r'\d+', x).group()))
    
    def extract_sort_key(filename):
        # Extract all numbers in the filename as integers, or return the filename as-is for non-numeric parts
        return [int(num) if num.isdigit() else num.lower() 
                for num in re.split(r'(\d+)', filename)]
    
    return sorted(filenames, key=extract_sort_key)


def is_accessible(path : str, file_name : str) -> None:
    """Checks if input file is accessible.
       Stops program and displays:
        - FileNotFoundError: if file is not accessible/found.
        - TypeError: if path argument for input file is not a string.

    Args:
        path (str): path to file being checked.
        file_name (str): name of file being checked for display in error.

    Raises:
        FileNotFoundError: Not able to read input file.
        PermissionError: Not able to access {file_name}.')
    """
    if not os.path.exists(path):
        print(f'FileNotFoundError: The file {file_name} does not exist: {path}')
        raise FileNotFoundError(f'The file {file_name} does not exist at: {path}.')


    if not os.access(path, os.R_OK):
        print(f'PermissionError: Not able to access {file_name}.')
        raise PermissionError(f'Not able to access {file_name}.')


def get_first_rna_chain(pdb : str) -> str:
    """Identifies first RNA chain in .pdb/.cif file.

    Args:
        pdb (str): Path to input .pdb/.cif file.

    Returns:
        str: Chain identifier of the first chain containing RNA.
    """

    structure = prody.parsePDB(pdb)
        
    for ch in structure.getHierView():
        chid = ch.getChid()
        if structure.select(f'chain {chid} and nucleotide'):
            print(f'\nAutomatically selecting chain: {chid}\n'
                  'Use --chain to manually specify an RNA chain.\n')
            return chid
        
    print(f'No RNA chain found in {pdb}.\n'
                    'Verify that your pdb structure contains RNA.\n'
                    'Or manually set an RNA chain with (--chain) option.')
    raise ValueError(f'No RNA chain found in {pdb}.\n'
                    'Verify that your pdb structure contains RNA.\n'
                    'Or manually set an RNA chain with (--chain) option.')


def is_rna_chain(pdb : str, chain : str) -> None:
    """Identifies if a chain an a .pdb/.cif file contains rna.

    Args:
        pdb (str): Path to input .pdb/.cif file.
        chain (str): Chain identifier for RNA chain.

    Returns:
        None
    """
    structure = prody.parsePDB(pdb)
    chains = chain.split(',')
    hv = structure.getHierView()
    chids = list(hv)

    for chain in chains:
        if not structure.select(f'chain {chain}'):
            print(f'KeyError: The input {pdb} does not contain the chain: {chain}.\n'
            'Input a valid chain id using the (--chain) option.\n'
            'NOTE: Chain identifiers are case sensitive.\n'
            f'Chain(s) present in {pdb}:\n{chids}')
            raise KeyError(f'The input {pdb} does not contain the chain: {chain}.'
            'Input a valid chain id using the (--chain) option.'
            'NOTE: Chain identifiers are case sensitive.')
        elif not structure.select(f'chain {chain} nucleic'):
            print(f'KeyError: Chain {chain} of {pdb} does not contain rna.\n'
            'Input a valid RNA chain id using the (--chain) option.\n'
            'NOTE: Chain identifiers are case sensitive.'
            f'Chain(s) present in {pdb}:\n{chids}')
            raise KeyError(f'Chain {chain} of {pdb} does not contain rna.'
            'Input a valid RNA chain id using the (--chain) option.\n'
            'NOTE: Chain identifiers are case sensitive.')


def get_file_paths(
    analysis : str, name :str , pdb : str, state : int
    )-> tuple[str, str, str, str, list[str], str, str]:
    """Gets paths to required input files and check if they are accessible.

    Args:
        analysis (str): path to directory containing fpocket outputs.
        name (str): user specified filename prefix for analysis outputs.
        pdb (str): Path to input .pdb file.
        state (int): Structural state to analyze.

    Returns:
        str: valid path to *.pdb file
        str: valid path to fpocket generated *_out.pdb file
        str: valid path to fpocket generated *_out.pqr file
        str: valid path to fpocket generated *_info.txt file
        str: valid paths to pockets/pocket*_atm.pdb files
        str: 4 character pdb code
        str: filename prefix for analysis and figure output files
    """
    cwd = os.getcwd()
    analysis_basename = os.path.basename(analysis)[0:-4]

    if not pdb:
        pdb_basename = analysis_basename.replace('_clean', '')
        pdb = os.path.join(cwd, f'{pdb_basename}.pdb')
        is_accessible(pdb, 'pdb')

    pdb_out = os.path.join(cwd, analysis,
                           f'{analysis_basename}_out.pdb')
    is_accessible(pdb_out, 'pdb_out')

    pqr_out = os.path.join(cwd, analysis,
                           f'{analysis_basename}_pockets.pqr')
    is_accessible(pqr_out, 'pqr_out')

    info_txt = os.path.join(cwd, analysis,
                            f'{analysis_basename}_info.txt')
    is_accessible(info_txt, 'info_txt')

    pockets_dir = os.path.join(cwd, analysis, 'pockets')
    
    # Sort pockets_out filenames naturally (eg: 1,2,3...10,11....20,21)
    pockets_out = natsorted(glob(f'{pockets_dir}/pocket*_atm.pdb'))

    # Get PDB identifiers from the path to the fpocket out.pdb.
    pdb_basename = os.path.basename(pdb)
    pdb_code = pdb_basename[0:4]

    if state:
        name = f'{name}_state{state}'

    return pdb, pdb_out, pqr_out, info_txt, pockets_out, pdb_code, name


def calc_npr(atom_group: object) -> tuple[float, float]:
    """Generates the normalized PMI ratio for an input object.

    Args:
        atomgroup (object): prody atomgroup

    Returns:
        float: normalized principle ratio 1 = I1/I3
        float: normalized principle ratio 2 = I2/I3
    """
    # calculates inertia tensor for atomgroup
    inertia_tensor = calc_inertia_tensor(atom_group)

    # calculates the principle moments of inertia for the atomgroup
    pmi = calc_pmi(inertia_tensor)
    I1 = pmi[0]
    I2 = pmi[1]
    I3 = pmi[2]

    # calcualtes normalize PMI ratios
    npr1 = I1 / I3
    npr2 = I2 / I3

    return npr1, npr2


def calc_inertia_tensor(atom_group: object) -> np.array:
    """Generates an interia tensor for a prody atomgroup object.

    Args:
        atomgroup (object): prody atomgroup

    Returns:
        np.array: 3x3 inertia tensor
    """
    # Calculate center of mass
    totmass = 0.0
    x_com, y_com, z_com = 0, 0, 0

    for atom in atom_group:

        atom_coords = atom.getCoords()

        if atom.getResname() == 'STP':
            atom_mass = 1
        else:
            atom_mass = atom.getMass()

        x_com += atom_coords[0] * atom_mass
        y_com += atom_coords[1] * atom_mass
        z_com += atom_coords[2] * atom_mass
        totmass += atom_mass

    x_com /= totmass
    y_com /= totmass
    z_com /= totmass

    # Make and populate inertia tensor

    I = []
    for _ in range(9):
        I.append(0)

    for atom in atom_group:
        atom_coords = atom.getCoords()
        if atom.getResname() == 'STP':
            atom_mass = 1
        else:
            atom_mass = atom.getMass()

        temp_x, temp_y, temp_z = atom_coords[0], atom_coords[1], atom_coords[2]
        temp_x -= x_com
        temp_y -= y_com
        temp_z -= z_com

        I[0] += atom_mass * (temp_y ** 2 + temp_z ** 2)
        I[4] += atom_mass * (temp_x ** 2 + temp_z ** 2)
        I[8] += atom_mass * (temp_x ** 2 + temp_y ** 2)
        I[1] -= atom_mass * temp_x * temp_y
        I[3] -= atom_mass * temp_x * temp_y
        I[2] -= atom_mass * temp_x * temp_z
        I[6] -= atom_mass * temp_x * temp_z
        I[5] -= atom_mass * temp_y * temp_z
        I[7] -= atom_mass * temp_y * temp_z

    # generate inertia tensor
    inertia_tensor = np.array([(I[0:3]), (I[3:6]), (I[6:9])])

    return inertia_tensor


def calc_pmi(inertia_tensor : np.array) -> list:
    """Calaculates sorted eigen values (principle moments of inertia) 
    for an imput inertia tensor.

    Args:
        inertia_tensor (array): 3x3 inertia tensor

    Returns:
        list: sorted principle moments of inertia (PMI) [I1, I2, I3]
    """
    # calculate eigen values and eigen vectors
    eigvals, _ = np.linalg.eig(inertia_tensor)

    # sort eigen values > (I1, I2, I3)
    eig_ord = np.argsort(eigvals)

    sorted_eigvals = eigvals[eig_ord]

    return sorted_eigvals


def get_offset(pdb : str, chain : str, offset : int) -> int:
    """Calculates offset (or difference) between the nucleotide index (start at 1)
    and residue number for the first nucleotide in the first chain.

    Args:
        pdb (str): Path to input .pdb file.
        chain (str): Chain identifier for RNA chain.
        offset (int): Sequence offset between .pdb and .nsd file (default=None).

    Returns:
        int: Nucleotide offset of PDB chain.
    """

    if ',' in chain:
        chain = chain.split(',')[0]

    if pdb.endswith('.pdb'):
        structure = parsePDB(pdb)
    elif pdb.endswith('.cif'):
        structure = parsePDB(pdb)
    else:
        print(f'Could not parse structure: {pdb}')
        return None

    if structure:
        hv = structure.getHierView()
        ch = hv[chain]

    else:
        print(f'Could not parse structure: {pdb}')
        return None
    
    if ch:
        resn = ch.getResnums().min()
        offset = resn - 1
    
    else:
        print(f'Chain {chain} not found in structure: {pdb}')
        return None
    
    return offset

# Function to get the last processed state
def get_last_processed_state(state_tracker_filename):
    if os.path.exists(state_tracker_filename):
        with open(state_tracker_filename, 'r') as f:
            return int(f.read().strip())  # Read the last state
    else:
        return 0  # If no tracker file, start from state 1

# Function to update the state tracker
def update_last_processed_state(state_tracker_filename, state):
    with open(state_tracker_filename, 'w') as f:
        f.write(str(state))  # Write the current state to the file
        