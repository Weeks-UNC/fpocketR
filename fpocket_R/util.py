#!/usr/bin/env python3
# -----------------------------------------------------
# Utilities functions for fpocket-R
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2022
#
# Version 1.0.0
#
# -----------------------------------------------------
import os
from glob import glob
from prody import *
import numpy as np


def is_accessible(path, file_name):
    """Checks if input file is accessible.
       Stops program and displays:
        - FileNotFoundError: if file is not accessible/found.
        - TypeError: if path argument for input file is not a string.

    Args:
        path (str): path to file being checked.
        file_name (str): name of file being checked for display in error.

    Raises:
        FileNotFoundError: not able to read input file.
    """
    try:
        if not os.access(path, os.R_OK):
            raise FileNotFoundError(f'FileNotFoundError: Unable to access/find '
                                    f'a valid {file_name} file.\n{path}\n')

    except FileNotFoundError as error:
        print(f'{error}')
        if file_name in ['pdb', 'nsd']:
            print(f'Provide path to the input {file_name} file using the '
                  f'-{file_name} flag.')
        raise error

    except TypeError:
        print(f'TypeError: no path provided for the {file_name} file.\n')
        if file_name in ['pdb', 'nsd']:
            print(f'Provide path to the input {file_name} file using the '
                  f'-{file_name} flag.')
        exit()


def contains_structure(pdb, chain):
    if chain != "None":
        pdb_structure = parsePDB(pdb)
        chain = chain.split(',')[0]
        # if not pdb_structure.numAtoms('nucleic'):
        #     print(f'\nERROR: The input file {pdb} does not contain an RNA.\n'
        #           'Input a pdb file containing an RNA using the (-pdb) flag.')
        #     exit()
        if not pdb_structure.select(f'chain {chain}'):
            print(f'\nKeyError: The input file {pdb} does not contain a chain {chain}.\n'
                  'Input a valid chain id using the (-c) flag')
        #     exit()
        # if not pdb_structure.select(f'chain {chain} and nucleic'):
        #     print(f'\nError: The input file {pdb} does not contain RNA in chain {chain}.\n'
        #           'Input a chain id containing an RNA using the (-c) flag')
        #     exit()


def get_file_paths(analysis, name, pdb, state):
    """Gets paths to required input files and check if they are accessible.

    Args:
        analysis (str): path directory contianing fpocket outputs for analysis
        name (str): user specified filename prefix for analysis outputs
        pdb (str): Path to input .pdb file.

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
    pockets_out = glob(f'{pockets_dir}/pocket*_atm.pdb')

    # Get PDB identifiers from the path to the fpocket out.pdb.
    pdb_basename = os.path.basename(pdb)
    pdb_code = pdb_basename[0:4]
    pdb_name = pdb_basename[0:-4]

    # Creates default name varible for naming output files.
    if name is None:
        name = pdb_name

    if state is not None:
        name = f'{name}_state{state}'

    return pdb, pdb_out, pqr_out, info_txt, pockets_out, pdb_code, name


def calc_npr(atom_group: object):
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


def calc_inertia_tensor(atom_group: object):
    """Generates an interia tensor for a prody atomgroup object.

    Args:
        atomgroup (object): prody atomgroup

    Returns:
        np.array: inertia tensor
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


def calc_pmi(inertia_tensor):
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
