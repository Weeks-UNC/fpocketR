#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Module for parsing fpocket outputs and analyzing pocket characteristics
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2022
#
# Version 0.1.0
#
# -----------------------------------------------------------------------------
from re import findall
from prody import *
import numpy as np
import pandas as pd
from fpocket_R import util


def analyze_pockets(pdb, pdb_out, name, info_txt, pockets_out, pdb_code, chain, state, ligandchain, ligand,
                    m, M, i, D, A, p, distancefilter, qualityfilter, offset,
                    pqr_out, analysis):

    # Parses pdb files and returns prody structure objects.
    ligand_structure = parsePDB(pdb)
    out_structure = parsePDB(pdb_out)

    # Get offset of pdb_ligand file.
    if offset is None:
        ligand_polymer = parsePDBHeader(pdb, 'polymers')
        offset = get_offset(ligand_polymer, chain)

    # Sets ligand chain to first pdb chain by default.
    if ligandchain is None:
        ligandchain = chain[0]

    # Creates a dataframe with pocket characteristics from fpocket info.txt
    pc_df = get_characteristics(info_txt, pdb_code, m, M, i, D, state)

    # rna_coords = out_structure.nucleic.copy()
    rna_coords = out_structure.protein.copy()

    # If pockets are detected calculate features and add to pc_df.
    if out_structure.select('resname STP'):
        # Get atomgroup for the ligand and rna.
        ligand_coords = get_ligand_coords(
            ligand_structure, ligand, ligandchain, name)
        stp_coords = out_structure.select('resname STP').copy()

        add_basic_characteristics(stp_coords, pockets_out,
                                  qualityfilter, pc_df, chain)

        if ligand_coords:
            add_ligand_characteristics(
                stp_coords, ligand_coords, distancefilter, pc_df)

    return pc_df, rna_coords, offset


def get_offset(ligand_polymer, chain):
    """Locates PDB nucleotide offset as documented in the dbrefs section of
    the input .pdb file. NOTE: Can be reported incorrectly.

    Args:
        ligand_polymer (object): ProDy parsed PDB header with polymer info.
        chain (str): Desired chain identifier (default='A').

    Returns:
        int: Nucleotide offset of PDB chain.
    """

    for idx, ch in enumerate(ligand_polymer):
        if ch.chid == chain[0]:
            ch_idx = idx
    ligand_polymer_chain = ligand_polymer[ch_idx]
    dbref = ligand_polymer_chain.dbrefs[0]
    offset = dbref.first[0] - 1
    return offset


def get_ligand_coords(ligand_structure, ligand, ligandchain, name):
    """Gets coordinates for all atoms in a RNA binding ligand.

    Args:
        ligand_structure (str): Path to .pdb file with RNA-ligand complex.
        ligand (str): Three character ligand identifier (optional).
        ligandchain (str): Chain containing ligand.
        name (str): name of input pdb file

    Returns:
        object: ProDy atomgroup of all atoms in known RNA ligand.
    """
    if ligand == 'n' or ligand == 'no' or ligand == 'none':
        return None

    elif ligand and 2 <= len(ligand) <= 3:
        try:
            ligand_coords = ligand_structure.select(
                f'chain {ligandchain} and resname {ligand}').copy()
        except AttributeError:
            print(f'Chain {ligandchain} of {name}.pdb does not contain a ligand named: '
                  f'{ligand}\n'
                  'Please provide a valid ligand chain (-lc) '
                  'and ligand residue name (-l).')
            exit()

    elif ligand_structure.select(
            f'chain {ligandchain} and hetatm and not ion and not water') is None:
        return None

    else:
        ligand_sele = ligand_structure.select(
            f'chain {ligandchain} and hetatm and not ion and not water').copy()
        hetatm_resn = np.unique(ligand_sele.getResnames()).tolist()

        if len(hetatm_resn) == 1 \
                and len(hetatm_resn[0]) == 3 \
                and hetatm_resn[0] != 'GTP':
            ligand_coords = ligand_sele.select(
                f'resname {hetatm_resn[0]}').copy()
            print(f'Using {hetatm_resn[0]} as ligand for analysis.')

        elif len(hetatm_resn) == 1 and len(hetatm_resn[0]) < 3:
            print('No ligand detected.\n')
            return None

        elif len(hetatm_resn) > 0:
            input_resn = input(
                f'Detected heteroatoms: {hetatm_resn}.\n'
                'Input the three character name of your '
                'desired ligand or None for no ligand:')
            print()

            if 'None' in input_resn or 'none' in input_resn or 'no' in input_resn:
                return None
            elif input_resn in hetatm_resn:
                ligand_coords = ligand_sele.select(
                    f'resname {input_resn}').copy()
                print(f'Using {input_resn} as ligand for analysis.')
            else:
                print(f'NameError: {input_resn} is not a not a valid ligand.\n'
                      'Proceeding with no known ligand...\n')
                return None

        else:
            print('No ligand detected.\n')
            return None

    return ligand_coords


def get_characteristics(info_txt, pdb_code, m, M, i, D, state):
    """Creates a dataframe containing characteristics for
        all fpocket generates pockets.

    Args:
        info_txt (string): path to *info.txt file outputed by fpocket
        pdb_code (string): 4 digit identifier for the PDB structure
        m (float): Specifies the minimum radius for an a-sphere.
        M (float): Specifies the maximum radius for an a-sphere.
        i (int): Specifies the minimum number of a-spheres per pocket.
        D (float): Specifies the a-sphere clustering distance for pockets.
        state (int): Structural state to analyze.

    Returns:
        DataFrame: Displays the characteristics and properties most relvant
        to scoring each pocket.
    """

    pc_d = {'Parameters': [], 'PDB': [], 'State': [], 'Pocket': [], 'Score': [],
            'Drug score': [], 'a-sphere': [], 'SASA': [],
            'Volume': [], 'Hydrophobic density': [],
            'Apolar a-sphere proportion': [],
            'Hydrophobicity score': [], 'Polarity score': []}

    with open(info_txt, 'r') as f:

        poc_count = 1

        for row in f:
            if 'Pocket' in row:
                pc_d['Parameters'].append(f'-m {m} -M {M} -i {i} -D {D}')
                pc_d['PDB'].append(pdb_code)
                pc_d['State'].append(state)
                pc_d['Pocket'].append(int(row.split(' ')[1].strip()))
                pocket = int(row.split(' ')[1].strip())
            if 'Score :' in row:
                if 'Druggability' not in row:
                    pc_d['Score'].append(float(row.split(':')[1].strip()))
            if 'Druggability Score :' in row:
                pc_d['Drug score'].append(float(row.split(':')[1].strip()))
            if 'Number of Alpha Spheres :' in row:
                pc_d['a-sphere'].append(int(row.split(':')[1].strip()))
            if 'Total SASA :' in row:
                pc_d['SASA'].append(float(row.split(':')[1].strip()))
            if 'Volume :' in row:
                pc_d['Volume'].append(float(row.split(':')[1].strip()))
            if 'Mean local hydrophobic density' in row:
                pc_d['Hydrophobic density'].append(
                    float(row.split(':')[1].strip()))
            if 'Apolar alpha sphere proportion' in row:
                pc_d['Apolar a-sphere proportion'].append(
                    float(row.split(':')[1].strip()))
            if 'Hydrophobicity score' in row:
                pc_d['Hydrophobicity score'].append(
                    float(row.split(':')[1].strip()))
            if 'Polarity score:' in row:
                pc_d['Polarity score'].append(float(row.split(':')[1].strip()))
            if 'Flexibility' in row and poc_count == pocket:
                poc_count += 1

        pc_df = pd.DataFrame.from_dict(pc_d)
        pc_df.insert(loc=4, column='Filter', value='Fail')
        pc_df.insert(loc=4, column='Type', value='Novel')

    return pc_df


def add_basic_characteristics(stp_coords, pockets_out,
                              qualityfilter, pc_df, chain):
    """Adds characteristics to the pocket characteristics DataFrame that do
        not require a ligand to calculate.
    PocketNT: nucleotides in contact with pocket,
    Geometric Center: X,Y,Z coordinates of pocket center,
    Filter: Pocket quality filter (Pass or Fail)

    Args:
        stp_coords (object): ProDy atom group of a-sphere coordinates.
        rna_coords (object): ProDy atom group of rna coordinates.
        distancefilter (float): Distances in angstroms used to determine
                                RNA nucleotide contacts with pockets.
        pc_df (DataFrame): Characteristics and properities for each pocket.
        c (str): Chain identifier for desired RNA chain (default='A').

    """
    pocket_center_l = []
    pocketNT_l = []

    for i, residue in enumerate(stp_coords.iterResidues()):
        pocket = residue.getCoords()

        # Calculate the geometric center of each pocket.
        pocket_center = calcCenter(pocket)
        pocket_center_l.append(list(pocket_center.round(decimals=2)))

        # Calculate pocketNT (nucleotides near each pocket)
        x = pockets_out[i]
        with open(x, 'r') as f:
            tmp0 = f.read()
            pc_dict = {'Pocket': [(int(x[1]))] for x in findall(
                r'(\w+\s\w+\s*)(\d)(:)(\n)', tmp0)}

            pc_dict.update({x[1].strip(): float(x[2]) for x in findall(
                r'(\w+\s\w+\s*-\s*)(.+):\s*([\d.-]+)(\n)', tmp0)})

        pdb = parsePDB(x)

        if ',' in chain:
            chain = chain.replace(',', ' & ')
        selection = pdb.select(f'chain {chain}')
        nt = selection.getResnums().tolist()
        pocketNT_l.append(np.unique(nt).tolist())

    # Add pocketNT and geometric center to pc dataframe.
    pc_df['PocketNT'] = pocketNT_l
    pc_df['Geometric center'] = pocket_center_l

    # Add pocket filter (Pass or Fail) to pc dataframe.
    pc_df.loc[pc_df['Score'] > qualityfilter, 'Filter'] = 'Pass'


def add_ligand_characteristics(stp_coords, ligand_coords,
                               distancefilter, pc_df):
    """Adds characteristics to the pocket characteristics DataFrame that
       require the presence of a ligand to calculate.
    Pocket overlap: Ratio of a-spheres in contact with ligand.
    Ligand overlap: Ratio of ligand atoms in contact with a-spheres.
    Center criteria: Check to verify that at least one atom in the
                    ligand is within the specified distance filter
                    to the geometric center of the pocket.
    Type: Indicates if pocket is known or novel based on its
          overlap and offset with the specified ligand.

    Args:
        stp_coords (object): ProDy atom group of a-sphere coordinates.
        ligand_coords (object): ProDy atom group of ligand coordinates.
        distancefilter (float): Distances in angstroms used to determine
                                RNA nucleotide contacts with pockets.
        pc_df (DataFrame): Characteristics and properities for each pocket.
    """

    pocket_overlap_l = []
    ligand_overlap_l = []
    center_criteria_l = []

    ligand_npr1_l = []
    ligand_npr2_l = []

    # Iterates through each pocket (residue).
    for _, residue in enumerate(stp_coords.iterResidues()):
        pocket = residue.getCoords()

        # Calculate the normalize PMI ratios for each pocket
        npr1, npr2 = util.calc_npr(ligand_coords)
        ligand_npr1_l.append(npr1)
        ligand_npr2_l.append(npr2)

        # Calculate ratio of a-spheres in each pocket that is near the ligand.
        pocket_overlap_sele = residue.select(
            f'within {distancefilter} of ligand', ligand=ligand_coords)
        if pocket_overlap_sele is None:
            pocket_overlap_l.append(0)
        else:
            # Calculates total number of a-spheres in pocket.
            _, pocket_stp_count = np.unique(
                residue.getResnums(), return_counts=True)

            # Calculates number of overlapped a-spheres in pocket.
            _, pocket_overlap_stp_count = np.unique(
                pocket_overlap_sele.getResnums(), return_counts=True)

            pocket_overlap_l.append(float(
                pocket_overlap_stp_count / pocket_stp_count))

        # Calculate ratio of a-spheres in each pocket that is near the ligand.
        ligand_overlap_sele = ligand_coords.select(
            f'within {distancefilter} of pocket', pocket=residue)
        if ligand_overlap_sele is None:
            ligand_overlap_l.append(0)
        else:
            # Calculates total number of atoms in ligand.
            ligand_atom_count = ligand_coords.numAtoms()

            # Calculates number of overlapped atoms in ligand.
            ligand_overlap_atom_count = ligand_overlap_sele.numAtoms()

            ligand_overlap_l.append(
                ligand_overlap_atom_count / ligand_atom_count)

        # Check center criteria between pocket and ligand.

        pocket_center = calcCenter(pocket)

        center_criteria_sele = ligand_coords.select(
            f'within {distancefilter} of center', center=pocket_center)

        if center_criteria_sele is None:
            center_criteria_l.append(0)
        elif center_criteria_sele.numAtoms() is False:
            center_criteria_l.append(0)
        else:
            center_criteria_l.append(1)

    # Add Ligand overlap and Center criteria to pc dataframe.
    pc_df['Pocket overlap'] = pocket_overlap_l
    pc_df['Ligand overlap'] = ligand_overlap_l
    pc_df['Center criteria'] = center_criteria_l
    pc_df['Ligand npr1'] = ligand_npr1_l
    pc_df['Ligand npr2'] = ligand_npr2_l

    # Add pocket type (Known or Novel) to pc_df.
    pc_df.loc[(pc_df['Ligand overlap'] >= 0.25) &
              (pc_df['Pocket overlap'] >= 0.25) &
              (pc_df['Center criteria'] >= 1), 'Type'] = 'Known'
