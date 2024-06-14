#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Module for parsing fpocket outputs and analyzing pocket characteristics
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2022
#
# Version 1.0.1
#
# -----------------------------------------------------------------------------
import os
import math
from re import findall, sub
from prody import *
import numpy as np
import pandas as pd
import requests
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import QED
import trimesh
from fpocketR import util


def analyze_pockets(
    pdb : str, pqr_out : str, pdb_out : str, analysis : str, name : str,
    info_txt : str, pockets_out : list[str], pdb_code : str, chain : str,
    state, ligandchain : str, ligand : str, m : float, M : float, i : int,
    D : float, A : int, p : float, qualityfilter : float
    ) -> tuple[pd.DataFrame, prody.AtomGroup]:

    # Parses pdb files and returns prody structure objects.
    ligand_structure = parsePDB(pdb)
    out_structure = parsePDB(pdb_out)

    # Sets ligand chain to first pdb chain by default.
    if ligandchain is None:
        ligandchain = chain[0]

    # Create real_sphere.pdb ouput be combinding the pqr_out and pdb_out.
    get_real_sphere(pqr_out, pdb_out, analysis, name)

    # Creates a dataframe with pocket characteristics from fpocket info.txt
    pc_df = get_characteristics(info_txt, pdb_code, name, m, M, i, D, A, p, state)

    rna_coords = out_structure.copy()

    # If pockets are detected calculate features and add to pc_df.
    if out_structure.select('resname STP'):

        # Get atomgroup for rna and add pocket characteristics.
        stp_coords = out_structure.select('resname STP').copy()

        add_basic_characteristics(stp_coords, pockets_out, qualityfilter,
                                  pc_df, chain, analysis, name)
        
        # Get atomgroup for ligand and add ligand characteristics.
        if ligand in ('n', 'N', 'no', 'No', 'none', 'None'):
            ligand_coords = None
            ligand = None
        else:
            ligand_coords, ligand = get_ligand_coords(
                ligand_structure, ligand, ligandchain, analysis, name)
            if ligand_coords:
                add_ligand_characteristics(analysis, stp_coords, ligand_coords,
                                           ligand, pc_df)

    return pc_df, rna_coords


# -----------------------------------------------------------------------------


def get_real_sphere(
    pqr_file : str, pdb_file : str, analysis : str, name : str
    ) -> None:
    """Gets fpocket pocket a-sphere radii form a pqr file and encodes it into
    into the B factor column of a *real_sphere.pdb output file.

    Args:
        pqr_file (string): Path to fpocket *out.pqr file
                           containing a-spheres radii.
        pdb_file (string): Path to fpocket *out.pdb file.
        analysis (str): Path directory contianing fpocket outputs for analysis.
        name (str): User specified filename prefix for analysis outputs.
    """
    structure = parsePDB(pdb_file)
    pockets = parsePQR(pqr_file)

    for residue in structure.iterResidues():
        resnum = residue.getResnum()
        resname = residue.getResname()
        if resname == 'STP':
            for atom in residue.iterAtoms():
                coords = atom.getCoords()
                for pqr_atom in pockets.iterAtoms():
                    pqr_resnum = pqr_atom.getResnum()
                    pqr_coords = pqr_atom.getCoords()
                    if resnum == pqr_resnum and (coords == pqr_coords).all():
                        radius = pqr_atom.getRadius()
                        atom.setBeta(radius)
    writePDB(f'{analysis}/{name}_out_real_sphere.pdb', structure)


def get_ligand_coords(
    ligand_structure : prody.AtomGroup,
    ligand : str,
    ligandchain : str,
    analysis : str,
    name : str
    ) -> tuple[prody.AtomGroup, str]:
    """Gets coordinates and residue name for RNA-binding ligand.

    Args:
        ligand_structure (object): Prody structure of RNA-ligand complex.
        ligand (str): Ligand residue name (usually a 3-letter code).
        ligandchain (str): Chain identifier for desired ligand.
        analysis (str): Path directory contianing fpocket outputs for analysis.
        name (str): Name of input pdb file.

    Returns:
        object: ProDy atomgroup of all atoms in known RNA ligand.
        str: Ligand residue name (usually a 3-letter code).
    """
    if ligand and 2 <= len(ligand) <= 3:
        try:
            ligand_coords = ligand_structure.select(
                f'resname {ligand}').copy()
        except AttributeError:
            print(f'Chain {ligandchain} of {name}.pdb does not contain a '
                  f'ligand named: {ligand}\n'
                  'Please provide a valid ligand chain (--ligandchain) '
                  'and ligand residue name (--ligand).')
            exit()

    elif ligand_structure.select(f'chain {ligandchain} and hetatm '
                                 'and not ion and not water') is None:
        return (None, None)

    else:
        ligand_sele = ligand_structure.select(
            f'chain {ligandchain} and hetatm and not ion and not water').copy()
        hetatm_resn : list = np.unique(ligand_sele.getResnames()).tolist()
        
        if len(hetatm_resn) > 1:
                resnames_qeds : dict = {}
                for resname in hetatm_resn:
                    try:
                        response = requests.get(
                            f'https://files.rcsb.org/ligands/download/{resname}_model.sdf')
                        
                        with open(f'{analysis}/{resname}_model.sdf', 'wb') as f:
                            f.write(response.content)
                        mol = Chem.MolFromMolFile(f'{analysis}/{resname}_model.sdf')
                        qed = QED.default(mol)
                        mw : float = Chem.rdMolDescriptors.CalcExactMolWt(mol)
                        pat = Chem.MolFromSmarts("[#6]")
                        num_of_carbon : int =len(mol.GetSubstructMatches(pat))
                        print(f'# of C: {num_of_carbon}')
                        if mw < 100.0 or num_of_carbon < 2:
                            qed = 0
                            continue
                         
                    except:
                        print(f'Error: Not able to calculate QED score for {resname}.\n')
                        qed = 0
                    
                    resnames_qeds[resname] = qed

                # Remove heteroatoms with nan qed scores
                # filtered = {k: v for k, v in resnames_qeds.items() if v is not np.nan}
                # resnames_qeds.clear()
                # resnames_qeds.update(filtered)

                hetatm_resn = [max(resnames_qeds, key=resnames_qeds.get)]

        if len(hetatm_resn) == 1 \
                and 2 <= len(hetatm_resn[0]) <= 3:
            ligand_coords = ligand_sele.select(
                f'chain {ligandchain} and resname {hetatm_resn[0]}').copy()
            print(f'Using {hetatm_resn[0]} as ligand for analysis.')
            ligand = hetatm_resn[0]
            return (ligand_coords, ligand)

        elif len(hetatm_resn) == 1 and 2 > len(hetatm_resn[0]) > 3:
            print('No ligand detected.\n')
            return (None, None)

        elif len(hetatm_resn) > 0:
            while True:
                input_resn = input(
                    f'Detected heteroatoms: {hetatm_resn}.\n\n'
                    'Input the target ligand ID (case-sensitive; "none" for no ligand): ')
                print('\n')

                if input_resn in ('n', 'N', 'no', 'No', 'none', 'None'):
                    return (None, None)
                elif input_resn in hetatm_resn:
                    ligand_coords = ligand_sele.select(
                        f'chain {ligandchain} and resname {input_resn}').copy()
                    print(f'Using {input_resn} as ligand for analysis.')
                    ligand = input_resn
                    break
                else:
                    print(f'NameError: {input_resn} is not a not a valid '
                          'ligand ID.\n'
                          'Input "None" to proceed without a ligand.\n')

        else:
            print('No ligand detected.\n')
            return (None, None)
    if ligand_coords is None:
        return (None, None)
    else:
        return (ligand_coords, ligand)


def get_characteristics(
    info_txt : str, pdb_code :str, name : str, m : float, M : float,
    i : int, D : float, A : int, p : float, state : int
    ) -> pd.DataFrame:
    """Creates a dataframe containing characteristics for
        all fpocket generates pockets.

    Args:
        info_txt (string): path to *info.txt file outputed by fpocket.
        pdb_code (string): 4 digit identifier for the PDB structure.
        m (float): Specifies the minimum radius for an a-sphere.
        M (float): Specifies the maximum radius for an a-sphere.
        i (int): Specifies the minimum number of a-spheres per pocket.
        D (float): Specifies the a-sphere clustering distance for pockets.
        A (int): # of electroneg. atoms to define a polar a-sphere.
        p (float): Max. ratio of apolar a-spheres in a pocket.
        state (int): Structural state to analyze.

    Returns:
        DataFrame: Displays the characteristics and properties most relvant
        to scoring each pocket.
    """

    pc_d = {'Parameters': [], 'PDB': [], 'State': [], 'Pocket': [],
            'Score': [], 'Drug score': [], 'a-sphere': [],
            'SASA': [], 'Volume': [], 'Hydrophobic density': [],
            'Apolar a-sphere proportion': [],
            'Hydrophobicity score': [], 'Polarity score': []}

    with open(info_txt, 'r') as f:

        poc_count = 1

        for row in f:
            if 'Pocket' in row:
                pc_d['Parameters'].append(f'-m {m} -M {M} -i {i} -D {D} '
                                          f'-A {A} -p {p}')
                pc_d['PDB'].append(name)
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


def add_basic_characteristics(
    stp_coords : prody.AtomGroup, pockets_out : list[str], 
    qualityfilter : float, pc_df : pd.DataFrame, chain :str,
    analysis : str, name :str
    ) -> None:
    """Adds characteristics to the pocket characteristics DataFrame that do
        not require a ligand to calculate.
    PocketNT: nucleotides in contact with pocket,
    Geometric Center: X,Y,Z coordinates of pocket center,
    Filter: Pocket quality filter (Pass or Fail)
    
    Args:
        stp_coords (object): ProDy atom group of a-sphere coordinates.
        pockets_out list[str]: Valid paths to pockets/pocket*_atm.pdb files.
        qualityfilter (float): Minimum fpocket score filter for pockets.
        pc_df (DataFrame): Characteristics and properities for each pocket.
        chain (str): Chain identifier for desired RNA chain (default='A').
        analysis (str): path directory contianing fpocket outputs for analysis.
        name (str): Name of input pdb file.
    """
    pocketNT_l = []
    pocket_npr1_l = []
    pocket_npr2_l = []

    for i, _ in enumerate(stp_coords.iterResidues()):
        # Export surface obj files for each pocket.
        cmd.load(f'{analysis}/{name}_out_real_sphere.pdb')
        cmd.remove(f'not resn STP or not resi {i+1}')
        cmd.hide('spheres')
        cmd.show('surface')
        cmd.set('surface_quality', 1)
        cmd.save(f'{analysis}/pockets/pocket{i+1}_surf.obj',
                 f'pocket_{i+1}_surface')
        cmd.reinitialize()

        # Create mesh object for each pocket.
        mesh = trimesh.load(f'{analysis}/pockets/pocket{i+1}_surf.obj')

        # Calc inertia tensor for each pocket.
        inertia = mesh.moment_inertia

        # calculate eigen values and eigen vectors.
        eigvals, _ = np.linalg.eig(inertia)

        # sort eigen values > (I1, I2, I3)
        eig_ord = np.argsort(eigvals)

        # Calculate the normalize PMI ratios for each pocket.
        sorted_eigvals = eigvals[eig_ord]
        I1 = sorted_eigvals[0]
        I2 = sorted_eigvals[1]
        I3 = sorted_eigvals[2]
        npr1 = I1/I3
        npr2 = I2/I3
        pocket_npr1_l.append(npr1)
        pocket_npr2_l.append(npr2)

        # Calculate pocketNT (nucleotides near each pocket).
        x = pockets_out[i]
        with open(x, 'r') as f:
            tmp0 = f.read()
            pc_dict = {'Pocket': [(int(x[1]))] for x in findall(
                r'(\w+\s\w+\s*)(\d)(:)(\n)', tmp0)}

            pc_dict.update({x[1].strip(): float(x[2]) for x in findall(
                r'(\w+\s\w+\s*-\s*)(.+):\s*([\d.-]+)(\n)', tmp0)})

        structure = parsePDB(x)

        if ',' in chain:
            chain = chain.replace(',', ' & ')
        selection = structure.select(f'chain {chain}')
        nt = selection.getResnums().tolist()
        pocketNT_l.append(np.unique(nt).tolist())

    # Add pocketNT and pocket npr data to pc dataframe.
    pc_df['PocketNT'] = pocketNT_l
    pc_df['Pocket npr1'] = pocket_npr1_l
    pc_df['Pocket npr2'] = pocket_npr2_l

    # Add pocket filter (Pass or Fail) to pc dataframe.
    pc_df.loc[pc_df['Score'] > qualityfilter, 'Filter'] = 'Pass'


def add_ligand_characteristics(
    analysis : str, stp_coords : prody.AtomGroup,
    ligand_coords : prody.AtomGroup, ligand : str, pc_df : pd.DataFrame
    ) -> None:
    """Adds characteristics to the pocket characteristics DataFrame that
       require the presence of a ligand to calculate.
    Pocket overlap:  Ratio of a-spheres in contact with ligand.
    Ligand overlap:  Ratio of ligand atoms in contact with a-spheres.
    Center criteria: Distance from ligand to the geometric center of pocket.
    Type: Indicates if pocket is known or novel based on:
          pocket overlap, ligand overlap, and center criteria.

    Args:
        analysis (str): path directory contianing fpocket outputs for analysis.
        stp_coords (object): ProDy atom group of a-sphere coordinates.
        ligand_coords (object): ProDy atom group of ligand coordinates.
        ligand (str): Ligand residue name (usually a 3-letter code).
        pc_df (DataFrame): Characteristics and properities for each pocket.
    """

    pocket_overlap_l = []
    ligand_overlap_l = []
    center_criteria_l = []

    # Calculate the normalize PMI ratios for the ligand.
    ligand_npr1, ligand_npr2 = util.calc_npr(ligand_coords)

    #Calculate ligand QED score.
    if not os.path.isfile(f'{analysis}/{ligand}_model.sdf'):
        try:
            response = requests.get(
                f'https://files.rcsb.org/ligands/download/{ligand}_model.sdf')
                
            with open(f'{analysis}/{ligand}_model.sdf', 'wb') as f:
                f.write(response.content)
        except:
            print(f'Error: Not able downlaod ligand model for {ligand}.\n')
            qed = np.nan
    
    try:        
        mol = Chem.MolFromMolFile(f'{analysis}/{ligand}_model.sdf')
        qed = QED.default(mol)
    except:
        print(f'Error: Not able calculate qed score for ligand {ligand}.\n')
        qed = np.nan

    # Iterates through each pocket (residue).
    for _, residue in enumerate(stp_coords.iterResidues()):
        pocket = residue.getCoords()

        # Calculate ratio of a-spheres in each pocket that is near the ligand.
        pocket_overlap_sele = residue.select(
            f'within 3 of ligand', ligand=ligand_coords)
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
            f'within 3 of pocket', pocket=residue)
        if ligand_overlap_sele is None:
            ligand_overlap_l.append(0)
        else:
            # Calculates total number of atoms in ligand.
            ligand_atom_count = ligand_coords.numAtoms()

            # Calculates number of overlapped atoms in ligand.
            ligand_overlap_atom_count = ligand_overlap_sele.numAtoms()

            ligand_overlap_l.append(
                ligand_overlap_atom_count / ligand_atom_count)
        
        # Calculate minimum distance between pocket center and ligand atoms.
        pocket_center = calcCenter(pocket)
        center_distance_l = calcDistance(pocket_center, ligand_coords)
        center_distance = min(center_distance_l)
        center_criteria_l.append(center_distance)

    # Add Ligand overlap and Center criteria to pc dataframe.
    pc_df['Pocket overlap'] = pocket_overlap_l
    pc_df['Ligand overlap'] = ligand_overlap_l
    pc_df['Center criteria'] = center_criteria_l

    # Add pocket type (Known or Novel) to pc_df.
    pc_df.loc[(pc_df['Ligand overlap'] >= 0.33) &
              (pc_df['Pocket overlap'] >= 0.33) &
              (pc_df['Center criteria'] <= 4), 'Type'] = 'Known'

    # Add ligand NPR and QED score to pc_df.
    pc_df['NPR1'] = np.nan
    pc_df['NPR2'] = np.nan
    pc_df['QED score'] = np.nan
    pc_df.loc[(pc_df['Type'] == 'Known'), 'NPR1'] = ligand_npr1
    pc_df.loc[(pc_df['Type'] == 'Known'), 'NPR2'] = ligand_npr2
    pc_df.loc[(pc_df['Type'] == 'Known'), 'Geometry'] = 'Balanced'
    pc_df.loc[pc_df.eval('NPR1 - NPR2 + 0.5 < 0'), 'Geometry'] = 'Rod-like'
    pc_df.loc[pc_df.eval('- NPR1 - NPR2 + 1.5 < 0'), 'Geometry'] = 'Sphere-like'
    pc_df.loc[pc_df.eval('NPR2 - 0.75 < 0'), 'Geometry'] = 'Disc-like'
    pc_df.loc[(pc_df['Type'] == 'Known'), 'QED score'] = qed


