#!/usr/bin/env python3
# -----------------------------------------------------
# fpocket analysis code
# Seth Veenbaas, Patrick Irving
# Weeks Lab, UNC-CH
# 2-9-2023
#
# Version 1.0.0
#
# -----------------------------------------------------
import argparse
import pandas as pd
from pymol import cmd
from prody import *
from fpocket_R import (analyze, pocket, figures, util)
confProDy(verbosity='none')
# -----------------------------------------------------


def pipline(pdb, nsd, chain, state, ligandchain, offset, qualityfilter,
            m, M, i, D, A, p, out, name, dpi, ligand, yes, zoom,
            connectpocket, alignligand):
    """Runs pocket finding pipeline

    Args:
        pdb (str): Path to input .pdb file.
        nsd (str)): Path to input secondary structure drawing.
        chain (str): Chain identifier for desired RNA chain (default='A').
        state (int): Structural state to analyze.
        ligandchain (str): Chain identifier for desired ligand (default=chain).
        offset (int): Sequence offset between .pdb and .nsd file (default=0).
        qualityfilter (float): Minimum fpocket score filter for pockets.
        m (float): Min. a-sphere radius in angstroms (default=3.0).
        M (float): Max. a-sphere radius in angstroms (default=5.7).
        i (int): Min. number of a-spheres per pocket (default=42).
        D (float): a-sphere clustering distance in angstroms (default=1.65).
        A (int): # of electroneg. atoms to define a polar a-sphere (default=3).
        p (float): Max. ratio of apolar a-spheres in a pocket (default=0.0).
        out (str): name of fpocket output parent directory name.
        name (str): Output file name prefix (default={pdb_name}).
        dpi (int): Figure resolution in dpi (default=300).
        ligand (str): Ligand residue name (usually a 3-letter code).
        yes (boolean): Overwrite output files and directories with same name.
        zoom (float): Zoom buffer distance (Å) for creating 3D figures.
        connectpocket (boolean): Connects pockets in 2D figure (Default=False).
        alignligand (boolean): Align ligand to pymol output (Default=True).

    Returns:
        str: Path to clean .pdb input file.
        str: Path directory contianing fpocket outputs for analysis.
        object: Pandas Dataframe with characteristics for each pocket.

    """

    # Check if pdb contains a file extension.
    if len(pdb.split('.')) < 2:
        # Fetch PDB ID from the PDB. Fetchs .pdb files then .cif.
        pdb_id_lower = pdb.lower()
        filename = fetchPDBs(f'{pdb_id_lower}', compressed=False, quiet=True)
        pdb = filename[0]

    # Checks if required input files are accessible/exist.
    print('Checking input files.')
    util.is_accessible(pdb, 'pdb')
    util.contains_structure(pdb, chain)

    # Runs fpocket on input pdb file and manages output files.
    analysis = pocket.find_pockets(
        pdb, chain, state, m, M, i, D, A, p, out, yes).strip('/')

    # Checks if the analysis directory is accessible.
    util.is_accessible(analysis, 'analysis directory')

    # Get paths to fpocket input and output file.
    (pdb, pdb_out, pqr_out, info_txt, pockets_out,
        pdb_code, name) = util.get_file_paths(analysis, name, pdb, state)

    # Analyze fpocket data and create pocket characteristics dataframe.
    (pc_df, rna_coords, offset) = analyze.analyze_pockets(
                                  pdb, pqr_out, pdb_out, analysis, name,
                                  info_txt, pockets_out, pdb_code, chain,
                                  state, ligandchain, ligand, m, M, i, D, A, p,
                                  qualityfilter, offset)

    # Generates 1D (.csv), 2D (.png, .svg), and 3D (.pdb, .pse, .png)
    figures.make_figures(pdb_code, pdb, state, pc_df, rna_coords, nsd, 
                         analysis, name, chain, dpi, zoom, offset, 
                         connectpocket, alignligand)

    return pc_df, out, pdb_code, name


# -----------------------------------------------------
def parseArgs():
    prs = argparse.ArgumentParser()

    prs.add_argument('-pdb', '--pdb', type=str, required=True,
                     help='Specify a .pdb file to run fpocket on before '
                     'starting analysis.')
    prs.add_argument('-nsd', '--nsd', type=str, required=False,
                     help='Specify an .nsd file for generating secondary '
                     'structure figures.')
    prs.add_argument('-c', '--chain', type=str, required=False, default='A',
                     help='Specify a chain from the input .pdb file ("A").')
    prs.add_argument('-s', '--state', type=int, required=False, default=None,
                     help='Specify which NMR states/model '
                     'you would like to analyze. Set to 0 for all (None).')
    prs.add_argument('-lc', '--ligandchain', type=str, required=False,
                     help='Chain containing ligand the from the '
                     'input .pdb file (--chain input).')
    prs.add_argument('-off', '--offset', type=int, required=False,
                     help='Specify an offset between the rna sequence and '
                     'starting nucleotide of the PDB structure.')
    prs.add_argument('-qf', '--qualityfilter', type=float, default=0.0, 
                     required=False, help='Specify minimum fpocket score for '
                     'pocket to pass the quality filter (0.0).')
    prs.add_argument('-m', type=float, default=3.00, required=False,
                     help='fpocket -m flag. Specifies the minimum radius '
                     'for an a-sphere (3.0).')
    prs.add_argument('-M', '--M', type=float, default=5.70, required=False,
                     help='fpocket -M flag. Specifies the maximium radius '
                     'for an a-sphere (5.70).')
    prs.add_argument('-i', '--i', type=int, default=42, required=False,
                     help='fpocket -i flag. Specifies the minimum number '
                     'of a-spheres per pocket (42).')
    prs.add_argument('-D', '--D', type=float, default=1.65, required=False,
                     help='fpocket -D flag. Specifies the a-sphere '
                     'clustering distance for forming pockets (1.65).')
    prs.add_argument('-A', '--A', type=int, default=3, required=False,
                     help='fpocket -A flag. Number of electronegative atoms '
                     'required to define a polar a-sphere (3).')
    prs.add_argument('-p', '--p', type=float, default=0, required=False,
                     help='fpocket -p flag. Maximum ratio of apolar a-spheres '
                     'in a pocket (0.0).')
    prs.add_argument('-o', '--out', type=str, required=False,
                     help='Specify name of fpocket output '
                     'parent directory name.')
    prs.add_argument('-n', '--name', type=str, required=False,
                     help='Specify name prefix for fpocket_out and '
                     'analysis_out subdirectories.')
    prs.add_argument('-dpi', '--dpi', type=int, default=300, required=False,
                     help='Sets figure resolution in dpi (300).')
    prs.add_argument('-df', '--distancefilter', type=float, default=4.5,
                     help='Distance filter (in Angstroms) for identifying '
                     'ligands close in space to pockets (4.5).')
    prs.add_argument('-l', '--ligand', type=str,
                     help='Three character residue name of desired ligand.')
    prs.add_argument('-y', '--yes', action='store_true', default=False,
                     help='Answers yes to user prompts for overwriting files.')
    prs.add_argument('-z', '--zoom', type=float, default=5.0, required=False,
                     help='Zoom buffer (Å) for creating 3D figures (5.0).')
    prs.add_argument('-cp', '--connectpocket', action='store_true',
                     default=False, help='Visually connects pockets in 2D '
                     'figures (False).')
    prs.add_argument('-al', '--alignligand', action='store_false',
                     default=True, help='Aligns ligand to the output structure '
                     'containing pocket predictions (True).')

    args = prs.parse_args()
    return args


def main(pdb, nsd, chain, state, ligandchain, offset, qualityfilter,
         m, M, i, D, A, p, out, name, dpi, distancefilter,
         ligand, yes, zoom, connectpocket, alignligand):
    """Runs the fpocket analysis pipeline.
    Pipeline runs once: by default or if provided a user specified state
    is specified using the -s flag.

    Pipeline runs multiple times: if the -s flag is set to 0 (all).
    This feature is intended for analyzing NMR structures
    with several modeled states.
    """
    # Runs pipeline for a single state of the input structure.
    if state != 0:
        if out is None:
            out = f'fpocket-R_out-m_{m}-M_{M}-i_{i}-D_{D}-A_{A}-p_{p}'
        (pc_df, out_all, pdb_code_csv, name_all) = pipline(
            pdb, nsd, chain, state, ligandchain, offset, qualityfilter,
            m, M, i, D, A, p, out, name, dpi, ligand, yes, zoom,
            connectpocket, alignligand)

    # Runs pipeline for multiple states of the input structure.
    elif state == 0:
        if out is None:
            out = f'Multistate_{pdb.split(".")[0]}'

        try:
            num_states = parsePDBHeader(pdb, 'n_models')+1
        except:
            print('ERROR: Unable to perform multisate analysis.\n'
                  f'The header for {pdb} does not contain state information.\n')
            exit() 

        pc_all_states_df = pd.DataFrame()
        for state in range(1, num_states):
            print(f'\nAnalyzing state {state}/{num_states-1}...\n')

            (pc_df, out_all, pdb_code_csv, name_all) = pipline(
                pdb, nsd, chain, state, ligandchain, offset, qualityfilter,
                m, M, i, D, A, p, out, name, dpi, ligand, yes, zoom, 
                connectpocket, alignligand)

            pc_all_states_df = pd.concat([pc_all_states_df, pc_df])

        # Generates csv output containing pocket characteristics for all states.
        pc_all_states_df.to_csv(
            f'{out_all}/{pdb_code_csv}_all_states_pocket_characteristics.csv',
            index=True, float_format='%.2g')

        # Generates a 3D for pockets in all states.
        figures.get_all_states_3D_figure(out_all, name_all, dpi, chain, zoom)

    # Close pymol session.
    cmd.quit()


if __name__ == "__main__":
    main(**vars(parseArgs()))
