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
import os
import pandas as pd
from pymol import cmd
from prody import *
from fpocketR import (analyze, pocket, figures, util)
confProDy(verbosity='none')
# -----------------------------------------------------


def pipeline(
    pdb : str,
    nsd : str,
    chain : str,
    state : int,
    ligand : str,
    ligandchain : str,
    knownnt : list[int],
    offset : int,
    qualityfilter : float,
    m : float,
    M : float,
    i : int,
    D : float,
    A : float,
    p : float,
    out : str,
    name : str,
    dpi : int,
    yes : bool,
    zoom : float,
    connectpocket : bool,
    alignligand : str,
):   
    """Runs pocket finding pipeline

    Args:
        pdb (str): Path to input .pdb file.
        nsd (str)): Path to input secondary structure drawing.
        chain (str): Chain identifier for desired RNA chain (default='A').
        state (int): Structural state to analyze.
        ligandchain (str): Chain identifier for desired ligand (default=chain).
        offset (int): Sequence offset between .pdb and .nsd file (default=None).
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
        alignligand (str): Align ligand to pymol output (Default=True).

    Returns:
        str: Path to clean .pdb input file.
        str: Path directory contianing fpocket outputs for analysis.
        object: Pandas Dataframe with characteristics for each pocket.

    """

    # Checks if required input files are accessible/exist.
    print('Checking input files.')
    util.is_accessible(pdb, 'pdb')

    if chain is None:
        chain = util.get_first_rna_chain(pdb)
    else:
        util.is_rna_chain(pdb, chain)

    # Runs fpocket on input pdb file and manages output files.
    analysis, yes = pocket.find_pockets(
        pdb,
        chain,
        state,
        m,
        M,
        i,
        D,
        A,
        p,
        out,
        name,
        yes,
    )

    # Checks if the analysis directory is accessible.
    util.is_accessible(analysis, 'analysis directory')

    # Get paths to fpocket input and output file.
    (pdb, pdb_out, pqr_out, info_txt, pockets_out,
        pdb_code, name) = util.get_file_paths(analysis, name, pdb, state)

    # Analyze fpocket data and create pocket characteristics dataframe.
    (pc_df, rna_coords) = analyze.analyze_pockets(
        pdb,
        pqr_out,
        pdb_out,
        analysis,
        name,
        info_txt,
        pockets_out,
        pdb_code,
        chain,
        state,
        ligandchain,
        ligand,
        m,
        M,
        i,
        D,
        A,
        p,
        qualityfilter,
        knownnt,
    )
    
    offset = util.get_offset(pdb, chain, offset) if offset is None else offset

    # Generates 1D (.csv), 2D (.png, .svg), and 3D (.pdb, .pse, .png)
    pocket_cmap = figures.make_figures(
        pdb,
        state,
        pc_df,
        rna_coords,
        nsd,
        analysis,
        name,
        chain,
        dpi,
        zoom,
        offset,
        connectpocket,
        alignligand
    )

    return pc_df, out, pocket_cmap, chain, yes


# -----------------------------------------------------
def parseArgs():
    prs = argparse.ArgumentParser()
    pocket_type = prs.add_mutually_exclusive_group()

# Input options
    prs.add_argument(
        '-pdb',
        '--pdb',
        type=str,
        required=True,
        default=None,
        help='Path to a .pdb file, .cif file, or 4 charater PDB '
        'indentification code.',
    )
    prs.add_argument(
        '-nsd',
        '--nsd',
        type=str,
        required=False,
        default=None,
        help='Path to an .nsd or other secondary structure file '
        'for generating secondary structure figures.',
    )
    
# fpocket parameter options
    prs.add_argument(
        '-m',
        type=float,
        required=False,
        default=3.00,
        help='fpocket -m flag. Specifies the minimum radius '
        'for an a-sphere (3.0).',
    )
    prs.add_argument(
        '-M',
        '--M',
        type=float,
        required=False,
        default=5.70,
        help='fpocket -M flag. Specifies the maximium radius '
        'for an a-sphere (5.70).',
    )
    prs.add_argument(
        '-i',
        '--i',
        type=int,
        required=False,
        default=42,
        help='fpocket -i flag. Specifies the minimum number '
        'of a-spheres per pocket (42).',
    )
    prs.add_argument(
        '-D',
        '--D',
        type=float,
        required=False,
        default=1.65,
        help='fpocket -D flag. Specifies the a-sphere '
        'clustering distance for forming pockets (1.65).'
    )
    prs.add_argument(
        '-A',
        '--A',
        type=int,
        required=False,
        default=3,
        help='fpocket -A flag. Number of electronegative atoms '
        'required to define a polar a-sphere (3).',
    )
    prs.add_argument(
        '-p',
        '--p',
        type=float,
        required=False,
        default=0.0,
        help='fpocket -p flag. Maximum ratio of apolar a-spheres '
        'in a pocket (0.0).',
    )
    
# Output options
    prs.add_argument(
        '-o',
        '--out',
        type=str,
        required=False,
        default=None,
        help='Specify name of fpocket output parent directory name.',
    )
    prs.add_argument(
        '-n',
        '--name',
        type=str,
        required=False,
        default=None,
        help='Specify name prefix for fpocket_out and analysis_out subdirectories.',
    )
    prs.add_argument(
        '-y',
        '--yes',
        required=False,
        action='store_true',
        help='Answers yes to user prompts for overwriting files (False).',
    )

# Analysis optiopns
    prs.add_argument(
        '-s',
        '--state',
        type=int,
        required=False,
        default=None,
        help='Specify the NMR states/model you would like to analyze. 0 for all (None).',
    )
    prs.add_argument(
        '-c',
        '--chain',
        type=str,
        required=False,
        default=None,
        help='Specify a chain from the input .pdb file ("A").',
    )
    pocket_type.add_argument(
        '-l',
        '--ligand',
        type=str,
        required=False,
        default=None,
        help='PDB ligand identification code (≤ 3 characters).',
    )
    prs.add_argument(
        '-lc',
        '--ligandchain',
        type=str,
        required=False,
        default=None,
        help='Chain containing ligand the from the input .pdb file (--chain input).',
    )
    pocket_type.add_argument(
        '-nt',
        '--knownnt',
        type=parse_int,
        required=False,
        default=None,
        help='List residue ID of nts in known pocket (e.g. 1,2,3) (None).',
    )
    prs.add_argument(
        '-off',
        '--offset',
        type=int,
        required=False,
        default=None,
        help='Specify an offset between the '
        'starting nucleotide of the rna sequence and '
        'starting nucleotide of the PDB structure (usually = 0).\n'
        'You must manually enter a nucleotide offset if:\n'
        '1) Your PDB does not have a DBREF header.\n'
        '2) You are inputting a .cif file.\n'
        'offset (usually) = starting index of the PDB sequence - 1\n'
        'offset is typically: 0',
    )
    prs.add_argument(
        '-qf',
        '--qualityfilter',
        type=float,
        required=False,
        default=0.0,
        help='Specify minimum fpocket score for pocket (0.0).',
    )
    
# Figure options
    prs.add_argument(
        '-dpi',
        '--dpi',
        type=int,
        required=False,
        default=300,
        help='Sets figure resolution in dpi (300).',
    )
    prs.add_argument(
        '-z',
        '--zoom',
        type=float,
        required=False,
        default=5.0,
        help='Zoom buffer (Å) for creating 3D figures (5.0).',
    )
    prs.add_argument(
        '-cp',
        '--connectpocket',
        required=False,
        action='store_true',
        help='Visually connects pockets in 2D figures (False).',
    )
    prs.add_argument(
        '-al',
        '--alignligand',
        type=str,
        required=False,
        # action='store_false',
        help='Aligns ligand to the output structure with pockets (input pdb).',
    )

    args = prs.parse_args()
    return args

def parse_int(string : str) -> list[int]:
    try:
        # Split the string by commas and convert each part to an integer
        return [int(item) for item in string.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError("List of integers expected. Example: '1,2,3,4'")


def main(
    pdb : str,
    nsd : str,
    chain : str,
    state : int,
    ligand : str,
    ligandchain : str,
    knownnt : list[int],
    offset : int,
    qualityfilter : float,
    m : float,
    M : float,
    i : int,
    D : float,
    A : float,
    p : float,
    out : str,
    name : str,
    dpi : int,
    yes : bool,
    zoom : float,
    connectpocket : bool,
    alignligand : str,
):
    """Runs the fpocket analysis pipeline.
    Pipeline runs once: by default or if provided a user specified state
    is specified using the -s flag.

    Pipeline runs multiple times: if the -s flag is set to 0 (all).
    This feature is intended for analyzing NMR structures
    with several modeled states.
    """

    # Check if pdb contains a file extension.
    if len(pdb.split('.')) < 2:
        pdb = util.fetch_pdb(pdb)
    
    # Set pdb name
    if name is None:
        name = os.path.basename(pdb)[0:-4]

    if alignligand is None or True:
        alignligand = pdb
    elif len(alignligand.split('.')) < 2:
        alignligand = util.fetch_pdb(pdb)
    elif not os.path.isfile(alignligand):
        alignligand = None

    # Runs pipeline for a single state of the input structure.
    if state != 0:
        if out is None:
            out = f'fpocketR_out-m_{m}-M_{M}-i_{i}-D_{D}-A_{A}-p_{p}'
        (_, _, _, _, _) = pipeline(
            pdb,
            nsd,
            chain,
            state,
            ligand,
            ligandchain,
            knownnt,
            offset,
            qualityfilter,
            m,
            M,
            i,
            D,
            A,
            p,
            out,
            name,
            dpi,
            yes,
            zoom,
            connectpocket,
            alignligand,
        )

    # Runs pipeline for multiple states of the input structure.
    elif state == 0:
        if out is None:
            out = f'Multistate_{pdb.split(".")[0]}'

        try:
            structure = parsePDB(pdb)
            num_states = structure.numCoordsets()
        except:
            print('ERROR: Unable to perform multisate analysis.\n'
                  f'The header for {pdb} does not contain state information.\n')
            exit() 

        pc_all_states_df = pd.DataFrame()
        multistate_pocket_cmap = {}
        for state in range(1, num_states+1):
            print(f'\nAnalyzing state {state}/{num_states}...\n')

            (pc_df, pdb_out, pocket_cmap, chain, yes) = pipeline(
                pdb,
                nsd,
                chain,
                state,
                ligand,
                ligandchain,
                knownnt,
                offset,
                qualityfilter,
                m,
                M,
                i,
                D,
                A,
                p,
                out,
                name,
                dpi,
                yes,
                zoom,
                connectpocket,
                alignligand,
            )
            yes = yes
            pc_all_states_df = pd.concat([pc_all_states_df, pc_df])
            multistate_pocket_cmap[state]=pocket_cmap


        # Generates csv output containing pocket characteristics for all states.
        pc_all_states_df.to_csv(
            f'{pdb_out}/{name}_all_states_pocket_characteristics.csv',
            index=True, float_format='%.2g')

        # Generates a 3D for pockets in all states.
        print(f'Making all states 3D figure...')
        figures.get_all_states_3D_figure(
            num_states,
            pdb_out,
            name,
            multistate_pocket_cmap,
            alignligand,
            dpi,
            chain,
            zoom
        )

    # Close pymol session.
    cmd.quit()


if __name__ == "__main__":
    main(**vars(parseArgs()))
