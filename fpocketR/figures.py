#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Module for generating secondary and tertiary structure figures
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2022
#
# Version 1.0.0
#
# -----------------------------------------------------------------------------
import os
from glob import glob
from pylab import *
from prody import *
import rnavigate as rnav
import numpy as np
import seaborn as sns
from fpocketR import make3D


def make_figures(pdb_code, pdb, state, pc_df, rna_coords, nsd, analysis,
                 name, chain, dpi, zoom, offset, connectpocket, alignligand):

    # Get the rna sequnece length from the .pdb file.
    pdb_seq_len = rna_coords.numResidues()

    # Gets secondary structure input file required for making figures.

    try:
        if os.path.isfile(nsd) is False:
            raise Exception

    except TypeError:
        nsd = False

        # If not making 2D figures than get rna sequence length from .pdb file.
        rna_seq_len = pdb_seq_len
        pass

    except Exception:
        print('An valid .nsd file is required to make 2D figures.\n'
              'We suggest generating secondary structure files by:\n'
              '1) Converting .pdb to .ct files at:\n'
              'http://rnapdbee.cs.put.poznan.pl\n'
              '2) Making .nsd from .ct with StructureEditor avalible at:\n'
              'https://rna.urmc.rochester.edu/RNAstructureDownload.html\n')
        nsd = False

        # If not making 2D figures than get rna sequence length from .pdb file.
        rna_seq_len = pdb_seq_len
        pass

    else:
        # Makes RNAvigate object for rna secondary structure.
        rna_map = rnav.Sample(sample=name, ss=nsd)

        # Gets length of RNA. nsd_seq_len needed to make figures.
        nsd_seq_len = rna_map.data["ss"].length

        if pdb_seq_len != nsd_seq_len:
            print(f'WARNING: The PDB sequence length ({pdb_seq_len}) is different '
                  f'than the NSD seqence length ({nsd_seq_len}).\n'
                  f'Make sure these files are correct.\n')

        # Define RNA sequence length based on nsd file.
        rna_seq_len = nsd_seq_len

    # Generates csv output containing pocket characteristics.
    if 'Pocket' in pc_df.columns:
        print(f'\nNumber of pockets detected: {pc_df["Pocket"].max()}')
    if 'Known' in pc_df:
        print(
            f'Number of known pockets: {pc_df["Type"].value_counts()["Known"]}\n')
    pc_df.to_csv(f'{analysis}/{name}_out_pocket_characteristics.csv',
                 index=True, float_format='%.2g')

    # Get color maps to color 2D and 3D figures
    seq_cmap, pocket_cmap, nt_group_color_list = get_colorNT(
        pc_df, rna_seq_len, offset, rna_coords, nsd, chain)

    # Get 2D figures based on secondary structure in .nsd file.
    if nsd:
        get_2D_figure(nsd, seq_cmap, nt_group_color_list,
                      analysis, name, connectpocket)

    # Get 3D PyMol figures based on the real_sphere.pdb file.
    get_3D_figure(pdb, state, analysis, name, dpi,
                  chain, zoom, pocket_cmap, alignligand)
    
    return pocket_cmap

# -----------------------------------------------------------------------------

def get_pdb_offset(rna_coords, nt, offset, chain):
    """Used for for calculating nucleotide offsets for RNA with 2 chains.

    Args:
        rna_coords (object): Prody atom group of RNA coordinates.
        nt (int): Index of the nucleotide of interest for calculating an offset
        offset (int): Sequence offset between .pdb and .nsd file (1st chain).
        chain (str): Chain identifier for RNA chain.

    Returns:
        int: 2D structure index for the specified nucleotide
    """
    second_chain = str(chain.split(",")[1])
    offset2NT = rna_coords.select(
        f'chain {second_chain}').copy().getResnums()[0]
    first_chain = str(chain.split(",")[0])
    offsetNT = rna_coords.select(
        f'chain {first_chain}').copy().getResnums()[-1]
    offset2 = offset2NT - offsetNT - 1

    if nt < offset2NT:
        cmap_idx = nt - 1 - offset
    elif nt >= offset2NT:
        cmap_idx = nt - 1 - offset - offset2

    return cmap_idx


def get_colorNT(pc_df, rna_seq_len, offset, rna_coords, nsd, chain):
    """Generates a color map of pocket loations throughout the RNA.

    Args:
        pc_df (dataframe): Contains all the characteristics for each pocket.
        rna_seq_len (int): Number of nucleotides in RNA sequence.
        offset (int): Sequence offset between .pdb and .nsd file (1st chain).
        rna_coords (object): Prody atom group of RNA coordinates.
        nsd (str): Path to input secondary structure drawing.
        chain (str): Chain identifier for RNA chain.

    Returns:
        list(tuple): Per nucleotide color map. Index = NT sequence. Value = color.
        list(tuple): Per pocket color map. Index = pocket residue number. Value = color.
    """

    seq_cmap = [(1, 1, 1, 1)] * rna_seq_len
    pocket_cmap = {}
    nt_group_color_list = []

    # Cubehelix color map for known nucleotides/pockets.
    known_cmap = sns.color_palette(
        'ch: 1.1, rot=0.1, gamma=0.6, light=0.65, dark=0.35, hue=3, reverse=1',
        as_cmap=True)

    # Cubehelix color map for novel nucleotides/pockets.
    novel_cmap = sns.color_palette(
        'ch: 0.09, rot=2.55, gamma=0.85, light=0.45, dark=0.25, hue=1.25, reverse=1',
        as_cmap=True)

    # Creates DataFrame containing only pockets that Pass the quality filter.
    pocket_df = pc_df[pc_df['Filter'] == 'Pass']
    if len(pocket_df) is None:
        return None, None

    if pocket_df[pocket_df['Type'] == 'Known'] is None:
        known_len = 0
    else:
        known_len = len(pocket_df[pocket_df['Type'] == 'Known'])

    if pocket_df[pocket_df['Type'] == 'Novel'] is None:
        novel_len = 0
    else:
        novel_len = len(pocket_df[pocket_df['Type'] == 'Novel'])

    cmaps = {'Known': known_cmap, 'Novel': novel_cmap}
    lengths = {'Known': known_len - 1, 'Novel': novel_len - 1}
    counts = {'Known': 0, 'Novel': 0}

    colored_nt_list = []
    known_colored_nt_list = []

    for idx, row in pocket_df.iterrows():
        pocketNT = row['PocketNT']
        poc_type = str(row['Type'])
        poc_num = int(row['Pocket'])

        # Gets equally spaced colors for each pocket from color map.
        if lengths[poc_type] == 0:
            color = cmaps[poc_type](0)
        else:
            color = cmaps[poc_type](counts[poc_type] / lengths[poc_type])

        # Assigns color to each pocket.
        pocket_cmap[poc_num] = color
        counts[poc_type] += 1

        # Colors nucleotides with matching color to their associated pocket.
        # If a single nucleotide is in contact with multiple pockets:
        # -Priority is given to higher ranked pockets.
        # -Color of a known pocket will not be overwritten.
        cmap_idx_per_pocket = []
        for nt in pocketNT:
            # Offset used to match index of .pdb input to .nsd output.
            if len(chain) == 1:
                cmap_idx = nt - 1 - offset

            # If the RNA contains 2 chains get the second offset.
            elif len(chain.split(',')) == 2:
                cmap_idx = get_pdb_offset(rna_coords, nt, offset, chain)
            else:
                print('WARNING: 2D figures can only be made from RNA strctures \n'
                      'composed of a maxiumum of 2 chains (-c).')
                nsd = None
                break

            cmap_idx_per_pocket.append(cmap_idx)

            if nsd:
                if nt not in colored_nt_list and \
                        nt not in known_colored_nt_list:
                    seq_cmap[cmap_idx] = color
                elif poc_type == 'Known' and \
                        nt not in known_colored_nt_list:
                    seq_cmap[cmap_idx] = color

            colored_nt_list.append(nt)
            if poc_type == 'Known':
                known_colored_nt_list.append(nt)

        nt_group_color_list.append(
            {'sites': cmap_idx_per_pocket, 'color': color})

    return seq_cmap, pocket_cmap, nt_group_color_list


def get_2D_figure(nsd, seq_cmap, nt_group_color_list, analysis, name, connectpocket):
    """Uses RNAvigate to plot pocketNT onto the RNA secondary structure and
       saves resulting 2D figure as png and svg.

    Args:
        nsd (str): Path to input secondary structure drawing.
        pdb_ligand (str): Path to input .pdb file.
        chain (str): Indicated desired PDB chain (default='A').
        offset (int): Sequence offset between .pdb and .nsd file (default=0).
        seq_cmap list(tuple): RGB color codes for each nucleotide.
        analysis (str): Path to the analysis output directory.
        name (str): Output file name prefix (default=pdb_name).
        connectpocket (boolean): Connects pockets in 2D figure (Default=False).
    """
    print('Making 2D figure.\n')

    rna_map = rnav.Sample(sample=name, ss=nsd)

    # Makes 2D figures with transparent colored lines
    # which connect all the nucleotides associated with a pocket.
    if connectpocket:
        rna_map.set_data(name="fpocket", filepath=None,
                         instantiator=rnav.data.Annotation, seq_source="ss",
                         groups=nt_group_color_list)

        plot = rnav.plot_ss([rna_map], annotations=["fpocket"], colors=seq_cmap,
                            apply_color_to='background', sequence=True,
                            bp_style="line")

    # Makes 2D figures without connecting nucelotides.
    else:
        plot = rnav.plot_ss([rna_map], colors=seq_cmap, apply_color_to='background',
                            sequence=True, bp_style="line")

    ax = plot.axes[0, 0]
    x0, x1, y0, y1 = ax.axis()
    ax.axis((x0-1, x1+1, y0-1, y1+1))
    plot.set_figure_size()

    # Adds colored backgroud to nucleotides associated with the same pockets.
    colored_nts = [c != 'w' for c in seq_cmap]
    x = rna_map.data["ss"].xcoordinates[colored_nts]
    y = rna_map.data["ss"].ycoordinates[colored_nts]
    seq_cmap_2 = np.array(seq_cmap, dtype=object)[colored_nts]
    pts = 20
    plot.axes[0, 0].scatter(x, y, c=seq_cmap_2, s=pts**2, zorder=19)

    # Saves 2D figure was .png file.
    plt.savefig(f'{analysis}/{name}_2D.png', dpi=300, format=None,
                metadata=None, bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto', backend=None)
    # Saves 2D figure as .svg (editable) file.
    plt.savefig(f'{analysis}/{name}_2D.svg', dpi=300, format=None,
                metadata=None, bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto', backend=None)


def get_3D_figure(pdb, state, analysis, name, dpi,
                  chain, zoom, pocket_cmap, alignligand):
    print('Making 3D figure.\n')
    real_sphere_name = f'{name}_out_real_sphere'
    real_sphere_pdb = os.path.join(analysis, f'{real_sphere_name}.pdb')
    make3D.load_pdb(real_sphere_pdb)
    if alignligand:
        make3D.alignligand(pdb, real_sphere_name, state)
    make3D.set_default()
    if pocket_cmap:
        make3D.color_pockets(pocket_cmap)
    make3D.save_3D_figure(analysis, name, dpi, chain, zoom)


def get_all_states_3D_figure(pdb_out, pdb_code, pocket_cmap, dpi, chain, zoom):
    """Generates a 3D figure and pymol session file with all states.

    Args:
        pdb_out (str): File name for .pdb file.
        name_out (str): PDB identification code for output .pdb file.
        pocket_cmap list(tuple): Per pocket color map. Index = pocket residue number. Value = color.
        dpi (int): Figure resolution in dpi (dots per linear inch).
        chain (str): Chain identifier for desired RNA chain.
        zoom (float): Zoom buffer distance (Å) for creating 3D figures.
    """
    real_sphere_list = glob(f'{pdb_out}/*out/*out_real_sphere.pdb')
    for real_sphere in real_sphere_list:
        object_name = os.path.basename(real_sphere)[:-4]
        make3D.load_pdb(real_sphere)
        make3D.transparent_pocket(object_name, pocket_cmap)

    make3D.set_default()
    make3D.make_multistate()
    make3D.save_3D_figure(pdb_out, f'{pdb_code}_all', dpi, chain, zoom)
