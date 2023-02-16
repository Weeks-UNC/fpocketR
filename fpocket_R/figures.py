#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Module for generating secondary and tertiary structure figures
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2022
#
# Version 0.1.0
#
# -----------------------------------------------------------------------------


import rnavigate as MaP
import numpy as np
import seaborn as sns
from pylab import *
from prody import *
import make3D
import os


def make_figures(pdb_code, pdb, s, pc_df, rna_coords, nsd, analysis,
                 pqr_out, pdb_out, name, c, dpi, zoom, offset,
                 connectpocket, alignligand):

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
        rna_map = MaP.Sample(sample=name, ss=nsd)

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
        pc_df, rna_seq_len, offset, rna_coords, nsd, c)

    # Get 2D figures based on secondary structure in .nsd file.
    if nsd:
        get_2D_figure(nsd, seq_cmap, nt_group_color_list,
                      analysis, name, connectpocket)

    # Creates real_sphere.pdb.
    # Encodes a-sphere radii into pdb b-factor column.
    pqr_to_pdb(pqr_out, pdb_out, analysis, name)

    # Get 3D PyMol figures based on the real_sphere.pdb file.
    get_3D_figure(pdb_code, pdb, s, analysis, name, dpi,
                  c, zoom, pocket_cmap, alignligand)


def pqr_to_pdb(pqr_file, pdb_file, analysis, name):
    """Encodes fpocket pocket a-sphere radii into a single .pdb output file.
    a-sphere radii incoded into the B factor column of the output .pdb file.

    Args:
        pqr_file (string): Path to fpocket *out.pqr file
                           containing a-spheres radii.
        pdb_file (string): Path to fpocket *out.pdb file.
    """

    # Opens .pqr and .pdb inputs and .pdb output.
    with open(pqr_file, 'r') as pqr_file, open(pdb_file, 'r') as pdb_file, \
            open(f'{analysis}/{name}_out_real_sphere.pdb', 'a') as out:

        pqr_line = pqr_file.readline()
        pdb_lines = pdb_file.readlines()

        # Iterates through Header lines in pqr file
        while pqr_line.startswith("HEADER"):
            pqr_line = pqr_file.readline()

        # Iternates through each line of the pdb file.
        for pdb_line in pdb_lines:
            if pdb_line.startswith("HEADER") or pdb_line.startswith("ATOM"):
                line_str = pdb_line[0:80]
            elif pdb_line.startswith("HETATM") and 'APOL' not in pdb_line:
                line_str = pdb_line[0:80]
            elif pdb_line.startswith("HETATM") and \
                    pdb_line[24:55] == pqr_line[24:55]:
                line_str = pdb_line[0:60] + \
                    pqr_line[65:71] + pdb_line[66:80]
                # Advance to next line of pqr file.
                pqr_line = pqr_file.readline()
            else:
                print(f"ERROR: pqr_to_pdb. \n"
                      f"PDB input file: {pdb_file} \n"
                      f"PQR input file: {pqr_file}")
                break

            # Write line to out_real_sphere.pdb file
            out.write(line_str)
            out.write('\n')


def get_pdb_offset(rna_coords, nt, offset, c):
    """Used for for calculating nucleotide offsets for RNA with 2 chains.

    Args:
        rna_coords (object): prody atom group of RNA coordinates
        nt (int): index of the nucleotide of interest for calculating an offset
        offset (int): idex offset of the first chain
        c (str): specified chains to analyze

    Returns:
        int: 2D structure index for the specified nucleotide
    """
    second_chain = str(c.split(",")[1])
    offset2NT = rna_coords.select(
        f'chain {second_chain}').copy().getResnums()[0]
    first_chain = str(c.split(",")[0])
    offsetNT = rna_coords.select(
        f'chain {first_chain}').copy().getResnums()[-1]
    offset2 = offset2NT - offsetNT - 1

    if nt < offset2NT:
        cmap_idx = nt - 1 - offset
    elif nt >= offset2NT:
        cmap_idx = nt - 1 - offset - offset2

    return cmap_idx


def get_colorNT(pc_df, rna_seq_len, offset, rna_coords, nsd, c):
    """Generates a color map of pocket loations throughout the RNA.

    Args:
        pc_df (dataframe): Contains all the characteristics for each pocket.
        nsd_seq_len (int): Sequence length of the .nsd 2D structure file.

    Returns:
        list: Per nucleotide color map. Index = NT sequence. Value = color.
    """

    seq_cmap = [(1,1,1,1)] * rna_seq_len
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
            if len(c) == 1:
                cmap_idx = nt - 1 - offset

            # If the RNA contains 2 chains get the second offset.
            elif len(c.split(',')) == 2:
                cmap_idx = get_pdb_offset(rna_coords, nt, offset, c)
            else:
                print('WARNING: 2D figures can only be made from RNA strctures \n'
                      'composed of a maxiumum of 2 chains (-c).')
                nsd = None

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
        nsd (str)): Path to input secondary structure drawing.
        pdb_ligand (str): Path to input .pdb file.
        c (str): Indicated desired PDB chain (default='A').
        offset (int): Sequence offset between .pdb and .nsd file (default=0).
        seq_cmap (list of tuples): RGB color codes for each nucleotide.
        analysis (str): Path to the analysis output directory.
        name (str): Output file name prefix (default=pdb_name).
    """
    print('Making 2D figure.\n')

    rna_map = MaP.Sample(sample=name, ss=nsd)

    # Makes 2D figures with transparent colored lines
    # which connect all the nucleotides associated with a pocket.
    if connectpocket:
        rna_map.set_data(name="fpocket", filepath=None,
                         instantiator=MaP.data.Annotation, seq_source="ss",
                         groups=nt_group_color_list)

        plot = rna_map.plot_ss(annotations=["fpocket"], colors=seq_cmap,
                               apply_color_to='background', sequence=True,
                               bp_style="line")

    # Makes 2D figures without connecting nucelotides.
    else:
        plot = rna_map.plot_ss(colors=seq_cmap, apply_color_to='background',
                               sequence=True, bp_style="line")

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


def get_3D_figure(pdb_code, pdb, s, analysis, name, dpi,
                  c, zoom, pocket_cmap, alignligand):
    print('Making 3D figure.\n')
    real_sphere_name = f'{name}_out_real_sphere'
    real_sphere_pdb = os.path.join(analysis, f'{real_sphere_name}.pdb')
    # cmd.load(f'{real_sphere_pdb}')
    make3D.load_pdb(real_sphere_pdb)
    if alignligand:
        make3D.alignligand(pdb, real_sphere_name, s)
    make3D.set_default()
    if pocket_cmap:
        make3D.color_pockets(pocket_cmap)
    make3D.save_3D_figure(analysis, name, dpi, c, zoom)
