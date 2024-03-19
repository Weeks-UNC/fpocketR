#!/usr/bin/env python3
# -----------------------------------------------------
# pymol figure creation code
# Seth Veenbaas
# Weeks Lab, UNC-CH
# 2022
#
# Version 1.0.1
#
# -----------------------------------------------------
import os
from pymol import util
from pymol import cmd


def load_pdb(pdb: str) -> None:
    """load the input pdb file into a pymol session

    Args:
        pdb (str): path to pdb file
    """
    cmd.load(f'{pdb}', partial=1)


def alignligand(pdb: str, real_sphere_name: str, state: int) -> None:
    """aligns state specific structure in the output real_sphere.pdb file
    to the structure in the input pdb file.

    Args:
        pdb (str): Path to RNA-ligand complex .pdb file (target_object).
        real_sphere_name (str): Path to output real_sphere.pdb (mobile_object).
        state (int): State of structure to align.
    """
    pdb_base = os.path.basename(pdb)
    pdb_selection_name = os.path.splitext(pdb_base)[0]
    if not state:
        cmd.load(pdb)

        # aligns mobile object to target object
        cmd.align(f'{real_sphere_name}', f'{pdb_selection_name}')
        cmd.disable(f'{pdb_selection_name}')

    # Align to the state from the input pdb file that matches
    # the state of the output real_sphere.pdb file.
    else:
        cmd.load(pdb, multiplex='1')
        state_str = str(state)
        state_id = state_str.zfill(4)
        state_object_name = f'{pdb_selection_name}_{state_id}'
        cmd.delete(f'not {state_object_name} and not {real_sphere_name}')

        # aligns mobile object to target object
        cmd.align(f'{real_sphere_name}', f'{state_object_name}')
        cmd.disable(f'{state_object_name}')


def set_default() -> None:
    """Set default pymol settings in current pymol session:
    -raytracing performance
    -structure style

    Settings are from the pymolrc.pml file at:
    https://github.com/Weeks-UNC/small-scripts/tree/master/Pymol
    """
    # enable multi-thread processing
    cmd.set('max_threads', 16)

    # increase raytracing memory allowance
    cmd.set('hash_max', 2048)

    # change backround color
    cmd.bg_color('white')

    # sets specular reflection intensity
    cmd.set('specular', '0.1')

    # controls appearence of shadows for ray-traced images
    cmd.set('ray_shadows', 'off')

    # controls antiliasing/edge smoothing for ray-traced images
    cmd.set('antialias', '2')

    # orthoscopic turns on and off the perspective handling
    cmd.set('orthoscopic', 'off')

    # set trasparent background
    cmd.set('ray_opaque_background', '0')

    # raytraces full color without outlines
    cmd.set('ray_trace_mode', '0')

    # settings related to surface features
    cmd.set('surface_quality', '1')
    cmd.set('solvent_radius', '1.6')
    cmd.set('transparency', '0.6')
    cmd.set('surface_color', 'grey80')

    # settings related to rendering meshes
    cmd.set('mesh_quality', '2')
    cmd.set('mesh_type', '0')
    cmd.set('mesh_width', '0.5')
    cmd.set('mesh_radius', '0.015')

    # RNA style settings
    cmd.hide('sticks', 'polymer')
    cmd.set('cartoon_ring_mode', '3')
    cmd.set('cartoon_ring_finder', '1')
    cmd.remove('resn hoh')
    cmd.remove('inorganic and not resn STP')
    cmd.cartoon('oval', 'polymer')
    cmd.set('cartoon_oval_length', '0.75')
    cmd.set('cartoon_oval_width', '0.25')
    cmd.set_color('greyish', [0.625, 0.7, 0.7])
    cmd.set_color('novel', [0.0, 0.4314322351168238, 0.1118361643280874])
    cmd.set_color('known', [0.908605075491407, 0.3955005147576708, 0.0])
    cmd.color('greyish', 'polymer')
    util.cbao('organic')
    cmd.color('lightpink', '(byres polymer & name CA)')
    cmd.cartoon('automatic', '(byres polymer & name CA)')


def color_pockets(pocket_cmap: dict) -> None:
    """settings appearance of the a-spheres (STP) that compose pockets:
    -sets a-sphere radius
    -colors a-spheres based on input color map

    Args:
        pocket_cmap (Dict[int,tuple(int,int,int)]): 
            key: pocket index
            value: rgb color value (tuple)
    """
    cmd.hide('everything', 'resn STP')

    # Alters sphere radius to actual size of a-core.
    cmd.alter('resn STP', 'vdw = b - 1.65')

    # colors each pocket based on color map
    for pocket_num in pocket_cmap:
        cmd.show('spheres', f'resn STP and resi {str(pocket_num)}')
        rgb_tuple = pocket_cmap[pocket_num]
        rgb_list = list(rgb_tuple[0:3])
        cmd.set_color(f'pocket{str(pocket_num)}_color', rgb_list)
        cmd.color(f'pocket{str(pocket_num)}_color',
                  f'resn STP and resi {str(pocket_num)}')


def transparent_pocket(object_name : str, pocket_cmap : dict) -> None:
    cmd.alter(f'{object_name} and resn STP', 'vdw = b - 1.65')
    for pocket_num in pocket_cmap:
        rgb_tuple = pocket_cmap[pocket_num]
        rgb_list = list(rgb_tuple[0:3])
        cmd.set_color(f'pocket{str(pocket_num)}_color', rgb_list)
        cmd.set('surface_color', f'pocket{str(pocket_num)}_color', 
                f'{object_name} and resn stp and resi {pocket_num}')
    
    cmd.extract(f'{object_name}_pockets',
                f'{object_name} and resn STP')
    cmd.hide('sticks', f'{object_name}_pockets')
    cmd.hide('spheres', f'{object_name}_pockets')
    cmd.show('surface', f'{object_name}_pockets')
    cmd.set('transparency', '0.80', f'{object_name}_pockets')


def make_multistate() -> None:

    cmd.alter('resn STP', 'vdw = b - 1.65')
    cmd.set('cartoon_ring_finder', '0')
    cmd.set('cartoon_ring_mode', '1')
    cmd.set('cartoon_transparency', '0.90')
    cmd.set('transparency_mode', '3')


def save_3D_figure(
    path: str, name: str, dpi: int, chain: str, zoom: int
    ) -> None:
    """Saves pymol session file and raytraced png file.

    Args:
        analysis (str): path directory contianing fpocket outputs for analysis
        name (str): name of output pdb structure
        dpi (int): Figure resolution in dpi (dots per linear inch).
        c (str): chain of anzylzed structure
        zoom (float): Sets zoom buffer distance (Å) for creating 3D figures.
    """

    if ',' in chain:
        chain = chain.replace(',', '+')
    cmd.remove('hydrogens')
    cmd.set('ray_trace_fog', '0')
    cmd.orient()
    cmd.rotate('z', '90')
    cmd.zoom(f'chain {chain} and (byres polymer & name O2)', f'{zoom}', complete=1)
    cmd.save(f'{path}/{name}_out_real_sphere.pse')
    dimension = dpi * 8
    print(f'Ray tracing: {name}...')
    cmd.png(f'{path}/{name}_3D_{dpi}.png',
            width=dimension, height=dimension, dpi=dpi, ray=1)
    cmd.reinitialize()
