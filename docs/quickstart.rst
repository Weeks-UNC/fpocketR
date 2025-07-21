Quickstart Guide
===============

Run fpocketR from the command line to analyze RNA structures and visualize ligand binding pockets.

Activate the fpocketR Environment
--------------------------------

.. code-block:: bash

   conda activate fpocketR

Tip: For a full list of options, run:

.. code-block:: bash

   python -m fpocketR --help

Basic Structure Analysis
-----------------------

Analyze a local PDB file or fetch by PDB ID using the ``-pdb`` argument:

.. code-block:: bash

   python -m fpocketR -pdb 3e5c.pdb
   # or
   python -m fpocketR -pdb 3e5c

Outputs pockets in tertiary structure:
"""""""""""""""""""


.. figure:: /_static/images/3e5c_3D.png
    :width: 50%
    :align: center
    :alt: Tertiary structure


Pocket color legend:
"""""""""""""""""""

.. figure:: /_static/images/fpocketR_pocket_color_legend.png
   :width: 40%
   :align: center
   :alt: Pocket color legend

Secondary Structure Visualization
---------------------------------

Add a secondary structure diagram using the ``-ss`` argument:

.. code-block:: bash

   python -m fpocketR -pdb 2l1v.pdb -ss 2l1v.nsd

Outputs pockets in secondary structure:
"""""""""""""""""""

.. figure:: /_static/images/2l1v_2D.png
    :width: 30%
    :align: center
    :alt: Secondary structure


Multistate Analysis
-------------------

Analyze all NMR or Cryo-EM states using the ``--state 0`` argument:

.. code-block:: bash

   python -m fpocketR -pdb 2l1v.pdb -ss 2l1v.nsd --state 0

Outputs pocket density in 3D and 2D:
"""""""""""""""""""

.. figure:: /_static/images/2l1v_all_states_3D.png
    :width: 50%
    :align: center
    :alt: Tertiary structure (pocket density)

.. figure:: /_static/images/2l1v_2D_pocket_density.png
    :width: 30%
    :align: center
    :alt: Secondary structure (pocket density)

Apo/Holo Analysis
-----------------

Align ligand-bound (holo) and ligand-free (apo) structures for direct comparison using the ``--alignligand`` argument:

.. code-block:: bash

   python -m fpocketR -pdb 8f4o_apo.pdb --alignligand 2gdi_holo.pdb --knownnt 19,20,42,43

Outputs apo structure and aligned apo/holo structures:
"""""""""""""""""""

.. figure:: /_static/images/8f4o_apo_3D.png
    :width: 50%
    :align: center
    :alt: Apo structure and pocket

.. figure:: /_static/images/8f4o_apo_holo.png
    :width: 50%
    :align: center
    :alt: Apo and holo structures aligned

Additional Arguments
-------------------

Customize analysis with optional arguments:

- Select RNA chain: ``-c (--chain)``
- Select ligand: ``-l (--ligand)``
- Set raytracing resolution (lower = faster): ``-dpi (--dpi)``
- Specify output path: ``-o (--out)``

.. code-block:: bash

   python -m fpocketR -pdb 2gdi_holo.pdb --chain Y --ligand TPP --dpi 10 --out ./TPP_RS

* Output files and figures add to custom directory: ``./TPP_RS/2gdi_holo_clean_out/``.

Outputs low resolution (fast) tertiary structure:
"""""""""""""""""""

.. figure:: /_static/images/2gdi_holo_3D_10.png
   :width: 60%
   :align: center
   :alt: Tertiary structure (low resolution)
