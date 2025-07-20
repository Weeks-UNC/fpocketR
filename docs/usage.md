# Usage

Full list of arguments for fpocketR CLI.

| Option / Argument             | Type        | Description |
| :---------------------------- | :---------- | :---------- |
| `-pdb`, `--pdb` (Required)    | str         | Path to a .pdb file, .cif file, or 4 character PDB identification code. |
| `-ss`, `--ss`                 | str         | Path to an .ss or other secondary structure file for generating secondary structure figures. |
| `-m`                          | float       | Minimum radius for an a-sphere (Default: 3.0). |
| `-M`                          | float       | Maximum radius for an a-sphere (Default: 5.70). |
| `-i`                          | int         | Minimum number of a-spheres per pocket (Default: 42). |
| `-D`                          | float       | A-sphere clustering distance for forming pockets (Default: 1.65). |
| `-A`                          | int         | Number of electronegative atoms required to define a polar a-sphere (Default: 3). |
| `-p`                          | float       | Maximum ratio of apolar a-spheres in a pocket (Default: 0.0). |
| `-o`, `--out`                 | str         | Path to the output parent directory (Default: "./fpocketR_out"). |
| `-n`, `--name`                | str         | Output filename prefix and output subdirectory name (Default: "{PDB}_clean_out"). |
| `-y`, `--yes`                 | bool        | Answers yes to user prompts for overwriting files (Default: False). |
| `-s`, `--state`               | int         | Specify the NMR states/model to analyze. 0 for all (Default: None). |
| `-c`, `--chain`               | str         | Specify a chain from the input .pdb file (Default: <first_rna_chain>). |
| `-l`, `--ligand`              | str         | PDB ligand identification code (2-3 characters). |
| `-lc`, `--ligandchain`        | str         | Chain containing ligand from the input .pdb file (Default: <--chain input>). |
| `-nt`, `--knownnt`            | list[int]   | List residue IDs of nucleotides in known pocket (e.g. 1,2,3) (Default: None). |
| `-off`, `--offset`            | int         | Offset between starting nucleotide of RNA sequence and starting nucleotide of PDB structure (automatic). |
| `-qf`, `--qualityfilter`      | float       | Minimum fpocket score for pocket (Default: 0.0). |
| `-dpi`, `--dpi`               | int         | Figure resolution in dpi (Default: 300). |
| `-z`, `--zoom`                | float       | Zoom buffer (Å) for creating 3D figures (Default: 5.0). |
| `-cp`, `--connectpocket`      | bool        | Visually connects pockets in 2D figures (Default: False). |
| `-al`, `--alignligand`        | str | bool  | Align structure with pocket prediction (target structure) to an RNA structure with a ligand (mobile structure). |

**TIP:** To see all these options in your terminal, run:

```bash
python -m fpocketR --help
```
