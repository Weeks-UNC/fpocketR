"""The PDB object to represent tertiary structures with atomic coordinates.

This data can be used to filter interactions by 3D distance, and to visualize
profile and interactions data on interactive 3D structures.
"""

import Bio.PDB
from rnavigate import data
import numpy as np


class PDB(data.Sequence):
    """A class to represent RNA tertiary structures with atomic coordinates.

    This data can be used to filter interactions by 3D distance, and to visualize
    profile and interactions data on interactive 3D structures.

    Parameters
    ----------
    input_data : str
        path to a PDB or CIF file
    chain : str
        chain identifier of RNA of interest
    sequence : rnavigate.Sequence or str, optional
        A sequence to use as the reference sequence.
        This is required if the sequence cannot be found in the header
        Defaults to None.
    name : str, optional
        A name for the data set. Defaults to None.

    Attributes
    ----------
    sequence : str
        The RNA sequence
    length : int
        The length of the RNA sequence
    name : str
        A name for the data set
    path : str
        The path to the PDB or CIF file
    chain : str
        The chain identifier of the RNA of interest
    offset : int
        The offset between the sequence positions and the PDB residue indices
    pdb : Bio.PDB.Structure.Structure
        The PDB structure
    pdb_idx : np.array
        The PDB indices of the RNA
    pdb_seq : np.array
        The PDB sequence of the RNA
    distance_matrix : dict
        A dictionary of distance matrices for each atom type
    """

    def __init__(self, input_data, chain, sequence=None, name=None):
        """Construct PDB object based on an input PDB or CIF file."""
        self.offset = 0
        self.chain = chain
        if sequence is None:
            sequence = self.get_sequence(input_data)
        super().__init__(sequence, name=name)
        self.read_pdb(input_data)
        self.path = input_data
        self.distance_matrix = {}

    def get_sequence(self, pdb):
        """Find the sequence in the provided CIF or PDB file.

        Parameters
        ----------
        pdb : str
            path to a PDB or CIF file

        Returns
        -------
        sequence : string
            The RNA sequence
        """
        seqres = []
        with open(pdb) as file:
            for line in file.readlines():
                line = [field.strip() for field in line.split()]
                if line[0] == "SEQRES" and line[2] == self.chain:
                    seqres.extend(line[4:])
        return self.get_sequence_from_seqres(seqres)

    def get_sequence_from_seqres(self, seqres):
        """Used by get_sequence to parse the SEQRES entries.

        Parameters
        ----------
        seqres : list
            A list of SEQRES entries for the RNA chain of interest

        Returns
        -------
        sequence : string
            The RNA sequence
        """
        sequence = []
        for nt in seqres:
            valid = nt[0].upper() in "GUACT"
            if valid and (len(nt) == 3) and nt.endswith("TP"):
                sequence.append(nt[0])
            elif valid and (len(nt) == 1):
                sequence.append(nt[0])
            else:
                raise ValueError("invalid nt in SEQRES: " + nt)
        return "".join(sequence)

    def read_pdb(self, pdb):
        """Read a PDB or CIF file into the data structure.

        Parameters
        ----------
        pdb : str
            path to a PDB or CIF file
        """
        if pdb.split(".")[-1] == "pdb":
            parser = Bio.PDB.PDBParser(QUIET=True)
        elif pdb.split(".")[-1] == "cif":
            parser = Bio.PDB.MMCIFParser(QUIET=True)
        self.pdb = parser.get_structure("RNA", pdb)
        self.pdb_idx = []
        self.pdb_seq = []
        for res in self.pdb[0][self.chain].get_residues():
            res_id = res.get_id()
            res_seq = res.get_resname()
            if res_id[0] == " ":
                self.pdb_idx.append(res_id[1])
                self.pdb_seq.append(res_seq.strip())
        self.pdb_idx = np.array(self.pdb_idx)
        self.set_indices()

    def set_indices(self):
        """Uses self.data and self.sequence to set self.offset"""
        for i in range(len(self.sequence)):
            correct_offset = True
            for pdb_nt, pdb_idx in zip(self.pdb_seq, self.pdb_idx):
                if self.sequence[pdb_idx - i - 1] != pdb_nt:
                    correct_offset = False
                    break
            if correct_offset:
                self.offset = i
                break
        if not correct_offset:
            print("PDB entries could not be matched to sequence.")

    def get_pdb_idx(self, seq_idx):
        """Return the PDB index given the sequence index (0-indexed)."""
        return seq_idx + self.offset

    def get_seq_idx(self, pdb_idx):
        """Return the sequence index given the PDB index."""
        return pdb_idx - self.offset

    def is_valid_idx(self, pdb_idx=None, seq_idx=None):
        """Determines if a PDB or sequence index is in the PDB structure.

        Parameters
        ----------
        pdb_idx : int, optional
            A PDB index (1-indexed). Defaults to None.
        seq_idx : int, optional
            A sequence index (1-indexed). Defaults to None.

        Returns
        -------
        bool
            True if the index is in the PDB structure, False otherwise.
        """
        if pdb_idx is not None and pdb_idx in self.pdb_idx:
            return True
        elif seq_idx is not None and (seq_idx in (self.pdb_idx - self.offset)):
            return True
        else:
            return False

    def get_xyz_coord(self, nt, atom):
        """Return the x, y, and z coordinates for a given residue and atom.

        Parameters
        ----------
        nt : int
            The nucleotide of interest (1-indexed)
        atom : string or dict, defaults to "O2'"
            The atom to use for distance calculations. If a string, the same atom
            will be used for all residues. If a dict, the atom will be chosen based
            on the nucleotide type. If "DMS", the N1 atom will be used for A and G,
            and the N3 atom will be used for U and C.

        Returns
        -------
        xyz : list
            A list of x, y, and z coordinates
        """
        pdb_idx = self.get_pdb_idx(nt)
        if atom == "DMS":
            atom = {nt: "N1" for nt in "AG"} | {nt: "N3" for nt in "UC"}
        elif isinstance(atom, str):
            atom = {nt: atom for nt in "AUCG"}
        seq = self.sequence[nt - 1].upper().replace("T", "U")
        atom = atom[seq]
        xyz = [
            float(c) for c in self.pdb[0][self.chain][int(pdb_idx)][atom].get_coord()
        ]
        return xyz

    def get_distance(self, i, j, atom="O2'"):
        """Get the distance between given atom in nucleotides i and j (1-indexed).

        Parameters
        ----------
        i : int
            The first nucleotide
        j : int
            The second nucleotide
        atom : string or dict, defaults to "O2'"
            The atom to use for distance calculations. If a string, the same atom
            will be used for all residues. If a dict, the atom will be chosen based
            on the nucleotide type. If "DMS", the N1 atom will be used for A and G,
            and the N3 atom will be used for U and C.

        Returns
        -------
        distance : float
            The distance between the atoms
        """
        if atom in self.distance_matrix:
            return self.distance_matrix[atom][i - 1, j - 1]
        try:
            xi, yi, zi = self.get_xyz_coord(i, atom)
            xj, yj, zj = self.get_xyz_coord(j, atom)
            distance = ((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2) ** 0.5
        except KeyError:
            distance = np.nan
        return distance

    def get_distance_matrix(self, atom="O2'"):
        """Get the pairwise atomic distance matrix for all residues.

        Parameters
        ----------
        atom : string or dict, defaults to "O2'"
            The atom to use for distance calculations. If a string, the same atom
            will be used for all residues. If a dict, the atom will be chosen based
            on the nucleotide type. If "DMS", the N1 atom will be used for A and G,
            and the N3 atom will be used for U and C.

        Returns
        -------
        matrix : NxN numpy.ndarray
            A 2D array of pairwise distances. N is the length of the RNA.
        """
        if atom in self.distance_matrix:
            return self.distance_matrix[atom]
        x = np.full(self.length, np.nan)
        y = np.full(self.length, np.nan)
        z = np.full(self.length, np.nan)
        for i in self.pdb_idx - self.offset:
            try:
                x[i - 1], y[i - 1], z[i - 1] = self.get_xyz_coord(i, atom)
            except KeyError:
                pass
        a = x - x[:, np.newaxis]
        b = y - y[:, np.newaxis]
        c = z - z[:, np.newaxis]
        matrix = np.sqrt(a * a + b * b + c * c)
        # self.distance_matrix[atom] = np.nan_to_num(matrix, nan=1000)
        return matrix
