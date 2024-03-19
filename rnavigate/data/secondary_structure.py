#!/usr/bin/env python

###############################################################################
#  SecondaryStructure code originally based on RNAtools.py
#      Main contributors:
#           Gregg Rice
#           Anthony Mustoe
#           Tom Christy
#           Patrick Irving
###############################################################################

import sys
from os.path import isfile
import numpy as np
import pandas as pd
import math
import xml.etree.ElementTree as xmlet
import json
from rnavigate import data


class SecondaryStructure(data.Sequence):
    """Base class for secondary structures.

    Parent classes:
        rnav.data.Data(rnav.data.Sequence)

    Attributes:
        data (str | pandas.DataFrame): dataframe storing base-pairs
            One row for every nucleotide position
            Required columns:
                "Nucleotide"   the nucleotide position
                "Sequence"     the nucleotide letter
                "Pair"         the nucleotide position of the base-pair
                               0 indicates a single stranded nucleotide
            Optional columns:
                "X_coordinate" the x position of the nucleotide in the drawing
                "Y_coordinate" the y position of the nucleotide in the drawing
        sequence (str): sequence string
        nts (numpy.array): "Nucleotide" column of data
        pair_nts (numpy.array): "Pair" column of data
        header (str): header information from CT file
        xcoordinates (numpy array): "X_coordinate" column of data
        ycoordinates (numpy array): "X_coordinate" column of data
    """

    ###########################################################################
    # Initialization methods
    ###########################################################################

    def __init__(self, input_data, extension=None, autoscale=True, name=None, **kwargs):
        """Creates a SecondaryStructure object from a given file or dataframe.

        Args:
            input_data (str | pandas.DataFrame):
                A dataframe or filepath containing a secondary structure
                DataFrame should contain these columns:
                    ["Nucleotide", "Sequence", "Pair"]
                "Pair" column must be redundant.
                Filepath parsing is determined by file extension:
                    varna, xrna, nsd, cte, ct, dbn, bracket, json (R2DT), forna
        """
        if isinstance(input_data, pd.DataFrame):
            self.data = input_data
            self.filepath = "dataframe"
        elif isfile(input_data):
            if extension is None:
                extension = input_data.split(".")[-1].lower()
            read_file = {
                "varna": self.read_varna,
                "xrna": self.read_xrna,
                "nsd": self.read_nsd,
                "cte": self.read_cte,
                "ct": self.read_ct,
                "db": self.read_dotbracket,
                "dbn": self.read_dotbracket,
                "bracket": self.read_dotbracket,
                "json": self.read_r2dt,
                "forna": self.read_forna,
            }[extension]
            self.filepath = input_data
            self.data = read_file(**kwargs)
        super().__init__(input_data=self.data, name=name)
        if "X_coordinate" in self.data.columns and autoscale is True:
            self.transform_coordinates(scale=1, center=(0, 0))
        self.distance_matrix = None

    @classmethod
    def from_sequence(cls, input_data):
        seq = data.Sequence(input_data=input_data)
        df = pd.DataFrame(
            data={
                "Nucleotide": np.arange(seq.length) + 1,
                "Sequence": list(seq.sequence),
                "Pair": np.zeros(seq.length),
            }
        )
        df = df.astype(
            dtype={"Nucleotide": "Int32", "Sequence": "string", "Pair": "Int32"}
        )
        return cls(input_data=df)

    @classmethod
    def from_pairs_list(cls, input_data, sequence):
        """Creates a SecondaryStructure from a list of pairs and a sequence.

        Args:
            pairs (list): 1-indexed list of base pairs. e.g. [(1, 20), (2, 19)]
            sequence (str): sequence string. e.g., "AUCGUGUCAUGCUA"
        """
        structure = cls.from_sequence(sequence)
        structure.add_pairs(input_data)
        return structure

    ###########################################################################
    # properties
    ###########################################################################

    @property
    def nts(self):
        return self.data["Nucleotide"].values

    @property
    def pair_nts(self):
        return self.data["Pair"].values

    @property
    def ycoordinates(self):
        return self.data["Y_coordinate"].values

    @ycoordinates.setter
    def ycoordinates(self, values):
        self.data["Y_coordinate"] = values

    @property
    def xcoordinates(self):
        return self.data["X_coordinate"].values

    @xcoordinates.setter
    def xcoordinates(self, values):
        self.data["X_coordinate"] = values

    def __str__(self):
        """print the filepath and length of the RNA"""
        return f"Name = {self.filepath}, length = {len(self.pair_nts)}"

    ###########################################################################
    # Loading Data
    ###########################################################################

    def read_ct(self, structure_number=0):
        """Loads secondary structure information from a given ct file.

        Requires a properly formatted header!

        Args:
            structure_number (int, optional): If ct file contains multiple structures,
                uses the given structure. Defaults to 0.
        """
        fIN = self.filepath
        num, seq, bp, mask = [], "", [], []
        try:
            with open(fIN) as f:
                # process the first line (pops off the first element)
                spl = f.readline().split()
                nt_length = int(spl[0])
                header = " ".join(spl[1:])
                current_structure = 0
                for i, line in enumerate(f):
                    spl = line.split()
                    if (i + 1) % (nt_length + 1) == 0:
                        # catch the next header
                        current_structure += 1
                        if current_structure > structure_number:
                            break
                        header = " ".join(spl[1:])
                        continue
                    if current_structure == structure_number:
                        num.append(int(spl[0]))
                        seq += str(spl[1])
                        bp.append(int(spl[4]))
                        # check if there is masking info
                        if len(spl) == 7 and spl[6] == "1":
                            mask.append(1)
                        else:
                            mask.append(0)
        except Exception as e:
            raise IOError(f"{fIN} has invalid format or does not exist") from e
        if len(num) == 0:
            sys.exit(f"Structure {structure_number} was not found in the ct file")
        # check consistency!
        for idx, pair in enumerate(bp):
            if pair == 0:
                continue
            p1 = idx + 1
            p2 = bp[pair - 1]
            if p1 != p2:
                print(f"WARNING: inconsistency: {pair} paired to {p1} and {p2}")
        self.header = header
        return pd.DataFrame(
            {"Nucleotide": num, "Sequence": list(seq), "Pair": bp, "Mask": mask}
        )

    def read_varna(self):
        """Generates SecondaryStructure object data from a VARNA file,
        including nucleotide x and y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        tree = xmlet.parse(self.filepath)
        root = tree.getroot()
        # extract sequence, y and x coordinates
        sequence, xcoords, ycoords = "", [], []
        for nt in root.findall("./RNA/bases/nt"):
            base = nt.find("base").text
            sequence += base
            for i in nt:
                if i.get("r") == "pos":
                    xcoords.append(float(i.get("x")))
                    ycoords.append(float(i.get("y")))
        # extract pairing information
        basepairs = []
        for pair in root.findall("./RNA/BPs/bp"):
            i = int(pair.get("part5")) + 1
            j = int(pair.get("part3")) + 1
            basepairs.append((i, j))
        # store attributes
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i - 1] = j
            pairs[j - 1] = i
        return pd.DataFrame(
            {
                "Nucleotide": np.arange(len(sequence)) + 1,
                "Sequence": list(sequence),
                "Pair": pairs,
                "X_coordinate": xcoords,
                "Y_coordinate": ycoords,
            }
        )

    def read_xrna(self):
        """Generates SecondaryStructure object data from an XRNA file,
        including nucleotide x and y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        tree = xmlet.parse(self.filepath)
        root = tree.getroot()
        # extract sequence, x and y coordinates
        nucList = root.findall("./Complex/RNAMolecule/")
        nucLists = []
        for i in nucList:
            if i.tag == "NucListData":
                nucLists.append(i)
        sequence, xcoords, ycoords = "", [], []
        for nt in nucLists[0].text.split("\n"):
            if nt == "":
                continue
            line = nt.split()
            sequence += line[0]
            xcoords.append(float(line[1]))
            ycoords.append(float(line[2]))
        # extract pairing information
        basepairs = []
        for helix in root.findall("./Complex/RNAMolecule/BasePairs"):
            i_outter = int(helix.get("nucID"))
            j_outter = int(helix.get("bpNucID"))
            length = int(helix.get("length"))
            helix_list = [(i_outter + nt, j_outter - nt) for nt in range(length)]
            basepairs.extend(helix_list)
        # store attributes
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i - 1] = j
            pairs[j - 1] = i
        return pd.DataFrame(
            {
                "Nucleotide": np.arange(len(sequence)) + 1,
                "Sequence": list(sequence),
                "Pair": pairs,
                "X_coordinate": xcoords,
                "Y_coordinate": ycoords,
            }
        )

    def read_cte(self):
        """Generates SecondaryStructure object data from a CTE file, including
        nucleotide x and y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        return pd.read_table(
            self.filepath,
            sep=r"\s+",
            usecols=[0, 1, 4, 8, 10],
            names=["Nucleotide", "Sequence", "Pair", "X_coordinate", "Y_coordinate"],
            header=0,
        )

    def read_nsd(self):
        """Generates SecondaryStructure object data from an NSD file
        (format for RNAStructure StructureEditor), including nucleotide x and
        y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        basepairs, sequence, xcoords, ycoords = [], "", [], []
        with open(self.filepath, "r") as file:
            item = ""
            for line in file.readlines():
                line = line.strip().split(" ")
                if "Strand:[" in line:
                    item = "strand"
                    continue
                elif "]" in line:
                    item = ""
                elif "Pairs:[" in line:
                    item = "pairs"
                    continue
                if item == "strand":
                    nt = {}
                    for field in line:
                        if ":" in field:
                            key, value = field.split(":")
                            nt[key] = value
                    sequence += nt["Base"]
                    xcoords.append(float(nt["X"]))
                    ycoords.append(-float(nt["Y"]))
                elif item == "pairs":
                    for field in line:
                        if field.startswith("Pair"):
                            field = field.strip('Pair:"').split(":")
                            basepairs.append([int(nuc) for nuc in field])
        # store attributes
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i - 1] = j
            pairs[j - 1] = i
        return pd.DataFrame(
            {
                "Nucleotide": np.arange(len(sequence)) + 1,
                "Sequence": list(sequence),
                "Pair": pairs,
                "X_coordinate": xcoords,
                "Y_coordinate": ycoords,
            }
        )

    def read_dotbracket(self):
        """Generates SecondaryStructure object data from a dot-bracket notation
        file, including nucleotide x and y coordinates.
        """
        header, seq, bp_str, = (
            "",
            "",
            "",
        )
        with open(self.filepath) as f:
            for line in f:
                if line[0] in ">#":
                    header += line
                    continue
                if all(nt.upper() in "GUACT" for nt in line.strip()):
                    seq += line.strip()
                elif all(s in ".][)(}{><" for s in line.strip()):
                    bp_str += line.strip()
        assert len(bp_str) == len(seq), "db: seq and bp strings are mismatched"
        num = list(range(1, len(seq) + 1))
        bp = [0 for _ in num]
        opens = {"[": [], "{": [], "(": [], "<": []}
        sym_pair = {"]": "[", "}": "{", ")": "(", ">": "<"}
        alphabet = "abcdefghijklmnopqrstuvwxyz"
        sym_pair |= {c: c.upper() for c in alphabet}
        for i, sym in enumerate(bp_str):
            if sym == ".":
                continue
            elif sym in "[{(<" + alphabet.upper():
                opens[sym].append(i)
            elif sym in "]})>" + alphabet:
                j = opens[sym_pair[sym]].pop()
                bp[i] = j + 1
                bp[j] = i + 1
        # store attributes
        return pd.DataFrame(
            {"Nucleotide": np.arange(len(seq)) + 1, "Sequence": list(seq), "Pair": bp}
        )

    def read_r2dt(self):
        """Generates SecondaryStructure object data from an R2DT JSON file,
        including nucleotide x and y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        basepairs, sequence, xcoords, ycoords = [], "", [], []
        with open(self.filepath, "r") as file:
            file = json.load(file)
        for nt in file["rnaComplexes"][0]["rnaMolecules"][0]["sequence"]:
            if nt["residueName"] not in ["5'", "3'"]:
                sequence += nt["residueName"]
                xcoords.append(nt["x"])
                ycoords.append(nt["y"])
        sequence = sequence.replace("5'", "").replace("3'", "")
        for bp in file["rnaComplexes"][0]["rnaMolecules"][0]["basePairs"]:
            i = bp["residueIndex1"]
            j = bp["residueIndex2"]
            if [i, j] not in basepairs:
                basepairs.append([i, j])
        # store attributes
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i - 1] = j
            pairs[j - 1] = i
        return pd.DataFrame(
            {
                "Nucleotide": np.arange(len(sequence)) + 1,
                "Sequence": list(sequence),
                "Pair": pairs,
                "X_coordinate": xcoords,
                "Y_coordinate": ycoords,
            }
        )

    def read_forna(self):
        """Generates SecondaryStructure object data from a FORNA JSON file,
        including nucleotide x and y coordinates.
        """
        # Parse file and get sequence, xcoords, ycoords, and list of pairs.
        basepairs, sequence, xcoords, ycoords = [], [], [], []
        with open(self.filepath, "r") as file:
            file = json.load(file)
        rna_key = list(file["rnas"].keys())[0]
        rna = file["rnas"][rna_key]
        rna_length = rna["rnaLength"]
        for nt in rna["nodes"][:rna_length]:
            sequence.append(nt["name"])
            xcoords.append(nt["x"])
            ycoords.append(nt["y"])
        for pair in rna["pseudoknotPairs"]:
            basepairs.append(pair)
        for i, j in enumerate(rna["pairtable"][-rna_length:]):
            if (j != 0) and ([j, i + 1] not in basepairs):
                basepairs.append([i + 1, j])
        # store attributes
        pairs = [0] * len(sequence)
        for i, j in basepairs:
            pairs[i - 1] = j
            pairs[j - 1] = i
        return pd.DataFrame(
            {
                "Nucleotide": np.arange(len(sequence)) + 1,
                "Sequence": list(sequence),
                "Pair": pairs,
                "X_coordinate": xcoords,
                "Y_coordinate": ycoords,
            }
        )

    ###########################################################################
    # writing files
    ###########################################################################

    def write_sto(self, outfile, name="seq"):
        """ "write structure file out into Stockholm (STO) file format
        for use for infernal searches"""
        with open(outfile, "w") as out:
            # write header
            out.write("# STOCKHOLM 1.0\n\n")
            namelen = max(len(name), 12)
            out.write(name.ljust(namelen + 1))
            # write nts
            for nt in self.sequence:
                out.write(nt)
            out.write("\n")
            out.write("#=GC SS_cons".ljust(namelen + 1))
            # write dot-brackets
            for nt, pair in enumerate(self.pair_nts):
                if pair == 0:
                    out.write(".")
                elif pair > nt:
                    out.write("(")
                else:
                    out.write(")")
            out.write("\n//\n")

    def write_ct(self, out_file):
        """Writes a ct file from the SecondaryStructure object."""
        with open(out_file, "w") as write_file:
            write_file.write(f"{len(self.nts):6d} {self.filepath}\n")
            for nt, seq, pair in zip(self.nts, self.sequence, self.pair_nts):
                prev = nt - 1
                next = nt + 1 % (self.length + 1)  # last nt goes back to zero
                write_file.write(
                    f"{nt:5d} {seq} {prev:5d} {next:5d} {pair:5d} {nt:5d}\n"
                )

    def write_cte(self, outputPath):
        """writes the current structure out to CTE format for Structure Editor.

        Args:
            outputPath (string): path to output cte file to be created
        """
        # set scaling factors based on data source.
        xscale = {"xrna": 1.525 * 20, "varna": 0.469 * 65, "nsd": 30.5, "cte": 30.5}[
            self.ss_type
        ]
        yscale = {"xrna": -1.525 * 20, "varna": 0.469 * 65, "nsd": 30.5, "cte": 30.5}[
            self.ss_type
        ]
        w = open(outputPath, "w")
        w.write("{0:6d} {1}\n".format(self.length, self.filepath))
        line = "{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d} " + ";! X: {5:1f} Y: {6:1f}\n"
        for i in range(self.length):
            xcoord = round(xscale * self.xcoordinates[i])
            ycoord = round(yscale * self.ycoordinates[i])
            num, seq, cti = self.nts[i], self.sequence[i], self.pair_nts[i]
            # no nums after self.length, resets to zero
            nextnum = (num + 1) % (self.length + 1)
            cols = [num, seq, num - 1, nextnum, cti, xcoord, ycoord]
            w.write(line.format(*cols))
        w.close()

    def write_dbn(self, filename, rna_name):
        """Writes dot-bracket notation of current structure to a new file.

        Args:
            filename (str): name of new file
            rna_name (str): name of this RNA for the dbn header
        """
        with open(filename, "w") as write_file:
            write_file.write(f">{rna_name}\n")
            write_file.write(self.sequence + "\n")
            write_file.write(self.get_dotbracket())

    ###########################################################################
    # retrieve structure components or alternative representations
    ###########################################################################

    def get_pairs(self):
        """
        Returns a non-redundant list of base pairs i < j as a array of tuples.
        e.g., [(19,50),(20,49)....]
        """
        out = []
        for left, right in zip(self.nts, self.pair_nts):
            if right != 0 and right > left:
                out.append((left, right))
        return out

    def get_paired_nts(self):
        """Returns a list of residues that are paired."""
        paired_idx = self.pair_nts != 0
        paired_nts = self.data.loc[paired_idx, "Nucleotide"].values.tolist()
        return paired_nts

    def get_unpaired_nts(self):
        """Returns a list of residues that are paired."""
        unpaired_idx = self.pair_nts == 0
        unpaired_nts = self.data.loc[unpaired_idx, "Nucleotide"].values.tolist()
        return unpaired_nts

    def get_junction_nts(self):
        """
        Returns a list of residues at junctions (paired, but adjacent to
        an unpaired residue or the end of a chain)
        """
        junction_residues = []
        for i, nt in enumerate(self.pair_nts):
            if nt == 0:  # if its not paired
                continue
            pi = np.where(self.nts == nt)[0][0]
            neighbors = [i - 1, i + 1, pi - 1, pi + 1]
            neighbors_paired = any(self.pair_nts[neighbors] == 0)
            neighbors_not_in_chain = any([x not in self.nts for x in neighbors])
            if neighbors_paired or neighbors_not_in_chain:
                junction_residues.append(self.nts[i])
        return junction_residues

    def get_nonredundant_ct(self):
        """Returns the ct attribute in non-redundant form - only pairs i<j"""
        pairs = self.get_pairs()
        halfPlexCT = np.zeros_like(self.pair_nts)
        for i, j in pairs:
            halfPlexCT[i - 1] = j
        return halfPlexCT

    def get_helices(self, fill_mismatches=True, split_bulge=True, keep_singles=False):
        """Returns a dictionary of helices from the secondary structure.
        Keys are equivalent to list indices. Values are lists of paired
        nucleotides (1-indexed) in that helix. e.g. {0:[(1,50),(2,49),(3,48)}

        Args:
            fill_mismatches (bool, optional):
                Whether 1-1 and 2-2 bulges are replaced with base pairs.
                Defaults to True.
            split_bulge (bool, optional):
                Whether to split helices on bulges.
                Defaults to True.
            keep_singles (bool, optional):
                Whether to return helices that contain only 1 base-pair.
                Defaults to False.
        """
        # first step is to find all the helices
        rna = self.copy()
        if fill_mismatches:
            rna = rna.fill_mismatches()
        helices = {}
        nt = 0
        heNum = 0
        while nt < len(rna.pair_nts):
            if rna.pair_nts[nt] == 0:
                nt += 1
                continue
            else:
                # skip dups
                if nt > rna.pair_nts[nt]:
                    nt += 1
                    continue
                tempPairs = []
                stillHelix = True
                previous = rna.pair_nts[nt]
                # this is the default behavior which splits helices that are
                # separated by a single nucleotide bulge
                if split_bulge:
                    while stillHelix:
                        # see if the helix juts up against another one
                        if abs(rna.pair_nts[nt] - previous) > 1:
                            break
                        # add the pairs
                        elif rna.pair_nts[nt] != 0:
                            tempPairs.append((nt + 1, rna.pair_nts[nt]))
                            previous = rna.pair_nts[nt]
                            nt += 1
                        else:
                            break
                else:
                    while stillHelix:
                        # Handle bulges on the 3' end of the helix
                        pair = rna.pair_nts[nt]
                        if pair != 0 and not (abs(pair - previous) > 2):
                            tempPairs.append((nt + 1, pair))
                            previous = pair
                            nt += 1
                        # Handle bulges on the 5' end of the helix
                        elif rna.pair_nts[nt + 1] == rna.pair_nts[nt - 1] - 1:
                            nt += 1
                        else:
                            break
                # remove single bp helices
                if len(tempPairs) <= 1 and not keep_singles:
                    continue
                helices[heNum] = tempPairs
                heNum += 1
        return helices

    def extractPK(self, fill_mismatches=True):
        """Returns the pk1 and pk2 pairs from the secondary structure.

        Ignores single base pairs. PK1 is defined as the helix crossing the
        most other bps. if there is a tie, the most 5' helix is called pk1
        returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
        """

        def check_overlap(h1, h2):
            # only need to check one set of pairs from each
            # of the helices. Test is to see if they form
            # a cross hatching pattern
            h1_max = max(h1[0])
            h1_min = min(h1[0])
            h2_max = max(h2[0])
            h2_min = min(h2[0])
            if h1_max > h2_min > h1_min and h2_max > h1_max:
                return True
            if h2_max > h1_min > h2_min and h1_max > h2_max:
                return True
            return False

        # make a copy so we don't destroy the original object
        rna = self.copy()
        # get the helices by calling the extract helix function
        helices = rna.get_helices(fill_mismatches=fill_mismatches)
        heNum = len(helices)
        # do the helices overlap? Check for a crosshatching pattern
        # append them to a list if they have it.
        overlaps = []  # stores the helix number
        for i in range(heNum):
            for j in range(i + 1, heNum):
                if check_overlap(helices[i], helices[j]):
                    overlaps.append((i, j))
        # if there are no overlapping bps, return none
        if len(overlaps) == 0:
            return None, None
        # iterate through all pairs of overlapping helices and count
        # the number of pairs in each to determine helix crossing most bps
        crossbps = {}
        for i, j in overlaps:
            if i not in crossbps:
                crossbps[i] = 0
            if j not in crossbps:
                crossbps[j] = 0
            crossbps[i] += len(helices[j])
            crossbps[j] += len(helices[i])
        # convert to list and sort
        crossbps = list(crossbps.items())
        # first sort list by helix number -- i.e. 5' to 3'
        crossbps.sort(key=lambda x: x[0])
        # now sort by which helix crosses the most BPs. This two-step procedure
        # ensures that in the case of ties, we always choose 5'-most helix
        crossbps.sort(key=lambda x: x[1], reverse=True)
        pk1 = helices[crossbps[0][0]]
        pk2 = []
        for k, v in crossbps[1:]:
            pk2.extend(helices[k])
        return pk1, pk2

    def get_dotbracket(self):
        """Returns a dotbracket notation representation of the secondary
        structure.

        Pseudoknot levels:
            1: ()
            2: []
            3: {}
            4: <>
            5: Aa
            6: Bb
            7: Cc
            etc...
        """
        dbn = ["."] * self.length
        pair_list = self.get_pairs()
        pair_levels = np.zeros(len(pair_list), dtype=int)
        levels = ["()", "[]", "{}", "<>"]

        # check for pseudoknots
        level = 0
        while level in pair_levels:
            indexes = np.where(pair_levels == level)[0]
            for i in indexes:
                for j in indexes[indexes > i]:
                    left1, right1 = pair_list[i]
                    left2, right2 = pair_list[j]
                    if pair_levels[i] != pair_levels[j]:
                        continue
                    elif (left1 < left2 < right1) and (right1 < right2):
                        pair_levels[j] += 1
            level += 1

        for pair, level in zip(pair_list, pair_levels):
            dbn[pair[0] - 1] = levels[level][0]
            dbn[pair[1] - 1] = levels[level][1]
        dbn = "".join(dbn)
        return dbn

    def get_human_dotbracket(self):
        """Returns dotbracket notation string representing SecondaryStructure
        object. This is an experimental format designed to be more human
        readable, i.e. no counting of brackets required.


        1)  Letters, instead of brackets, are used to denote nested base pairs.
        2)  Each helix is assigned a letter, which is incremented one letter
            alphabetically from the nearest enclosing stem.
        3)  Non-nested helices (pseudoknots) are assigned canonical brackets.

        From this canonical dbn string:
            how many bases are in the base stem?
            how many nested helices are there?
            ((((....(((.[[..)))))(((...(((..]].))))))))
        Same question, new format:
            AABB....CCC.[[..cccbbBBB...CCC..]].cccbbbaa

        Pseudoknot levels:
            1: Aa, Bb, Cc, etc.
            2: [], 3: {}, 4: <>
        """
        dbn = list(self.get_dotbracket())
        alphabet = "abcdefghijklmnopqrstuvwxyz"
        closed = {n: 0 for n in alphabet}
        helices = self.get_helices(split_bulge=False, keep_singles=True)
        for pairs in helices.values():
            left_most, right_most = pairs[0]
            if dbn[left_most - 1] not in "()":
                continue
            for n in alphabet:
                if left_most > closed[n]:
                    letter = n
                    closed[n] = right_most
                    break
            for left, right in pairs:
                dbn[left - 1] = letter.capitalize()
                dbn[right - 1] = letter
        dbn = "".join(dbn)
        return dbn

    def get_structure_elements(self):
        """This code is not yet implemented.

        Returns a string with a character for each nucleotide, indicating
        what kind of structure element it is a part of.

        Characters:
            Dangling Ends (E)
            Stems (S)
            Hairpin Loops (H)
            Bulges (B)
            Internal Loops (I)
            MultiLoops (M)
            External Loops (X)
            Pseudoknot (P)
        """

    ###########################################################################
    # edit base pairs
    ###########################################################################

    def add_pairs(self, pairs, break_conflicting_pairs=False):
        """Add base pairs to current secondary structure

        Args:
            pairs (list): 1-indexed list of paired residues.
                e.g. [(1, 20), (2, 19)]
        """
        # get a non-redundant list of non-zero pairs
        pairs = list(set([tuple(sorted(p)) for p in pairs if 0 not in p]))
        # find conflicts within provided pairs
        unique, counts = np.unique(np.ravel(pairs), return_counts=True)
        conflicts = unique[counts > 1]
        conflicts = [[pair for pair in pairs if i in pair] for i in conflicts]
        if len(conflicts) > 0:
            print(
                "There are conflicting pairs:\n\t"
                + "\n\t".join([str(conflict) for conflict in conflicts])
            )
            raise ValueError("Input pairs list cannot contain conflicts")
        # find conflicts in structure and either break existing pair or raise
        conflicts = []
        for i, j in pairs:
            current_i = self.pair_nts[i - 1]
            current_j = self.pair_nts[j - 1]
            if current_i not in [0, j]:
                if break_conflicting_pairs:
                    idx = self.data["Nucleotide"] == current_i
                    self.data.loc[idx, "Pair"] = 0
                else:
                    conflicts.append([(i, j), (i, current_i)])
            if current_j not in [0, i]:
                if break_conflicting_pairs:
                    idx = self.data["Nucleotide"] == current_j
                    self.data.loc[idx, "Pair"] = 0
                else:
                    conflicts.append([(i, j), (current_j, j)])
        if len(conflicts) > 0:
            print(
                "There are conflicting pairs:"
                + "\n\t".join([str(conflict) for conflict in conflicts])
                + "use break_conflicting_pairs=True to remove existing pairs."
            )
            raise ValueError(
                "Conflicting pairs found in structure, use "
                "break_conflicting_pairs=True to remove existing pairs."
            )
        # if all is good, update the current structure
        for i, j in pairs:
            self.pair_nts[i - 1] = j
            self.pair_nts[j - 1] = i

    def break_pairs_region(self, start, end, break_crossing=True, inverse=False):
        """Removes pairs from the specified region (1-indexed, inclusive)

        WARNING: this deletes information

        Args:
            start (int): start position
            end (int): end position
            break_crossing (bool, optional): whether to keep pairs that cross
                over the specified region. Defaults to True.
            inverse (bool, optional): invert the behavior. Defaults to False.
        """
        region = pd.Interval(left=start, right=end, closed="both")
        for nt1, nt2 in self.get_pairs():
            pair_idx = self.data["Nucleotide"].isin([nt1, nt2])
            if inverse and (nt1 not in region or nt2 not in region):
                self.data.loc[pair_idx, "Pair"] = 0
            elif not inverse and (nt1 in region or nt2 in region):
                self.data.loc[pair_idx, "Pair"] = 0
            if break_crossing and nt1 < start < end < nt2:
                self.data.loc[pair_idx, "Pair"] = 0

    def fill_mismatches(self, mismatch=1):
        """Adds base pairs to fill 1,1 and optionally 2,2 mismatches.

        mismatch = 1 will fill only 1,1 mismatches
        mismatch = 2 will fill 1,1 and 2,2
        """
        for i in range(len(self.pair_nts) - 3):
            if self.pair_nts[i + 1] == 0:
                if self.pair_nts[i] - self.pair_nts[i + 2] == 2:
                    self.pair_nts[i + 1] = self.pair_nts[i] - 1
            if mismatch == 2 and self.pair_nts[i + 1] + self.pair_nts[i + 2] == 0:
                if self.pair_nts[i] - self.pair_nts[i + 3] == 3:
                    self.pair_nts[i + 1] = self.pair_nts[i] - 1
                    self.pair_nts[i + 2] = self.pair_nts[i] - 2
        return self

    def break_pairs_nts(self, nt_positions):
        """break base pairs at the given list of positions.

        Args:
            nt_positions (list of int): 1-indexed positions
        """
        nt_idx = self.data["Nucleotide"].isin(nt_positions)
        pair_nts = self.data.loc[nt_idx, "Pair"].to_numpy()
        nt_idx &= self.data["Nucleotide"].isin(pair_nts)
        self.data.loc[nt_idx, "Pair"] = 0

    def break_noncanonical_pairs(self):
        """Removes non-canonical basepairs from the secondary structure.

        WARNING: this deletes information.
        """
        seq = self.sequence.upper().replace("T", "U")
        for nt1, nt2 in self.get_pairs():
            bp = "".join([seq[nt1 - 1], seq[nt2 - 1]])
            if bp not in ["AU", "UA", "GC", "CG", "GU", "UG"]:
                mask = self.data["Nucleotide"].isin([nt1, nt2])
                self.data.loc[mask, "Pair"] = 0

    def break_singleton_pairs(self):
        """Removes singleton basepairs from the secondary structure.

        WARNING: This deletes information.
        """
        pair_nts = self.pair_nts
        for nt1, nt2 in self.get_pairs():
            if pair_nts[nt1 - 2] != nt2 - 1 and pair_nts[nt1] != nt2 + 1:
                mask = self.data["Nucleotide"].isin([nt1, nt2])
                self.data.loc[mask, "Pair"] = 0

    ###########################################################################
    # Make calculations based on structures
    ###########################################################################

    def get_distance_matrix(self, recalculate=False):
        """Based on Tom's contact_distance function, but instead returns
        the all pairs shortest paths matrix, and stores it as an attribute. If
        the attribute has already been set, it returns the attribute. This is
        faster than calling contact_distance pairwise to fill the matrix.

        Args:
            recalculate (bool, optional): Set to true to recalculate the matrix
                even if the attribute is set. In case changes to the structure
                have been made.
        """
        if (self.distance_matrix is not None) and not recalculate:
            return self.distance_matrix
        # this method will be used later to make sure a nucleotide hasn't
        # been visited and is within the bounds of the RNA
        def viable(nt):
            not_empty = nt != 0
            valid_so_far = not_empty
            if valid_so_far:
                in_range = min(self.nts) <= nt <= max(self.nts)
                valid_so_far = in_range
            if valid_so_far:
                not_visited = level[np.where(self.nts == nt)[0][0]] == -1
                valid_so_far = not_visited
            return valid_so_far

        # The default value of -1 means that a nt hasn't been visited.
        distance_matrix = np.full((len(self.nts), len(self.nts)), -1)
        for i, nt in enumerate(self.nts):
            # Level keeps track of how far each other nt is from nt
            level = distance_matrix[i, :]
            # create the queue and add nt, set its level as 0
            queue = list()
            queue.append(nt)
            level[i] = 0
            # while there are nucleotides that haven't been visited
            while len(queue) > 0:
                current_nt = queue.pop(0)
                current_i = np.where(self.nts == current_nt)[0][0]
                # Find all the neighbors of the current NT
                NTBack = current_nt - 1
                NTUp = current_nt + 1
                NTPair = self.pair_nts[current_i]
                # check all the neighbors to see if they should be added
                # and if so, record their level(ie their Contact Distance)
                for neighbor in [NTBack, NTUp, NTPair]:
                    if viable(neighbor):
                        queue.append(neighbor)
                        level[np.where(self.nts == neighbor)[0][0]] = (
                            level[current_i] + 1
                        )
            # set row of distance matrix for nt
            distance_matrix[i, :] = level
        # store the distance matrix from the search and return
        return distance_matrix.copy()

    def contact_distance(self, i, j):
        """Returns the contact distance between positions i and j"""
        if self.distance_matrix is None:
            self.distance_matrix = self.get_distance_matrix()
        return self.distance_matrix[i - 1, j - 1]

    def compute_ppv_sens(self, structure2, exact=True):
        """Compute the PPV and sensitivity between self and another
        SecondaryStructure object.

        Args:
            structure2 (SecondaryStructure): The SecondaryStructure to compare to.
            exact (bool, optional): True requires BPs to be exactly correct.
                                    False allows +/-1 bp slippage.
                                    Defaults to True.

        Returns:
            float, float, tuple: sensitivity, PPV, (TP, TP+FP, TP+FN)
        """
        if len(self.pair_nts) != len(structure2.pair_nts):
            raise ValueError("sequence lengths must be the same")
        trues = sum(self.pair_nts != 0) / 2
        positives = sum(structure2.pair_nts != 0) / 2
        true_positives = 0
        for i, v in enumerate(self.pair_nts):
            if v == 0 or v < i:
                continue
            try:
                if v == structure2.pair_nts[i]:
                    true_positives += 1
                elif not exact:
                    neighbors = [
                        structure2.pair_nts[i] - 1,
                        structure2.pair_nts[i] + 1,
                        structure2.pair_nts[i - 1],
                        structure2.pair_nts[i + 1],
                    ]
                    if v in neighbors:
                        true_positives += 1
            except IndexError:
                pass
        try:
            ppv = true_positives / float(positives)
            sens = true_positives / float(trues)
        except ZeroDivisionError:
            ppv, sens = 0.0, 0.0
        return sens, ppv, (true_positives, positives, trues)

    ###########################################################################
    # get other data objects based on this one
    ###########################################################################

    def copy(self):
        return self.get_aligned_data(self.null_alignment)

    def get_aligned_data(self, alignment):
        """Returns a new SecondaryStructure matching the alignment target.

        Args:
            alignment (data.Alignment): an alignment object used to map values
        """
        df = alignment.map_nucleotide_dataframe(self.data)
        df["Pair"] = df["Pair"].fillna(0)
        mask = df["Pair"] != 0
        df.loc[mask, "Pair"] = alignment.map_positions(df.loc[mask, "Pair"].values)
        return SecondaryStructure(
            input_data=df,
            autoscale=False,
            name=self.name,
        )

    def get_interactions_df(self):
        mask = self.data.eval("Pair != 0 & Pair > Nucleotide")
        df = self.data.loc[mask, ["Nucleotide", "Pair"]].copy()
        df["Structure"] = 1
        df = (
            df.rename(columns={"Nucleotide": "i", "Pair": "j"})
            .sort_values(by=["i", "j"])
            .astype({"Structure": "Int32"})
            .reset_index(drop=True)
        )
        return df

    def as_interactions(self, structure2=None):
        """returns list of i, j basepairs as rnavigate.Interactions data.

        Args:
            structure2 (rnavigate.data.SecondaryStructure, optional):
                If provided, basepairs from both structures are included and
                colored by structure (left, right, or both)
                defaults to None.
        """
        if structure2 is None:
            return data.StructureAsInteractions(
                input_data=self,
                sequence=self.sequence,
            )
        elif isinstance(structure2, SecondaryStructure):
            return data.StructureCompareTwo(
                input_data=[self, structure2],
                sequence=self.sequence,
            )
        elif isinstance(structure2, list):
            return data.StructureCompareMany(
                input_data=[self] + structure2,
                sequence=self.sequence,
            )
        raise ValueError(
            "structure2 must be a SecondaryStructure or list of " "SecondaryStructures"
        )

    ###########################################################################
    # edit other attributes
    ###########################################################################

    def normalize_sequence(self, t_or_u="U", uppercase=True):
        super().normalize_sequence(t_or_u=t_or_u, uppercase=uppercase)
        self.data["Sequence"] = list(self.sequence)

    def transform_coordinates(
        self, flip=None, scale=None, center=None, rotate_degrees=None
    ):
        """Perform transformation of structure coordinates

        Args:
            flip (str, optional):
                "horizontal" or "vertical".
                Defaults to None.
            scale (float, optional):
                new median distance of basepairs.
                Defaults to None.
            center (tuple of floats, optional):
                new center x and y coordinate.
                Defaults to None.
            rotate_degrees (float, optional):
                number of degrees to rotate structure.
                Defaults to None.
        """
        structure = StructureCoordinates(
            self.xcoordinates, self.ycoordinates, self.get_pairs()
        )
        if flip is not None:
            structure.flip(flip == "horizontal")
        if scale is not None:
            structure.scale(scale)
        if center is not None:
            structure.center(*center)
        if rotate_degrees is not None:
            structure.rotate(rotate_degrees)
        self.xcoordinates = structure.x
        self.ycoordinates = structure.y


###############################################################################
# end of SecondaryStructure class
###############################################################################


class StructureCoordinates:
    """Helper class to perform structure coordinate transformations"""

    def __init__(self, x, y, pairs=None):
        """initialize structure with arrays holding x and y coordinates. A list
        of base-paired positions is required only if scaling coordinates.

        Args:
            x (numpy.array): x coordinates
            y (numpy.array): y coordinates
            pairs (list of pairs, optional): list of base-paired positions.
        """
        self.x = x
        self.y = y
        self.pairs = pairs

    def get_center_point(self):
        """return the x, y coordinates for the center of structure

        Returns:
            tuple: x and y coordinates of structure center as floats
        """
        x_center = (max(self.x) + min(self.x)) / 2
        y_center = (max(self.y) + min(self.y)) / 2
        return (x_center, y_center)

    def scale(self, median_bp_distance=1):
        """Scale structure such that median base-pair distance is constant.

        Args:
            median_bp_distance (int, optional):
                New median distance between all base-paired nucleotides.
                Defaults to 1.
        """
        bp_distances = []
        for bp in self.pairs:
            y_dist = self.y[bp[0] - 1] - self.y[bp[1] - 1]
            x_dist = self.x[bp[0] - 1] - self.x[bp[1] - 1]
            bp_distances.append((y_dist**2 + x_dist**2) ** 0.5)
        scale_factor = np.median(bp_distances)
        self.y *= median_bp_distance / scale_factor
        self.x *= median_bp_distance / scale_factor

    def flip(self, horizontal=True):
        """Flip structure vertically or horizontally

        Args:
            horizontal (bool, optional):
                whether to flip structure horizontally, otherwise vertically
                Defaults to True.
        """
        if horizontal:
            self.x *= -1
        else:
            self.y *= -1

    def center(self, x=0, y=0):
        """Center structure on given x, y coordinate

        Args:
            x (int, optional): x coordinate of structure center. Defaults to 0.
            y (int, optional): y coordinate of structure center. Defaults to 0.
        """
        x_center, y_center = self.get_center_point()
        self.x -= x_center + x
        self.y -= y_center + y

    def rotate(self, degrees):
        """Rotate structure on current center point

        Args:
            degrees (float): number of degrees to rotate structure
        """
        radians = math.radians(degrees)
        center_x, center_y = self.get_center_point()
        # Translate the coordinates to the origin
        translated_x = self.x - center_x
        translated_y = self.y - center_y
        # Perform rotation
        self.x = translated_x * math.cos(radians) - translated_y * math.sin(radians)
        self.y = translated_x * math.sin(radians) + translated_y * math.cos(radians)
        # Translate the coordinates back to the original center
        self.x += center_x
        self.y += center_y


class SequenceCircle(SecondaryStructure):
    def __init__(self, input_data, gap=30, name=None, **kwargs):
        length = input_data.length
        df = pd.DataFrame(
            {
                "Sequence": list(input_data.sequence),
                "Nucleotide": np.arange(length) + 1,
                "Pair": np.zeros(length),
            }
        )
        gap = 2 * np.pi * gap / 360
        self.gap = gap
        theta_between = (2 * np.pi - self.gap) * 1 / (length - 1)
        self.radius = (length + gap / theta_between) / np.pi
        theta_start = gap / 2
        theta = theta_between * np.arange(length) + theta_start
        df["Theta"] = theta
        super().__init__(
            input_data=df, extension=None, autoscale=False, name=name, **kwargs
        )
