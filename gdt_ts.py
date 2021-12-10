# !/usr/bin/env python3
import os
import sys

import numpy as np
from pymol import cmd

# Author: Martin Kim (martinkim0)

"""
Script that computes the global distance test total score (GDT-TS) used in
the CASP assessment for two protein structures - the ground truth and the 
prediction.

Usage:
- To compare a new pair of protein structures, create a new folder inside the 'data/'
directory and add the two PDB files, one ending with '_predicted' and the other with
'_true'. 
- To compute the scores for all pairs in 'data/':
>>> python gdt_ts.py
- To compute the score for one pair:
>>> python gdt_ts.py {pdb id}
"""

def load_pair(dir: str, root: str='data/'):
    path = root + dir
    if not os.path.exists(path):
        raise Exception('Path does not exist')

    contents = os.listdir(path)
    if len(contents) != 2:
        print(contents)
        raise Exception('Each data subdirectory must have exactly 2 files')

    pair = []
    for f in contents:
        if not f.endswith('.pdb'):
            raise Exception('Input files must be PDB format')

        name = f.split(".")[0]
        cmd.load(root + dir + '/' + f, name)
        pair.append(name)
    
    return tuple(pair)

@cmd.extend
def gdt(comp: str, ref: str, name: str):
    """
    Parameters
    ref: reference structure
    comp: predicted structure
    """
    # PyMOL syntax to select alpha-carbons
    ca_suffix = ' and name CA'
    ref_ca = ref + "_ca"
    comp_ca = comp + "_ca"
    cmd.select(ref_ca, selection=ref + ca_suffix)
    cmd.select(comp_ca, selection=comp + ca_suffix)

    # Default cutoff distances for GDT-TS
    cutoffs = [1., 2., 4., 8.]

    # Superimpose the structures - more robust than PyMOL's `align`
    cmd.super(comp_ca, ref_ca, object=name)
    mappings = cmd.get_raw_alignment(name)

    # Retrieve the distances between superimposed atoms
    distances = []
    for m in mappings:
        atom1 = f'{m[0][0]} and id {m[0][1]}'
        atom2 = f'{m[1][0]} and id {m[1][1]}'
        dist = cmd.get_distance(atom1, atom2)
        cmd.alter(atom1, f'b = {dist}')
        distances.append(dist)

    # Compute cutoff distances
    distances = np.array(distances)
    n = len(distances)
    gdts = [np.sum(distances <= c) / n for c in cutoffs]

    return np.mean(gdts)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise Exception("Invalid number of arguments")

    pairs = []
    if len(sys.argv) == 1:
        contents = os.listdir('data')
    else:
        contents = [sys.argv[1]]

    for dir in contents:
        if os.path.isdir('data/' + dir):
            pair = load_pair(dir)
            pair = sorted(pair)
            pairs.append(pair)

    for pair in pairs:
        name = pair[0].split('_')[0]
        score = gdt(*pair, name)
        print(f"GDT for {name}: {score}")
