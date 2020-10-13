"""
Top level for running GlauconiteFilter.
Runs all population generation (operations).
Runs plotting at end.
"""
import __future__

import sys

import glauconite.operators.operations as operations

def main_execute(vars):
    """
    This function takes the user variables and runs Glauconite

    Inputs:
    :param dict vars: dict of user variables which will govern how the
        programs runs
    """
    # This will run operations which will:
    # 1) filter ligands
    # 2) optionally convert from 1D smiles to 3D (mol2/PDB)

    sys.stdout.flush()

    smile_file_new_gen, new_gen_ligands_list = operations.populate_generation(vars)
    sys.stdout.flush()
#
