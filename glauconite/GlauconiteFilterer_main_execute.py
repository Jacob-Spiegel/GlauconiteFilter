"""
Top level for running GlauconiteFilterer.
Runs all population generation (operations).
Runs plotting at end.
"""
import __future__

import os
import glob
import sys
import shutil

import glauconite.operators.operations as operations

def main_execute(vars):
    """
    This function takes the user variables and runs Glauconite

    Inputs:
    :param dict vars: dict of user variables which will govern how the
        programs runs
    """

    # Unpack necessary variables
    # output_directory is the root output folder for the run
    output_directory = vars["output_directory"]

    # This will run operations which will:
    # 1) filter ligands
    # 2) optionally convert from 1D smiles to 3D (mol2/PDB)

    sys.stdout.flush()

    smile_file_new_gen, new_gen_ligands_list = operations.populate_generation(vars)
    sys.stdout.flush()
#
