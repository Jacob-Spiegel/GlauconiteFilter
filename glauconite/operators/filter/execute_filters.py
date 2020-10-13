"""
Top level for running filters.
"""
import __future__

import copy

import rdkit
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog("rdApp.*")

from glauconite.operators.filter.filter_classes.parent_filter_class import ParentFilter
from glauconite.operators.filter.filter_classes.get_child_filter_class import get_all_subclasses

import glauconite.operators.convert_files.gypsum_dl.gypsum_dl.MolObjectHandling as MOH
from glauconite.operators.filter.filter_classes.filter_children_classes import *


def make_run_class_dict(filters_to_use):
    """
    This will retrieve all the names of every child class of the parent class
    ParentFilter

    Inputs:
    :param list filters_to_use: list of filters to be used.
            defined in vars["chosen_ligand_filters"]

    Returns:
    :returns: dict child_dict: This dictionary contains all the names of the
        chosen filters as keys and the the filter objects as the items. returns
        None if no filters are specified by user.
    """

    if filters_to_use is None:
        # if the user turned off filters
        return None

    children = get_all_subclasses(ParentFilter)

    child_dict = {}
    for child in children:
        child_object = child()
        child_name = child_object.get_name()

        if child_name in filters_to_use:
            child_dict[child_name] = child_object

    return child_dict


def run_filter(vars, list_of_new_ligands, verbose=True):
    """
    This will run a filter of the Users choosing.

    This will take a list of lists of ligands to filter. list_of_new_ligands =
    [["CCC","Zinc123],["CCCC","Zinc1234]]

    Inputs:
    :param dict vars: User variables which will govern how the programs runs
    :param list list_of_new_ligands: list of lists containing all the newly
        generated ligands and their names
    :param bool verbose: print message if True; don't if False

    Returns:
    :returns: list ligands_which_passed_filter: a list of only the molecules
        which passed the filter. Excludes all molecules which failed.
    """

    # Get the already generated dictionary of filter objects
    filter_object_dict = vars["filter_object_dict"]
    start_num = len(list_of_new_ligands)
    # make a list of tuples for multi-processing Filter
    job_input = []
    for smiles_info in list_of_new_ligands:
        temp_tuple = tuple([smiles_info, filter_object_dict])
        job_input.append(temp_tuple)
    job_input = tuple(job_input)

    results = vars["parallelizer"].run(job_input, run_filter_mol)

    # remove mols which fail the filter
    ligands_failed_to_sanitize = [x[0] for x in results if x[1] == "Sanitize_fail"]
    ligands_which_failed_filter = [x[0] for x in results if x[1] == "Filter_fail"]
    ligands_which_passed_filter = [x[0] for x in results if x[1] == "Filter_Passed"]

    if verbose is True:
        print("######################")
        print("")
        print("Pass/fail Stats")
        print("")
        print("Total number of ligs starting: ", len(list_of_new_ligands))
        print("Number of ligs failed sanitization: ", len(ligands_failed_to_sanitize))
        print("Number of ligs failed the filters: ", len(ligands_which_failed_filter))
        print("Number of ligs which PASSED: ", len(ligands_which_passed_filter))
        print("")
        print("{}% PASSED".format(str(100*len(ligands_which_passed_filter)/start_num)))
        print("")
        print("######################")

    return ligands_which_passed_filter


def run_filter_mol(smiles_info, child_dict):
    """
    This takes a smiles_string and the selected filter list (child_dict) and
    runs it through the selected filters.

    Inputs:
    :param list smiles_info: A list with info about a ligand, the SMILES string
        is idx=0 and the name/ID is idx=1. example: smiles_info
        ["CCCCCCC","zinc123"]
    :param dict child_dict: This dictionary contains all the names of the
        chosen filters as keys and the the filter objects as the items Or None if
        User specifies no filters

    Returns:
    :returns: list smiles_info: list of the smiles_info if it passed the filter and
        "Filter_Passed".
        returns smiles_info and "Sanitize_fail" if the mol fails to sanitize.
        returns smiles_info and "Sanitize_fail" if the mol fails one or more filters.
    """

    smiles_string = smiles_info[0]

    mol = Chem.MolFromSmiles(smiles_string, sanitize=False)
    # try sanitizing, which is necessary later
    mol = MOH.check_sanitization(mol)
    if mol is None:
        return [smiles_info, "Sanitize_fail"]

    mol = MOH.try_deprotanation(mol)
    if mol is None:
        return [smiles_info, "Sanitize_fail"]

    mol = MOH.check_sanitization(mol)
    if mol is None:
        return [smiles_info, "Sanitize_fail"]

    # remove charge from mol objects. This affects some properties
    # such as: logP, Mol refractivity, and polar surface area
    # which can impact filters such as Ghose and VandeWaterbeemd
    # This is done because logP is traditionally applied to neutral molecules
    uncharger_obj = rdMolStandardize.Uncharger()
    mol = uncharger_obj.uncharge(mol)
    if mol is None:
        return [smiles_info, "Sanitize_fail"]

    if child_dict is not None:
        # run through the filters
        filter_result = run_all_selected_filters(mol, child_dict)

        # see if passed
        if filter_result is False:
            return [smiles_info, "Filter_fail"]
        # it passed return the smiles_info
        return [smiles_info, "Filter_Passed"]

    # This will return None
    return [smiles_info, "Filter_Passed"]


def run_filter_on_just_smiles(smile_string, child_dict):
    """
    This takes a smiles_string and the selected filter list (child_dict) and
    runs it through the selected filters.

    Inputs:
    :param str smile_string: A smiles_string. example: smiles_info
        ["CCCCCCC","zinc123"]
    :param dict child_dict: This dictionary contains all the names of the
        chosen filters as keys and the the filter objects as the items Or None if
        User specifies no filters

    Returns:
    :returns: str smile_string: smile_string if it passed the filter. returns
        False If the mol fails a filter.
    """

    mol = Chem.MolFromSmiles(smile_string, sanitize=False)
    # try sanitizing, which is necessary later
    mol = MOH.check_sanitization(mol)
    if mol is None:
        return False

    mol = MOH.try_deprotanation(mol)
    if mol is None:
        return False

    if child_dict is not None:
        # run through the filters
        filter_result = run_all_selected_filters(mol, child_dict)

        # see if passed
        if filter_result is False:
            return False
        # it passed return the smiles_info
        return smile_string
    # return the smile string
    return smile_string


def run_all_selected_filters(mol, child_dict):
    """
    Iterate through all of the filters specified by the user for a single
    molecule. returns True if the mol passes all the chosen filters. returns
    False if the mol fails any of the filters.

    Inputs:
    :param rdkit.Chem.rdchem.Mol object mol: An rdkit mol object to be tested
        if it passes the filters
    :param dict child_dict: This dictionary contains all the names of the
        chosen filters as keys and the the filter objects as the items

    Returns:
    returns bol bol: True if the mol passes all the filters. False if the mol
        fails any filters.
    """

    filters_failed = 0
    mol = MOH.check_sanitization(mol)
    if mol is None:
        return False
    for child in list(child_dict.keys()):
        mol_copy = copy.deepcopy(mol)
        filter_function = child_dict[child].run_filter
        if filter_function(mol_copy) is False:
            filters_failed = filters_failed + 1

    if filters_failed == 0:
        return True

    # failed one or more filters
    return False
