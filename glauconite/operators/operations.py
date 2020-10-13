"""
Also filters and converts SMILES to 3d SDFS.
"""
import __future__

import os
import sys
import glob

import glauconite.operators.filter.execute_filters as Filter
import glauconite.operators.convert_files.conversion_to_3d as conversion_to_3d


def get_usable_format(infile):
    """
    This code takes a string for an file which is formatted as an .smi file. It
    opens the file and reads in the components into a usable list.

    The .smi must follow the following format for each line:
        MANDATORY INFO
            part 1 is the SMILES string
            part 2 is the SMILES name/ID

        Optional info
            part -1 (the last piece of info) is the SMILES diversity score
                relative to its population
            part -2 (the second to last piece of info) is the fitness metric
                for evaluating
                - For default setting this is the Docking score
                - If you add a unique scoring function Docking score should be
                    -3 and that score function should be -2

            Any other information MUST be between part 2 and part -2 (this
            allows for the expansion of features without disrupting the rest of the code)

    Inputs:
    :param str infile: the string of the PATHname of a formatted .smi file to
        be read into the program

    Returns:
    :returns: list usable_list_of_smiles: list of SMILES and their associated
        information formatted into a list which is usable by the rest of Autogrow
    """

    # IMPORT SMILES FROM THE PREVIOUS GENERATION
    usable_list_of_smiles = []

    if os.path.exists(infile) is False:
        print("\nFile of Source compounds does not exist: {}\n".format(infile))
        raise Exception("File of Source compounds does not exist")

    with open(infile) as smiles_file:
        for line in smiles_file:
            line = line.replace("\n", "")
            parts = line.split("\t")  # split line into parts separated by 4-spaces
            if len(parts) == 1:
                parts = line.split(
                    "    "
                )  # split line into parts separated by 4-spaces

            choice_list = []
            for i in range(0, len(parts)):
                choice_list.append(parts[i])
            usable_list_of_smiles.append(choice_list)

    return usable_list_of_smiles

#############
# Main run GlauconiteFilterer operators to make a generation
#############
def populate_generation(vars):
    """
    This will run all of the filters and if chosen convert to 3D.

    Inputs:
    :param dict vars: a dictionary of all user variables

    Returns:
    :returns: str full_generation_smiles_file: the name of the .smi file
        containing the new population
    :returns: list full_generation_smiles_list: list with the new population
        of ligands
    :returns: bool None: returns None twice if any step failed. This will
        result in the program ending
    """

    # Get the Source compound list. This list is the full population from
    # either the previous generations or if its Generation 1 than the its the
    # entire User specified Source compound list If either has a SMILES that
    # does not sanitize in RDKit it will be excluded and a printout of its
    # Name and SMILES string will be printed.

    seed_list = get_usable_format(vars["source_compound_file"])

    # Save source compounds list
    save_ligand_list(
        vars["output_directory"],
        seed_list,
        "Initial_SMILES",
    )

    # Run Glauconite
    passed_ligands = Filter.run_filter(vars, seed_list, vars["verbose"])
    passed_lig_names = [x[1] for x in passed_ligands]
    failed_filters = [x for x in seed_list if x[1] not in passed_lig_names]

    # save those which passed and those that failed
    save_ligand_list(
        vars["output_directory"],
        passed_ligands,
        "SMILES_Passed_All_Filters",
    )
    save_ligand_list(
        vars["output_directory"],
        failed_filters,
        "SMILES_Failed_Filter",
    )
    sys.stdout.flush()

    # CONVERT SMILES TO .sdf USING GYPSUM and convert .sdf to .pdb with rdkit
    # This will output sdf files into a folder. The .smi.0.sdf file is not a
    # valid mol, but all the others will be valid the 1st Smiles in the
    # original .smi file is saved as .smi.1.sdf and 2nd file is saved as
    # .smi.2.sdf
    if vars["convert_to_3D"] is True:
        smiles_to_convert_file = vars["output_directory"] + "SMILES_Passed_All_Filters.smi"
        conversion_to_3d.convert_to_3d(vars, smiles_to_convert_file, vars["output_directory"])
        get_list_of_3D_SMILES(vars, passed_ligands)

    sys.stdout.flush()

    return smiles_to_convert_file, passed_ligands

def get_list_of_3D_SMILES(vars, new_generation_smiles_list):
    """
    This will obtain and save the list of SMILES in the same order as
    found in NEW_SMILES.smi but with 3D variant information from PDBS.

    Will save to vars["output_directory"] +  "New_SMILES_After_3D_Conversion.smi"

    Inputs:
    :param dict vars: a dictionary of all user variables
    :param list new_generation_smiles_list: list of all 1/2D SMILES
    """
    list_of_3D_SMILES = []
    PDBs_dir = vars["output_directory"] + os.sep + "PDBs" + os.sep
    for mol_info in new_generation_smiles_list:
        short_id = mol_info[1].split(")")[-1]
        pdb_files = glob.glob(PDBs_dir + short_id + "__*.pdb")
        pdb_files.sort()
        for pdb_pose in pdb_files:
            base_info = os.path.basename(pdb_pose).replace(".pdb", "")
            with open(pdb_pose) as f:
                SMILES_string = f.readline().replace("\n", "")
                SMILES_string = SMILES_string.replace("REMARK Final SMILES string: ", "")
            list_of_3D_SMILES.append("\t".join([SMILES_string, mol_info[1], base_info]))

    # Save all info to a .smi file
    list_of_3D_SMILES = "\n".join(list_of_3D_SMILES)
    with open(vars["output_directory"] +  "New_SMILES_After_3D_Conversion.smi", "w") as f:
        f.write(list_of_3D_SMILES)

#############
# Saving Output files for generations and seeds
#############

def save_ligand_list(output_directory, list_of_chosen_ligands, nomenclature_tag):
    """
    Save the list of ligands. nomenclature_tag is a string such as "Mutation"
    or "Crossover" describing what this data is used
    for.

    Inputs:
    :param dict output_directory: the directory of the run to save the
    :param list list_of_chosen_ligands: The formatted list of ligands to seed
    :param str nomenclature_tag: The str describing the ligand list
        -ie seeding_mutations is the list that seeded the mutations while
            mutations would be the list of mutations generated from the
            seeding_mutations list
        -ie. mutation, crossover
    """

    # make a folder for the Seed files
    seed_folder_path = output_directory + os.sep

    # check if folders exist, if not make them
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    if not os.path.isdir(seed_folder_path):
        os.makedirs(seed_folder_path)

    output_file_name = "{}{}.smi".format(
        seed_folder_path, nomenclature_tag
    )

    # save to a new output smiles file. ie. save to ranked_smiles_file
    with open(output_file_name, "w") as output:
        for line in list_of_chosen_ligands:
            output_line = "\t".join(line) + "\n"
            output.write(output_line)

    sys.stdout.flush()
