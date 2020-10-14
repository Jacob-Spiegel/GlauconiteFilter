# GlauconiteFilter

[![DOI](https://zenodo.org/badge/303535253.svg)](https://zenodo.org/badge/latestdoi/303535253)

GlauconiteFilter is a free, open-source program for parallelized ADME-PK chemical filtration. GlauconiteFilter uses the AutoGrow4 filtration algorithm to SMILES strings. It is a highly parallelized algorithm that has been test on both SMP and MPI high performance computing environments. GlauconiteFilter provides nine predefined chemical filters, however, its intuitive plugin-style architecture allows for the easy expansion of filter options. GlauconiteFilter also provides the option to automate 1/2D SMILES to 3D conversion using Gypsum-DL.

GlauconiteFilter accepts a list of SMILES (in .smi format) and will output a list of compounds
that pass and fail a given filter. GlauconiteFilter also allows users to apply chemical property filters. Additionally, GlauconiteFilter options the automated conversion of the children from 1/2D SMILES to 3D representations (PDB and 3D-SDF) using Gypsum-DL.

![Image of Predefined_Filters](https://github.com/Jacob-Spiegel/GlauconiteFilter/blob/main/Figures/Predefined_Filters.png)

# Using GlauconiteFilter

GlauconiteFilter is primarily a python based program which can be run through command-line and GUI interface. More details can be found within `/PATH/GlauconiteFilter/tutorial/tutorial.md` or by running:

`python /PATH/GlauconiteFilter/RunGlauconiteFilter.py -h`

## Command-Line Interface

In a terminal with a modern Python environment run either:

`python /PATH/GlauconiteFilter/RunGlauconiteFilter.py -j /PATH/To/JSON_Paramerter_File.json`

or

`python /PATH/GlauconiteFilter/RunGlauconiteFilter.py --source_compound_file  /PATH/To/SMILES_File.smi`

Where `/PATH/To/JSON_Paramerter_File.json` is a JSON file containing all user parameters and `/PATH/To/SMILES_File.smi` is a file of all SMILES to be crossed


## GUI Interface

The GUI interface requires the additional dependency of GOOEY (https://github.com/chriskiehl/Gooey). To use the GUI, please run the following command from a terminal with a modern Python environment:

`python /PATH/GlauconiteFilter/RunGlauconiteFilter_GUI.py `


# Dependencies

GlauconiteFilter requires the following dependencies:
- RDKit
- func_timeout
- scipy
- numpy
- mpi4py (optional but required for MPI multiprocessing)
- Gooey (optional but required for GUI interface)

GlauconiteFilter was tested on a Python 3.7.7 environment on a Ubuntu OS. A docker image is provided for users with unsupported operating systems, such as Windows.

# Acknowledgment and Citing

Much of this code is take directly and/or adapted from AutoGrow4. This program also relies on Gypsum-DL and Dimorphite-DL for ligand handling and multiprocessing. Please remember to cite the following papers:

## GlauconiteFilter Citation:

- Spiegel, J.O. GlauconiteFilter: an open-source program for automated ADME-PK filtering. (2020) https://doi.org/10.5281/zenodo.4087648
[![DOI](https://zenodo.org/badge/303535253.svg)](https://zenodo.org/badge/latestdoi/303535253)

## AutoGrow4 Citation:

- Spiegel, J.O., Durrant, J.D. AutoGrow4: an open-source genetic algorithm for de novo drug design and lead optimization. J Cheminform 12, 25 (2020). https://doi.org/10.1186/s13321-020-00429-4

## Gypsum-DL Citation:

- Ropp, P.J., Spiegel, J.O., Walker, J.L. et al. Gypsum-DL: an open-source program for preparing small-molecule libraries for structure-based virtual screening. J Cheminform 11, 34 (2019). https://doi.org/10.1186/s13321-019-0358-3

## AutoGrow4 Dimorphite-DL:

- Ropp, P.J., Kaminsky, J.C., Yablonski, S. et al. Dimorphite-DL: an open-source program for enumerating the ionization states of drug-like small molecules. J Cheminform 11, 14 (2019). https://doi.org/10.1186/s13321-019-0336-9

