# CALYPSO - Structural Analysis and Hierarchical Clustering of Solved Amyloid Structures

### Author: Jack P Connor
### Contact: bsjpc@leeds.ac.uk

### Last Updated 13-11-2025

<p align="center">
  <img src="Figures/Calypso.png" alt="img" width="600"/>
</p>

## ABSTRACT

Over 600 amyloid structures have been solved to date to near-atomic resolution. This has highlighted an enormous diversity of fibril structures conforming to the canonical cross-β amyloid fold. Using α-synuclein and tau amyloid structures as models, we show that they can be clustered into topologically distinct fold families. In these classes, the same, or similar, regions pair in different ways to generate topologies that can be hierarchically clustered. Despite their different topologies, the fibrils have similar stability, as determined by FoldX. The results provide a framework to classify newly solved fibril structures as belonging to an existing class or forming a new topological cross-β fold. Furthermore, this enables comparisons between fibrils found in disease and those formed in vitro. The workflow has been automated, enabling users to interrogate new amyloid structures as they emerge using this pipeline.

For more information please see the associated paper: https://doi.org/10.1016/j.str.2025.07.005

## Features
* Web scrape the [Amyloid Atlas](https://people.mbi.ucla.edu/sawaya/amyloidatlas/) to access known published structures
* Identify unique chains from each PDB file (intra-PDB variation)
* Calculate structural similarity between each unique chain using RMSD
* Hierarchical clustering of amyloid structures into defined groups based on their structural properties
* Compare stabilities for each amyloid structure using FoldX
* Define stabilising regions
* Analyse how the stabilising regions interact for each structural group identified by RMSD clustering
* Calculate the propensity for each residue in the ordered fibril core to be in a B-sheet

<p align="center">
  <img src="Figures/graphical_abstract.png" alt="img" width="600"/>
</p>

## INSTALLATION

### FoldX
Ensure FoldX is installed on your system and add 'export FOLDX_BINARY=/your/path/to/foldx.exe' to your .bashrc file.
For more information, follow the instuctions given by PyFoldx: https://github.com/leandroradusky/pyfoldx.

Important note - Under the academic licence, the FoldX license agreement ends on the 31st December of each year. To continue using FoldX the FoldX.exe needs to be redownloaded and added to your path as described in PyFoldX's installation.

### Conda Virtual Environment
conda env create -f environment.yaml
conda activate calypso

## Running the analysis
The analysis can be run from a GUI by running calypso.py where you can adjust run parameters.
Alternatively the analysis can be ran without a GUI using run_analyis.sh. In this script you can change the run parameters as desired.

## Using local/unpublished PDBs
Move any PDBs that are cannot be retrieved online using PyMol's fetch command into the Local folder for them to be included in the analysis.



