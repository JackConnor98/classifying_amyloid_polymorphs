# CALYPSO - Structural Analysis and Hierarchical Clustering of Solved Amyloid Structures

### Author: Jack P Connor
### Contact: bsjpc@leeds.ac.uk

### Last Updated 12-10-2025

<p align="center">
  <img src="Figures/Calypso.png" alt="img" width="600"/>
</p>

## ABSTRACT

More than 500 amyloid structures have been solved to date to near-atomic resolution. This has highlighted an enormous diversity of fibril structures conforming to the canonical cross-β amyloid fold. Using α-synuclein and tau amyloid structures as models, we show that they can be clustered into topologically distinct fold families. In these classes, the same, or similar, regions pair in different ways to generate topologies that can be hierarchically clustered. Despite their different topologies, the fibrils have similar stability, as determined by FoldX. The results provide a framework to classify newly solved fibril structures as belonging to an existing class or forming a new topological cross-β fold. Furthermore, this enables comparisons between fibrils found in disease and those formed in vitro. The workflow has been automated, enabling users to interrogate new amyloid structures as they emerge using this pipeline.

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
For more information, follow the instuctions given by PyFoldx: https://github.com/leandroradusky/pyfoldx

### Python Packages
pymol-open-source, os, pandas, itertools, PIL, math, numpy, scipy, string, matplotlib, sys, csv, bs4, requests, PyPDF2, re, gzip, shutil, pyfoldx 

### R Packages
dplyr, stringr, tidyr, readr, ggplot2, dendextend, RColourBrewer, tibble, bio3d, purrr, zoo, pastecs, ggpubr, ggtext, igraph, ggrepel, FSA, XML, xml2

## Running the analysis
The analysis pipeline is run in a Linux terminal using the script run_analysis.sh. In this script, you can change the run parameters as desired.



