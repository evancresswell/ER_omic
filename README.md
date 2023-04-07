# ER_omic
Omic analysis using Expectation Reflection

### Repo Table of Contents
- [Anaconda Environment Setup](#Anaconda-Environment-Setup)
	- Setting up conda environment for ER and PYDCA simulations
- [Tutorial Directories](#Tutorial-Directories)
	- Outline of tutorials for ER omic analysis


## Protein Structure Prediction with DCA-Expectation Reflection
=======================

This repository provides tutorial examples for using Expectation Reflection for omic analysis. In the different jupyter notebooks you can run the different steps and applications of ER in Omic analysis including (see description of notebooks for topic specific tutorials):
* Data acuisition and preprocessing 
* Parallel computation of Direct Information and position coevolution
* Method result plotting
* Comparison with other popular methods
* Generation of Energy Landscape

Feel free to contact <evancresswell@gmail.com> or <vipulp@niddk.nih.gov > regarding questions on implementation and execution of repo

#### Package Requirements
- Anaconda/Miniconda - to run the proper enviornment for code. If using another python package manager see the following files for enviornment requirements 
    - Expectation Reflection Environment: DCA_ER_requirments.txt
    - PYDCA Enviornment: PDYCA_requirements.txt
    - Plotting Environment: plotting_requirements.txt

### Anaconda Environment Setup
[Back to Top](#Repo-Table-of-Contents)

#### Tutorial Directories
This repository shows users how to use Expectation Reflection for omic analysis of protein and genome sequences. The different applications are outlined in the 3 main directories.
- Protein Contact Prediction
	- PDB2MSA.ipynb -- Generates Direct Information (used to predict contacts protein tertiary structure) from coevolutionary information generated by ER
	- Method_Comparison.ipynb -- Generates plots to show efficacy of PCP by ER and other popular methods.
- SARS-CoV2 Genome-Wide Interaction
	- Outlines the process for large-scale alignment of SARS-CoV2 genomes sequences using Biowulf swarm simulations
- ER Energy Landscape
	- LINE1_landscape.ipynb -- Generates ER-defined energy landscapes for LINE1 evolution
		- analysis of pre-defined clusters
		- good visualization of clusters in PCA space
	- spike_landscape.ipynb -- Generates ER-defined energy landscapes for SARS-CoV2 spike protein
		- analysis of variants and within-variant spectral clusters
		- generation and visualization of ER-mutated sequences
		- entrenchment
- Utilities
	- [codon_mapping.py](https://github.com/evancresswell/ER_omic/blob/main/utilities/codon_mapping.py) -- Python module for mapping from genome sequence to resulting protein sequence of amino acids
	- [data_processing.py](https://github.com/evancresswell/ER_omic/blob/main/utilities/data_processing.py) -- Python module for processing Multiple Sequence Alignments (MSAs) in preparation of inferring coevolutionary couplings
	- [direct_info.py](https://github.com/evancresswell/ER_omic/blob/main/utilities/direct_info.py) -- Python module for calculating direct information from w/b coupling and local bias generated by Expectation Reflection (ER)
	- [er_energy.py](https://github.com/evancresswell/ER_omic/blob/main/utilities/er_energy.py) -- Python module for calculating sequence energy given sequence background (ie w/b from specific sequence group), also contains functions for mutating sequences using ER-inferred coupling (w)
	- [expectation_reflection.py](https://github.com/evancresswell/ER_omic/blob/main/utilities/expectation_reflection.py) -- Python module for calculatino the coevolutionary coupling matrix (w) and position bias (b) for a given group of processed sequences.
	- [tools.py](https://github.com/evancresswell/ER_omic/blob/main/utilities/tools.py) -- Python Module with miscellaneous functions for processing, computation and visualization of ER analysis
	
[Back to Top](#Repo-Table-of-Contents)

