# ER_omic
Omic analysis using Expectation Reflection

## Protein Structure Prediction with DCA-Expectation Reflection
=======================

This repository is supplementarty to the publication (PUBLICATION LINK). In the different jupyter notebooks you can run the different steps of our method including:
* Data acuisition and preprocessing 
* Parallel computation of Direct Information
* Method result plotting
* Comparison with other popular methods

Feel free to contact <evancresswell@gmail.com> or <vipulp@niddk.nih.gov > regarding questions on implementation and execution of repo

#### Package Requirements
- Anaconda/Miniconda - to run the proper enviornment for code. If using another python package manager see the following files for enviornment requirements 
    - Expectation Reflection Environment: DCA_ER_requirments.txt
    - PYDCA Enviornment: PDYCA_requirements.txt

### Table of Contents
- [Anaconda Environment Setup](#Anaconda-Environment-Setup)
	- Setting up conda environment for ER and PYDCA simulations
- [PDB to MSA Mapping](#PDB-to-MSA-Mapping)
	- Given a PDB structure we find the best matching MSA to infer connections
- [Expectation Reflection](#Expectation-Reflection)
- [Other Methods](#Other-Methods)
	- Introduction of interfacing with pydca (https://github.com/KIT-MBS/pydca)
	- Presentation/Plotting of MF-DCA and PLM-DCA methods
	- Comparison of methods (De
- [Results](#Results)
	- Result Drivers used to generate figures for papers.

#### PDB to MSA Mapping
[Back to Top](#Table-of-Contents)

