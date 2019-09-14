# CAMELOT.py

##### Last updated: September 2019
##### Version: 0.1.2

## Overview
`CAMELOT.py` is a Python package that automates the process of learning bonded and non-bonded interaction parameters for coarse-grained simulations. 

`CAMELOT.py` was written by [Alex Holehouse](www.holehouse.wustl.edu) during his time in the [Pappu lab](www.pappu.wustl.edu) and is based on ideas developed by Dr. Kiersten Ruff which were applied to systematically coarse-grain various systems in the original CAMELOT manuscript [1]. This version represents a Python re-write of the original code, and we plan to re-structure this code into a standard Python package.

In particular, CAMELOT provides a general, reference framework for applying GPBO to learn parameters for coarse-grained simulations and can be integrated to work 

## Application
Specifically, ``CAMELOT.py`` contains two independent stages.

1. Stage 1 analyzes all-atom simulations and extracts out bonded-terms using an inverse-Boltzmann approach to parameterize bond lengths, angles, and dihedrals. For more details on this see the accompanying manuscript by Ruff _et al._ [1]. 

2. Stage 2 integrates Gaussian Process Bayesian Optimization into the procedure of learning non-bonded interactions. The code utilizes the [Gaussian Process for Machine Learning](http://www.gaussianprocess.org/gpml/code/matlab/doc/) software package, although this is opaque to the users as CAMELOT.py autogenerates all of the underlying code needed to interact with and run the associated parameterization.

## Usage
``CAMELOT.py`` is fully developed in terms of the underlying models, although remains less user-friendly than it could/should be. We are working on this - if you need access to   `CAMELOT.py` now the safest thing is to contact Alex directly (alex dot holehouse at wustl dot edu). 


## References
[1] Ruff, K. M., Harmon, T. S. & Pappu, R. V. CAMELOT: A machine learning approach for coarse-grained simulations of aggregation of block-copolymeric protein sequences. J. Chem. Phys. 143, 243123 (2015).