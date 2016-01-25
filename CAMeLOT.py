## CAMeLOT V 0.1.1
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
## Copyright 2015 - 2016
##

#-----------------------------------------------------------------------------------
# 
#  ________  ________  _____ ______   _______   ___       ________  _________   
# |\   ____\|\   __  \|\   _ \  _   \|\  ___ \ |\  \     |\   __  \|\___   ___\ 
# \ \  \___|\ \  \|\  \ \  \\\__\ \  \ \   __/|\ \  \    \ \  \|\  \|___ \  \_| 
#  \ \  \    \ \   __  \ \  \\|__| \  \ \  \_|/_\ \  \    \ \  \\\  \   \ \  \  
#   \ \  \____\ \  \ \  \ \  \    \ \  \ \  \_|\ \ \  \____\ \  \\\  \   \ \  \ 
#    \ \_______\ \__\ \__\ \__\    \ \__\ \_______\ \_______\ \_______\   \ \__\
#     \|_______|\|__|\|__|\|__|     \|__|\|_______|\|_______|\|_______|    \|__|
#
# 
#
#-----------------------------------------------------------------------------------
#
# VERSION 0.1.1
# Copyright 2015 
# 
# Coarse grain approach by Kiersten Ruff, Tyler Harmon, and Rohit Pappu
# Original code by Kiersten Ruff
# Current implementation by Alex Holehouse
# 
#
# Ruff, K. M., Harmon, T. S. & Pappu, R. V. 
# CAMELOT: A machine learning approach for coarse-grained simulations of aggregation of 
# block-copolymeric protein sequences. J. Chem. Phys. 143, 243123 (2015).
#
# 
# We report the development and deployment of a coarse-graining method that is well suited for computer
# simulations of aggregation and phase separation of protein sequences with block-copolymeric
# architectures. Our algorithm, named CAMELOT for Coarse-grained simulations Aided by MachinE
# Learning Optimization and Training, leverages information from converged all atom simulations that
# is used to determine a suitable resolution and parameterize the coarse-grained model. To parameterize
# a system-specific coarse-grained model, we use a combination of Boltzmann inversion, non-linear
# regression, and a Gaussian process Bayesian optimization approach. The accuracy of the coarsegrained
# model is demonstrated through direct comparisons to results from all atom simulations.
# We demonstrate the utility of our coarse-graining approach using the block-copolymeric sequence
# from the exon 1 encoded sequence of the huntingtin protein. This sequence comprises of 17 residues
# from the N-terminal end of huntingtin (N17) followed by a polyglutamine (polyQ) tract. Simulations
# based on the CAMELOT approach are used to show that the adsorption and unfolding of the wild
# type N17 and its sequence variants on the surface of polyQ tracts engender a patchy colloid like
# architecture that promotes the formation of linear aggregates. These results provide a plausible
# explanation for experimental observations, which show that N17 accelerates the formation of linear
# aggregates in block-copolymeric N17-polyQ sequences. The CAMELOT approach is versatile and
# is generalizable for simulating the aggregation and phase behavior of a range of block-copolymeric
# protein sequences.
#
#-----------------------------------------------------------------------------------
#


# CAMELOT imports
from configs              import CAMELOT_VERSION
from extract_parameters   import Simulations
from keyfileParser        import KeyfileParser
from optimization         import Optimization
from utils                import print_logo


# we have to import argparse to use it
import argparse 


print_logo(CAMELOT_VERSION)

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-keyfile", "-k", help="Filename of keyfile") 
    args = parser.parse_args()

    
    if args.keyfile:
        keyfilename = args.keyfile
    else:
        print "No keyfile passed"
        print ""
        print "Please pass a keyfile using"
        print "-k <keyfilename>"
        exit(1)



KeyFile        = KeyfileParser(keyfilename)

SimObj         = Simulations(KeyFile)            
(L,K, DH)      = SimObj.build_bond_lengths()
params_by_res  = SimObj.build_dihedral_angles()
damps          = SimObj.build_damping_parameters()
masses         = SimObj.build_mass_parameters()
ISP            = SimObj.build_initial_starting_parameters()
(K,L)          = SimObj.build_bond_angles()
SimObj.build_interresidue_distances()

OptObj = Optimization(KeyFile, SimObj.get_sequence_vector(), SimObj.get_charge_vector(), SimObj.get_residue_rgs())

OptObj.write_initial_sigma_eps()
OptObj.build_moltemplate_input()
OptObj.build_code()



