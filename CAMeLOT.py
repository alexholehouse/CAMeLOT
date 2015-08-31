## CAMeLOT V 0.1.0
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
## Copyright 2015
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
#
# VERSION 0.1.0
# Copyright 2015 
#
# By Kiersten Ruff, Tyler Harmon, Alex Holehouse, and Rohit Pappu
#
# Paper Abstract:
#
#
#
# 
#-----------------------------------------------------------------------------------
#



from configs              import CAMELOT_VERSION
from extract_parameters   import Simulations
from keyfileParser        import KeyfileParser
from optimization         import Optimization
from utils                import print_logo




print_logo(CAMELOT_VERSION)

KeyFile        = KeyfileParser('ok')

SimObj         = Simulations(KeyFile)            
#(L,K, DH)      = SimObj.build_bond_lengths()
#params_by_res  = SimObj.build_dihedral_angles()
#damps          = SimObj.build_damping_parameters()
#masses         = SimObj.build_mass_parameters()
#masses         = SimObj.build_initial_starting_parameters()
#K,L            = SimObj.build_bond_angles()

OptObj = Optimization(KeyFile, SimObj.get_sequence_vector(), SimObj.get_charge_vector(), SimObj.get_residue_rgs())

OptObj.write_initial_sigma_eps()
OptObj.build_moltemplate_input()
OptObj.build_GaussianProcess_runner_code()
OptObj.build_GaussianProcess_function_code()
OptObj.build_GaussianProcess_update_code()



