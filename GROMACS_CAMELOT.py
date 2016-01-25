import mdtraj as md
import numpy as np

## GROMACS_CAMELOT creates an interface identical to CAMPARITraj for 
## dealing with GROMACS simulations instead of CAMPARI simulations. Rather
## Than creating a single analysis engine which can deal with both, MDTraj
## is ~99% of the way for vanilla GROMACS trajectories so we've just reproduced
## the CAMPARITraj interface here

class GMXTrajException(Exception):
    pass

class GMXTraj:

    def __init__(self, xtc_filename, pdb_filename):
        

        # read into a native MDTraj object
        self.MDObj = md.load(xtc_filename, top=pdb_filename)
    
        
