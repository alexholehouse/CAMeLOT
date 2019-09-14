## CAMeLOT V 0.1.2
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
## Copyright 2015 - 2017
##


import os
import sys
import itertools
from utils import *

import numpy as np
from numpy.core.umath_tests import inner1d
from scipy.stats import norm
import scipy.ndimage
import matplotlib.pyplot as plt
import mdtraj as md

import CTraj.CTTrajectory as CT
import CTraj.CTExceptions as CTException

from CAMeLOTExceptions import SimulationsException
from keyfileParser import KeyfileParser

# Comments on coding style
#
# Function names
#   I've tried to use a consistent naming scheme for function names. Specifically
#
#   -)  get_*     Carry out some logic and return an appropriate value.
#
#   -)  build_*   Carry out some logic and return an appropriate value as well as
#                 writing an appropriately formatted file. These are (generally)
#                 the functions which build the parameter files.
#
#
#

class Simulations:
        

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def __init__(self, KeyFileObj):
        """
        This is the initialization function for the parameter extraction process. This function sets all the relevant variables
        for parameter extraction, reads in the trajectories, determines what the sequences is and how long it is, and initializes
        a whole bunch of stuff.

        """


        

        # set key constants
        self.kB   = 0.0019872041               # kcal/mol/K
        
        self.eta  = (0.629e-3)*(1e-10)*(1e-15) # visocity of water
        self.amu2kg = 1.66053892e-27           # conversion from amu to kg

        # Read the keyfile info
        self.TEMP                       = KeyFileObj.TEMP                        # temperature
        self.SIMROOT                    = KeyFileObj.SIMROOT
        self.BOOTSTRAP_NUM_ITER         = KeyFileObj.BOOTSTRAP_NUM_ITER 
        self.BOOTSTRAP_SIZE             = KeyFileObj.BOOTSTRAP_SIZE
        

        # define filenames
        self.BOND_DEFINITION_FILE       = KeyFileObj.BOND_DEFINITION_FILE
        self.BOND_PARAMETER_FILE        = KeyFileObj.BOND_PARAMETER_FILE
        self.DIHEDRAL_DEFINITION_FILE   = KeyFileObj.DIHEDRAL_DEFINITION_FILE
        self.DIHEDRAL_PARAMETER_FILE    = KeyFileObj.DIHEDRAL_PARAMETER_FILE
        self.ANGLE_DEFINITION_FILE      = KeyFileObj.ANGLE_DEFINITION_FILE
        self.ANGLE_PARAMETER_FILE       = KeyFileObj.ANGLE_PARAMETER_FILE
        self.MBT_PARAMETER_FILE         = KeyFileObj.MBT_PARAMETER_FILE
        self.EBT_PARAMETER_FILE         = KeyFileObj.EBT_PARAMETER_FILE
        self.AAT_PARAMETER_FILE         = KeyFileObj.AAT_PARAMETER_FILE 
        self.AT_PARAMETER_FILE          = KeyFileObj.AT_PARAMETER_FILE
        self.BB13_PARAMETER_FILE        = KeyFileObj.BB13_PARAMETER_FILE  
        self.DAMPENING_PARAMETER_FILE   = KeyFileObj.DAMPENING_PARAMETER_FILE 
        self.MASS_PARAMETER_FILE        = KeyFileObj.MASS_PARAMETER_FILE 
        self.INITIAL_XYZ_FILE           = KeyFileObj.INITIAL_XYZ_FILE
        self.RES_RES_DISTANCE_FILE      = KeyFileObj.RES_RES_DISTANCE_FILE
        self.RES_RES_BINS               = KeyFileObj.RES_RES_BINS
        
        self.PLOT_BOND_HISTOGRAMS       = KeyFileObj.PLOT_BOND_HISTOGRAMS
        self.PLOT_ANG_HISTOGRAMS        = KeyFileObj.PLOT_ANG_HISTOGRAMS
        self.PLOT_DIHEDRAL_HISTOGRAMS   = KeyFileObj.PLOT_DIHEDRAL_HISTOGRAMS

        # create the 'plots' directory if it doesn't exist
        if self.PLOT_ANG_HISTOGRAMS or self.PLOT_BOND_HISTOGRAMS or self.PLOT_DIHEDRAL_HISTOGRAMS:
            if not os.path.exists('plots'):
                os.makedirs('plots')


        self.TRAJECTORY_FILE            = KeyFileObj.TRAJECTORY_FILE
        self.PDB_FILE                   = KeyFileObj.PDB_FILE

        # >> Determine the number of replica based on the number of directories (1...n)
        # =================================================================================================================
        self.n_replica = self.get_number_of_replicas()
        
        
        # >> Read in all replicas, and ensure we have at least 1!
        # =================================================================================================================
        # iterate through each replica and read in the trajectory
        self.replica_vector = []
        
        self.STDMessage("Reading in trajectories",msgType='PHASE')
        for i in xrange(1,self.n_replica+1):
            
            self.replica_vector.append(CT.CTTrajectory('%s/%i/%s'%(self.SIMROOT,i, self.TRAJECTORY_FILE),'%s/%i/%s'%(self.SIMROOT,i, self.PDB_FILE)))
                    
        if len(self.replica_vector) == 0:
            raise SimulationsException('Failed to find any simulations at %s. Ensure simulations are located in monotonically incremented directories and __START.pdb and __traj.xtc files are found in both directories' % self.SIMROOT)
        # =================================================================================================================



        # >> Check that all replicas have the same number of residues
        # =================================================================================================================
        # get the number of residues in the first trajectory (this must exist as we explicitly checked it does)
        exampleProtein = self.replica_vector[0].proteinTrajectoryList[0]

        nRes=len(exampleProtein.get_aminoAcidSequence())
            
        r=0                
        for replica in self.replica_vector:
            r=r+1
            self.STDMessage("Checking trajectory %i " % r, subMessage=True)
            if not nRes == len(replica.proteinTrajectoryList[0].get_aminoAcidSequence()):
                raise SimulationsException('Replica %i does not have the same number of residues as replica 1' % r)
        # =================================================================================================================

        # >> Determine if sequences are capped or not
        # =================================================================================================================

        # We DO NOT count CAPS as residues, but they are in the trajectory, so this is essentially defining a global offset
        # for the residue ID indices

        # if res 1 is ACE
        if exampleProtein.get_aminoAcidSequence()[0].split('-')[0] == 'ACE':
            self.start_resid = 1
        else:
            self.start_resid = 0

        # if final res is NME (the end_resid is INCLUSIVE - so (NME 
        if exampleProtein.get_aminoAcidSequence()[-1].split('-')[0] == 'NME' or exampleProtein.get_aminoAcidSequence()[-1].split('-')[0] == 'NAC':
            self.end_resid   = len(exampleProtein.get_aminoAcidSequence()) - 2
        else:
            self.end_resid   = len(exampleProtein.get_aminoAcidSequence()) - 1 

        # actual number of residues being examined
        self.n_residues = (self.end_resid - self.start_resid)+1

        # .......................
        # | build cycle vectors |
        # .......................
        #
        # Having estabilished if the sequence is capped or not we now define a bunch of cyle vectors.
        # Cycle vectors are useful because we constantly find ourself looping over the same two sets of vectors (i.e.
        # the resids from the perspective of the simulation and a 0 based index for Simulation internal values).
        # By constructing the cyclevectors here we have a pre-defined vector to cycle over depending on what we want
        # to do. This probably seems like overkill, but it makes avoiding things like off-by-one errors and other stupid
        # issues much easier
        #
        
        self.resVector = range(self.start_resid, self.end_resid+1)
        self.idxVector = range(0, len(self.resVector))
        self.doubleVec = zip(self.resVector,self.idxVector)

        # > Build some final datastructures **which depend on the preceding initialization**
        # =================================================================================================================
        self.charge_vector   = self.get_charge_vector()
        self.sequence_vector = self.get_sequence_vector()


        # > Finally 
        # =================================================================================================================
        self.logwriter = build_logwriter(KeyFileObj.LOGFILE)
            
        
        with open(KeyFileObj.LOGFILE, 'w') as fh:
            fh.write('')

        self.logwriter('CAMELOT initializing...')

        
       


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def get_charge_vector(self):
        """
        Returns a vector as long as the true amino acid sequence (i.e. without caps if they're present) where
        ASP and GLU are -1, ARG and LYS and 1, and every other residue is 0

        PRECONDITIONS: Must be run after the rest of the initialzation has been run!
        
        """

        charge_vector = []         
        AAs = self.replica_vector[0].proteinTrajectoryList[0].get_aminoAcidSequence()

        for res in self.resVector:
            if AAs[res].find('ASP') > -1 or AAs[res].find('GLU') > -1:
                charge_vector.append(-1)
            elif AAs[res].find('LYS') > -1 or AAs[res].find('ARG') > -1:
                charge_vector.append(1)
            else:
                charge_vector.append(0)
        return charge_vector



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def get_sequence_vector(self):
        """
        Returns a vector as long as the true amino acid sequence (i.e. without caps if they're present) each position
        is the 1 letter code for that amino acid.

        PRECONDITIONS: Must be run after the rest of the initialzation has been run!
        
        """

        sequence_vector = []

        AAs = self.replica_vector[0].proteinTrajectoryList[0].get_aminoAcidSequence()

        for res in self.resVector:
            sequence_vector.append(threeToOne(AAs[res].split('-')[0]))

        return sequence_vector



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def get_number_of_replicas(self):
        """
        Returns the number of directories that monotonically/linear increment from 1. Specifically, this assumes that your
        trajectory directory structure is as follows;

        self.SIMROOT --|
                       +-- 1
                       +-- 2
                       +-- 3
                       +-- .
                       +-- .
                       +-- n

        So this function iterates through the subdirectories of self.SIMROOT and until it cannot find the next directory
        from increasing the increment counter.

        """
        subDirs=next(os.walk(self.SIMROOT))[1]
        
        # iterate from 0 to +inf
        for i in itertools.count():

            # if i+1 is not in subD
            if str(i+1) not in subDirs:
                return i



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_damping_parameters(self):
        """
        This function determines and writes dampening parameters for the Brownian Dynamic simulations. Specifically, it calculates
        the water visocity at the simulation temperature, and then calculates a residue-by-residue dampening parameter based on
        the mass and the visocity.

        """

        self.STDMessage("Building custom dampening parameters...",msgType='PHASE')
        
        # get water visocity
        visocity       = self.get_water_visocity(self.TEMP)

        # get a vector of the masses of each residue
        residue_masses = self.build_mass_parameters(writeFile=False)

        # get a vector of the Rgs of each residue
        residue_rgs    = self.get_residue_rgs()

        # for each residue calculate the dampening parameter as
        # mass / (6*pi*eta*r) where
        # - mass is residue mass
        # - eta is solution visocity
        # - r is the residue radius of gyration ('size')
        
        raw_damp = []
        for idx in self.idxVector:
            raw_damp.append((residue_masses[idx] * self.amu2kg) / (6 * np.pi * visocity * residue_rgs[idx]))

        # get the most common damping value to use as a normalization constant
        # to elaborate - this code will find the residue most frequently observed, and use that value as a normalization
        # constant such that it becomes 1 and everything else becomes in reference to that residue. This is basically an 
        # implementation detail - could normalize relative to glycine or something else.
        norm = max(set(raw_damp), key=raw_damp.count)

        # recale based on the normalization parameter
        norm_damp = np.array(raw_damp)/norm

        # finally, we write to disk
        with open(self.DAMPENING_PARAMETER_FILE,'w') as fh:
            self.STDMessage("Building custom dampening parameters...", msgType='WRITING',subMessage=True)

            # the no newline is deliberate!
            for idx in self.idxVector:
                fh.write('scale %i %4.2f ' % (idx+1, norm_damp[idx]))

        self.STDMessage("Custom dampening parameters generated", msgType='SUCCESS', subMessage=True)

        return norm_damp
            

        
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def get_water_visocity(self, T):
        """
        Function to calculate the visocity of water in the correct units based on the provided temperature

        """

        # Returns the temperature dependent viscoity of water in kg per (ang * fs)
        #
        # Some math
        #
        # Visocity is reported in mPa * s ==> g / (m * s)
        #
        # We want it in our units, which are KG / (angstroms * fs) --> multiply the SI units of visocity y
        # (1e-3)*(1e-10)*(1e-15)

        units_prefactor = (1e-3)*(1e-10)*(1e-15)

        def viscocity_vs_temp(T):
            """ Function which takes in a temperature and returns visocity in standard units """
            #
            # Functional form which describes classic visocity behaviour of water
            # (Pa.S) = (2.414 * 10-5) * 10^(247.8/(T-140))

            # So we do this and *1000 to get the value in mPa * S
            #
            
            return 1000 * (2.414*10**(-5) * 10**(247.8/(T-140)))

        return viscocity_vs_temp(T)*units_prefactor
    


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_mass_parameters(self, writeFile=True):
        """
        Function to calculate the mass of each residue. Writes the masses to the MASS_PARAMETER_FILE 
        if writeFile is set to True (default).

        """

        if writeFile:
            self.STDMessage("Building mass parameters...")

        # get the CTraj protein object
        prot = self.replica_vector[0].proteinTrajectoryList[0]

        mass_vector = []
        
        # get the masses from the residues in that protein
        for res in self.resVector:
            mass_vector.append(prot.get_residue_mass(res))

        if writeFile:
            self.STDMessage("Writing mass parameters...", msgType='WRITING')

            # for each mass identified write to the mass definition file
            with open(self.MASS_PARAMETER_FILE,'w') as fh:

                for idx in self.idxVector:
                    fh.write('@atom:Res%i %4.1f\n' % (idx+1, mass_vector[idx]))

            self.STDMessage("Mass parameters generated", msgType='SUCCESS')
            
        
        return mass_vector



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_initial_starting_parameters(self):
        """
        Function to deterine the initial starting positios of residues for the CG simulation. This is done by taking the last
        frame of the first trajectory, and using the center of mass positions for each of the residues. 

        This is a convenient way to ensure we (hopefully!) have a non-classhing initial conformation, though we may want 
        to add a different way of doing this in the future?
        
        This initial starting conformation is written to file

        """

        self.STDMessage("Building initial starting parameters...", msgType='PHASE')

        # grab the first trajectory. Note we may want to allow for a specific trajectory
        # to be chosen to initialze? Possible update later on.
        prot = self.replica_vector[0].proteinTrajectoryList[0]

        # construct an empty center of mass (COM) vector
        COM_vector = np.zeros((self.n_residues, 3))
        
        # for each residue
        for (res, idx) in self.doubleVec:
            # set the COM position for each residue to the last frame from the
            # trajectory simulation
            COM_vector[idx] = prot.get_residue_COM(res)[-1]

        # recalibrate (COM positions come in Angstroms - need nanonmeters)
        COM_vector = COM_vector*10

        # write starting positions to file
        with open(self.INITIAL_XYZ_FILE,'w') as fh:
            self.STDMessage("Writing initial starting parameters...", msgType='WRITING',subMessage=True )

            for idx in self.idxVector:
                fh.write('$atom:Res%i $mol:... @atom:Res%i %1.1f %4.2f %4.2f %4.2f\n' %(idx+1, idx+1, self.charge_vector[idx], COM_vector[idx][0], COM_vector[idx][1], COM_vector[idx][2],))

        self.STDMessage("Initial starting parameters built", msgType='SUCCESS', subMessage=True)

        return COM_vector



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def get_residue_rgs(self):
        """
        Function which returns the a vector with the average radius of gyration for each residue
        """
        # get the CTraj protein object

        rg_vector = []
        
        # for each residue
        for res in self.resVector:
            tmp = []

            # in each replica
            for replica in self.replica_vector:                
                prot = replica.proteinTrajectoryList[0]             

                # note the /10.0 converts from Angstroms to nm
                tmp.append(prot.get_radius_of_gyration(res,res)/10.0)

            # calculate mean associated with the residue from all replicas
            rg_vector.append(np.mean(tmp))


        # return rg_vector indexes from 0
        return 10*np.array(rg_vector)
        


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_bond_lengths(self):
        """
        Function which determines the optimal bond length and stiffness based on a Boltzmann Inversion approach.

        Writes values out to file as well as return L and K
        L = bond lengths
        K = stiffness of bond
        
        """

        self.STDMessage("Extracting center of mass distances from atomistic simulation replicas...")
        

        i_to_i1_distance      = []
        i_to_i1_distance_sqrd = []

        # We cycle through to the penultimate residue (no bond from ultimate to something)
        for res in self.resVector[0:-1]:
            
            sys.stdout.write('.')
            sys.stdout.flush()
            tmp_distances = np.array([])
                
            # for each replica get the distance between the res-th and the subsequent residue
            for replica in self.replica_vector:
                
                prot = replica.proteinTrajectoryList[0]

                # /10.0 to convert into nm
                tmp_distances = np.concatenate((tmp_distances, prot.get_interResidueCOMDistance(res,res+1)/10.0))

            
                                                                        
            i_to_i1_distance.append(tmp_distances)
            i_to_i1_distance_sqrd.append(np.square(tmp_distances))

        # > Now carry out bootstrapped subsampling
        # ==============================================================================
        sample_size = i_to_i1_distance[0].shape[0]

        # lengths
        L = np.zeros((self.BOOTSTRAP_NUM_ITER, self.n_residues-1))

        # stiffness
        K = np.zeros((self.BOOTSTRAP_NUM_ITER, self.n_residues-1))

        # distance histogram
        DH = np.zeros((self.BOOTSTRAP_NUM_ITER, self.n_residues, 100))
    
        # ALL DATA doesn't actually get used, but is returned and may be
        # useful for a comparison of how bootstrapping vs. straight at averaging 
        # does
        ALL_DATA = [np.array([]) for _ in range(self.n_residues-1)]        
        
        print ""
        self.STDMessage('Creating summary histogram fits...')
        for bts_iteration in xrange(0,self.BOOTSTRAP_NUM_ITER):

            # generate a vector of random indicies
            random_index = np.random.randint(0, sample_size,(self.BOOTSTRAP_SIZE))

            # note -1 because we're looking at *bonds* not residues!
            tmpdata = []
            for idx in self.idxVector[0:-1]:
                
                # extract mean distance between centers of mass for this residue based on a random
                # selection of frames
                L[bts_iteration][idx] = np.mean(i_to_i1_distance[idx][random_index])

                # extract apparent stiffness based on width of standard devaition
                K[bts_iteration][idx] = (self.kB*self.TEMP)/(2*(np.mean(i_to_i1_distance_sqrd[idx][random_index]) - np.square(L[bts_iteration][idx])))

                # extract histogram - NOTE that the histogram boundaries are hardcoded as 0 to 1 nm, which is a pretty
                # good range and for single bead per residue shouldn't ever be too small, BUT if we go to higher resolution
                # in later iterations may need to increase the max
                (DH[bts_iteration][idx], tmp) = np.histogram(i_to_i1_distance[idx][random_index],np.arange(0,1.01,0.01))
                
                # add the results from this bootstrap sample to the ALL_DATA vector, which is collecting all the data
                # used for bootstrapping calculations for reside idx
                ALL_DATA[idx] = np.concatenate((ALL_DATA[idx], i_to_i1_distance[idx][random_index]))
                
        # if we're generating plots...
        if self.PLOT_BOND_HISTOGRAMS:

            # write summary histograms fitted to normal distribution for manual inspection
            self.STDMessage('Creating summary histogram fits...',msgType='WRITING')
        
            for idx in self.idxVector[0:-1]:
                # histogram from 0 to 10 angstroms (0 to 1 nm) with 0.1 angstrom bins
                # tmp is just the bins so gets ignored and overwritten every time
                
                # fit to a normal distribution
                bin_cs  = np.arange(0,1.01,0.01)
                mu, std = norm.fit(ALL_DATA[idx])

                # plot histogrammed data
                plt.hist(ALL_DATA[idx], bins=bin_cs, normed=True,alpha=0.3, color='b')

                # plot normal fit
                p = norm.pdf(bin_cs, mu, std)
                plt.plot(bin_cs, p, 'k', linewidth=2)                

                # add nice titles and save!
                plt.title("COM bond fit results: mu = %.2f,  std = %.2f" % (mu, std))
                plt.xlabel('Distance (nm)')
                plt.ylabel('Normalized Count')
                plt.savefig('plots/BONDS_res_%i.png'%(idx),dpi=150)
                plt.close()

        # calculate mean and standard error on parameters
        L_mean = np.mean(L,0)
        K_mean = np.mean(K,0)

        # yeah right now we don't actually use the standard error for anything - perhaps this
        # but written out for analysis?
        L_STDER = np.std(L,0)/np.sqrt(self.BOOTSTRAP_NUM_ITER)
        K_STDER = np.std(K,0)/np.sqrt(self.BOOTSTRAP_NUM_ITER)

        # >> Write bonded parameters out to file
        self.STDMessage('Writing bond definition to [%s]...'%self.BOND_DEFINITION_FILE,msgType='WRITING')

        with open(self.BOND_DEFINITION_FILE,'w') as fh:
            for idx in self.idxVector[0:-1]:
                i=idx+1
                j=idx+2
                fh.write('$bond:Res%iRes%i @bond:Res%iRes%i $atom:Res%i $atom:Res%i\n' %(i,j,i,j,i,j) )

        self.STDMessage('Writing bond parameters to [%s]...'%self.BOND_PARAMETER_FILE,msgType='WRITING')
        with open(self.BOND_PARAMETER_FILE,'w') as fh:
            for idx in self.idxVector[0:-1]:
                i=idx+1
                j=idx+2

                # note we /100 and * 10 so the values come out as angstroms and angstrom appropriate (100 because K comes
                # from an L squared term so have to account for squared units). Note this assumes the distances were in nm
                # before (so we convert from A->nm->A - yes this is silly, legacy code is a pain!)
                fh.write('bond_coeff @bond:Res%iRes%i %4.2f %4.2f \n' %(i,j,K_mean[idx]/100,L_mean[idx]*10) )

        self.STDMessage("Bonded parameters generated", msgType='SUCCESS')
                        
        return (L, K, ALL_DATA)
        


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_bond_angles(self):
        """
        Function which determines the optimal bond angles and stiffness based on a Boltzmann Inversion approach to
        give an ideal angle value and stifness.

        Writes values out to file as well as return L and K
        L = bond lengths
        K = stiffness of bond
        
        """


        ## ----------------------------------------------------------------------------
        ## STAGE 1 - EXTRACT THE BOND ANGLES FROM THE ALL ATOM TRAJECTORIES
        ##

        angle_by_res= []
        self.STDMessage('Extracting bond angles from atomistic simulation replicas...', msgType='STATUS')

        # cycle through each residue extracting the appropriate theta vector on the
        # COM of groups of 3 residues
        L=[]
        K=[]
        for res in self.resVector[0:-2]:
            tmp_angles = np.array([])
            sys.stdout.write('.')
            sys.stdout.flush()

            vals_temp = np.array([])
            for replica in self.replica_vector:                            
                prot = replica.proteinTrajectoryList[0]
                i  = res
                j  = i + 1
                k  = j + 1

                # extract COM vector between the 3 residues
                # (so get 2 vectors)
                # *10 so we move into Angstroms now...
                v1 = prot.get_interResidueCOMVector(i,j)*10
                v2 = prot.get_interResidueCOMVector(k,j)*10

                ang_denominator = np.linalg.norm(v1,axis=1)*np.linalg.norm(v2,axis=1)

                ang_vector  = (180/np.pi)*np.arccos(inner1d(v1,v2) / ang_denominator)

                vals_temp = np.concatenate((vals_temp, ang_vector))

            # calculate the mean value from all replicas
            # NOTE we're not bootstrapping here - as far as I can see it doesn't actually buy us but
            # this might not be true. Implementing a bootstrapping step here would be easy though
            vals_mean = np.mean(vals_temp)
            
            L.append(vals_mean)
            
            # get angle stiffness
            K.append( (np.square(180/np.pi))*(self.kB*self.TEMP)/(2*(np.mean(np.square(vals_temp)) - (vals_mean**2))))

            # if we're generating plots...
            if self.PLOT_ANG_HISTOGRAMS:
                # fit to a normal distribution (not 0 180 is always gonna be correct)
                bin_cs=np.arange(0,181,1)
                mu, std = norm.fit(vals_temp)

                # plot histogrammed data
                plt.hist(vals_temp, bins=bin_cs, normed=True,alpha=0.3, color='b')
            
                # plot normal fit
                p = norm.pdf(bin_cs, mu, std)
                plt.plot(bin_cs, p, 'k', linewidth=2)                
            
                # add nice titles and save!
                plt.title("Bond angle fit results: mu = %.2f,  std = %.2f" % (mu, std))
                plt.xlabel('Degrees')
                plt.ylabel('Probability')
                plt.savefig('plots/ANG_res_%i.png'%(res),dpi=150)
                plt.close()

        ## ----------------------------------------------------------------------------
        ## STAGE 2 - Write the derived parameters out to file
        ##
        
        self.STDMessage('Writing bond angle definition to [%s]...'%self.ANGLE_DEFINITION_FILE,msgType='WRITING')

        with open(self.ANGLE_DEFINITION_FILE,'w') as fh:
            for idx in self.idxVector[0:-2]:
                i=idx+1
                j=idx+2
                k=idx+3
                fh.write('@angle:Res%iRes%iRes%i @atom:Res%i @atom:Res%i @atom:Res%i @bond:Res%iRes%i @bond:Res%iRes%i\n'  % (i,j,k,i,j,k,i,j,j,k,))

        self.STDMessage('Writing dihedral parameters to [%s]...'%self.ANGLE_PARAMETER_FILE,msgType='WRITING')

        with open(self.ANGLE_PARAMETER_FILE,'w') as fh:
            for idx in self.idxVector[0:-2]:
                i=idx+1
                j=idx+2
                k=idx+3
                fh.write('angle_coeff  @angle:Res%iRes%iRes%i %4.2f %4.2f\n'  % (i,j,k, K[idx], L[idx]))
        
        self.STDMessage("Dihedral parameters generated", msgType='SUCCESS')

        return (K,L)
            




    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_dihedral_angles(self):
        """
        Function which determines the optimal dihedral angles by fitting them to a Fourrier Sum. This is useful because one
        of the bond functional forms in LAMMPS takes uses 3 (?) parameter Fourier Sum.

        Bit of a monster function, but works pretty well. Basically, I implemented an inteligent iterative approach to the Fourier
        sum fitting because we need the parameters to be within some limits, but by default scipy doesn't allow constrains
        to be placed on the fitting algrorithm, so I had to implement my own cost function to get around this. Doesn't give *exactly*
        the same answer as the MATLAB implementation, but does pretty well and is clearly close enough!


        """
        
        
        ## ----------------------------------------------------------------------------
        ## STAGE 1 - EXTRACT THE DIHEDRAL ANGLES FROM THE ALL ATOM TRAJECTORIES
        ##

        theta_by_res= []
        self.STDMessage('Extracting dihedral angles from atomistic simulation replicas...', msgType='STATUS')

        # cycle through each residue extracting the appropriate theta vector on the
        # COM of groups of 4 residues
        for res in self.resVector[0:-3]:
            tmp_theta = np.array([])
            sys.stdout.write('.')
            sys.stdout.flush()

            for replica in self.replica_vector:                            
                prot = replica.proteinTrajectoryList[0]
                i  = res
                j  = i + 1
                k  = j + 1
                l  = k + 1                

                # extract COM vector between the 4 residues
                # (so get 3 vectors)
                # x10 so we move into Angstroms now...
                b1 = prot.get_interResidueCOMVector(j,i)*10
                b2 = prot.get_interResidueCOMVector(k,j)*10
                b3 = prot.get_interResidueCOMVector(l,k)*10
                
                n1_numerator = np.cross(b1,b2)
                n1_denominator = np.linalg.norm(np.cross(b1,b2),axis=1)

                n2_numerator = np.cross(b2,b3)
                n2_denominator = np.linalg.norm(np.cross(b2,b3),axis=1)

                n1 = n1_numerator / np.transpose([n1_denominator,n1_denominator,n1_denominator])
                n2 = n2_numerator / np.transpose([n2_denominator,n2_denominator,n2_denominator])
                
                b2n = b2 / np.transpose([np.linalg.norm(b2,axis=1),np.linalg.norm(b2,axis=1),np.linalg.norm(b2,axis=1)])
                m1 = np.cross(n1,b2n)

                # avoids calculating the full dot product and then getting the diagonal - this is a 
                # really efficient vectorized way to get this (seriously this shaves multiple seconds per replica)
                
                x = inner1d(n1,n2)
                y = inner1d(m1,n2)

                tmp_theta = np.concatenate((tmp_theta, 180/np.pi*np.arctan2(y,x)))

            # NOTE!
            # We multiply by -1 because this gives a dihedral form where the chirality of the final CG model
            # matches the all atom model (note that this is a relective operations around 0). This, in itself
            # is kind of interesting!
            theta_by_res.append(-1*tmp_theta)

        # this to add a new line
        print ""
        ## ----------------------------------------------------------------------------
        ## STAGE 2 - Fit angle histograms to 3 term Fourier Series equation 
        ##
        
        # For each residue we
        # 1) Histogram the data and generate an empyrical PDF from that histogram 
        
        # ugh...
        bins_edges  = np.arange(-180,181,1)
        bin_centers = np.arange(-179.5,180.5,1)

        n_angles = len(theta_by_res)
        print ""
        self.STDMessage('Fitting dihedral angles to Fourier Series function...', msgType='STATUS')
        params_by_res = []
        for i in xrange(0, n_angles):  
            self.STDMessage('Fitting angle along %i-%i-%i-%i vector.'%(i, i+1, i+2, i+3), msgType='STATUS')
            sys.stdout.write('.')
            sys.stdout.flush()
            
            # histogram the data (density means our histogram Y axis is now in 
            # density units (0 < density < 1) rather than count)
            (vals,b)=np.histogram(theta_by_res[i],bins_edges,density=True)

            # Convert probability into energy and normalize onto some scale such 
            # that the lowest minima is 0
            vals_EN = -np.log(vals)*self.kB*self.TEMP - np.min(-np.log(vals)*self.kB*self.TEMP)

            # set any values which have become inft to the max observed energy
            vals_EN[ vals_EN == np.inf] = max(vals_EN[np.isfinite(vals_EN)])

            # 2) Fit that data to a 3 term Fourier series
        
            # lower and upper bounds
            LB = [0.0, -360.0,  0.0, -360.0,  0.0, -540.0]
            UB = [10.0, 360.0, 10.0,  360.0, 10.0,  540.0]
        
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            # Internal function which returns a custom optimizer function
            #
            def make_optimization_function(LB_setter, UB_setter, penalty_multiplier_setter):
            
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                def funct(x, P1, P2, P3, P4, P5, P6):
                    # define our 3 term Fourier Series function (T1/2/3 are the three terms)
                    T1 = P1*(1 - np.cos(np.deg2rad(1*x - P2)))
                    T2 = P3*(1 - np.cos(np.deg2rad(2*x - P4)))
                    T3 = P5*(1 - np.cos(np.deg2rad(3*x - P6)))
                
                    # when make_optimization_function returns 'funct' the following
                    # variables have been defined BY the make_optimization_function
                    # - functional programming FTW!
                    LB = LB_setter
                    UB = UB_setter
                    penalty_multiplier = penalty_multiplier_setter                
                
                    # we have to define limits manually because SciPy does not allow 
                    # automatic parameter constraints. To do this we define a flat bottom
                    # function for each parameter which = 0 within the range we care about but grows
                    # rapidly outside that range. How rapidly depends on the penalty multiplier
                
                    penalty = 0
                
                    # IDX sets the index into the LB and UB vector, which specifies
                    # bounds for each parameter
                    IDX=0
                    for P in [P1, P2, P3, P4, P5, P6]:                
                    
                        if P > UB[IDX]:
                            penalty=(P-UB[IDX])*penalty_multiplier + penalty
                        
                        if P < LB[IDX]:
                            penalty=(LB[IDX]-P)*penalty_multiplier + penalty
                        
                        IDX=IDX+1
                    
                    # sum with penality
                    return T1 + T2 + T3 + penalty
                #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            

                return funct
            #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
            # now we iterativly do this with a more and more painful penality_multiplier 
            penalty_multiplier=10000
            boundary_fail = True
            while boundary_fail:
            
                # construct a cost function with a defined penality value
                threeParamFunct = make_optimization_function(LB, UB, penalty_multiplier)
            
                # set initial goodness and past parameter set
                goodness    = 1000000000000000 
                best_params = [-1,-1,-1,-1,-1,-1] 
                                        
                # set initial guesses for parameters
                np.random.seed()
                guess       = [0.5, 90, 0.5, 90, 0.5, 90]
                
                # try 20 different runs with random starting positions for the 0-10 parameters
                for run in xrange(0,20):
                    try:

                        # optimizie using our customized function with a set penality term, where we *only* optimize 
                        # the 6 parameters defined in guess
                        (params, tmp) = scipy.optimize.curve_fit(threeParamFunct, bin_centers, vals_EN, p0=guess,maxfev=100000)
                        
                        # using parameters see how well we did by determinig the goodness of fit to the emprical histogram
                        tmpgood = sum(abs((threeParamFunct(bin_centers, params[0],params[1],params[2],params[3],params[4],params[5])-vals_EN)))

                        # if this fit was better than the previous one update the current gold standard (goodness) and the current best
                        # parameter set (best_params)
                        if tmpgood < goodness:
                            goodness = tmpgood
                            best_params = params

                        # New random parameters to try again with
                        # (while the coefficients *can* go 0 to 10 we find better results come from lower intial guesses (i.e constraining
                        # the initial guess to the 0 to 3 interval)
                        guess = [np.random.rand()*3, 90, np.random.rand()*3, 90, np.random.rand()*3, 90]
                    except RuntimeError:
                        self.STDMessage('Runtime error fitting Fourier series to angle starting on r %i'%i, msgType='WARNING')
                        self.STDMessage("Don't worry - we'll try again with different initial parameters", msgType='WARNING')

                    
                # set paramst to the best parameters you saw
                params = best_params

                # check ALL parameters lie inside the boundaries boundaries
                IDX = 0
                old_pm = penalty_multiplier
                for P in params:
                    if P > UB[IDX]:
                        penalty_multiplier=penalty_multiplier*10
                        break
                    if P < LB[IDX]:
                        penalty_multiplier=penalty_multiplier*10
                        break
                    IDX=IDX+1
                
                # if we didn't update the penality multiplier then all parameters looked good
                if old_pm == penalty_multiplier:
                    boundary_fail=False

            # if shit has hit the fan (can't imagine this happening)
            if best_params[0] == -1:
                raise SimulationsException("FUNDEMENTAL ERROR: Despite our best effors we could not fit the dihedral associated with the %i-%i-%i-%i stretch to a 3-term Fourier Series. This is either indicative of a bug in how we do the fitting, or a more fundemental issue. Note this doesn't mean we could't get a  *good* fit, it means we literally couldn't fit it at all. Please contact alex.holehouse@wustl.edu because this is bad news!")


            # save the parameters!
            params_by_res.append(params)

            # plot fit (note need to make this optional!)
            if self.PLOT_DIHEDRAL_HISTOGRAMS:
                plt.plot(bin_centers,threeParamFunct(bin_centers, params[0],params[1],params[2],params[3],params[4],params[5],))
                plt.plot(bin_centers, vals_EN)
                plt.title("Theta angle fit (goodness = %4.2f" % goodness)
                plt.xlabel('Angle (degrees)')
                plt.ylabel('Probability')
                plt.savefig('plots/DIHEDRAL_res_%i.png'%(i),dpi=150)
                plt.close()
        
        ## ----------------------------------------------------------------------------
        ## STAGE 3 - Write the derived parameters out to file
        ##
        
        # write summary histograms fitted to normal distribution for manual inspection
        self.STDMessage('Writing dihedral definition to [%s]...'%self.DIHEDRAL_DEFINITION_FILE,msgType='WRITING')

        with open(self.DIHEDRAL_DEFINITION_FILE,'w') as fh:
            for idx in self.idxVector[0:-3]:
                i=idx+1
                j=idx+2
                k=idx+3
                l=idx+4
                fh.write('@dihedral:Res%iRes%iRes%iRes%i @atom:Res%i @atom:Res%i @atom:Res%i @atom:Res%i @bond:Res%iRes%i @bond:Res%iRes%i @bond:Res%iRes%i\n'  % (i,j,k,l,i,j,k,l,i,j,j,k,k,l))

        self.STDMessage('Writing dihedral parameters to [%s]...'%self.DIHEDRAL_PARAMETER_FILE,msgType='WRITING')
        with open(self.DIHEDRAL_PARAMETER_FILE,'w') as fh:
            for idx in self.idxVector[0:-3]:
                i=idx+1
                j=idx+2
                k=idx+3
                l=idx+4
                fh.write('dihedral_coeff  @dihedral:Res%iRes%iRes%iRes%i %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n'  % (i,j,k,l, params_by_res[idx][0],params_by_res[idx][1],params_by_res[idx][2],params_by_res[idx][3],params_by_res[idx][4],params_by_res[idx][5]))
        
        self.write_other_dihedral_angles(n_angles)

        self.STDMessage("Dihedral parameters generated", msgType='SUCCESS')
                                
        return params_by_res


        
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def write_other_dihedral_angles(self, n_angles):
        """
        Function which writes out a bunch of empty files which are required by LAMMPS to use the 3-term Fourier series

        """

        ## Middle bond torsional coefficients 
        with open(self.MBT_PARAMETER_FILE,'w') as fh:
            for idx in self.idxVector[0:-3]:
                i=idx+1
                j=idx+2
                k=idx+3
                l=idx+4
                fh.write('dihedral_coeff  @dihedral:Res%iRes%iRes%iRes%i mbt %4.1f %4.1f %4.1f %4.1f \n'  % (i,j,k,l, 0,0,0,0))

        ## end bond torsional coefficients 
        with open(self.EBT_PARAMETER_FILE,'w') as fh:
            for idx in self.idxVector[0:-3]:
                i=idx+1
                j=idx+2
                k=idx+3
                l=idx+4
                fh.write('dihedral_coeff  @dihedral:Res%iRes%iRes%iRes%i ebt %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n'  % (i,j,k,l, 0,0,0,0,0,0,0,0))

        ## angle-angle torsion coefficients
        with open(self.AAT_PARAMETER_FILE,'w') as fh:
            for idx in self.idxVector[0:-3]:
                i=idx+1
                j=idx+2
                k=idx+3
                l=idx+4
                fh.write('dihedral_coeff  @dihedral:Res%iRes%iRes%iRes%i aat %4.1f %4.1f %4.1f\n'  % (i,j,k,l, 0,0,0,))

        ## angle torsion coefficients
        with open(self.AT_PARAMETER_FILE,'w') as fh:
            for idx in self.idxVector[0:-3]:
                i=idx+1
                j=idx+2
                k=idx+3
                l=idx+4
                fh.write('dihedral_coeff  @dihedral:Res%iRes%iRes%iRes%i at %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f %4.1f\n'  % (i,j,k,l, 0,0,0,0,0,0,0,0))

        ## bond-bond 1-3 torsion coefficients
        with open(self.BB13_PARAMETER_FILE,'w') as fh:
            for idx in self.idxVector[0:-3]:
                i=idx+1
                j=idx+2
                k=idx+3
                l=idx+4
                fh.write('dihedral_coeff  @dihedral:Res%iRes%iRes%iRes%i bb13 %4.1f %4.1f %4.1f\n'  % (i,j,k,l, 0,0,0,))


    def build_interresidue_distances(self):
        """
        Function which generates a CSV file with the full inter-residue distances (i.e. a redundant set of data). Note
        this calculates ALL the distances...
        
        """
        # to avoid building up a huge memory burden we write each residue by residue 

        self.STDMessage("Generating inter-residue distances",msgType='PHASE')

        # remove any content which existed before
        with open(self.RES_RES_DISTANCE_FILE, 'w') as fh:
            fh.write('')

        
            
        self.logwriter('inter-residue distance bins:\n%s'%(str(self.RES_RES_BINS/10.0)))
        
            
        # cycle through each residue
        for res in self.resVector[0:-1]:
            print "On res %i"%res
            
            # and the cycle through the same residues
            for res2 in xrange(res+1, self.resVector[-1]+1):
                print "----> On res2 %i"%res2

                # and each over tmp
                tmp_distances = np.array([])
                for replica in self.replica_vector:

                    prot = replica.proteinTrajectoryList[0]
                    
                    # note the /10.0 is so we convert from Angstroms to nm
                    tmp_distances = np.concatenate((tmp_distances, prot.get_interResidueCOMDistance(res, res2)/10.0))

                # and build a histogram with the defined binsize
            
                # we we divide the bins by 10 because the distances are in nanometers...            
                (vals, b) = np.histogram(tmp_distances, self.RES_RES_BINS/10.0)


                # convert into pmf
                vals = vals / float(np.sum(vals))                

                with open(self.RES_RES_DISTANCE_FILE, 'a') as fh:
                    fh.write("%i, %i, " % (res, res2))
                    np.savetxt(fh, vals, delimiter=',', newline=', ', fmt="%4.4f")
                    fh.write("\n")


                



        
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    def STDMessage(self, msg,msgType='STATUS', subMessage=False):
        if subMessage:
            print ">    [%s] : %s" %(msgType, msg)
        else:
            print "> [%s] : %s" %(msgType, msg)

                

                



