## CAMeLOT V 0.1.0
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
##

from CGParameterGroup import CGParameterGroup
from configs import CAMELOT_VERSION

class KeyfileParser:

    def __init__(self, keyFileName):
        
        keywordDict = self.parse_default(keyFileName)

        self.CAMELOT_VERSION             = CAMELOT_VERSION
        

        self.SIMROOT                     = keywordDict['SIMROOT']
        self.PYTHON_BIN                  = keywordDict['PYTHON_BIN']
        self.LAMMPS_BIN                  = keywordDict['LAMMPS_BIN']


        self.BOOTSTRAP_NUM_ITER          = keywordDict['BOOTSTRAP_NUM_ITER']
        self.BOOTSTRAP_SIZE              = keywordDict['BOOTSTRAP_SIZE']

        self.PLOT_BOND_HISTOGRAMS        = keywordDict['PLOT_BOND_HISTOGRAMS']
        self.PLOT_ANG_HISTOGRAMS         = keywordDict['PLOT_ANG_HISTOGRAMS']
        
        self.INITIAL_XYZ_FILE            = keywordDict['INITIAL_XYZ_FILE']
        self.MASS_PARAMETER_FILE         = keywordDict['MASS_PARAMETER_FILE']

        self.BOND_DEFINITION_FILE        = keywordDict['BOND_DEFINITION_FILE']
        self.BOND_PARAMETER_FILE         = keywordDict['BOND_PARAMETER_FILE']
        self.DIHEDRAL_DEFINITION_FILE    = keywordDict['DIHEDRAL_DEFINITION_FILE']
        self.DIHEDRAL_PARAMETER_FILE     = keywordDict['DIHEDRAL_PARAMETER_FILE']
        self.ANGLE_DEFINITION_FILE       = keywordDict['ANGLE_DEFINITION_FILE']
        self.ANGLE_PARAMETER_FILE        = keywordDict['ANGLE_PARAMETER_FILE']
        self.TRAJECTORY_FILE             = keywordDict['TRAJECTORY_FILE']
        self.PDB_FILE                    = keywordDict['PBD_FILE']
        self.MBT_PARAMETER_FILE          = keywordDict['MBT_PARAMETER_FILE']
        self.EBT_PARAMETER_FILE          = keywordDict['EBT_PARAMETER_FILE']
        self.AAT_PARAMETER_FILE          = keywordDict['AAT_PARAMETER_FILE']
        self.AT_PARAMETER_FILE           = keywordDict['AT_PARAMETER_FILE']
        self.BB13_PARAMETER_FILE         = keywordDict['BB13_PARAMETER_FILE']
        self.DAMPENING_PARAMETER_FILE    = keywordDict['DAMPENING_PARAMETER_FILE']
        
        

        # optimization keywords
        self.OPT_KAPPA                   = keywordDict['OPT_KAPPA']     
        self.OPT_LJ_CUTOFF               = keywordDict['OPT_LJ_CUTOFF']   
        self.OPT_COUL_CUTOFF             = keywordDict['OPT_COUL_CUTOFF']
        self.OPT_MOLTEMPLATE_FILE        = keywordDict['OPT_MOLTEMPLATE_FILE']
        self.OPT_DIELECT                 = keywordDict['OPT_DIELECT']
        self.OPT_MIXING                  = keywordDict['OPT_MIXING']    
        self.INITIAL_LJ_PARAMS_FILE      = keywordDict['INITIAL_LJ_PARAMS_FILE']
        self.OPT_ITERATIONS              = keywordDict['OPT_ITERATIONS']

        # coarse grained alphabet, bounds, names and which parameters
        # are to be optimized
        self.CG_GROUPS                   = keywordDict['CG_GROUPS']
        


    def parse_default(self, keyFileName):
        keywordDict = {}
        #keywordDict['SIMROOT'] = '/work/alex/SPA/SPA_100_jedd/SPA100'
        keywordDict['PYTHON_BIN'] = '/opt/pappuPython/bin/python'
        keywordDict['LAMMPS_BIN'] = '/home/kruff/LAMMPS/lammps-16Dec13/src/lmp_openmpi5'
        
        keywordDict['SIMROOT'] = '/work/kruff/N17permutants/09_2014_new/Nt17Qn/monomer/Q40'
        
        keywordDict['BOOTSTRAP_NUM_ITER'] = 20
        keywordDict['BOOTSTRAP_SIZE'] = 1000
        keywordDict['BOND_DEFINITION_FILE'] = 'BONDED_DEF.txt'
        keywordDict['BOND_PARAMETER_FILE'] = 'BONDED_PARAMS.txt'
        keywordDict['DIHEDRAL_DEFINITION_FILE'] = 'DIHEDRAL_DEF.txt'
        keywordDict['DIHEDRAL_PARAMETER_FILE'] = 'DIHEDRAL_PARAMS.txt'
        keywordDict['ANGLE_DEFINITION_FILE'] = 'ANGLE_DEF.txt'
        keywordDict['ANGLE_PARAMETER_FILE'] = 'ANGLE_PARAMS.txt'

        keywordDict['PLOT_BOND_HISTOGRAMS'] = True
        keywordDict['PLOT_ANG_HISTOGRAMS'] = True

        keywordDict['MBT_PARAMETER_FILE']  = 'MBT_PARAMS.txt'
        keywordDict['EBT_PARAMETER_FILE']  = 'EBT_PARAMS.txt'
        keywordDict['AAT_PARAMETER_FILE']  = 'AAT_PARAMS.txt'
        keywordDict['AT_PARAMETER_FILE']   = 'AT_PARAMS.txt'
        keywordDict['BB13_PARAMETER_FILE'] = 'BB13_PARAMS.txt'
        keywordDict['DAMPENING_PARAMETER_FILE'] = 'DAMP_PARAMS.txt'
        keywordDict['MASS_PARAMETER_FILE'] = 'MASS_PARAMS.txt'
        keywordDict['INITIAL_XYZ_FILE']    = 'INITIAL_POS.txt'

        keywordDict['TRAJECTORY_FILE'] = "N_005___traj.xtc"
        #keywordDict['TRAJECTORY_FILE'] = "__traj.xtc"
        keywordDict['PBD_FILE'] = "N_005___START.pdb"
        #keywordDict['PBD_FILE'] = "__START.pdb"

        keywordDict['INITIAL_LJ_PARAMS_FILE'] = 'INITIAL_LJ_PARAMS.txt'

        ## Moltemplate keywords below here - these will probably be things
        ## one should think more deelpy about and indeed detailed guidance
        ## on what these mean is going to be critical!
        keywordDict['OPT_KAPPA']             = float(0.1)
        keywordDict['OPT_LJ_CUTOFF']         = float(21.1)
        keywordDict['OPT_COUL_CUTOFF']       = float(15.0)
        keywordDict['OPT_MOLTEMPLATE_FILE']  = 'Moltemplate_config.lt'
        keywordDict['OPT_DIELECT']           = float(80)
        keywordDict['OPT_MIXING']            = 'arithmetic'
        keywordDict['OPT_ITERATIONS']        = 500
        

        # CG alphabet defined here
        # Basically in the actual keyfile we're going to have each grouping defined
        # as one group per line, where the line defines the residues, the group name, and upper and lower bounds
        # for sigma and epsilon - I'm thinking format like this
        #
        # CG_GROUP residues:KRED, name:chg, eps_lb:0.2, eps_ub:2.0, sig_lb:8.5, sig_ub:10.5
        #
        # This would be nice because if you wanted to use default values for sig you just don't
        # include a sigma value OR your specificy a value with;
        #
        # sig:<val>
        #
        # For eps you can specify a range or a default as 
        # eps:<val> 
        #
        # but you MUST provide one of the two (i.e. cannot work out what Epislon should be from
        # scratch!

        keywordDict['CG_GROUPS'] = []

        line="residues:KRED, name:chg, eps_lb:0.2, eps_ub:2.0, sig_lb:8.5, sig_ub:10.5"
        keywordDict['CG_GROUPS'].append(CGParameterGroup(line))

        line="residues:VFMILYQNWH, name:strong, eps_lb:0.5, eps_ub:1.5"
        keywordDict['CG_GROUPS'].append(CGParameterGroup(line))

        line="residues:APGSTC, name:weak, eps_lb:0.01, eps_ub:0.3"
        keywordDict['CG_GROUPS'].append(CGParameterGroup(line))

        return keywordDict
