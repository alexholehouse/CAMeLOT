## CAMeLOT V 0.1.2
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
## Code written by Alex Holehouse based on original principles developed by Kiersten Ruff
## Copyright 2015 - 2020
##
## Guide to adding new keywords
## *********************************
## 
## Adding new keywords to CAMELOT allows additional functionality to be added 
##
## 1. Add the keyword to the expected_keywords list
##
## 2. If the keyword cannot have a default value add to the required_keywords list. ELSE define a default
##    value in 
##
## 3. Add to the appropriate place in the parse_keyfile method (i.e. if its a string add to the string)

import numpy as np

from CAMeLOTExceptions import KeyFileException
from CGParameterGroup import CGParameterGroup
from configs import CAMELOT_VERSION

def is_comment_line(line):
    """
    A comment line is any line where the first readable character is a '#'
    
    """
    if len(line[0:line.find('#')].strip()) == 0:
        return True
    else:
        return False


def remove_comments(line):
    """
    Removes all the the content after a comment character
    and returns the whitespace/newline stripped string
    left of the comment charater ('#')

    """
    return line.split('#')[0].strip()


class KeyfileParser:

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    def __init__(self, keyfile_name):

        # expected keywords contains a list of possible keywords which can be read form the keyfile. Importantly if these
        # keywords are *not* included then they are set to their default values
        self.expected_keywords = ['LOGFILE', 'TEMP', 'SIMROOT', 'PYTHON_BIN', 'LAMMPS_BIN', 'MATLAB_GP_CODE_DIR', 'BOOTSTRAP_NUM_ITER', 'BOOTSTRAP_SIZE', 
                                  'PLOT_BOND_HISTOGRAMS', 'PLOT_ANG_HISTOGRAMS','PLOT_DIHEDRAL_HISTOGRAMS', 'INITIAL_XYZ_FILE', 'MASS_PARAMETER_FILE',
                                  'BOND_DEFINITION_FILE', 'BOND_PARAMETER_FILE', 'DIHEDRAL_DEFINITION_FILE', 'DIHEDRAL_PARAMETER_FILE', 'ANGLE_DEFINITION_FILE', 
                                  'ANGLE_PARAMETER_FILE', 'MBT_PARAMETER_FILE', 'EBT_PARAMETER_FILE', 'AAT_PARAMETER_FILE', 
                                  'AT_PARAMETER_FILE', 'BB13_PARAMETER_FILE', 'DAMPENING_PARAMETER_FILE', 'RES_RES_DISTANCE_FILE', 'RES_RES_BIN_START', 'RES_RES_BIN_END', 
                                  'RES_RES_BIN_SIZE', 'OPT_KAPPA', 'OPT_LJ_CUTOFF', 'OPT_COUL_CUTOFF', 'OPT_MOLTEMPLATE_FILE', 'OPT_DIELECT', 'OPT_MIXING', 
                                  'INITIAL_LJ_PARAMS_FILE', 'OPT_ITERATIONS', 'SIM_NSTEPS', 'SIM_DCD_OUT', 'SIM_THERMO_OUT', 'SIM_EQUIL_FRACTION' , 'CG_GROUPS',
                                'DISABLE_SANITY_CHECKS', 'SIMULATION_ENGINE', 'PDB_FILE','TRAJECTORY_FILE']
        
        # required keywords are those which MUST be included in the keyfile - i.e. CAMELOT can't set default values for these
        # keywords. Note that required_keywords is a subset of expected_keywords
        self.required_keywords = ['TEMP', 'SIMROOT', 'PYTHON_BIN', 'LAMMPS_BIN', 'MATLAB_GP_CODE_DIR', 'BOOTSTRAP_NUM_ITER', 'BOOTSTRAP_SIZE',                                 
                                  'RES_RES_BIN_START', 'RES_RES_BIN_END', 'RES_RES_BIN_SIZE', 'OPT_KAPPA', 'OPT_COUL_CUTOFF', 'OPT_DIELECT', 
                                  'OPT_ITERATIONS', 'SIM_NSTEPS', 'SIM_DCD_OUT', 'SIM_THERMO_OUT', 'SIM_EQUIL_FRACTION', 'CG_GROUPS', 'SIMULATION_ENGINE',
                                  'PDB_FILE','TRAJECTORY_FILE']
    
        # zero out a couple of key dictionaries
        self.keyword_lookup = {}
        self.DEFAULTS = {}
        
        # prepare default values
        default_arguments = self.parse_default()

        # set keyfile keywords
        self.parse(keyfile_name)

        # checks that all our ducks are in a row re keywords
        self.sanity_check_keywords_and_set_defaults(default_arguments)

        
        # Assign keywords to object variables        
        self.CAMELOT_VERSION             = CAMELOT_VERSION
        self.LOGFILE                     = self.keyword_lookup['LOGFILE']
        self.SIMULATION_ENGINE           = self.keyword_lookup['SIMULATION_ENGINE']
        
        self.TEMP                        = self.keyword_lookup['TEMP']
        
        self.SIMROOT                     = self.keyword_lookup['SIMROOT']
        self.PYTHON_BIN                  = self.keyword_lookup['PYTHON_BIN']
        self.LAMMPS_BIN                  = self.keyword_lookup['LAMMPS_BIN']
        self.MATLAB_GP                   = self.keyword_lookup['MATLAB_GP_CODE_DIR']

        self.BOOTSTRAP_NUM_ITER          = self.keyword_lookup['BOOTSTRAP_NUM_ITER']
        self.BOOTSTRAP_SIZE              = self.keyword_lookup['BOOTSTRAP_SIZE']

        self.PLOT_BOND_HISTOGRAMS        = self.keyword_lookup['PLOT_BOND_HISTOGRAMS']
        self.PLOT_ANG_HISTOGRAMS         = self.keyword_lookup['PLOT_ANG_HISTOGRAMS']
        self.PLOT_DIHEDRAL_HISTOGRAMS    = self.keyword_lookup['PLOT_DIHEDRAL_HISTOGRAMS']
        
        self.INITIAL_XYZ_FILE            = self.keyword_lookup['INITIAL_XYZ_FILE']
        self.MASS_PARAMETER_FILE         = self.keyword_lookup['MASS_PARAMETER_FILE']

        self.BOND_DEFINITION_FILE        = self.keyword_lookup['BOND_DEFINITION_FILE']
        self.BOND_PARAMETER_FILE         = self.keyword_lookup['BOND_PARAMETER_FILE']
        self.DIHEDRAL_DEFINITION_FILE    = self.keyword_lookup['DIHEDRAL_DEFINITION_FILE']
        self.DIHEDRAL_PARAMETER_FILE     = self.keyword_lookup['DIHEDRAL_PARAMETER_FILE']
        self.ANGLE_DEFINITION_FILE       = self.keyword_lookup['ANGLE_DEFINITION_FILE']
        self.ANGLE_PARAMETER_FILE        = self.keyword_lookup['ANGLE_PARAMETER_FILE']
        self.TRAJECTORY_FILE             = self.keyword_lookup['TRAJECTORY_FILE']
        self.PDB_FILE                    = self.keyword_lookup['PDB_FILE']
        self.MBT_PARAMETER_FILE          = self.keyword_lookup['MBT_PARAMETER_FILE']
        self.EBT_PARAMETER_FILE          = self.keyword_lookup['EBT_PARAMETER_FILE']
        self.AAT_PARAMETER_FILE          = self.keyword_lookup['AAT_PARAMETER_FILE']
        self.AT_PARAMETER_FILE           = self.keyword_lookup['AT_PARAMETER_FILE']
        self.BB13_PARAMETER_FILE         = self.keyword_lookup['BB13_PARAMETER_FILE']
        self.DAMPENING_PARAMETER_FILE    = self.keyword_lookup['DAMPENING_PARAMETER_FILE']     

        # inter-residue stuff
        self.RES_RES_DISTANCE_FILE       = self.keyword_lookup['RES_RES_DISTANCE_FILE']

        self.RES_RES_BIN_START           = self.keyword_lookup['RES_RES_BIN_START']
        self.RES_RES_BIN_END             = self.keyword_lookup['RES_RES_BIN_END']
        self.RES_RES_BIN_SIZE            = self.keyword_lookup['RES_RES_BIN_SIZE']
        self.RES_RES_BINS                = np.arange(self.RES_RES_BIN_START, self.RES_RES_BIN_END, self.RES_RES_BIN_SIZE)

        # optimization keywords
        self.OPT_KAPPA                   = self.keyword_lookup['OPT_KAPPA']     
        self.OPT_LJ_CUTOFF               = self.keyword_lookup['OPT_LJ_CUTOFF']   
        self.OPT_COUL_CUTOFF             = self.keyword_lookup['OPT_COUL_CUTOFF']
        self.OPT_MOLTEMPLATE_FILE        = self.keyword_lookup['OPT_MOLTEMPLATE_FILE']
        self.OPT_DIELECT                 = self.keyword_lookup['OPT_DIELECT']
        self.OPT_MIXING                  = self.keyword_lookup['OPT_MIXING']    
        self.INITIAL_LJ_PARAMS_FILE      = self.keyword_lookup['INITIAL_LJ_PARAMS_FILE']
        self.OPT_ITERATIONS              = self.keyword_lookup['OPT_ITERATIONS']

        # simulation keyword        
        self.SIM_NSTEPS                  = self.keyword_lookup['SIM_NSTEPS'] 
        self.SIM_DCD_OUT                 = self.keyword_lookup['SIM_DCD_OUT'] 
        self.SIM_THERMO_OUT              = self.keyword_lookup['SIM_THERMO_OUT']   
        self.SIM_EQUIL_FRACTION          = self.keyword_lookup['SIM_EQUIL_FRACTION'] 
        
        # coarse grained alphabet, bounds, names and which parameters
        # are to be optimized
        self.CG_GROUPS                   = self.keyword_lookup['CG_GROUPS']

        # Having intialized all the object variables we run any self consistency
        # checks which may be important. One can imagine adding a large set of 
        # sanitation checks here such that we know all the input parameters are
        # valid before we start anything...
        
        # perform any/all sanity checks needed (i.e. checking files exits,
        # directories are writable, basically everything we might need to ensure the
        # CG procedure will work
        self.run_sanity_checks()

        
        

    def run_sanity_checks(self):
        
        # check CG groups
        self.check_CGGroups_are_unique()

        # check specific keyword variables
        if self.SIMULATION_ENGINE not in ['GROMACS', 'CAMPARI']:
            raise KeyFileException('Unexpected atomistic simulation engine used - expected GROMACS or CAMPARI but got %s'%self.SIMULATION_ENGINE)
        
    

    def sanity_check_keywords_and_set_defaults(self, default_arguments):
        """
        Simple check to make sure that we're not missing any key keywords. Having
        performed the relevant check we then assign default keywords to the remaining
        undefined keywords in this specific run.

        Does not return anything but throws an exception if anything seems off.
        

        """
        # check we have all the required keywords
        for keyword in self.required_keywords:
            if keyword not in self.keyword_lookup:
                raise KeyFileException('The required keyword [%s] was not found in the keyfile and must be included. Please correct your keyfile and retry' % keyword)

        
        # check all the expected keywords are defined in either the default or parsed keyword list
        for keyword in self.expected_keywords:
            if (keyword not in self.keyword_lookup) and (keyword not in default_arguments):
                raise KeyFileException('The expected keyword [%s] was not found in either keyfile and there was no pre-defined default value. Please correct your keyfile and retry' % keyword)
                
        # finally safe in the knoweldge that everything is accounted for set the default values
        for keyword in self.expected_keywords:
            if (keyword not in self.keyword_lookup):
                self.keyword_lookup[keyword] = default_arguments[keyword]

        # final check to make sure we match the expected and actual number
        if not len(self.keyword_lookup) == len(self.expected_keywords):
            raise KeyFileException('The total read keywords does not match the expected number (%i expected, but %i found)' %(len(self.expected_keywords),len(self.keyword_lookup)))
            




    #-----------------------------------------------------------------
    #    
    def parse(self, filename):
        """
        Main function reads in a keyfile and extracts out the relevant details based on the keywords. Keywoerds are assigned to 
        the self.keywords_lookup

        Keywords must be defined as 

        KEYWORD : VALUE # comment 

        All read variables are correctly cast and assigned to the self.keyword_lookup variable

        """

        self.keyword_lookup['CG_GROUPS'] = []
        
        
        # read the keyfile lines
        with open(filename, 'r') as fh:
            contents = fh.readlines()

        # for each line
        for line in contents:

            # if it's a comment line skip the whole line
            if is_comment_line(line):
                continue
                
            # remove comment section (at end of line)
            un_comment = remove_comments(line)

            # split based on the keyword/value separator
            splitline = un_comment.split(':')

            if len(splitline) > 2:
                raise KeyFileException('On keyword %s - found multiple keyword-value separators...' % splitline[0].strip())

            if len(splitline) < 2:
                raise KeyFileException('On line [%s] - no keyword separator found...' % line)

            # if get here must have keyword/keyvalue                
            putative_keyword = splitline[0].strip().upper()
            putative_value   = splitline[1].strip()
            
            if putative_keyword in self.expected_keywords:
                
                ## ** CHECK TO ENSURE WE DON'T OVERWRITE KEYWORDS **
                # check if we've seen this keyword before - if we're trying to overwrite raise an exception
                if putative_keyword in self.keyword_lookup.keys():

                    if putative_keyword == "CG_GROUPS":
                        # this is OK - we can have multiple CG_groups
                        pass

                    else:
                        raise KeyFileException('Found a second occurence of the [%s] keyword. Please correct your keyfile and retry' % putative_keyword)
                        
                # special case for CG parameter definitions
                if putative_keyword == "CG_GROUPS":
                    self.keyword_lookup['CG_GROUPS'].append(CGParameterGroup(putative_value))
                
                # integer keywords
                elif putative_keyword in ['BOOTSTRAP_NUM_ITER', 'BOOTSTRAP_SIZE', 'SIM_NSTEPS', 'SIM_DCD_OUT', 'SIM_THERMO_OUT','OPT_ITERATIONS']:
                    self.keyword_lookup[putative_keyword] = int(putative_value)

                # float keywords
                elif putative_keyword in ['TEMP', 'RES_RES_BIN_START', 'RES_RES_BIN_END', 'RES_RES_BIN_SIZE', 'OPT_KAPPA', 'OPT_LJ_CUTOFF', 'OPT_COUL_CUTOFF', 
                                          'OPT_DIELECT', 'SIM_EQUIL_FRACTION']:            
                    self.keyword_lookup[putative_keyword] = float(putative_value)

                # string keywords
                elif putative_keyword in ['SIMROOT','PYTHON_BIN','LAMMPS_BIN', 'LOGFILE','MATLAB_GP_CODE_DIR', 'BOND_DEFINITION_FILE',
                                          'BOND_PARAMETER_FILE', 'DIHEDRAL_DEFINITION_FILE', 'DIHEDRAL_PARAMETER_FILE', 'ANGLE_DEFINITION_FILE',
                                          'ANGLE_PARAMETER_FILE','MBT_PARAMETER_FILE', 'EBT_PARAMETER_FILE', 'AAT_PARAMETER_FILE', 'AT_PARAMETER_FILE',
                                          'BB13_PARAMETER_FILE', 'DAMPENING_PARAMETER_FILE', 'MASS_PARAMETER_FILE', 'INITIAL_XYZ_FILE','RES_RES_DISTANCE_FILE',
                                          'TRAJECTORY_FILE', 'PDB_FILE', 'INITIAL_LJ_PARAMS_FILE','OPT_MOLTEMPLATE_FILE', 'OPT_MIXING', 'SIMULATION_ENGINE']:

                    # for case insensitive keywords set value to upper
                    if putative_keyword in ['SIMULATION_ENGINE']:
                        putative_value = putative_value.upper()
                        
                    self.keyword_lookup[putative_keyword] = str(putative_value)

                # boolean keywords
                elif putative_keyword in ['PLOT_BOND_HISTOGRAMS', 'PLOT_ANG_HISTOGRAMS', 'PLOT_DIHEDRAL_HISTOGRAMS', 'DISABLE_SANITY_CHECKS']:            

                    if putative_value.upper() == "TRUE":
                        self.keyword_lookup[putative_keyword] = True
                    elif putative_value.upper() == "FALSE":
                        self.keyword_lookup[putative_keyword] = False
                    else:
                        raise KeyFileException('Expected True or False for [%s] keyword but got [%s]. Please correct your keyfile and retry' % (putative_keyword, putative_value))
                                                               
            else:
                raise KeyFileException('Found an unsupported keyword - [%s] ' % putative_keyword)




    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    def check_CGGroups_are_unique(self):
        """
        Function which ensures the names associated with all the CGGroups are
        unique. This is important because we use these names to define placerholders
        for parameters during the optimization process, so their uniqueness is
        kind of critical.

        """
        
        names=[]
        
        for group in self.CG_GROUPS:
            if group.name in names:
                raise KeyFileException('One of the CGGroup names defined in the key file [%s] was not unique. Please ensure all names are unique or you will FSU' % group.name)
            else:
                names.append(group.name)


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                   
    def assign_missing_values_to_defaults(defaults):
        """
        Function which will scan through the self.keyword_lookup dictionary and determine if non-required keywords have not yet
        been set - those which have not been set 


        """
        pass

    

    def parse_default(self):
        """
        Builds a dictionary where non-required keywords are set to their default value

        """

        # initialize
        keywordDict = {}

        # LOGFILE
        keywordDict['LOGFILE']                  = 'camelot.log'        

        # Molecular Mechanics filenames
        keywordDict['BOND_DEFINITION_FILE']     = 'BONDED_DEF.txt'
        keywordDict['BOND_PARAMETER_FILE']      = 'BONDED_PARAMS.txt'
        keywordDict['DIHEDRAL_DEFINITION_FILE'] = 'DIHEDRAL_DEF.txt'
        keywordDict['DIHEDRAL_PARAMETER_FILE']  = 'DIHEDRAL_PARAMS.txt'
        keywordDict['ANGLE_DEFINITION_FILE']    = 'ANGLE_DEF.txt'
        keywordDict['ANGLE_PARAMETER_FILE']     = 'ANGLE_PARAMS.txt'
        keywordDict['MBT_PARAMETER_FILE']       = 'MBT_PARAMS.txt'
        keywordDict['EBT_PARAMETER_FILE']       = 'EBT_PARAMS.txt'
        keywordDict['AAT_PARAMETER_FILE']       = 'AAT_PARAMS.txt'
        keywordDict['AT_PARAMETER_FILE']        = 'AT_PARAMS.txt'
        keywordDict['BB13_PARAMETER_FILE']      = 'BB13_PARAMS.txt'
        
        # mass/LD parameter files
        keywordDict['DAMPENING_PARAMETER_FILE'] = 'DAMP_PARAMS.txt'
        keywordDict['MASS_PARAMETER_FILE']      = 'MASS_PARAMS.txt'
        keywordDict['INITIAL_XYZ_FILE']         = 'INITIAL_POS.txt'
        keywordDict['RES_RES_DISTANCE_FILE']    = 'RES_RES_DIS.txt'
        
        # Starting LJ parameters
        keywordDict['INITIAL_LJ_PARAMS_FILE']   = 'INITIAL_LJ_PARAMS.txt'

        # moltempate filename
        keywordDict['OPT_MOLTEMPLATE_FILE']     = 'Moltemplate_config.lt'
        keywordDict['OPT_LJ_CUTOFF']            = float(25.0) # note this is a big LJ cutoff but actually doesn't seem to play any role whatsoever...
       
        # should we plot histograms (default no!) 
        keywordDict['PLOT_BOND_HISTOGRAMS']     = False
        keywordDict['PLOT_ANG_HISTOGRAMS']      = False
        keywordDict['PLOT_DIHEDRAL_HISTOGRAMS'] = False

        # mixing rule 
        keywordDict['OPT_MIXING']              = 'arithmetic'

        # disable sanity checks
        keywordDict['DISABLE_SANITY_CHECKS']   = False

                                   #

        return keywordDict
        


        
                
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    
