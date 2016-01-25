## CAMeLOT V 0.1.1
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
##
##
#
# GaussianOptTools is a somewhat bizarre file. It represents a set of stand alone
# commandline utilites which are used by CAMeLOT's MATLAB Bayesian Optimization
# framework to do file and trajectory manipulation. MATLAB is really horrible
# for 'real' coding, so we write all the update/file manipulation/cost function
# evaluation code in Python taking advantage of it's excellent fileystem integration
# and our custom trajectory manipulation code built on top of MDTraj (CGTraj) to
# evaluate the cost function.
#
#


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def set_fail_flag(msg=None):
    """
    Function which sets the FAIL_STATUS.dat to 1. FAIL_STATUS.dat
    is read by the MATLAB code on every iteration and if it's found
    to be 1 the MATLAB code aborts. This provides a convenient and
    modular way to communicate failure between the Python portion
    of CAMeLOT and the MATLAB portion.

    The fact that FAIL_STATUS.dat is checked on every iteration
    also reflects the fact that we now have a continous evaluation
    of how things are going, rather than relying on some kind
    cross-language interupt like shenannigans which is likely
    to go balls up.

    We also write a message to the FAIL_MSG.dat file if required

    """
    with open('FAIL_STATUS.dat',w) as fh:
        fh.write("1\n")

    if msg:
        with open('FAIL_MSG.dat',w) as fh:
            fh.write("%s\n"%msg)

    exit(1)

        

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def update_parameter(name, value, parameter_file):
    """
    Function which takes a parameter name, value, and filename
    and updates the parameter in the file according to value.

    The parameter file should have the format
    
    ........................
    name \t value \n

    ........................
    
    Where name is a string and value is a float

    """

    # parse the parameter file
    parameter_dict = read_parameter_file(parameter_file)

    # if we cannot find the parameter in question in the parsed
    # parameter file
    if name not in parameter_dict:
        msg = 'ERROR: unable to find parameter [%s] in file [%s]\nThis will trigger a safe ABORT of CAMeLOT' % (name, parameter_file)
        raise GaussianOptToolsException(msg)
        set_fail_flag(msg)

    # update the parameter dictionary
    parameter_dict[name] = float(value)

    # write the update set of parameters back to file
    write_parameter_file(parameter_file, parameter_dict)
    
        

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def write_parameter_file(parameter_file, parameter_dict):
    """
    Function which generates a new parameter file using
    the correct format, containing the parameters
    defined in parameter_dict

    """
    with open(parameter_file,'w') as fh:
        for key in parameter_dict:
            fh.write('%s\t%4.3f\n'%(key, parameter_dict[key]))



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def read_parameter_file(parameter_file):
    """
    Function which reads in the parameter file and
    returns key-value dictionary where the key is
    the parameter name and the value is the associated
    value

    """
    parameter_dict = {}

    try:
        with open(parameter_file,'r') as fh:
            content=fh.readlines()
    except IOError:
        msg = 'ERROR: unable to find parameter file [%s]\nThis will trigger a safe ABORT of CAMeLOT' % parameter_file
        raise GaussianOptToolsException(msg)
        set_fail_flag(msg)
        
    for line in content:
        stripped=line.strip()
        line_split = stripped.split()
        parameter_dict[line_split[0]] = float(line_split[1])

    return parameter_dict


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def update_MolTemplate_file(moltemplate_parent, moltemplate_output, parameter_file):
    """
    Function which takes a moltemplate template (i.e. a moltemplate file where parameters to be optimized are missing)
    and 
    """
    
        
        
