## CAMeLOT V 0.1.0
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
##

from CAMeLOTExceptions import KeyFileException
from configs import MAX_NAME_LEN


class CGParameterGroup:
    """
    CGParameterGroup is an inteligently initialized object which
    takes a formatted paarmaeter line from a keyfile and converts 
    it into a set of rules for a set of residues.

    Specifically, the line being parsed defines how some group
    of residues should be defined in terms of their coarse grained
    parameters. This can be either by providing default values or
    by assigning a value to be optimized by the CAMeLOT optimization
    procedure.

    Each entry in a CGParameter line is of the format

    key:value, key:value, key:value

    IMPORTANTLYLY the characters ',' and ':' are used to define
    seperate entries in the line, and to separate out key/value
    pairing, so should NOT be used in the line (!!!)



    """

    def __init__(self, line):


        # initialize (we sanity check for semantics after all
        # items are accounted for)
        self.name      = None # parameter group name (i.e. what we're calling this
                              # parameter group - will be used as a prefix for LAMMPS
                              # parameter names. MUST BE UNIQUE 

        self.residues  = None # string of residues defined by this set - i.e. what are 
                              # the one letter amino acid codes associated with the
                              # residues to be included in this CGGroup

        self.eps       = None # epislon value
        self.sig       = None # sigma value
        self.eps_ub    = None # upper bound on epislon for optimization
        self.eps_lb    = None # lower bound on epsilon for optimization
        self.sig_ub    = None # upper bound on sigma for optimization
        self.sig_lb    = None # lower bound on sigma for optimization


        # remove trailing newline and split into seperate entries
        line=line.strip()
        linelist = line.split(",")

        ##
        ## Parse each item as separated by ','
        for item in linelist:
            
            try:
                key   = item.split(":")[0].lower().strip()
                value = item.split(":")[1].strip()
            except:
                raise KeyFileException("ERROR: Unable to the keyfile parameters associated with %s" % item)
            

            # residues string
            if key == "residues":
                self.residues = value.upper().strip()

            # CG group name
            elif key == "name":
                self.name = value.lower().strip()

            # default epsilon value
            elif key == "eps":
                try:
                    self.eps = float(value)
                except ValueError:
                    raise KeyFileExceptions("ERROR: Unable to interpret %s as an epislon value" % value)

            # default sigma value
            elif key == "sig":
                try:
                    self.sig = float(value)
                except ValueError:
                    raise KeyFileExceptions("ERROR: Unable to interpret %s as a sigma value" % value)

            # epsilon lower bound for optimization
            elif key == "eps_lb":
                try:
                    self.eps_lb = float(value)
                except ValueError:
                    raise KeyFileExceptions("ERROR: Unable to interpret %s as an epislon lower bound value" % value)

            # sigma lower bound for optimization
            elif key == "sig_lb":
                try:
                    self.sig_lb = float(value)
                except ValueError:
                    raise KeyFileExceptions("ERROR: Unable to interpret %s as a sigma lower bound value" % value)

            # epsilon upper bound for optimization
            elif key == "eps_ub":
                try:
                    self.eps_ub = float(value)
                except ValueError:
                    raise KeyFileExceptions("ERROR: Unable to interpret %s as an epislon upper bound value" % value)

            #  sigma upper bound for optimization
            elif key == "sig_ub":
                try:
                    self.sig_ub = float(value)
                except ValueError:
                    raise KeyFileExceptions("ERROR: Unable to interpret %s as a sigma upper bound value" % value)

            else:
                raise KeyFileExceptions("ERROR: Unexpected key in coarse-grained parameter definition [%s]" % key)

        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        #                                      Done parsing entries in the line
        # <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

        ########
        # Sanity checks


        # if we didn't define any residues or a name for the CGGroup
        if (self.residues is None) or (self.name is None):
            raise KeyFileExceptions("The coarse grained parameters defined by line [%s] failed to define a name and/or residues" % line)
            
        
        # 1) Convert the residue string to a set (in list form)
        self.residues = list(set(list(self.residues)))

        # 2) check the name is not too long
        if len(self.name) > MAX_NAME_LEN:
            raise KeyFileExceptions("ERROR: Name [%s] for coarse-grained parameter definition is too long, please ensure names are less than %i" % (self.name, MAX_NAME_LEN))

        # 3) if we defined an upper bound need a lower bound (cannot be bounded on one side only)
        if bool(self.eps_ub is None) != bool(self.eps_lb is None):
            # this is an XOR
            raise KeyFileExceptions("ERROR: Coarse-grained parameters %s must have both an upper and lower bound for epsilon"%self.name)

        if bool(self.sig_ub is None) != bool(self.sig_lb is None):
            # this is an XOR
            raise KeyFileExceptions("ERROR: Coarse-grained parameters %s must have both an upper and lower bound for sigma"%self.name)
            
        # 4) if we've defined an upper bound must have a lower bound and the lower bound must actually be lower!!!
        if self.sig_ub:
            if self.sig_ub < self.sig_lb:
                raise KeyFileExceptions("ERROR: Coarse-grained parameters %s must a sigma upper bound which is larger than the lower bound"%self.name)

        if self.eps_ub:
            if self.eps_ub < self.eps_lb:
                raise KeyFileExceptions("ERROR: Coarse-grained parameters %s must a epsilon upper bound which is larger than the lower bound"%self.name)

        # 5) epislon must have an upper/lower bound or a specific value (note we've already checked ub/lb consistency so we know if an ub is present
        #    there must be an lb
        if (self.eps_ub is None) and (self.eps is None):
            raise KeyFileExceptions("ERROR: Coarse-grained parameters %s requires a range for epsilon optimization OR a set value"%self.name)

        # 6) We don't actually have to check for sigma because if no ub/lb is provided for sigma then the code automatically uses the derived
        #    RG value for the sigma value (distance of closest approach)

        # finally set optimization flags (note we've checked for ub/lb consistency so only need to check ub here to set flags)
        if self.eps_ub:
            self.opt_eps=True
        else:
            self.opt_eps=False
        
        if self.sig_ub:
            self.opt_sig=True
        else:
            self.opt_sig=False



        


            
                


        
