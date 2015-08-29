## CAMeLOT V 0.1.0
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
##

import time
from CAMeLOTExceptions import OptimizationException
from MolTemplate import MolTemplate

class Optimization:
    
    def __init__(self, KeyFileObj, sequence_vector, charge_vector, rg_vector):
    
        # extract the amino acid sequence
        self.sequence_vector    = sequence_vector
        self.charge_vector      = charge_vector
        self.rg_vector          = rg_vector
        self.idxVector          = range(0,len(self.rg_vector))
        self.KeyFileObj         = KeyFileObj

        if not (len(self.rg_vector) == len(self.sequence_vector) == len(self.charge_vector)):
            raise OptimizationException('Unable to match Rg, charge, and sequence vector')

        # parse the keyFile for all relevant parameters
        self.parseKeyfile(KeyFileObj)


        # convert that amino acid sequence into a CG sequence
        # as defined by the keyfile

        #self.charge_seq  = self.build_charge_seq()

        (self.initial_sigmas, self.initial_sigma_coef) = self.build_initial_sigmas()
        self.initial_epsilons = self.build_initial_epsilons()



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def parseKeyfile(self, keyFileObj):
        """
        Function to extract the relevant keywords from the keyFileObj - we extract these out
        to help simplify further use of keywords 

        """

        self.GP_MATLAB_PATH = "PATH FOR MATLAB CODE"

        # build a mapping residue to group
        self.CG_parameters_list = keyFileObj.CG_GROUPS
            
        self.INITIAL_LJ_PARAMS_FILE = keyFileObj.INITIAL_LJ_PARAMS_FILE

        # get the relevant keyfile parameters into a dictionary which
        # we pass to the MolTemplate object
        relevant_keywords = ['OPT_KAPPA', 'OPT_LJ_CUTOFF', 'OPT_LJ_CUTOFF', 'OPT_COUL_CUTOFF', 'OPT_MOLTEMPLATE_FILE', 'OPT_DIELECT', 'OPT_MIXING']

        self.MolTemplateParameters = {}
        for keyword in relevant_keywords:
            self.MolTemplateParameters[keyword] = getattr(keyFileObj,keyword)            
        self.MolTemplateParameters['CAMELOT_VERSION'] = keyFileObj.CAMELOT_VERSION



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def get_residue_CGParameterGroup(self, res):
        """
        Returns the CGParameterGroup object which the specified residue comes from

        """
        for group in self.CG_parameters_list:
            if res in group.residues:
                return group



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    def build_initial_sigmas(self):
        """
        Sigma value can be defined in one of three ways

        1) Sigma values can be earmarked as optimization targets
        
        2) Sigma values can be defined with EXPLICIT values from the keyfile

        3) Sigma values can be derived from the atomistic simulation based on 
           the residue radius of gyration (recommended) (2*Rg)

        In all cases the coefficient is a totally dependent parameter which
        is equal to 2.5* the sigma value
        """

        sigma_vector      = []
        sigma_coef_vector = []
        for (idx,res) in zip(self.idxVector,self.sequence_vector):

            # extact radius of gyration
            rg = self.rg_vector[idx]

            pg = self.get_residue_CGParameterGroup(res)

            # if we're optimizing this sigma value create a placeholder
            # using the cg group name 
            if pg.opt_sig:
                sigma_vector.append("sig_opt_%s"%pg.name)
                sigma_coef_vector.append("coff_opt_%s"%pg.name)

            # if a specific sigma value was defined
            elif pg.sig:
                sigma_vector.append('%4.2f'%(pg.sig))
                sigma_coef_vector.append('%4.2f'%(pg.sig*2.5))
                
            # if we're using the radius of gyration
            else:
                sigma_vector.append('%4.2f'%(rg*2))
                sigma_coef_vector.append('%4.2f'%(rg*5))

        return (sigma_vector, sigma_coef_vector)



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_initial_epsilons(self):
        """
        Epsilon values can be defined in one of two ways

        1) Epsilon values can be earmarked as optimization targets

        2) Epsilon values can be explicitly set

        NOTE that unlike sigma values there is **NO WAY** to derive the epsilon
        value without it being optimized or explictly set from the keyfile
   
        """
        
        epsilon_vector = []
        for res in self.sequence_vector:

            pg = self.get_residue_CGParameterGroup(res)

            # if residue has optimized eps set
            if pg.opt_eps:
                epsilon_vector.append("eps_opt_%s"%pg.name)
            else:
                epsilon_vector.append('%4.2f'%(pg.eps))

        return epsilon_vector



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def write_initial_sigma_eps(self):
        """ 
        Function which writes out the initial LJ parameters to file. Note the 
        INITIAL_LJ_PARAMS_FILE basically defines our degrees of freedom for optimization
        """

        with open(self.INITIAL_LJ_PARAMS_FILE,'w') as fh:
            for idx in self.idxVector:
                res=idx+1
                fh.write('pair_coeff  @atom:Res%i @atom:Res%i lj/cut/coul/debye %s %s %s\n' %(res, res, self.initial_epsilons[idx], self.initial_sigmas[idx], self.initial_sigma_coef[idx]))
            
                

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def determine_number_of_parameters_to_optimize(self, display=False):
        """
        Function whic returns the total number of (independent) parameters to be
        optimized - essentially this defines the dimensionality of the space the
        Gaussian process optimization has to navigate.

        """                
        params_to_opt = []

        for group in self.CG_parameters_list:
            if group.opt_eps:
                params_to_opt.append("eps_opt_%s" % group.name)
            if group.opt_sig:
                params_to_opt.append("sig_opt_%s" % group.name)

        params_to_opt = list(set(params_to_opt))

        return (len(params_to_opt), params_to_opt)



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def determine_parameter_min_and_max(self):
        """
        Function which returns a dictoinary, where the keys give a parameter name and the values are 
        a two position tuple with the min and max value

        """                

        params_min_max = {}

        for group in self.CG_parameters_list:
            if group.opt_eps:
                params_min_max["eps_opt_%s" % group.name] = (group.eps_lb, group.eps_ub)
            if group.opt_sig:
                params_min_max["sig_opt_%s" % group.name] = (group.sig_lb, group.sig_ub)

        return params_min_max
       


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_moltemplate_input(self):
        ## Function which constructs
        ## an initial file with all parameters filled in
        ## and placeholders for the parameters were trying to estimate
                
        # initialize the object
        MolTemplateObject = MolTemplate(self.MolTemplateParameters)

        # build the moltemplate file
        MolTemplateObject.build_moltemplate_input(self.KeyFileObj.INITIAL_XYZ_FILE,
                                                  self.KeyFileObj.MASS_PARAMETER_FILE ,
                                                  self.KeyFileObj.INITIAL_LJ_PARAMS_FILE,
                                                  self.KeyFileObj.BOND_DEFINITION_FILE,
                                                  self.KeyFileObj.BOND_PARAMETER_FILE,
                                                  self.KeyFileObj.ANGLE_DEFINITION_FILE,
                                                  self.KeyFileObj.ANGLE_PARAMETER_FILE,
                                                  self.KeyFileObj.DIHEDRAL_DEFINITION_FILE,
                                                  self.KeyFileObj.DIHEDRAL_PARAMETER_FILE,
                                                  self.KeyFileObj.MBT_PARAMETER_FILE ,
                                                  self.KeyFileObj.EBT_PARAMETER_FILE,
                                                  self.KeyFileObj.AT_PARAMETER_FILE ,
                                                  self.KeyFileObj.AAT_PARAMETER_FILE ,
                                                  self.KeyFileObj.BB13_PARAMETER_FILE)
        
        
        
        
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def build_GaussianProcess_runner_code(self):
        """
        The Gaussiang process runner code (found in cg_opt.m) is the primary optimization code which is called


        """

        
        (n_params,param_list) = self.determine_number_of_parameters_to_optimize()
        parameter_min_max     = self.determine_parameter_min_and_max()

        # this whole function generates a bunch of runner code - specifically it creates
        # 1) A matlab parent script
        # 2) A matlab function called by the parent script to execute the Gaussian process
        
        with open('cg_opt.m','w') as fh:

            fh.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
            fh.write('%%   \n')
            fh.write('%%   \n')
            fh.write('%% THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    \n')
            fh.write('%% Code generated at ' + str(time.strftime("%c")) + "\n")
            fh.write('%% Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'] + "\n")
            fh.write('%%   \n')
            fh.write('%%   \n')
            fh.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
            
            fh.write('\n')
            fh.write('\n')

            fh.write('% Initialize Gaussian process\n')
            fh.write("addpath %s\n" % self.GP_MATLAB_PATH)
            fh.write("startup\n")
            fh.write('\n')
            fh.write("num_opt = %i\n" %  n_params)
            fh.write("type_opt = {'")

            # dynamically construct a matlab cell with the parameter types to be optimized (will match the names
            # used throughout the process)
            num_params = len(param_list)

            fh.write('% Set parameters to be optimized\n')
            for parameter_index in xrange(0,num_params):

                # if we're on the last parameter
                if parameter_index+1 == num_params:
                    fh.write("%s'" % param_list[parameter_index])
                else:
                    fh.write("%s', '" % param_list[parameter_index])
                    
            fh.write("}\n")

            # write the actual anonymous function call which is to be used
            fh.write('F = @(x)cg_opt_func(x, num_opt, type_opt);\n');

            # write the option s
            fh.write('\n')
            fh.write('opt = bo_default_opt();\n')
            fh.write('opt.dims = num_opt; \n')
            fh.write('\n')

            # dynamically construct matlab code defining the minimum and maximum values for each
            # of the parameter
            fh.write('% Set min and max for optimization procedure\n')
            for parameter_index in xrange(0,num_params):
                
                parameter_name = param_list[parameter_index]
                fh.write('opt.mins(%i)  = %3.3f;\n'  % (parameter_index+1,parameter_min_max[parameter_name][0]))
                fh.write('opt.maxes(%i) = %3.3f;\n' % (parameter_index+1,parameter_min_max[parameter_name][1]))
                
            fh.write('\n')

            # set the number of iterations
            fh.write('% Set number of optimization iterations [OPT_ITERATIONS]\n')
            fh.write('opt.max_iters = %i;\n'%(self.KeyFileObj.OPT_ITERATIONS))
            fh.write('opt.save_trace = true;\n')

            fh.write('\n')
            
            # write the actual optimization call out!
            fh.write('% run optimization\n')
            fh.write("[ms,mv,T] = bayesopt(F,opt);\n")
            fh.write('\n')
            fh.write("fprintf('\\n**Optimization finished.\\n');\n")
            fh.write("fprintf('Minimum function value found in grid was \%f\\n',mv);\n")
            fh.write("fprintf('Point with minimum value was...\\n');\n")
            fh.write("disp(ms);\n")

  
