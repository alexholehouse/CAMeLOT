## CAMeLOT V 0.1.0
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
##

import time
from CAMeLOTExceptions import OptimizationException
from MolTemplate import MolTemplate

class Optimization:
    
    def __init__(self, KeyFileObj, sequence_vector, charge_vector, rg_vector):
        """
        Initialization of the Optimization Object. This initializes the object by reeading
        in all the relevant information from the keyfile, and then parsing the derived values
        to build an initial set of sigma and epsilon values for each residue based on how
        we've set up our coarse graining.

        """
    
        # extract the amino acid sequence
        self.sequence_vector    = sequence_vector
        self.charge_vector      = charge_vector
        self.rg_vector          = rg_vector
        self.idxVector          = range(0,len(self.rg_vector))
        self.KeyFileObj         = KeyFileObj

        # check all the vector lengths are the same (1-bead-per-residue)
        if not (len(self.rg_vector) == len(self.sequence_vector) == len(self.charge_vector)):
            raise OptimizationException('Unable to match Rg, charge, and sequence vector')

        # parse the keyFile for all relevant parameters
        self.parseKeyfile(KeyFileObj)

        # extract and define the initial sigma and epsilon parameters which the optimization
        # procedure is going to use
        (self.initial_sigmas, self.initial_sigma_coef) = self.build_initial_sigmas()
        self.initial_epsilons = self.build_initial_epsilons()



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def parseKeyfile(self, keyFileObj):
        """
        Function to extract the relevant keywords from the keyFileObj - we extract these out
        to help simplify further use of keywords 

        """

        # get master file information
        self.GP_MATLAB_PATH = "PATH FOR MATLAB CODE"
        self.PYTHON_BIN     = keyFileObj.PYTHON_BIN
        self.LAMMPS_BIN     = keyFileObj.LAMMPS_BIN

        # build a mapping residue to group
        self.CG_parameters_list = keyFileObj.CG_GROUPS
            
        # the initial_LJ_PARAMS contains the filename to which initial LJ parameters
        # will be written to (by the write_initial_sigma_eps() function)
        self.INITIAL_LJ_PARAMS_FILE = keyFileObj.INITIAL_LJ_PARAMS_FILE

        # get the relevant keyfile parameters into a dictionary whichb
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
        Returns the CGParameterGroup object which the specified residue comes from.

        """
        for group in self.CG_parameters_list:
            if res in group.residues:
                return group



    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
    def build_initial_sigmas(self):
        """
        build_initial_sigmas returns two lists - list one defines the sigma values for each residuem and list
        two defines the sigma coefficienct values for each residue. How those sigma values are defined
        is described below.


        Sigma value can be defined in one of three ways

        1) Sigma values can be earmarked as optimization targets
        
        2) Sigma values can be defined with EXPLICIT values from the keyfile (i.e. we have *a prior* information
           about what we want the sigma value to be

        3) Sigma values can be derived from the atomistic simulation based on 
           the residue radius of gyration (recommended) (2*Rg)

        In all cases the coefficient is a totally dependent parameter which
        is equal to 2.5* the sigma value.

        
        """

        sigma_vector      = []
        sigma_coef_vector = []
        for (idx,res) in zip(self.idxVector,self.sequence_vector):

            # get the radius of gyration
            rg = self.rg_vector[idx]

            # get the associated CGGroup object which defines
            # this residue
            pg = self.get_residue_CGParameterGroup(res)

            # if we're optimizing this sigma value create a placeholder
            # using the cg group name 
            if pg.opt_sig:
                sigma_vector.append("sig_opt_%s"%pg.name)
                sigma_coef_vector.append("coff_sig_opt_%s"%pg.name)

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
        build_initial_epsilons returns a list containing a list of the epsilone values associated with each residue
        in the protein.

        Epsilon values can be defined in one of two ways

        1) Epsilon values can be earmarked as optimization targets

        2) Epsilon values can be explicitly set

        NOTE that unlike sigma values there is **NO WAY** to derive the epsilon
        value without it being optimized or explictly set from the keyfile
   
        """
        
        epsilon_vector = []

        # for each residue as defined by the sequence vector
        for res in self.sequence_vector:

            # get the associated CGGroup object which defines
            # this residue
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
        Function which returns the total number of (independent) parameters to be
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
        The Gaussiang process runner code (found in cg_opt.m) is the primary optimization code which is runs in MATLAB
        and acts as the framework around which the Gaussian process optimization occurs.

        This code itself calls another function which is custom generated by CAMeLOT (cg_opt_func.m) 


        """

        
        
        (num_params,param_list) = self.determine_number_of_parameters_to_optimize()
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
            fh.write('%% cg_opt - Optimization framework script which runs the BayesOpt package\n')
            fh.write('%% to perform Gaussian Process Bayesian optimization\n')
            fh.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
            
            fh.write('\n')
            fh.write('\n')

            fh.write('% Initialize Gaussian process\n')
            fh.write("addpath %s\n" % self.GP_MATLAB_PATH)
            fh.write("startup % this initializes the GP framework\n")
            fh.write('\n')
            fh.write("num_opt = %i\n" %  num_params)
            fh.write('% Set parameters to be optimized\n')


            # dynamically construct a matlab cell with the parameter types to be optimized (will match the names
            # used throughout the process)
            fh.write("type_opt = {'")
            for parameter_index in xrange(0,num_params):

                # if we're on the last parameter
                if parameter_index+1 == num_params:
                    fh.write("%s'" % param_list[parameter_index])
                else:
                    fh.write("%s', '" % param_list[parameter_index])
                    
            fh.write("}\n")

            # write the actual anonymous function call which is to be used. Note that we're essentially using a closure where
            # we're defining num_opt and type_opt HERE as we define the function, but when F is called it only
            # accepts 'x' (vector of values to be optimized) while num_opt and type_opt are already defined            
            fh.write('F = @(x)cg_opt_func(x, num_opt, type_opt);\n');

            # write the options
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



    def build_GaussianProcess_function_code(self):

        
        with open('cg_opt_func.m','w') as fh:
            
            fh.write('function value =  cg_opt_func(x, num_opt, type_opt)\n')
            fh.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
            fh.write('%%   \n')
            fh.write('%%   \n')
            fh.write('%% THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    \n')
            fh.write('%% Code generated at ' + str(time.strftime("%c")) + "\n")
            fh.write('%% Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'] + "\n")
            fh.write('%%   \n')
            fh.write('%%   \n')
            fh.write('%% cg_opt_func - Optimization function used by the BayesOpt to perform\n')
            fh.write('%% optimization\n')
            fh.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
            
            fh.write('\n')

            fh.write('  % Get current iteration\n')
            fh.write("  iteration = load('number_of_iterations.txt');\n")

            fh.write("\n")

            fh.write("  % Add to iteration number\n")
            fh.write("  iteration = iteration + 1\n")                    
            fh.write("  % Print iteration number to file\n")
            fh.write("  fileID=fopen('number_of_iterations.txt','w');\n")
            fh.write("  fprintf(fileID, '%4.0f', iteration);\n")
            fh.write("  fclose(fileID);\n")
            
            fh.write("\n")

            fh.write("  % cycle through each parameter estimate and write to file\n")
            fh.write("  for n = 1:num_opt\n")
            fh.write("     tempname=sprintf('%s_guess.txt', type_opt{n});\n")                     
            fh.write("     fileID=fopen(tempname, 'w');\n")
            fh.write("     fprintf(fileID, '%4.2f\m', roundn(x(n),-2));\n")
            fh.write("     fclose(fileID);\n")
            fh.write("     %if we're dealing with a sigma value also write the coefficient\n")
            fh.write("\n")
            fh.write("     if strcmp(tempname(1:3), 'sig') == 1\n")
            fh.write("        tempname=sprintf('coff_%s_guess.txt', type_opt{n});\n")
            fh.write("        fileID=fopen(tempname, 'w');\n")
            fh.write("        fprintf(fileID, '%4.2f\m', 2.5*roundn(x(n),-2));\n")
            fh.write("        fclose(fileID);\n")

            
            fh.write("\n")
            fh.write("\n")

            # this line defines the creation and updating of the OPT_MOLTEMPLATE_FILE to include the new parameters we just
            # wrote to file. update_moltemplate.py is a stand alone script which is 
            fh.write("  % this line will run a custom generated script which will update \n")
            fh.write("  % the moltemplate file with the current guesses \n")
            fh.write("  system('" + self.PYTHON_BIN + " update_moltemplate.py');\n")
            
            fh.write("\n")

            fh.write("  % this line will launch the simulation and then poll until the \n")
            fh.write("  % simulation has finished \n")
            fh.write("  system('" + self.PYTHON_BIN + " run_simulation.py');\n")
            
            fh.write("\n")

            fh.write("  % this line will run custom generated script which will compares\n")
            fh.write("  % the CG simulation against the all atom simulation \n")
            fh.write("  system('" + self.PYTHON_BIN + " evaluate_goodness.py');\n")
            fh.write("  value = load('current_goodness.txt');\n")

            

    def build_GaussianProcess_update_code(self):
        """
        Creates a file called update_moltemplate.py which is a customized script which updates the 
        moltemplate file we created previously with the current interation of guesses
        

        """


        (num_params,param_list) = self.determine_number_of_parameters_to_optimize()

        with open('update_moltemplate.py','w') as fh:

            fh.write('import os\n')
            fh.write('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n')
            fh.write('#   \n')
            fh.write('#   \n')
            fh.write('# THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    \n')
            fh.write('# Code generated at ' + str(time.strftime("%c")) + "\n")
            fh.write('# Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'] + "\n")
            fh.write('#   \n')
            fh.write('#   \n')
            fh.write('# update_moltemplate.py\n')
            fh.write('# \n\n')
            fh.write('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n')

            # for each parameter 
            fh.write("\n")
            fh.write("try:\n")
            fh.write("   with open('number_of_iterations.txt','r') as rfh:\n")
            fh.write("      iteration = int(rfh.read())\n")
            fh.write("except IOError:\n")
            fh.write("   iteration = 0\n")
            fh.write("\n")
                     


            # create a list of the parameters which are being optimized
            fh.write('parameters = [');
            for parameter_index in xrange(0,num_params):

                if parameter_index+1 == num_params:
                    fh.write("'%s'] "%param_list[parameter_index])
                else:
                    fh.write("'%s', "%param_list[parameter_index])
            fh.write('\n')
            fh.write('\n')
            fh.write('current_guess = {}\n')
            fh.write('\n')
            fh.write('# for each parameter read in the current best guess\n')
            fh.write('for i in parameters:\n')
            fh.write('   try:\n')
            fh.write("      with open(i+'_guess.txt','r') as rfh:\n")
            fh.write('         current_guess[i] = float(rfh.read())\n')
            fh.write('   except IOError:\n')
            fh.write("      current_guess[i] = 0.0\n")
            
            fh.write("\n")
            fh.write("   # sigma parmeters have a coefficient...\n")
            fh.write("   if i[0:3] == 'sig':\n")
            fh.write("      try:\n")
            fh.write("         with open('coff_' + i + '_guess.txt','r') as rfh:\n")
            fh.write("            current_guess['coff_'+i] = float(rfh.read())\n")
            fh.write('      except IOError:\n')
            fh.write("         current_guess['coff_'+i] = 0.0\n")

            fh.write("print current_guess\n")
            fh.write("\n")
            fh.write("# create next iteration directory\n")
            fh.write("print 'Creating framework for iteration %i' % iteration\n")
            fh.write("os.makedirs(str(iteration))\n")
            fh.write("with open('" + self.MolTemplateParameters['OPT_MOLTEMPLATE_FILE'] +"','r') as fh:\n")
            fh.write("  content = fh.readlines()\n")
            fh.write("\n")
            fh.write("newlines= []\n")
            fh.write("for line in content:\n")
            fh.write("   newline = line\n")
            fh.write("   for key in current_guess:\n")
            fh.write("      newline = newline.replace(key, str(current_guess[key]))\n")
            fh.write("   newlines.append(newline)\n")
            fh.write("\n")
            fh.write("with open('%i/sytem.lt'%iteration,'w') as fh:\n")
            fh.write("   for line in newlines:\n")
            fh.write("      fh.write(line)\n")
            fh.write("\n")
            fh.write("\n")
            
                     
            
                     
                     
                     
                     

                     
            
                
            
            
        

    def build_GaussianProcess_simulation_code(self):
        """
        Creates a file which runs the actual LAMMPS simulation

        """
        pass

    def build_GaussianProcess_target_code(self):\

        pass


            
                     
