## CAMeLOT V 0.1.1
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
##

import time
from CAMeLOTExceptions import OptimizationException
from MolTemplate import MolTemplate
from utils import writeline

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
        self.GP_MATLAB_PATH = keyFileObj.MATLAB_GP
        self.PYTHON_BIN     = keyFileObj.PYTHON_BIN
        self.LAMMPS_BIN     = keyFileObj.LAMMPS_BIN
        self.TEMP           = keyFileObj.TEMP
        self.DAMP_PARAMS    = keyFileObj.DAMPENING_PARAMETER_FILE

        # get some simulation parameters
        self.SIM_DCD_OUT        = keyFileObj.SIM_DCD_OUT
        self.SIM_NSTEPS         = keyFileObj.SIM_NSTEPS
        self.SIM_THERMO_OUT     = keyFileObj.SIM_THERMO_OUT
        self.SIM_EQUIL_FRACTION = keyFileObj.SIM_EQUIL_FRACTION

        # some other relevant information
        self.RES_RES_BIN_START     = keyFileObj.RES_RES_BIN_START
        self.RES_RES_BIN_END     = keyFileObj.RES_RES_BIN_END
        self.RES_RES_BIN_SIZE      = keyFileObj.RES_RES_BIN_SIZE 
        self.RES_RES_DISTANCE_FILE = keyFileObj.RES_RES_DISTANCE_FILE

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
                # the approach below was true when using a hybrid potential, but apparently for a non-hybird potential you
                # don't need to define the potential type on *each* line. 
                #fh.write('pair_coeff  @atom:Res%i @atom:Res%i lj/cut/coul/debye %s %s %s\n' %(res, res, self.initial_epsilons[idx], self.initial_sigmas[idx], self.initial_sigma_coef[idx]))
                fh.write('pair_coeff  @atom:Res%i @atom:Res%i %s %s %s\n' %(res, res, self.initial_epsilons[idx], self.initial_sigmas[idx], self.initial_sigma_coef[idx]))
            
                

    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def determine_number_of_parameters_to_optimize(self, display=False):
        """
        Function which returns the total number of (independent) parameters to be
        optimized - essentially this defines the dimensionality of the space the
        Gaussian process optimization has to navigate.

        In all cases, we're only optimizing sigma or epsilon parameters, and in all cases

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
        
        
        
        
    
    def build_code(self):
        self.build_GaussianProcess_runner_code()
        self.build_GaussianProcess_function_code()
        self.build_GaussianProcess_update_code()
        self.build_GaussianProcess_simulation_code()
        self.build_GaussianProcess_target_code()
        self.build_system_setup_file()
        self.build_nvt_minimization_file()
        self.build_nvt_simulation_code()
        
        
        
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
            fh.write('%   \n')
            fh.write('%   \n')
            fh.write('% THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    \n')
            fh.write('% Code generated at ' + str(time.strftime("%c")) + "\n")
            fh.write('% Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'] + "\n")
            fh.write('%   \n')
            fh.write('% cg_opt - Optimization framework script which runs the BayesOpt package\n')
            fh.write('% to perform Gaussian Process Bayesian optimization\n')
            fh.write('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
     
            fh.write('\n')
            fh.write('\n')

            fh.write('% Initialize Gaussian process\n')
            fh.write("addpath %s\n" % self.GP_MATLAB_PATH)
            fh.write("startup % this initializes the GP framework\n")
            fh.write('\n')
            fh.write('delete goodness.txt\n % delete the goodness file if it exists')
            fh.write("num_opt = %i;\n" %  num_params)
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
                    
            fh.write("};\n")

            # write initial parameters
            writeline('for n = 1:num_opt',fh)
            writeline("   tempname=sprintf('%s_guess.txt', type_opt{n});", fh)
            writeline("   fileID=fopen(tempname, 'w');", fh)
            writeline("   fprintf(fileID, '%4.2f\n');"%0.0, fh)
            writeline("   % if we're dealing with a sigma value also write the coefficient", fh)
            writeline("   ", fh)
            writeline("   if strcmp(tempname(1:3), 'sig') == 1", fh)
            writeline("      tempname=sprintf('coff_%s_guess.txt', type_opt{n});", fh)
            writeline("      fileID=fopen(tempname, 'w');", fh)
            writeline("      fprintf(fileID, '%4.2f\n');" % 0.0 , fh)
            writeline("      fclose(fileID);", fh)
            writeline("   end", fh)
            writeline("end", fh)
            writeline("fileID=fopen('number_of_iterations.txt','w');", fh)
            writeline("fprintf(fileID, '%i');" % 0, fh)
            writeline("fclose(fileID);", fh)
            


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
        """


        """

        
        with open('cg_opt_func.m','w') as fh:
            
            writeline('function overlap =  cg_opt_func(x, num_opt, type_opt)', fh)
            writeline('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', fh)
            writeline('%   ', fh)
            writeline('%   ', fh)
            writeline('% THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    ', fh)
            writeline('% Code generated at ' + str(time.strftime("%c")), fh)
            writeline('% Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'], fh)
            writeline('%   ', fh)
            writeline('%   ', fh)
            writeline('% cg_opt_func - Optimization function used by the BayesOpt to perform', fh)
            writeline('% optimization', fh)
            writeline('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', fh)
            
            writeline('    ', fh)

            writeline('  % Get current iteration', fh)
            writeline("  iteration = load('number_of_iterations.txt');", fh)

            writeline("    ", fh)

            writeline("  % Add to iteration number", fh)
            writeline("  iteration = iteration + 1;", fh)                    
            writeline("    ", fh)
            writeline("  % Print iteration number to file", fh)
            writeline("  fileID=fopen('number_of_iterations.txt','w');", fh)
            writeline("  fprintf(fileID, '%d', iteration);", fh)
            writeline("  fclose(fileID);", fh)
            
            writeline("    ", fh)

            writeline("  % cycle through each parameter estimate and write to file", fh)
            writeline("  for n = 1:num_opt", fh)
            writeline("     tempname=sprintf('%s_guess.txt', type_opt{n});", fh)                     
            writeline("     fileID=fopen(tempname, 'w');", fh)
            writeline("     fprintf(fileID, '%4.2f\n', roundn(x(n),-2));", fh)
            writeline("     fclose(fileID);", fh)
            writeline("     %if we're dealing with a sigma value also write the coefficient", fh)
            writeline("    ", fh)
            writeline("     if strcmp(tempname(1:3), 'sig') == 1", fh)
            writeline("        tempname=sprintf('coff_%s_guess.txt', type_opt{n});", fh)
            writeline("        fileID=fopen(tempname, 'w');", fh)
            writeline("        fprintf(fileID, '%4.2f\n', 2.5*roundn(x(n),-2));", fh)
            writeline("        fclose(fileID);", fh)
            writeline("     end", fh)
            writeline("  end", fh)

            
            writeline("    ", fh)
            writeline("    ", fh)

            # this line defines the creation and updating of the OPT_MOLTEMPLATE_FILE to include the new parameters we just
            # wrote to file. update_moltemplate.py is a stand alone script which is 
            writeline("  % this line will run a custom generated script which will update ", fh)
            writeline("  % the moltemplate file with the current guesses ", fh)
            writeline("  system('" + self.PYTHON_BIN + " update_moltemplate.py');", fh)
            
            writeline("", fh)

            writeline("  % this line will launch the simulation and then poll until the ", fh)
            writeline("  % simulation has finished ", fh)
            writeline("  system('" + self.PYTHON_BIN + " run_simulation.py');", fh)
            
            writeline("", fh)

            writeline("  % this line will run custom generated script which will compares", fh)
            writeline("  % the CG simulation against the all atom simulation ", fh)
            writeline("  overlap = evaluate_simulation(iteration, %i)" % len(self.sequence_vector), fh)
            writeline("", fh)
            writeline("  % Save the current overlap score to a textfile for convenience ", fh)
            writeline("  fid = fopen('goodness.txt', 'at');", fh)
            writeline("  fprintf(fid, '%i   %4.2f\n', iteration, overlap);", fh)
            writeline("  fclose(fid);", fh)
            writeline('    ', fh)
            writeline('% END OF FILE --------------------------------', fh)
            
            

            

    def build_GaussianProcess_update_code(self):
        """
        Creates a file called update_moltemplate.py which is a customized script which updates the 
        moltemplate file we created previously with the current interation of guesses
        

        """


        (num_params,param_list) = self.determine_number_of_parameters_to_optimize()

        with open('update_moltemplate.py','w') as fh:

            writeline('import os', fh)
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)
            writeline('#   ', fh)
            writeline('#   ', fh)
            writeline('# THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    ', fh)
            writeline('# Code generated at ' + str(time.strftime("%c")) + "", fh)
            writeline('# Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'] + "", fh)
            writeline('#   ', fh)
            writeline('#   ', fh)
            writeline('# update_moltemplate.py', fh)
            writeline('# \n', fh)
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)

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
            fh.write("      current_guess[i] = 1.01\n")
            
            fh.write("\n")
            fh.write("   # sigma parmeters have a coefficient...\n")
            fh.write("   if i[0:3] == 'sig':\n")
            fh.write("      try:\n")
            fh.write("         with open('coff_' + i + '_guess.txt','r') as rfh:\n")
            fh.write("            current_guess['coff_'+i] = float(rfh.read())\n")
            fh.write('      except IOError:\n')
            fh.write("         current_guess['coff_'+i] = 1.01\n")
            fh.write("# The __TEMPLATE__ value updates a header iteration counter\n")
            fh.write("current_guess['__TEMPLATE__'] = str(iteration)\n")
            fh.write("\n")
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
            fh.write("      # note we add a leading space because we're looking for separate values\n")
            fh.write("      newline = newline.replace(' ' + key, ' ' + str(current_guess[key]))\n")
            fh.write("   newlines.append(newline)\n")
            fh.write("\n")
            fh.write("with open('%i/stMol.lt'%iteration,'w') as fh:\n")
            fh.write("   for line in newlines:\n")
            fh.write("      fh.write(line)\n")
            fh.write("\n")
            fh.write("\n")
            
                                                                                                
                     
    def build_nvt_minimization_file(self):
        
        with open("run.in.min", 'w') as fh:
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)
            writeline('#   ', fh)
            writeline('#   ', fh)
            writeline('# THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    ', fh)
            writeline('# Code generated at ' + str(time.strftime("%c")), fh)
            writeline('# Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'], fh)
            writeline('#   ', fh)
            writeline('#   ', fh)
            writeline('# run.in.min', fh)
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)

            writeline('# -------- CAMeLOT MINIMIZATION: ---------', fh)
            writeline('# This file defines a CAMeLOT minimization file', fh)
            writeline('  ', fh)
            writeline('# -- Init section --', fh)
            writeline('include system.in.init', fh)
            writeline('   ', fh)
            writeline('# -- Atom definition section -- ', fh)   
            writeline('read_data system.data', fh)
            writeline('   ', fh)
            writeline('# -- Settings Section --', fh)        
            writeline('include system.in.settings', fh)
            writeline('   ', fh)
            writeline('# -- Run section --', fh)
            writeline('dump     1 all custom 50 traj_min.lammpstrj id mol type x y z ix iy iz', fh)
            writeline('minimize 1.0e-5 1.0e-7 10000 10000', fh)
            writeline('write_restart  system_after_min.rst',fh)
            writeline('   ', fh)
            writeline('#------------- END OF MINIMIZATION FILE --------------', fh)


    def build_system_setup_file(self):
        
        with open("system.lt", 'w') as fh:
            writeline('#   ', fh)
            writeline('#   ', fh)
            writeline('# THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    ', fh)
            writeline('# Code generated at ' + str(time.strftime("%c")), fh)
            writeline('# Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'], fh)
            writeline('#   ', fh)
            writeline('#   ', fh)
            writeline('# system.lt', fh)
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)
            writeline('   ', fh)
            writeline('import "stMol.lt"', fh)
            writeline('write_once("Data Boundary") {', fh)
            writeline('  0 255.1 xlo xhi', fh)
            writeline('  0 255.1 ylo yhi', fh)
            writeline('  0 255.1 zlo zhi', fh)
            writeline('}', fh)
            writeline('   ',fh)
            writeline('# Create a single peptide', fh)
            writeline('peptide = new stMol [1].move(50.0, 100.0, 100.0) ', fh)
            writeline('    ', fh)
            writeline('#------------- END OF SYSTEM DEFINITION --------------', fh)
      

    def build_nvt_simulation_code(self):
        """
        Function which constructs custom code which itself constructs a LAMMPS configuration file
        needed to run a production NVT simulation.

        Right now (0.1.0) we don't allow customization of the NVT configuration file, but the whole
        point of creating a function which writes code which writes a configuration file is that
        you can configure that second codeset to include custom variables as defined during the
        initial setup. 

        The point is, the 'generate_nvt.py' file which ends up getting created is a stand alone
        custom file which will build a system-bespoke LAMMPS configuration file as a totally
        independent script.

        This means we can use it as part of the MATLAB drive GPBO procedure, where the MATLAB
        code calls these generate_nvt.py script, and keep the idea of creating a system specific
        simulation set for optimization.

        """



        with open('generate_nvt.py', 'w') as fh:

            writeline('import random', fh)
            writeline('import time', fh)
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)
            writeline('#   ', fh)
            writeline('#   ', fh)
            writeline('# THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    ', fh)
            writeline('# Code generated at ' + str(time.strftime("%c")), fh)
            writeline('# Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'] + "", fh)
            writeline('#', fh)
            writeline('#', fh)
            writeline('# generate_nvt.py', fh)
            writeline('# This code is to be executed in the main simulation root subdirectory, and will', fh)
            writeline('# generate an NVT production simulation LAMMPS configuration file', fh)
            writeline('# with a unique and random SEED.', fh)
            writeline('# That config. file should then be moved into an iteration subdirectory', fh)
            writeline('# where the simulation should be run', fh)
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)
            writeline('', fh)
            writeline('# >>>>>>>>>>>>>>>>>>>>>', fh)
            writeline('# First up we need to define the Lagevin string...', fh)
            writeline('', fh)
            writeline("# Seed the PRNG using the OS' entropy pool", fh)
            writeline("random.seed()", fh)
            writeline('', fh)

            writeline("# First read in the dampening parameters and construct", fh)
            writeline("# a custom scale line", fh)                     
            writeline("with open('%s', 'r') as fh:" % self.DAMP_PARAMS, fh)
            writeline("   DPs = fh.readlines()", fh)            
            writeline("DP = DPs[0].strip()", fh)
            writeline(' ', fh)
                            
            writeline("TEMP = %3.3f"%self.TEMP, fh)
            writeline("RAND = random.randint(0,10000)", fh)
            writeline("dampline = 'fix fxlan  all langevin ' + str(TEMP) + ' ' + str(TEMP) + ' 2000.0 ' + str(RAND) + ' ' + DP", fh)
            
            writeline("with open('run.in.nvt', 'w') as fh:", fh)

            writeline("   fh.write('# -------- CAMeLOT NVT simulation: ---------\n')"   , fh)            
            writeline("   fh.write('# This file defines a CAMeLOT simulation file,\n')", fh)
            writeline("   fh.write('# and was generated by the generate_nvt.py script\n')", fh)
            writeline("   fh.write('# at ' + str(time.strftime('%c')) +'\n')", fh)

            writeline("   fh.write('\n')", fh)
                
            writeline("   fh.write('# -- Init section --\n')", fh)
            writeline("   fh.write('include system.in.init\n')", fh)
                
            writeline("   fh.write('\n')", fh)

            writeline("   fh.write('# -- Atom definition section -- \n')", fh)   
            writeline("   fh.write('read_restart system_after_min.rst\n')", fh)

            writeline("   fh.write('\n')", fh)
                
            writeline("   fh.write('# -- Settings Section --\n')    ", fh)    
            writeline("   fh.write('include system.in.settings\n')", fh)
            
            writeline("   fh.write('\n')", fh)
                
            writeline("   fh.write('# -- Run section --\n')", fh)
            writeline("   fh.write('timestep        2.0\n')", fh)
            writeline("   fh.write('dump            1 all dcd %i optimization.dcd\n')" % self.SIM_DCD_OUT, fh)
            writeline("   fh.write('neigh_modify    delay 1\n')", fh)
            writeline("   fh.write('\n')", fh)
            writeline("   fh.write('# -- Langevin section --\n')", fh)
            writeline("   fh.write('# To use Langevin dynamics in LAMMPS you need both fix langevin and fix nve.\n')", fh)
            writeline("   fh.write('# (See http://lammps.sandia.gov/doc/fix_langevin.html for details.)\n')", fh)
            writeline("   fh.write('\n')", fh)
            writeline("   fh.write('%s\n'%dampline)", fh)                
            writeline("   fh.write('fix fxnve all nve\n')", fh)
            writeline("   fh.write('\n')", fh)
            writeline("   ", fh)                
            writeline("   fh.write('thermo_style    custom step temp pe etotal press vol epair ebond eangle edihed\n')", fh)
            writeline("   fh.write('thermo          %i  # time interval for printing out thermo data\n')"%self.SIM_THERMO_OUT, fh)
            writeline("   ", fh)                
            writeline("   fh.write('restart         %i  restart_nvt\n')"%self.SIM_NSTEPS, fh)
            writeline("   fh.write('run		%i\n')"%self.SIM_NSTEPS, fh)
            writeline("   fh.write('write_restart  system_after_nvt.rst\n')", fh)    
            writeline("   fh.write('\n')", fh)
            writeline("   fh.write('\n')", fh)
      


    
        

    def build_GaussianProcess_simulation_code(self):
        """
        Creates a file which runs the actual LAMMPS simulation

        """
        
        with open('run_simulation.py', 'w') as fh:
            writeline('from shutil import copy', fh)            
            writeline('import os', fh)          
            writeline('import subprocess', fh)          
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)
            writeline('#   ', fh)
            writeline('#   ', fh)
            writeline('# THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    ', fh)
            writeline('# Code generated at ' + str(time.strftime("%c")), fh)
            writeline('# Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'] + "", fh)
            writeline('#', fh)
            writeline('#', fh)
            writeline('# run_simulation.py', fh)
            writeline('#', fh)
            writeline('# This code should be called by the MATLAB code cg_opt_func.m', fh)
            writeline('# and does the following things', fh)
            writeline('#', fh)
            writeline('# 1) Creates a new run.nvt.in configuration file (new SEED)', fh)
            writeline('#', fh)
            writeline('# 2) Moves that run.nvt.in file into the next iteration directory, which', fh)
            writeline('#    should already exist having been created by update_moltemplate.py, which', fh)
            writeline('#    itself should have been called immedietly preceding this function', fh)
            writeline('#', fh)
            writeline('# 3) Copies the standard run.in.min into the next iteration directory', fh)
            writeline('#', fh)
            writeline('# 4) Copies the standard system.lt file into the next iteration directory', fh)
            writeline('#', fh)
            writeline('# 5) Runs moltemplate on the system', fh)
            writeline('#', fh)
            writeline('# 6) launch the LAMMPS simulation in that directory', fh)
            writeline('#', fh)
            writeline('# 7) blocks until completed, then returns', fh)          
            writeline('#', fh)
            writeline('#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><>', fh)
            writeline(' ', fh)
            writeline("# Before we do that we have to figure out what iteration we're on", fh)
            writeline("", fh)
            writeline("try:", fh)
            writeline("   with open('number_of_iterations.txt','r') as rfh:", fh)
            writeline("      iteration = int(rfh.read())", fh)
            writeline("except IOError:", fh)
            writeline("   iteration = 0", fh)
            writeline("", fh)
            writeline('# 1) make run.nvt.in', fh)
            writeline('os.system("%s generate_nvt.py")'%self.PYTHON_BIN, fh)
            
            writeline('# 2) move run.nvt.in', fh)
            writeline("os.rename('run.in.nvt', '%i/run.in.nvt'%iteration)", fh)

            writeline('# 3) copy run.min.in', fh)
            writeline("copy('run.in.min', '%i/run.in.min'%iteration)", fh)

            writeline('# 4) copy system.lt', fh)
            writeline("copy('system.lt', '%i/system.lt'%iteration)", fh)

            writeline('# 5) run moltemplate', fh)
            writeline('CWD = os.getcwd()', fh)
            writeline('simdir = "%s/%i"%(CWD, iteration)',fh)
            writeline("p = subprocess.Popen(['sh', '/work/alex/coarse_graining/CAMeLOT/moltemplate/moltemplate3.sh', 'system.lt'], cwd=simdir)", fh)
            writeline("p.wait()",fh)
            
            
            writeline("p = subprocess.Popen(['/home/kruff/LAMMPS/lammps-16Dec13/src/lmp_openmpi5', '-i', 'run.in.min'], cwd=simdir)", fh)
            writeline("p.wait()",fh)
            
            writeline("p = subprocess.Popen(['/home/kruff/LAMMPS/lammps-16Dec13/src/lmp_openmpi5', '-i', 'run.in.nvt'], cwd=simdir)", fh)
            writeline("p.wait()",fh)
            


            
        

    def build_GaussianProcess_target_code(self):

        with open('evaluate_simulation.m','w') as fh:

            writeline('function mean_overlap = evaluate_simulation(iteration, nres)', fh)
            writeline('   ', fh)
            writeline('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', fh)
            writeline('%   ', fh)
            writeline('%%   ', fh)
            writeline('%% THIS CODE WAS DYNAMICALLY GENERATED BY CAMeLOT    ', fh)
            writeline('%% Code generated at ' + str(time.strftime("%c")), fh)
            writeline('%% Generated by CAMeLOT version ' + self.MolTemplateParameters['CAMELOT_VERSION'], fh)
            writeline('%%   ', fh)
            writeline('%% cg_opt - Optimization framework script which runs the BayesOpt package', fh)
            writeline('%% to perform Gaussian Process Bayesian optimization', fh)
            writeline('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%', fh)
            writeline('   ', fh)
            writeline('% source the dcd reading functionality', fh)
            writeline('addpath /work/alex/tools/matlab_functions/matdcd-1.0/   ', fh)            
            writeline('% box length', fh)
            writeline('bl = [255.1, 255.1, 255.1];', fh)
            writeline('      ', fh)
            writeline('% discard the first 20% as equlibration', fh)
            writeline('equil = %i ;'%(self.SIM_EQUIL_FRACTION*(self.SIM_NSTEPS/self.SIM_DCD_OUT)), fh)
            writeline('binwidth=%3.3f;'%(self.RES_RES_BIN_SIZE), fh)
            writeline('binhalf=binwidth/2;', fh)
            writeline('      ', fh)
            writeline('% remove the first 1e5 steps', fh)
            writeline('equil = 101;', fh)                        
            writeline('   ', fh)
            writeline('% get bead positions', fh)
            writeline("da = readdcd([num2str(iteration) '/optimization.dcd'],1:nres);",fh)
            writeline('xyzmat=da(equil:end,:);', fh)     
            writeline('        ', fh)       
            writeline('for i=1:nres-1', fh)
            writeline('    for j=1+1:nres', fh)
            writeline('        positions_i = xyzmat(:, 1+3*(i-1):3+3*(i-1));', fh)
            writeline('        positions_j = xyzmat(:, 1+3*(j-1):3+3*(j-1));', fh)
            writeline('        ', fh)
            writeline('        dPos = positions_i - positions_j;', fh)
            writeline('        % correct for periodicity in x/y/z', fh)
            writeline('        xcorrect=find(abs(dPos(:,1))>bl(1)*0.5);', fh)
            writeline('        dPos(xcorrect,1)=bl(1)-abs(dPos(xcorrect,1));', fh)
            writeline('        ', fh)
            writeline('        ycorrect=find(abs(dPos(:,2))>bl(2)*0.5);', fh)
            writeline('        dPos(ycorrect,2)=bl(2)-abs(dPos(ycorrect,2));', fh)
            writeline('        ', fh)
            writeline('        zcorrect=find(abs(dPos(:,3))>bl(3)*0.5);', fh)
            writeline('        dPos(zcorrect,3)=bl(3)-abs(dPos(zcorrect,3));', fh)
            writeline('        ', fh)
            writeline('        % determine full set of distance value', fh)
            writeline('        raw_distance{i,j} = sqrt(sum(dPos.^2,2));', fh)
            writeline('        ', fh)
            writeline("        % note we have to offset the bin edges because numpy and matlab have different", fh)
            writeline("        % interepetations of how the edges should be used...", fh)
            writeline("        vals = histc(raw_distance{i,j},%3.3f:%3.3f:%3.3f-(2*binwidth))';" % (self.RES_RES_BIN_START, self.RES_RES_BIN_SIZE, self.RES_RES_BIN_END), fh)
            writeline('        vals = vals./sum(vals);', fh)
            writeline('        ', fh)
            writeline('        CG_histogram{i,j}  =  vals;', fh)
            writeline('    end', fh)
            writeline('end', fh)
            writeline('        ', fh)
            writeline('        ', fh)
            writeline('% extract pre-computed all atom distance histograms', fh)            
            writeline("AA_distance_distributions=load('%s');"%(self.RES_RES_DISTANCE_FILE), fh)
            writeline('        ', fh)
            writeline('for i = 1:size(AA_distance_distributions,1)', fh)
            writeline('    line = AA_distance_distributions(i,:);', fh)
            writeline('    res1 = line(1);', fh)
            writeline('    res2 = line(2);', fh)
            writeline('    AA_histogram{res1,res2} = line(3:end);', fh)
            writeline('end', fh)
            writeline('        ', fh)
            writeline('% compute all overlaps', fh)
            writeline('for i=1:nres-1', fh)
            writeline('    for j=i+1:nres', fh)
            writeline('       overlap(i,j) =  (2-sum(abs(AA_histogram{i,j} -  CG_histogram{i,j})))/2;', fh)
            writeline('    end', fh)
            writeline('end', fh)
            writeline('     ', fh)
            writeline('% create a vector of the relevant overlaps as defined by the minimium', fh)
            writeline('% separation of interest (minsep)', fh)
            writeline('minsep=5;', fh)
            writeline('relevant = [];', fh)
            writeline('for i=1:nres-minsep', fh)
            writeline('    for j=i+minsep:nres', fh)
            writeline('       relevant = [relevant; (2-sum(abs(AA_histogram{i,j} -  CG_histogram{i,j})))/2];', fh)
            writeline('    end', fh)
            writeline('end',fh)
            writeline('      ',fh)
            writeline('mean_overlap = 1 - mean(relevant);', fh)                      
            writeline('   ', fh)            
            writeline('   ', fh)            
            writeline('   ', fh)            
            writeline('   ', fh)
            

  
        



            
                     
