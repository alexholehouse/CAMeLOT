## CAMeLOT V 0.1.0
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
##

import time
from CAMeLOTExceptions import MolTemplateException

class MolTemplate:


    def __init__(self, parameter_dictionary):
        """
        MolTemplate must be initialized with a dictionary containing the following parameters
        which are used for the BD simulations
        
        kappa            - inverse of the Debye length (inverse distance units)
        lj-cutoff        - global cutoff for LJ term (default is 24.1 [angstroms])
        coulombic-cutoff - global cutoff for Coulombic term (default is 15.0 [angstroms])
        dielectric       - environment dielectric  (default is 80)
        mixing-rule      - LJ mixing rule (default is arithmetic)
        CAMELOT-version  - Version of CAMeLOT being run


        """

        # initialize the status flag
        self.MTFStatus = 'not ready'

        # note this could EASILY be extended to give more options
        parameter_keywords = ['OPT_KAPPA', 'OPT_LJ_CUTOFF', 'OPT_LJ_CUTOFF', 'OPT_COUL_CUTOFF', 'OPT_MOLTEMPLATE_FILE', 'OPT_DIELECT', 'OPT_MIXING', 'CAMELOT_VERSION']

        self.parameters = {}

        # cycle through the expected keywords as defined above and raise exception if an associated value
        # is not found
        for keyword in parameter_keywords:
            try:
                self.parameters[keyword] = parameter_dictionary[keyword]
            except KeyError:
                raise MolTemplateException('Did not provide MolTemplate with a %s parameter'%keyword)
            

    def build_moltemplate_input(self, intial_positions_file, 
                                mass_parameters_file, 
                                pairwise_potential_file, 
                                bonds_definition_file,
                                bonds_parameters_file, 
                                angle_definition_file,
                                angle_parameters_file,
                                dihedral_definition_file,
                                dihedral_parameters_file,
                                mbt_parameters_file,
                                ebt_parameters_file,
                                at_parameters_file,
                                aat_parameters_file,
                                bb13_parameters_file):

        """
        Moltemplate is a general way to define molecules/forcefields for LAMMPS simulations. It provides a simple templating
        language which can be used to build LAMMPS compatible input files. CAMeLOT makes extensive use of Moltemplate to 
        define the simulation system. For more information on Moltemplate see http://www.moltemplate.org/

        build_moltemplate_input is a single function which constructs a moltemplate input file, which can then be run
        using moltemplate to build a LAMMPS compatible file


        """


        # first initialize the moltemplate file and write all the relevant header information
        self.initialize_file()

        
        # next write all the sections required
        self.write_section(intial_positions_file, 'Data atoms', 'write', ['atomId molId   atomType   charge   x      y        z'])

        self.write_section(mass_parameters_file, 'Data masses', 'write_once', ['atomType  mass'])

        self.write_section(pairwise_potential_file, 'In Settings', 'write_once', ['As usual, force-field parameters ("coeffs") go in the "In Settings" section:','Pairwise (non-bonded) interactions:','          atomType1 atomType2   epsilon sigma cutoff_lj'])

        self.write_section(bonds_definition_file, 'Data Bonds', 'write_once', ['bond-id   bond-type      atom-id1  atom-id2'])

        self.write_section(bonds_parameters_file, 'In Settings', 'write_once', ['            bond-type        k     r0'])

        self.write_section(angle_definition_file, 'Data Angles By Type', 'write_once', ['angle-type      atomType1 atomType2 atomType3  bondType1 bondType2'])
        
        self.write_section(angle_parameters_file, 'In Settings', 'write_once', ['               angle-type         k    theta0'])

        self.write_section(dihedral_definition_file, 'Data Dihedrals By Type', 'write_once', ['dihedral-type AtomType1 AtomType2 AtomType3 AtomType4 bondType1 btyp2 btyp3'])

        self.write_section(dihedral_parameters_file, 'In Settings', 'write_once', ['dihedral_coeff dihedralType  K1  phi1   K2  phi2 K3 phi3'])

        self.write_section(mbt_parameters_file, 'In Settings', 'write_once', ['dihedral_coeff dihedralType  A1 A2 A3 r2', 'MiddleBondTorsion Coeffs'])
        
        self.write_section(ebt_parameters_file, 'In Settings', 'write_once', ['dihedral_coeff dihedralType  B1 B2 B3 C1 C2 C3 r1 r3','EndBondTorsion Coeffs'])

        self.write_section(at_parameters_file, 'In Settings', 'write_once', ['dihedral_coeff dihedralType  D1 D2 D3 E1 E2 E3 theta1 theta2','AngleTorsion Coeffs'])

        self.write_section(aat_parameters_file, 'In Settings', 'write_once', ['dihedral_coeff dihedralType  M theta1 theta2', 'AngleAngleTorsion Coeffs'])

        self.write_section(bb13_parameters_file, 'In Settings', 'write_once', ['dihedral_coeff dihedralType  N r1 r3','BondBond13 Coeffs'])

        self.finalize_file()
        
        
    def initialize_file(self):
        """
        Function which ensures we're writing sections to an empty file and then
        writes the moltemplate header file

        """
        
        # this just zero's out the file and will fail if we can't write to the directory
        with open(self.parameters['OPT_MOLTEMPLATE_FILE'], 'w') as fh:
            fh.write("")

        # update the staus flag
        self.MTFStatus='ready'

        now = time.strftime("%c")

        with open(self.parameters['OPT_MOLTEMPLATE_FILE'], 'a') as fh:
            fh.write("stMol {\n")
            fh.write("\n")
            fh.write(" # This file was generated by CAMeLOT version %s \n" % str(self.parameters['CAMELOT_VERSION']))
            fh.write(" # File generated at %s\n" % time.strftime("%c"))
            fh.write(" \n")
            fh.write(' write_once("In Init") {\n')
            fh.write('  # -- Default Styles for "2bead" --\n')
            fh.write('  units           real\n')
            fh.write('  atom_style      full\n')
            fh.write('  bond_style      harmonic\n')
            fh.write('  angle_style     harmonic\n')
            fh.write('  pair_style      lj/cut/coul/debye %3.1f %3.1f %3.1f \n' % (self.parameters['OPT_KAPPA'], self.parameters['OPT_LJ_CUTOFF'], self.parameters['OPT_COUL_CUTOFF']) )
            fh.write('  dihedral_style  class2 \n')
            fh.write('  pair_modify     mix %s\n' % self.parameters['OPT_MIXING'])
            fh.write('  dielectric      %3.1f \n' % self.parameters['OPT_DIELECT'])
            fh.write(' }\n')
            fh.write('\n')
            

    def finalize_file(self):
        
        # check we've initialized the file for writing
        if not (self.MTFStatus == "ready"):
            raise MolTemplateException("File not ready to be written to - write_section() called before initialize_file()")
        
        with open(self.parameters['OPT_MOLTEMPLATE_FILE'], 'a') as fh:
            
            # write the opening section definition
            fh.write(' #  Note: "$atom:CA" is a variable storing the atom ID for a "CA" atom in a\n') 
            fh.write(' #        molecule.  Atoms in different molecules have unique ID numbers.\n')
            fh.write(' #\n')
            fh.write(' #  Note: "@atom:CA" denotes an atom type number (corresponding to the "CA" \n')
            fh.write(' #        atom type in the 2bead molecule).  CA atoms in different molecules \n')
            fh.write(' #        share the same type number.\n')
            fh.write(' #\n')
            fh.write(' #  Note: The "..." in "$mol:..." tells moltemplate that this molecule may be a\n')
            fh.write(" #        part of a larger molecule, and (if so) to use the larger object's\n")
            fh.write(" #        molecule ID number as it's own.\n")

            fh.write("\n")
            fh.write("} # stMol\n")

        # update write status flag
        self.MTFStatus='not-ready'

        print "Finalized moltemplate file!"
    


    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    def write_section(self, filename, section_name, write_style, comments):
        """
        write_section writes a moltemplate section to an existing file


        filename is the parameter file which is going to be written to this
        section.

        section_name is a string giving the section name

        write_style is one of two options ('write' or 'write_once')

        comments is a list of comments which go before the main data


        """

        # define allowed write-styles and check we comply!
        allowed_write_styles = ["write","write_once"]        
        if write_style not in allowed_write_styles:
            raise MolTemplateException("Undefined write style provided when writing moltemplate file")
            
        # read in the parameter file
        with open(filename,'r') as fh:
            content = fh.readlines()

        # check we've initialized the file for writing
        if not self.MTFStatus == 'ready':
            raise MolTemplateException("File not ready to be written to - write_section() called before initialize_file()")
        
        with open(self.parameters['OPT_MOLTEMPLATE_FILE'], 'a') as fh:
            
            # write the opening section definition
            fh.write(' %s("%s") {\n' % (write_style,section_name))

            # write comments
            for line in comments:
                fh.write('   # %s\n' % line)
            
            for line in content:
                fh.write("   %s"%line)
            fh.write(" }\n\n")

        print "Wrote %s to moltemplate configuration file..." % section_name
                
                        
        


                               
