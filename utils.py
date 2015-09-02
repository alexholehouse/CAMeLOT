## CAMeLOT V 0.1.0
## Coarse-grained simulations Aided by Machine Learning, Optimization and Training
## Pappu lab, Washington University in St. Louis
## Copyright 2015
##

def threeToOne(res):
    
    THREE_TO_ONE = {'ALA':'A', 
                    'CYS':'C',
                    'ASP':'D',
                    'GLU':'E',
                    'PHE':'F',
                    'GLY':'G',
                    'HIS':'H', 
                    'HID':'H', 
                    'HIE':'H', 
                    'ILE':'I',
                    'LYS':'K',
                    'LEU':'L',
                    'MET':'M',
                    'ASN':'N',
                    'PRO':'P',
                    'GLN':'Q',
                    'ARG':'R',
                    'SER':'S',
                    'THR':'T',
                    'VAL':'V',
                    'TRP':'W',
                    'TYR':'Y'}
    return THREE_TO_ONE[res]


def print_logo(version):
    print ""
    print "  ________  ________  _____ ______   _______   ___       ________  _________    "
    print " |\   ____\|\   __  \|\   _ \  _   \|\  ___ \ |\  \     |\   __  \|\___   ___\  "
    print " \ \  \___|\ \  \|\  \ \  \ \__\ \  \ \   __/|\ \  \    \ \  \ \  \|___ \  \_|  "
    print "  \ \  \    \ \   __  \ \  \|__|  \  \ \  \_|/_\ \  \    \ \  \ \  \   \ \  \   "
    print "   \ \  \____\ \  \ \  \ \  \    \ \  \ \  \_|\ \ \  \____\ \  \_\  \   \ \  \  "
    print "    \ \_______\ \__\ \__\ \__\    \ \__\ \_______\ \_______\ \_______\   \ \__\ "
    print "     \|_______|\|__|\|__|\|__|     \|__|\|_______|\|_______|\|_______|    \|__| "
    print ""
    print "Version: " + str(version)

def writeline(line, fh):
    """
    writeline is a stupidly simple function which takes a string and writes it to the open 
    file associated with fh in exatcly the formatted it's provided (i.e. all custom characters
    are include) with a trailing newline.

    This makes writing code that generates code or config files much easier...
    """
    fh.write("%s\n" % repr(line)[1:-1])

