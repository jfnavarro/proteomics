#! /usr/bin/env python
# @Created by Jose Fernandez
# dec 10th, 2011

from lxml import etree
import os
import getopt
import sys
import percolatorUtils as per
import generalUtils as utils
    
def usage():
    print "this script generates a MAYUS csv input file from a percolator pout.xml file"
    print "Usage : generateMayuInput.py <pout.xml> [-o, --output <output.csv>] [-h, --help] [-v, --verbose]"
    print "--output : the name of the output file that will contain the PSMs in mayu csv format"
    
def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        mayuOutputFile = "output.csv"
        verbose = False
        try:
            opts, args = getopt.getopt(sys.argv[2:], "o:hv", ["output=", "help", "verbose"])
        except getopt.GetoptError, err:
            # print help information and exit:
            print str(err) # will print something like "option -a not recognized"
            usage()
            sys.exit(2)
        
        for o, a in opts:
            if o == "-v":
                verbose = True
            elif o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-o", "--output"):
                mayuOutputFile = a
            else:
                assert False, "unhandled option"
        
        if(os.path.isfile(argv[0])):
            infile = argv[0]
        else:
            sys.stderr.write("Error: XML file not found\n")
            sys.exit()
                
        parser = etree.XMLParser(ns_clean=False, huge_tree=True)
        try:
            tree = etree.parse(infile,parser)
        except Exception, inst:
            sys.stderr.write("Unexpected error opening %s: %s\n" % (infile, inst))
            sys.exit()
             
        if(verbose):
            print "Reading " + str(argv[0])   
                         
        elems = tree.getroot()
        percolatorPSMs = per.getPSMs(elems)
        
        if(len(percolatorPSMs) > 0):
            if(verbose):
                print "writing in " + str(mayuOutputFile)   
            utils.writeMayuInput(percolatorPSMs,mayuOutputFile)
        else:
            sys.stderr.write("the input file does not contain any psm")
            sys.exit()
            
        if(verbose):
            print "Mayu file generated"
        
if __name__ == "__main__":
    main(sys.argv[1:]) 
    