#! /usr/bin/env python
# @Created by Jose Fernandez
# dec 10th, 2011
# This scripts generates an input file compatible with Fido from a percolator output file


from lxml import etree
import os
import getopt
import sys
import percolatorUtils as per
import  generalUtils as utils

def usage():
    print "this script generates a psm-proteins and unique peptide-proteins graph files for FIDO as well as a file containing all the target and decoy proteins names from a input file from a percolator pout.xml file"
    print "Usage : generateFidoInput.py <pout.xml> [-o, --output <output.txt>] [-p, --pattern ] [-d, --database <database.txt>] [-h, --help] [-v, --verbose]"
    print "--output : the name of the file that will contain the PSMs, probabilities and proteins"
    print "--database : the name of the file that will contain the target and decoy protein names"
    print "--pattern : the pattern that identifies the decoy hits"
    
def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        fidoOutputFile = "output.txt"
        fidoDBfile = "db.txt"
        verbose = False
        decoy_prefix = "random"
        try:
            opts, args = getopt.getopt(sys.argv[2:], "o:p:d:hv", ["output=", "pattern=", "database=", "help", "verbose"])
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
                fidoOutputFile = a
            elif o in ("-p", "--pattern"):
                decoy_prefix = a
            elif o in ("-d", "--database"):
                fidoDBfile = a
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
        percolatorPeptides = per.getPeptides(elems)
        percolatorPSMs = per.getPSMs(elems)
        
        if(len(percolatorPSMs) > 0):
            if(verbose):
                print "writing in " + str(fidoOutputFile)   
            utils.writeFidoInput(percolatorPSMs,"psms_" + fidoOutputFile,fidoDBfile,decoy_prefix)
            if(len(percolatorPeptides) > 0):
                utils.writeFidoInput(percolatorPeptides,"peptides_" + fidoOutputFile,fidoDBfile,decoy_prefix)
        else:
            sys.stderr.write("the input file does not contain any peptide or psms\n")
            sys.exit()
            
        if(verbose):
            print "Fido file generated"
        
if __name__ == "__main__":
    main(sys.argv[1:]) 
    