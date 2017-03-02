#! /usr/bin/env python
# @Created by Jose Fernandez
# dec 10th, 2011

import sys
import os
import getopt
import twohitProteinEstimator as twohit
import parsimonyProteinEstimator as parsimony
import percolatorUtils as per
from lxml import etree

def usage():
    print "Takes as input a percolator output file and compute the prior belief a protein is present using parsimony"
    print "Usage : percolatorGammaEstimaator.py pout.xml  [v, --verbose] [h, --help] [f, --fdr]"
    print "--fdr threshold of the peptides"
    
def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        fdr = 0.01
        verbose = False
        try:
            opts, args = getopt.getopt(sys.argv[2:], "f:", ["fdr="])
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
            elif o in ("-f", "--fdr"):
                fdr = a
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
        
        proteins_percolator= per.getProteins(elems)
        peptides = per.getPeptides(elems)
        
        ntarget_percolator = len([prot for prot in [prot for prot in proteins_percolator if not prot.isdecoy]])
        ndecoy_percolator = len([prot for prot in [prot for prot in proteins_percolator if prot.isdecoy]])
        
        print "Found " + str(ntarget_percolator) + " proteins"
        
        print "doing parsimony at " + str(fdr)
        
        parsimony_tmp = [prot for prot,pep in parsimony.getParsimony(peptides, fdr, True)]
        ntarget_parsimony = len([prot for prot in parsimony_tmp if prot.find("random") == -1])
        ndecoy_parsimony = len([prot for prot in parsimony_tmp if prot.find("random") != -1])
        
        print "Found " + str(ntarget_parsimony) + " proteins at " + str(fdr) + " fdr doing parsimony"
        
        print "Computed Gamma is " + str(float(ntarget_parsimony) / float(ntarget_percolator))
        
if __name__ == "__main__":
    main(sys.argv[1:]) 
    
    