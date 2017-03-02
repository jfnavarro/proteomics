#! /usr/bin/env python
# @Created by Jose Fernandez

from lxml import etree
import os
import getopt
import sys
import percolatorUtils as per
import generalUtils as utils

def usage():
    print "this script generates a tab delimited output for PSMs, peptides and proteins from a percolator pout.xml file"
    print "Usage : parsePercolator.py <pout.xml> [-o, --output <output.txt>] [-h, --help] [-c, --cutoff] [-v, --verbose]"
    print "--output : the name of the file that will contain the tab delimited list of psms|peptides|proteins and probabilities"
    print "--cutoff : you want to report only the psms,peptides and proteins below a certain threshold?"
    print "--qvalue : uses the qvalues insted of the PEPs for the threshold"
    
def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        OutputFile = "output.txt"
        verbose = False
        fdr = 1.0
        try:
            opts, args = getopt.getopt(sys.argv[2:], "o:hvc:", ["output=", "help", "verbose", "cutoff="])
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
                OutputFile = a
            elif o in ("-f", "--fdr"):
                fdr = float(a)
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
        percolatorPeptides = per.getPeptides(elems)
        percolatorProteins = per.getProteins(elems)
        
        if(verbose):
            print "Read " + str(len(percolatorPSMs)) + " PSMs"
            print "Read " + str(len(percolatorPeptides)) + " Peptides"
            print "Read " + str(len(percolatorProteins)) + " Proteins"
            
        if(fdr < 1.0 and fdr > 0.0):
            percolatorPSMs = [psm for psm in percolatorPSMs if psm.qvalue <= fdr]
            percolatorPeptides = [pep for pep in percolatorPeptides if pep.qvalue <= fdr]
            percolatorProteins = [prot for prot in percolatorProteins if prot.qvalue <= fdr]
          
        if(len(percolatorPSMs) > 0):
            if(verbose):
                print "writing in " + str(OutputFile)   
            utils.writePsms(percolatorPSMs, "psms_" + OutputFile)
            if(len(percolatorPeptides) > 0):
                utils.writePeptides(percolatorPeptides,"peptides_" + OutputFile)
            if(len(percolatorProteins) > 0):
                utils.writeProteins(percolatorProteins,"proteins_" + OutputFile)
        else:
            sys.stderr.write("the input file does not contain any information\n")
            sys.exit()
        
        if(verbose):
            print "Percolator file parsed"
            
if __name__ == "__main__":
    main(sys.argv[1:]) 
    