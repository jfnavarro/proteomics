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
    print "this script generates a PROSOLVE txt input file from a percolator pout.xml file"
    print "Usage : generateProsolveInput.py <pout.xml> [-o, --output <output.txt>] [-h, --help] [-v, --verbose] [-e, --hidden] [-m, --maxpsms]"
    print "--hidden : the percolator file was run with hidden decoys in the target database (normal decoys will be discarded)"
    print "--maxpsms : the output file containing PSMs will have a maximum of maxpsms number of PSMs per peptide"
    print "--output : the name of the two output file containing PSMs and peptides"
    
def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        magnusOutputFile = "output.txt"
        verbose = False
        hidden = False
        maxpsms = 1000
        try:
            opts, args = getopt.getopt(sys.argv[2:], "o:hvem:", ["output=", "help", "verbose","hidden","maxpsms"])
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
                magnusOutputFile = a
            elif o in ("-e", "--hidden"):
                hidden = True
            elif o in ("-m", "--maxpsms"):
                maxpsms = a
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
        
        psms = dict()
        for psm in percolatorPSMs:
            peptide_clean = psm.peptide[:-2][2:]
            if psms.has_key(peptide_clean):
                if(psms[peptide_clean].score < psm.score):
                    psms[peptide_clean] = psm
            else:
                psms[peptide_clean] = psm
        
        percolatorPeptides = per.getPeptides(elems)
        
        if(len(percolatorPSMs) > 0 and len(percolatorPeptides) > 0):
            if(verbose):
                print "writing in " + str(magnusOutputFile) + " with max psms = " + str(maxpsms)
            utils.writeMagnusInput(percolatorPSMs,"psms_"+magnusOutputFile,hidden,maxpsms)
            utils.writeMagnusPeptides(percolatorPeptides,psms,"peptides_"+magnusOutputFile,hidden)
        else:
            sys.stderr.write("the input file does not contain any peptide or psms\n")
            sys.exit()
            
        if(verbose):
            print "Magnus files generated"
        
if __name__ == "__main__":
    main(sys.argv[1:]) 
    