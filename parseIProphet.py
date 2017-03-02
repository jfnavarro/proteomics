#! /usr/bin/env python
# @Created by Jose Fernandez
'''
Created on Jun 18, 2012

@author: jfn

This script parses files from iProphet or Protein prophet to generate tab delimited files containing all the PSMs, peptides and proteins
'''
#! /usr/bin/env python
# @Created by Jose Fernandez

from lxml import etree
import os
import getopt
import sys
import tppUtils as tpp
import generalUtils as utils

def usage():
    print "this script generates a tab delimited output for PSMs, peptides and proteins from iProphet and/or Protein Prophet output files"
    print "Usage : parseIProphet.py <interact.iproph.pep.xml> <interact.iproph.prot.xml> [-o, --output <output.txt>] [-d, --decoy <prefix>] [-h, --help] [-v, --verbose]"
    print "--output : the name of the file that will contain the tab delimited list of proteins/peptides/psms and probabilities"
    print "--prefix : prefix used to identify the decoy hits"
    
def main(argv):
    if( len(argv) < 2):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        OutputFile = "output.txt"
        decoy_prefix = "random"
        verbose = False
        try:
            opts, args = getopt.getopt(sys.argv[3:], "o:d:hv", ["output=","decoy=", "help", "verbose"])
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
            elif o in ("-d", "--decoy"):
                decoy_prefix = a
            else:
                assert False, "unhandled option"
        
        if(os.path.isfile(argv[0]) and os.path.isfile(argv[1])):
            infile = argv[0]
            infile2 = argv[1]
        else:
            sys.stderr.write("Error: XML file not found\n")
            sys.exit()
                
        parser = etree.XMLParser(ns_clean=False, huge_tree=True)
        
        try:
            tree = etree.parse(infile,parser)
            tree2 = etree.parse(infile2,parser)
            
        except Exception, inst:
            sys.stderr.write("Unexpected error opening %s or %s: %s\n" % (infile,infile2, inst))
            sys.exit()
             
        if(verbose):
            print "Reading " + str(argv[0])   
            print "Reading " + str(argv[1])  
                         
        elems = tree.getroot()
        elems2 = tree2.getroot()
        tppPSMs,tppPeptides = tpp.readIprophetPSMsPeptides(elems,decoy_prefix)
        tppProteins = tpp.readIprophetProteins(elems2,decoy_prefix)

        if(len(tppPSMs) > 0 and len(tppPeptides) > 0 and len(tppProteins) > 0):
            if(verbose):
                print "writing in " + str(OutputFile)   
            utils.writePsms(tppPSMs, "psms_" + OutputFile)
            utils.writePeptides(tppPeptides,"peptides_" + OutputFile)
            utils.writeProteins(tppProteins,"proteins_" + OutputFile)
        else:
            sys.stderr.write("\nThe input files does not contain any information\n")
            sys.exit()
        
        if(verbose):
            print "iProphet and Protein Prophet files parsed"
            
if __name__ == "__main__":
    main(sys.argv[1:]) 
    