#! /usr/bin/env python
# @Created by Jose Fernandez
# dec 10th, 2011

import sys
import os
import getopt
import percolatorUtils as per
from lxml import etree
import random
def usage():
    print "Takes as input a percolator output file and generates a list of proteins obtained by two peptide rule"
    print "Usage : twohitProteinEstimator.py pout.xml  [v, --verbose] [h, --help] [c, --cutoff] [q, --qvalue] [o, --outfile]"
    print "--outfile : name for the output file"
    print "--cutoff : score of qvalue threshold"
    print "--qvalue : uses the qvalues insted of the PEPs for the threshold"

def getProteins(peptides):
    proteins_dict = dict()
    protein_set = set()
    ##asuming peptides are unique, then I get the list of proteins that contains at least 2 peptides
    print "computing two peptide rule from " + str(len(peptides)) + " peptides"
    for peptide in peptides:
        for protein in peptide.proteins:
            if(proteins_dict.has_key(protein)):
                protein_set.add(protein)
            else:
                proteins_dict[protein] = ""
    return list(protein_set)

def getProteinsUnique(peptides):
    proteins_dict = dict()
    protein_set = set()
    ##asuming peptides are unique, then I get the list of proteins that contains at least 2 peptides
    print "computing unique(selecting a random protein when shared peptides) two peptide rule from " + str(len(peptides)) + " peptides"
    for peptide in peptides:
        random_protein = 0
        if(len(peptide.proteins) > 1):
            random_protein = random.randint(0,len(peptide.proteins)-1)
        protein = peptide.proteins[random_protein]
        if(proteins_dict.has_key(protein)):
            protein_set.add(protein)
        else:
            proteins_dict[protein] = ""
    return list(protein_set)

def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        outfile = "output.txt"
        fdr = 0.1
        verbose = False
        qvalue = False
        try:
            opts, args = getopt.getopt(sys.argv[2:], "c:o:q", ["cutoff=","output=,qvalue"])
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
            elif o in ("-c", "--cutoff"):
                fdr = a
            elif o in ("-o", "--output"):
                outfile = a
            elif o in ("-q", "--qvalue"):
                qvalue = True
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

        ##READING ELEMENTS##
        peptides = per.getPeptides(elems)
        print "Read " + str(len(peptides)) + " peptides"
        if(qvalue):
            peptides = [peptide for peptide in peptides if peptide.qvalue <= float(fdr)]
        else:
            peptides = [peptide for peptide in peptides if peptide.pep <= float(fdr)]
        print "There are " + str(len(peptides)) + " peptides at " +str(fdr) + " threshold"
        proteins_filter = getProteins(peptides)
        proteins_filter_unique = getProteinsUnique(peptides)
        print "\nWriting results to " + outfile
        outf = open(outfile, "w")
        outf2 = open("unique_" + outfile, "w")
        outf.write("Protein\tProbability\tPeptides\n")
        outf2.write("Protein\tProbability\tPeptides\n")
        
        for protein in proteins_filter:
            outf.write(protein + "\t" + str(1.0) + "\t\n")
        outf.close()
        
        for protein in proteins_filter_unique:
            outf2.write(protein + "\t" + str(1.0) + "\t\n")
        outf2.close()
        
        print "\n%d proteins were written" % (len(proteins_filter))
        print "\n%d proteins(unique shared peptides) were written" % (len(proteins_filter_unique))
        
if __name__ == "__main__":
    main(sys.argv[1:]) 