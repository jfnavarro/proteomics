#! /usr/bin/env python
# @Created by Jose Fernandez
# dec 10th, 2011

import sys
import os
from operator import itemgetter, attrgetter
import getopt
import random
import copy
from pygraph.classes.graph import graph
from pygraph.algorithms.accessibility import connected_components
import itertools
import percolatorUtils as per
import generalUtils as utils
from lxml import etree

def GetPeptides(protein_set, protein_dict):
    peptides = itertools.chain(*[protein_dict[p] for p in protein_set])
    return set(peptides)

def getExactSolution(subgraph, protein_dict):
    # Given a subgraph of connected proteins, get the minimum subset that explains all peptides 
    to_cover = GetPeptides(subgraph, protein_dict)
    min_subsets = []
    found = False
    for k in range(1, len(subgraph) + 1):
        # get all the protein subsets of length k 
        subsets = list(itertools.combinations(subgraph, k))
        # check if any of these subsets cover all the peptides 
        for sub in subsets:
            covered = GetPeptides(sub, protein_dict)
            if len(covered) == len(to_cover):
                found = True
                min_subsets.append(sub)
                # we stop here since all the other sets are larger than the current 
        if found:
            break 

    # if more than one minimal subset, choose a random one 
    idx = 0 
    if len(min_subsets) > 1:
        idx = random.randint(0,  len(min_subsets) - 1)
        #print "Warning: %d subsets of size %d. Choosing subset with index %d" % (len(min_subsets), len(min_subsets[0]),  idx)

    return min_subsets[idx]

def getParsimony(peptides,fdr,qvalue):
    peptides_filter,proteins_filter = filterOut(peptides,float(fdr),qvalue)
    parsimonylist = list()
    if(len(peptides_filter) > 0 or len(proteins_filter) > 0):
        backup_dict = copy.deepcopy(proteins_filter)
        parsimonylist = GetExactParsimonyList(peptides_filter,proteins_filter)
        parsimonylist = filter(lambda p: p[0] in parsimonylist, backup_dict.items())
    return parsimonylist

def GetExactParsimonyList(peptides, protein_dict):

    parsimony_list = list()
    explained = set()
    for pr, pe in protein_dict.iteritems():
        has_unique_pe = any(len(peptides[p]) == 1 for p in pe) 
        if has_unique_pe:
            parsimony_list.append(pr)
            explained = explained.union(pe)
    tmp = set(peptides.keys())
    unexp_peptides = tmp.difference(explained) 
    for pr in protein_dict.keys():
        protein_dict[pr] = protein_dict[pr].intersection(unexp_peptides)
    print "There are " + str(len(unexp_peptides)) + " unexplained peptides"
    
    gr = graph()
    for pe in unexp_peptides: 
        proteins = peptides[pe]
        if(len(proteins)  == len([node for node in proteins if not gr.has_node(node)])):
            gr.add_nodes(proteins)
        # add edges betwen all these proteins 
            for i in range(len(proteins) - 1):
                for j in range(i + 1, len(proteins)):
                    gr.add_edge((proteins[i], proteins[j]))
    print "%d unexplained peptides, %d proteins in the graph, with %d connections " % (len(unexp_peptides), len(gr.nodes()), len(gr.edges()))
    
    con_components = connected_components(gr)
    subgraphs_dict = {}
    for pr, comp in con_components.iteritems():
        if comp in subgraphs_dict:
            subgraphs_dict[comp] = subgraphs_dict[comp].union(set([pr]))
        else:
            subgraphs_dict[comp] = set([pr]) 
    subgraphs = subgraphs_dict.values()
    tmp = [len(s) for s in subgraphs]
    print "%d independent subgraphs, with %d to %d proteins" % (len(subgraphs), min(tmp), max(tmp))

    for sub in subgraphs:
        min_subset = getExactSolution(sub, protein_dict)
        parsimony_list += min_subset
    print "%d proteins in the final list" % (len(parsimony_list), )

    print "Done parsimony"
    return parsimony_list

def filterOut(peptides,threshold,qvalue):
    
    if(qvalue):
        print "Filtering only peptides with q value < " + str(threshold)  
    else:
        print "Filtering only peptides with PEP < " + str(threshold) 
        
    proteins_list = {}
    peptides_list = {}
    for peptide in peptides:
        if ( (qvalue and peptide.qvalue < float(threshold)) or (not qvalue and peptide.pep < float(threshold)) ):
            peptides_list[peptide.peptide] = peptide.proteins
            for p in peptide.proteins:
                if proteins_list.has_key(p):
                    proteins_list[p] = proteins_list[p].union(set([peptide.peptide]))
                else:
                    proteins_list[p] = set([peptide.peptide])
  
    print "%d peptides mapping to %d proteins were identified"  % (len(peptides_list), len(proteins_list))  
    return peptides_list, proteins_list

def usage():
    print "Takes as input a percolator output file and generates a list of proteins obtained by parsimony"
    print "Usage : parsimonyProteinEstimator.py pout.xml  [v, --verbose] [h, --help] [c, --cutoff] [q, --qvalue]  [o, --outfile]"
    print "--outfile : name for the output file"
    print "--cutoff : score of qvalue threshold"
    print "--qvalue : uses the qvalues insted of the PEPs for the threshold (PEP by default)"
    
def main(argv):
    if( len(argv) < 2):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        outfile = "output.txt"
        fdr = 0.1
        verbose = False
        qvalue = False
        try:
            opts, args = getopt.getopt(sys.argv[2:], "c:o:q", ["cutoff=","output=","qvalue"])
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
            elif o in ("-q", "--qvalue"):
                qvalue = True
            elif o in ("-o", "--output"):
                outfile = a
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
        
        protein_filter = getParsimony(peptides,fdr,qvalue)
        protein_filter.sort(key = lambda p: len(p[1]), reverse = True)
        
        print "\nWriting results to " + outfile
        outf = open(outfile, "w")
        outf.write("Protein\tProbability\tPeptides\n")
        for protein, peptides in protein_filter:
            outf.write(protein + "\t" + str(1.0) + "\t") 
            for peptide in peptides:
                outf.write(str(peptide) + "\t")
            outf.write("\n")
        outf.close()
        print "\n%d proteins were written" % (len(protein_filter))
        
if __name__ == "__main__":
    main(sys.argv[1:]) 
