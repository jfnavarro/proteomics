#! /usr/bin/env python
# @Created by Jose Fernandez
'''
Created on Jun 20, 2012

@author: jfn
'''
import os
from operator import itemgetter, attrgetter
import getopt
import generalUtils as util
import sys
import random

def estimate_pvalues(target, decoy, pattern):
    n = len(decoy)
    combined = target + decoy
    combined.sort(reverse=True) #Higher score is better
    pvalues = []
    r = 0.0
    for psm in combined:
        isDecoy = psm[-1].find(pattern) != -1
        r += isDecoy
        pvalue = (r+1)/(n+1)
        if isDecoy == 0:
            pvalues.append(pvalue)
    return sorted(pvalues)

##TAKEN FROM VIKTOR
def qvalue2pvalue(qvalues, ties):
    '''Converts qvalues to pvalues'''
    '''Accepts list of tuples (qvalue, label), and booleans for ties and order'''
    qvalues.sort(reverse=True)          #Start from highest qvalue
    pi0 = qvalues[0][0]
    total = len(qvalues)
    total_null = pi0*total
    pvalues = []
    prev_q = None
    for i, t in enumerate(qvalues):
        q, label = t
        if ties or q != prev_q:         # If ties is True, treat every qvalue as a new one (regardless of ties)
            accepted = total-i
            accepted_null = q*accepted  # Expected number of nulls over threshold
            p = accepted_null/total_null
            pvalues.append((p, label))
        elif q == prev_q:               # If ties is False, a tie is equal to the last one
            pvalues.append((pvalues[-1][0], label))  #If qvalue is same, pvalue becomes same
        prev_q = q
    pvalues.sort()  # Smallest first
    return pvalues

##TAKEN FROM VIKTOR
def pep2pvalue(peps, ties):
    '''Converts PEPs to pvalues'''
    '''Accepts list of tuples (PEP, label), and booleans for ties and order'''
    peps.sort(reverse=True)                         # Start from highest PEP
    TOTAL_NULL = sum([pep for pep, label in peps])
    accepted_null = TOTAL_NULL
    pvalues = []
    prev_pep = None
    for i, t in enumerate(peps):
        pep, label = t
        if ties or pep != prev_pep:     # If ties is True, treat every PEP as a new one (regardless of ties)
            accepted_null -= pep
            p = accepted_null/TOTAL_NULL
            pvalues.append((p, label))
        elif pep == prev_pep:                   # If ties is False, a tie is equal to the last one
            pvalues.append((pvalues[-1][0], label))     # If qvalue is same, pvalue becomes same
        prev_pep = pep
    pvalues.sort()                                  # Smallest first
    return pvalues


##TAKEN FROM VIKTOR
def remove_zeroes(pvalues):
    nonzeropvalues = [i for i in pvalues if i > 0.0]
    lowest = min(nonzeropvalues)
    num = len(pvalues) - len(nonzeropvalues)
    newlist = nonzeropvalues + num * [lowest]
    return newlist

    
def usage():
    print "this scripts plots pvalues and calibration for multiple set of probabilities and their respective protein name in a metaFile, the metaFile must be like this : filename title"
    print "Usage : plot_pvalues.py <metaFile.txt>  [d, --prefix <prefix>]  [e, --hidden <pattern>] [l, --label <text>] [t, --ties]  [v, --verbose] [h, --help]"
    print "--hidden : the target database contains hidden decoys with pattern = "
    print "--prefix : the decoy prefix used to identify decoys"
    print "--label : the main label to be used to name the plots"
    print "--ties : give ties same values"
    
def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        
        hidden_decoys = False
        hidden_pattern = ""
        decoy_pattern = "random"
        verbose = False
        ties = True
        label = ""
        
        try:
            opts, args = getopt.getopt(sys.argv[2:], "e:d:vhl:t", ["hidden=", "prefix=", "verbose", "help", "label", "ties"])
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
            elif o in ("-e", "--hidden"):
                hidden_decoys = True
                hidden_pattern = a
            elif o in ("-d", "--prefix"):
                decoy_pattern = a
            elif o in ("-l", "--label"):
                label = a
            elif o in ("-t", "--ties"):
                ties = a
            else:
                assert False, "unhandled option"
            
        if(os.path.isfile(argv[0])):
            infile = argv[0]
        else:
            sys.stderr.write("Error: file " + str(argv[0]) + " not found\n")
            sys.exit()
        
        colors = ["red","blue","yellow","black","brown","pink","cyan","darkblue","darkred"]
        hits = list()
        names = list()
        ##READING ELEMENTS##
        for line in open(infile).readlines():
            words = line.split()
            prob_file = words[0]
            names.append(words[1])
            if(verbose):
                print "reading file " +str(prob_file)
            if(os.path.isfile(prob_file)):
                hits.append(util.importer(prob_file))
            else:
                sys.stderr.write("Error: file " + str(prob_file) + " not found\n")
                sys.exit()  

        pvalues = list(xrange(len(hits))) 
        pvalues2 = list(xrange(len(hits)))
        
        for x in xrange(len(hits)):
            ##ESTIMATING PVALUES FROM PEP##
            if(hidden_decoys):
                pvalues[x] = pep2pvalue([(ele.pep,ele.protein) for ele in [ele for ele in hits[x] if ele.protein.find(decoy_pattern) == -1]],ties)
                pvalues[x] = [ele[-2] for ele in pvalues[x] if ele[-1].find(hidden_pattern) != -1]
                pvalues[x] = remove_zeroes(pvalues[x])
                pvalues[x] = sorted(pvalues[x],reverse=False)
                pvalues2[x] = estimate_pvalues([(ele.score,ele.protein) for ele in hits[x] if ele.protein.find(hidden_pattern) != -1],\
                                               [(ele.score,ele.protein) for ele in hits[x] if ele.protein.find(decoy_pattern) != -1],\
                                               decoy_pattern)
            else:
                pvalues[x] = pep2pvalue([(ele.pep,ele.protein) for ele in hits[x]],ties)
                pvalues[x] = [ele[-2] for ele in pvalues[x] if ele[-1].find(decoy_pattern) != -1]
                pvalues[x] = remove_zeroes(pvalues[x])
                pvalues[x] = sorted(pvalues[x],reverse=False)
                pvalues2[x] = estimate_pvalues([(ele.score,ele.protein) for ele in hits[x] if ele.protein.find(decoy_pattern) == -1],\
                                               [(ele.score,ele.protein) for ele in hits[x] if ele.protein.find(decoy_pattern) != -1],\
                                               decoy_pattern)
                      

        if(verbose):
            print "generatings plots"
        ##PLOTTING##
        util.plotHist(pvalues2,names,colors,"estimated $p$ value", "#target " + label, label + "_pvalue_estimated.png")
        util.plot_pvalues_calibration(pvalues,names,colors,label + "_pvalue_calibration_from_pep.png")     
        util.plot_pvalues_calibration(pvalues2,names,colors,label + "_pvalue_calibration.png")    

        if(verbose):
            print "plots generated"
            
if __name__ == "__main__":
    main(sys.argv[1:]) 
    

