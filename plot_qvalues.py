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


def getPValues(combined,order=True,decoy_prefix="random"):
    ##combined sorted in best hit first order
    ##score has to be unnormalized
    combined = sorted(combined,key=attrgetter('score'),reverse=order)
    nDecoys = 0
    posSame = 0
    negSame = 0
    p = []
    prevScore = -4711.4711
    for hit in combined:
        if(hit.score != prevScore):
            for i in xrange(1,posSame):
                p.append(float( nDecoys + ( negSame / (posSame + 1)) * (i + 1) ) ) 
            nDecoys += negSame
            negSame = 0
            posSame = 0
            prevScore = hit.score
        if(hit.protein.find(decoy_prefix) == -1):
            posSame += 1
        else:
            negSame += 1
    ##careful with the order here
    for pi in xrange(len(p)):
        if(nDecoys > 0):
            p[pi] = float(p[pi]) / float(nDecoys)
        else:
            p[pi] = 1.0
        
    return sorted(p)
    
def estimatePi0(pvalues,numBoot=100):
    pBoot = []
    lambdas = []
    pi0s = []
    n = len(pvalues)
    numLambda = 100
    maxLambda = 0.5
    pvalues = sorted(pvalues,reverse = False)
    if(n > 0):
        for ix in xrange(numLambda):
            lambdav = ((ix + 1) / float(numLambda)) * maxLambda
            start = lower_bound(pvalues,lambdav)
            Wl = float(n - start) 
            pi0 = float((Wl / n) / (1 - lambdav))
            if (pi0 > 0.0):
                lambdas.append(lambdav)
                pi0s.append(pi0)
            
        if(len(pi0s) == 0):
            print "Error calculating Pi0"
            sys.exit()
        
        minPi0 = min(pi0s)
        mse = []
        for i in xrange(len(pi0s)):
            mse.append(0.0)
        
        for boot in xrange(numBoot):
            pBoot = bootstrap(pvalues)
            n = len(pBoot)
            for ix in xrange(len(lambdas)):
                start = lower_bound(pBoot,lambdas[ix])
                Wl = float(n - start)   #float(distance(start,pBoot[n]))
                pi0Boot = float((Wl / n) / (1 - lambdas[ix]))
                mse[ix] +=  float( ((pi0Boot - minPi0) * (pi0Boot - minPi0)) ) 
    
        minIx = mse.index(min(mse))  #distance(mse.begin(), min(mse), mse);
        pi0 = max(min(pi0s[minIx], 1.0), 0.0);
        return pi0;
    else:
        return -1
            
def bootstrap(pvalues):
    max_size = 1000
    n = len(pvalues);
    num_draw = min(n, max_size);
    out = []
    for ix in xrange(num_draw):
        out.append(float(random.choice(pvalues)))
    return sorted(out,reverse=False)

def lower_bound(liste,element):
    ret = 0
    for ix in xrange(len(liste)):
        if(float(liste[ix]) >= element):
            return ix
    return ret



def estimateQvalues(combined,pi0=1.0,prefix="random",order=False,Ties=False,countdecoys = True,fdr=0.01, tdratio=0.0):
    ##assuming combined sorted in descending order
    nTargets = 0
    nDecoys = 0
    qemp = 0.0
    qest = 0.0
    prevQemp = 0.0
    prevQest = 0.0
    ndecoys_fdr = 0
    ntargets_fdr = 0
    ndecoys_fdr_emp = 0
    ntargets_fdr_emp = 0
    suma = 0.0
    prev_prob = -1000
    qvaluesEst = []
    qvaluesEmp = []
    combined = sorted(combined,key=attrgetter('pep'),reverse=order)

    for hit in combined:
        
        isdecoy = hit.protein.find(prefix) != -1
        
        if(isdecoy):
            nDecoys += 1
        else:
            nTargets += 1
            
        if( (Ties and hit.pep != prev_prob) or (not Ties)):
            if(nTargets > 0):
                qemp = float(nDecoys) / float(nTargets)
                
            if(countdecoys):
                suma = float(suma) + float(hit.pep)
                qest = float(suma) / float(nTargets + nDecoys)
            elif(nTargets > 0):
                if(not isdecoy):
                    suma = float(suma) + float(hit.pep)
                qest = float(suma) / float(nTargets)
                
            if(qemp < prevQemp):
                qemp = prevQemp
            else:
                prevQemp = qemp
            
            if(qest < prevQest):
                qest = prevQest
            else:
                prevQest = qest
        
        prev_prob = hit.pep
        
        if(qest <= fdr):
            if(isdecoy):
                ndecoys_fdr +=1
            else:
                ntargets_fdr +=1
        ##this is not totally correct cos I might modify emp qvalues later      
        if(qemp <= fdr):
            if(isdecoy):
                ndecoys_fdr_emp +=1
            else:
                ntargets_fdr_emp +=1
                
        qvaluesEst.append(qest)
        qvaluesEmp.append(qemp)
    
    print "Estimating qvalues with " + str(nTargets) + " targets and " + str(nDecoys) + " decoys"

    if(nTargets > 0 and nDecoys > 0): 
        if(tdratio == 0.0):   
            factor = float( pi0 * (float(nTargets) / float(nDecoys)) )
        else:
            factor = float( pi0 * (tdratio))
    else:
        factor = 1.0

    for i in xrange(0,len(qvaluesEmp)):
        qvaluesEmp[i] *= factor   
         
    return sorted(qvaluesEmp),sorted(qvaluesEst),ndecoys_fdr,ntargets_fdr,ndecoys_fdr_emp,ntargets_fdr_emp


def estimateFDR_hidden_decoys(target,pi0=1.0,prefix="entrapment",order=False,Ties=False,countdecoys = True,fdr=0.01, tdratio=0.0):
    
    qvaluesEmp = []
    qvaluesEst = []
    nDecoys = 0
    nTargets = 0
    ndecoys_fdr = 0
    ntargets_fdr = 0
    ndecoys_fdr_emp = 0
    ntargets_fdr_emp = 0
    qemp = 0.0
    qest = 0.0
    suma = 0.0
    prevQest = 0.0
    prevQemp = 0.0
    prev_prob = -1000
    target = sorted(target,key=attrgetter('pep'),reverse=order)
    
    for hit in target:
        isdecoy = hit.protein.find(prefix) != -1
        
        if(isdecoy):
            nDecoys += 1
        else:
            nTargets += 1
        
        if( (Ties and hit.pep != prev_prob) or (not Ties)):
            if(nTargets):
                qemp = float(nDecoys * 2) / float(nTargets)
        
            if(countdecoys):
                suma = float(suma) + float(hit.pep)
                qest = float(suma) / float(nTargets + nDecoys)
            elif(nTargets):
                if(not isdecoy):
                    suma = float(suma) + float(hit.pep)
                qest = float(suma) / float(nTargets)
                
            if(qemp < prevQemp):
                qemp = prevQemp
            else:
                prevQemp = qemp
            
            if(qest < prevQest):
                qest = prevQest
            else:
                prevQest = qest
        
        prev_prob = hit.pep
            
        if(qest <= fdr):
            if(hit.protein.find(prefix) != -1):
                ndecoys_fdr +=1
            else:
                ntargets_fdr +=1

        ##this is not totally correct cos I might modify emp qvalues later      
        if(qemp <= fdr):
            if(isdecoy):
                ndecoys_fdr_emp +=1
            else:
                ntargets_fdr_emp +=1
                
        qvaluesEst.append(qest)
        qvaluesEmp.append(qemp)
    
    print "Estimating qvalues with " + str(nTargets) + " targets and " + str(nDecoys) + " decoys"

    if(nTargets > 0 and nDecoys > 0):    
        if(tdratio == 0.0):   
            factor = float( pi0 * (float(nTargets) / float(nDecoys)) )
        else:
            factor = float( pi0 * (tdratio))
    else:
        factor = 1.0

    for i in xrange(0,len(qvaluesEmp)):
        qvaluesEmp[i] *= factor   
         
    return sorted(qvaluesEmp),sorted(qvaluesEst),ndecoys_fdr,ntargets_fdr,ndecoys_fdr_emp,ntargets_fdr_emp
    
def usage():
    print "this scripts plots estimated and empirical q values for multiple set of probabilities and their respective protein name in a metaFile, the metaFile must be like this : filename title"
    print "Usage : plot_qvalues.py <metaFile.txt>  [d, --prefix <prefix>]  [e, --hidden <pattern>] [l, --label <text>] [f, --fdr <number>] [i, --pi0] [r, --tdratio] [t, --ties] [c, --countdecoys] [v, --verbose] [h, --help]"
    print "--hidden : the target database contains hidden decoys with pattern = (q values will be estimated from target and hidden decoys)"
    print "--prefix : the decoy prefix used to identify decoys"
    print "--label : the main label to be used to name the plots"
    print "--fdr : the fdr threshold to be used to cut off"
    print "--pi0 : use pi0 to adjust the empirical qvalues"
    print "--tdratio : target decoy ratio to adjust empirical q values"
    print "--ties : count elements with same probability as a cluster"
    print "--countdecoys : no include decoys when computing estimated q values"
    
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
        ties = False
        countdecoys = True
        label = ""
        fdr = 0.01
        tdratio = 0.0
        pi0 = False
        try:
            opts, args = getopt.getopt(sys.argv[2:], "e:d:vhl:f:i:r:tc", ["hidden=", "prefix=", "verbose", "help", "label=", "fdr=", "pi0=", "tdratio=","ties","countdecoys"])
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
            elif o in ("-f", "--fdr"):
                fdr = float(a)
            elif o in ("-r", "--tdratio"):
                tdratio = float(a)
                print "Using Target Decoy Ratio = " + str(tdratio)
            elif o in ("-i", "--pi0"):
                pi0 = True
                print "Using pi0 to adjust empirical q values"
            elif o in ("-c", "--countdecoys"):
                countdecoys = True
                print "Not counting decoys when computing empirical q values"
            elif o in ("-t", "--ties"):
                ties = True
                print "Counting clusters as one"
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
            if(line != ""):
                prob_file = words[0]
                names.append(words[1])
                if(os.path.isfile(prob_file)):
                    hits.append(util.importer(prob_file))
                    if(verbose):
                        print "reading file " +str(prob_file)
                else:
                    sys.stderr.write("Error: file " + str(prob_file) + " not found\n")
                    sys.exit()  
        
        pi0s = list(xrange(len(hits))) 
        qvaluesEmps = list(xrange(len(hits))) 
        qvaluesEsts = list(xrange(len(hits))) 
        
        if(hidden_decoys):
            print "Estimating qvalues for hidden decoys mode..."
        
        for x in xrange(len(hits)):
            ##ESTIMATING PI0##
            if(hidden_decoys):
                hits[x] = [ele for ele in hits[x] if ele.protein.find(decoy_pattern) == -1]
                if(pi0):
                    pi0s[x] = estimatePi0(getPValues(hits[x],True,hidden_pattern))
                else:
                    pi0s[x] = 1.0
            else:
                if(pi0):
                    pi0s[x] = estimatePi0(getPValues(hits[x],True,decoy_pattern))
                else:
                    pi0s[x] = 1.0
                          
            if(pi0s[x] > 1.0 or pi0s[x] < 0.0):
                pi0s[x] = 1.0
                
            ##ESTIMATING QVALUEs##
            ndecoys = 0
            ntargets = 0
            ndecoysEmp = 0
            ntargetEmp = 0
            if(hidden_decoys):
                qvaluesEmps[x],qvaluesEsts[x],ndecoys,ntargets,ndecoysEmp,ntargetEmp = estimateFDR_hidden_decoys(hits[x],pi0s[x],hidden_pattern,False,ties,countdecoys,fdr,tdratio)
            else:
                qvaluesEmps[x],qvaluesEsts[x],ndecoys,ntargets,ndecoysEmp,ntargetEmp = estimateQvalues(hits[x],pi0s[x],decoy_pattern,False,ties,countdecoys,fdr,tdratio)
            if(pi0):
                print "pi0 file " + str(x+1) + " :" + str(pi0s[x])
            print "elements with estimated qvalue below " + str(fdr) + " in file " + str(x+1) + " :"  + str(ntargets) 
            print "elements with empirical qvalue below " + str(fdr) + " in file " + str(x+1) + " :"  + str(ntargetEmp) 
            print "false positive elements with estimated qvalue below " + str(fdr) + " in file " + str(x+1) + " :"  + str(ndecoys) 
            print "false positive elements with empirical qvalue below " + str(fdr) + " in file " + str(x+1) + " :"  + str(ndecoysEmp) 
            
        if(verbose):
            print "generatings plots"
        
        ##PLOTTING##
        util.plotHist(qvaluesEsts,names,colors,"estimated $q$ value", "#target " + label, label + "_qvalue_estimated.png",fdr)
        util.plotHist(qvaluesEmps,names,colors,"empirical $q$ value", "#target " + label, label + "_qvalue_empirical.png",fdr)
        util.plotCorrelation(qvaluesEmps,qvaluesEsts,names,colors,"empirical $q$ value","estimated $q$ value",label + "_qvalue_estimated_VS_qvalue_empirical_low_range.png",1.0)
        util.plotCorrelation(qvaluesEmps,qvaluesEsts,names,colors,"empirical $q$ value","estimated $q$ value",label + "_qvalue_estimated_VS_qvalue_empirical.png",fdr)        

        if(verbose):
            print "plots generated"
            
if __name__ == "__main__":
    main(sys.argv[1:]) 
    