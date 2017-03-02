#! /usr/bin/env python
# @Created by Jose Fernandez
# dec 10th, 2011
# This scripts generates an input file compatible with Fido from a percolator output file


from lxml import etree
import os
import getopt
import sys
import  generalUtils as utils
from operator import itemgetter, attrgetter

def usage():
    print "this script generates a psm-proteins and unique peptide-proteins graph files for FIDO as well as a file containing all the target and decoy proteins names from a input file from a sqt file"
    print "Usage : generateFidoInput.py <target.sqt> <decoy.sqt> [-o, --output <output.txt>] [-p, --pattern ] [-d, --database <database.txt>] [-m, --max] [-h, --help] [-v, --verbose]"
    print "--output : the name of the file that will contain the PSMs, probabilities and proteins"
    print "--database : the name of the file that will contain the target and decoy protein names"
    print "--pattern : the pattern that identifies the decoy hits"
    print "--max : max number of hits per spectrum"
    
def main(argv):
    if( len(argv) < 2):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        fidoOutputFile = "output.txt"
        fidoDBfile = "db.txt"
        verbose = False
        decoy_prefix = "random"
        m = 1
        try:
            opts, args = getopt.getopt(sys.argv[3:], "o:p:d:hvm:", ["output=", "pattern=", "database=", "help", "verbose", "max="])
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
            elif o in ("-m", "--max"):
                m = int(a)
            else:
                assert False, "unhandled option"
        
        if(os.path.isfile(argv[0]) and os.path.isfile(argv[1])):
            target = argv[0]
            decoy = argv[1]
        else:
            sys.stderr.write("Error: XML file not found\n")
            sys.exit()
             
        if(verbose):
            print "Reading " + str(argv[0])   
            print "Reading " + str(argv[1]) 

        targetPSMs = utils.parseSqt(target,m,False)
        decoyPSMs = utils.parseSqt(decoy,m,True)
        
        if(verbose):
            print "read " + str(len(targetPSMs)) + " target PSMs"
            print "read " + str(len(decoyPSMs)) + " decoy PSMs"
        
        ##run qvality
        tmp_file_target = "qvality_tmp_target.txt"
        tmp_file_decoy = "qvality_tmp_decoy.txt"
        tmp_qvality = "qvality_tmp.txt"
        
        targetPSMs = sorted(targetPSMs,key=attrgetter('score'),reverse=True)
        decoyPSMs = sorted(decoyPSMs,key=attrgetter('score'),reverse=True)
        PSMs = targetPSMs + decoyPSMs
        PSMs = sorted(PSMs,key=attrgetter('score'),reverse=True)
        
        target_file = open(tmp_file_target,"w")
        for target_psm in targetPSMs:
            target_file.write(str(target_psm.score) + "\n")
        target_file.close()
        
        decoy_file = open(tmp_file_decoy,"w")
        for decoy_psm in decoyPSMs:
            decoy_file.write(str(decoy_psm.score) + "\n")
        decoy_file.close()
        
        command = "qvality -d -o " + "\"" + str(tmp_qvality) + "\" "  + str(tmp_file_target) + " " + str(tmp_file_decoy) + " "
        if(verbose):
            print command
        retcode = os.system(command)            
        if(retcode != 0):
            sys.stderr.write("Error: executing qvality\n")
            sys.exit()
         
        if(not verbose):    
            os.remove(tmp_file_decoy)
            os.remove(tmp_file_target)
        
        results_qvality = open(tmp_qvality,"r").readlines()
        results_qvality.pop(0)
        if( len(results_qvality)  != len(PSMs) ):
            sys.stderr.write("the results from qvality contain less rows than the read PSMs\n")
            sys.exit()
        else:
            for i in xrange(0,len(results_qvality)):
                row = results_qvality[i]
                words = row.split()
                score = float(words[0])
                pep = float(words[1])
                qvalue = float(words[2])
                if(score != round(PSMs[i].score,len(str(score).split(".")[1])) ):
                    sys.stderr.write("The score of PSM : " + str(PSMs[i].name) + " is " + str(PSMs[i].score) + " and does not correspond to the score read from qvality that is " +str(score) + "\n")
                else:
                    PSMs[i].pep = pep
                    PSMs[i].qvalue = qvalue
                    
        if(not verbose):
            os.remove(tmp_qvality)
        
        utils.writeFidoInput(PSMs, fidoOutputFile, fidoDBfile, decoy_prefix)
        
        if(verbose):
            print "done"
        
if __name__ == "__main__":
    main(sys.argv[1:]) 
    