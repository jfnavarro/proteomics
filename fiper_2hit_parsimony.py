#! /usr/bin/env python
# @Created by Jose Fernandez
# dec 10th, 2011

import sys
import os
import getopt
import twohitProteinEstimator as twohit
import parsimonyProteinEstimator as parsimony
import percolatorUtils as per
import generalUtils as utils
from lxml import etree
from pylab import *
import shlex, subprocess
import fnmatch
def usage():
    print "Takes as input a percolator output file and generates a plot of the number of proteins identified at different q values for percolator, parsimony and two peptide rule"
    print "Usage : fiper_2hit_parsimony.py pout.xml  [v, --verbose] [h, --help] [o, --outfile] [q, --qvalue] [m, --mayus] [p, --prefix]"
    print "--outfile : name for the output file"
    print "--qvalue : uses the qvalues insted of the PEPs for the threshold"
    print "--mayus : include mayus analysis using the database given as argument"
    print "--prefix : prefix used to identify decoys"
   
def getProteins(peptides):
    proteins = set()
    for peptide in peptides:
        for protein in peptide.proteins:
            proteins.add(protein)
    return proteins


def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        outfile = "output.png"
        verbose = False
        qvalue = False
        db_file = ""
        mayus = False
        prefix = "random"
        try:
            opts, args = getopt.getopt(sys.argv[2:], "o:qm:p:v", ["output=,qvalue,mayus=,prefix="])
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
                outfile = a
            elif o in ("-q", "--qvalue"):
                qvalue = True
            elif o in ("-m", "--mayus"):
                mayus = True
                db_file = a
            elif o in ("-p", "--prefix"):
                prefix = a
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
 
        thresholds = [0.0015,0.0025,0.005,0.0075,0.01,0.015,0.020,0.030,0.050,0.060,0.075,0.090,0.1,0.15,0.2]
        
        proteins_percolator= per.getProteins(elems)
        peptides = per.getPeptides(elems)
        psms = per.getPSMs(elems)
        
        percolator_proteins_target = list()
        percolator_proteins_decoy = list()
        twohit_proteins_target = list()
        twohit_proteins_decoy = list()
        twohit_proteins_target_unique = list()
        twohit_proteins_decoy_unique = list()
        parsimony_proteins_target = list()
        parsimony_proteins_decoy = list()
        
        mayu_proteins_target = list()
        mayu_proteins_decoy = list()
        mayu_proteins_target_1psm = list()
        mayu_proteins_decoy_1psm = list()
        mayu_proteins_target_2psm = list()
        mayu_proteins_decoy_2psm = list()
        mayu_thresholds = list()
        mayu_thresholds_1psm = list()
        mayu_thresholds_2psm = list()
        two_hit_thresholds = list()
        two_hit_thresholds_unique = list()
        parsimony_thresholds = list()
        ##READING ELEMENTS##
        maxtarget = 0
        maxdecoy = 0
        for t in thresholds:
            print "obtaining proteins from percolator at " + str(t)
            if(qvalue):
                ntarget_percolator = len([prot for prot in [prot for prot in proteins_percolator if not prot.isdecoy] if prot.qvalue <= t])
                ndecoy_percolator = len([prot for prot in [prot for prot in proteins_percolator if prot.isdecoy] if prot.qvalue <= t])
            else:
                ntarget_percolator = len([prot for prot in [prot for prot in proteins_percolator if not prot.isdecoy] if prot.pep <= t])
                ndecoy_percolator = len([prot for prot in [prot for prot in proteins_percolator if prot.isdecoy] if prot.pep <= t])
                
            percolator_proteins_target.append(ntarget_percolator)
            percolator_proteins_decoy.append(ndecoy_percolator)
            print "percolator gives " + str(ntarget_percolator) + " target proteins and " + str(ndecoy_percolator) + " decoy proteins"
            if(maxtarget < ntarget_percolator):
                maxtarget = ntarget_percolator
            if(maxdecoy < ndecoy_percolator):
                maxdecoy = ndecoy_percolator
            
            if(qvalue):
                peptides_filtered = [pep for pep in peptides if pep.qvalue <= t]
                psms_filtered = [psm for psm in psms if pep.qvalue <= t]
            else:
                peptides_filtered = [pep for pep in peptides if pep.pep <= t]
                psms_filtered = [psm for psm in psms if pep.pep <= t]
            
            print "doing mayus analysis at " + str(t)
            psm_file = "mayu_tmp_" + str(t) + ".csv"
            
            #utils.writeMayuInput(psms,psm_file) 
            #utils.writeMayuInputUniqueProteins(psms,psm_file) 
            utils.writeMayuInputFromPeptides(peptides,psm_file)
            #utils.writeMayuInputFromPeptidesUniqueProteins(peptides,psm_file)    
            
            #utils.writeMayuInput(psms_filtered,psm_file)
            #utils.writeMayuInputUniqueProteins(psms_filtered,psm_file)  
            #utils.writeMayuInputFromPeptides(peptides_filtered,psm_file)
            #utils.writeMayuInputFromPeptidesUniqueProteins(peptides_filtered,psm_file)
            
            var = " -B " + psm_file + " -C " + str(db_file) + " -E \"" + str(prefix) +"\"" + " -G " + str(t) 
            #var = " -B " + psm_file + " -C " + str(db_file) + " -E \"" + str(prefix) +"\"" + " -G 1.0 "
            if(verbose):
                print var
            retcode = os.system("Mayu.pl " + var)            
            if(retcode != 0):
                sys.stderr.write("Error: executing MAYU\n")
                sys.exit()
            os.remove(psm_file)
            output_mayu = ""
            for file in os.listdir('.'):
                if fnmatch.fnmatch(file, '*_main_1.06.txt'):
                    output_mayu = file
                if fnmatch.fnmatch(file, '*_main_1.06.csv'):
                    os.remove(file)
            if(not os.path.isfile(output_mayu)):
                sys.stderr.write("Error: finding MAYU output\n")
                sys.exit()
            f = open(output_mayu, "r")
            lines = f.readlines()
            last_line = lines[len(lines)-1]
            columns = last_line.split()
            prot_fdr_mayus = float(columns[18])
            num_targets_mayus = int(columns[13])
            num_decoys_mayus = int(columns[14])
            prot_fdr_single = float(columns[23])
            num_targets_single_psm = int(columns[19])
            mum_decoys_single_psm = int(columns[20])
            num_targets_multi_psm = int(columns[24])
            mum_decoys_multi_psm = int(columns[25])
            prot_fdr_multi = float(columns[28])
            os.remove(output_mayu)
            print "computed a " +str(prot_fdr_mayus) + " global FDR, a " + str(prot_fdr_single) + " single FDR and a " + str(prot_fdr_multi) + " multi FDR"
            print "mayu gives  " + str(num_targets_mayus) + " target proteins and " + str(num_decoys_mayus) + " decoy proteins"
            print "mayu gives (1single PSM) " + str(num_targets_single_psm) + " target proteins and " + str(mum_decoys_single_psm) + " decoy proteins"
            print "mayu gives (2min PSM) " + str(num_targets_multi_psm) + " target proteins and " + str(mum_decoys_multi_psm) + " decoy proteins"
            mayu_proteins_target.append(num_targets_mayus)
            mayu_proteins_decoy.append(num_decoys_mayus)
            mayu_thresholds.append(prot_fdr_mayus)
            mayu_proteins_target_1psm.append(num_targets_single_psm)
            mayu_proteins_decoy_1psm.append(mum_decoys_single_psm)
            mayu_thresholds_1psm.append(prot_fdr_single)
            mayu_proteins_target_2psm.append(num_targets_multi_psm)
            mayu_proteins_decoy_2psm.append(mum_decoys_multi_psm)
            mayu_thresholds_2psm.append(prot_fdr_multi)
            f.close     
        
            if(maxtarget < int(num_targets_mayus)):
                maxtarget = int(num_targets_mayus)
            if(maxdecoy < int(num_decoys_mayus)):
                maxdecoy = int(num_decoys_mayus)
                
            print "doing two peptide rule at " + str(t)
            two_rule_tmp = twohit.getProteins(peptides_filtered) 
            ntarget_2peprule = len([prot for prot in two_rule_tmp if prot.find("random") == -1])
            ndecoy_2peprule = len([prot for prot in two_rule_tmp if prot.find("random") != -1])
            twohit_proteins_target.append(ntarget_2peprule)
            twohit_proteins_decoy.append(ndecoy_2peprule)
            
            if(int(num_targets_mayus) > 0 ):
                tmp_target = float(ntarget_2peprule) / float(num_targets_mayus)
            else:
                tmp_target = 0
            
            if(int(num_decoys_mayus) > 0 ):
                tmp_decoy = float(ndecoy_2peprule) / float(num_decoys_mayus)
            else:
                tmp_decoy = 0
            
            if(tmp_target > 0):
                fdr_2hit = float((tmp_decoy / tmp_target) * prot_fdr_mayus)
            else:
                fdr_2hit = 0

            two_hit_thresholds.append(fdr_2hit)
            
            print "two peptide rule gives " + str(ntarget_2peprule) + " target proteins and " + str(ndecoy_2peprule) + " decoy proteins and " +str(fdr_2hit) + " qvalue"
            if(maxtarget < ntarget_2peprule):
                maxtarget = ntarget_2peprule
            if(maxdecoy < ndecoy_2peprule):
                maxdecoy = ndecoy_2peprule
            
            print "doing unique two peptide rule at " + str(t)
            two_rule_tmp_unique = twohit.getProteinsUnique(peptides_filtered) 
            ntarget_2peprule_unique = len([prot for prot in two_rule_tmp_unique if prot.find("random") == -1])
            ndecoy_2peprule_unique = len([prot for prot in two_rule_tmp_unique if prot.find("random") != -1])
            twohit_proteins_target_unique.append(ntarget_2peprule_unique)
            twohit_proteins_decoy_unique.append(ndecoy_2peprule_unique)

            if(num_targets_mayus > 0 ):
                tmp_target = float(ntarget_2peprule_unique) / float(num_targets_mayus)
            else:
                tmp_target = 0
            
            if(num_decoys_mayus > 0 ):
                tmp_decoy = float(ndecoy_2peprule_unique) / float(num_decoys_mayus)
            else:
                tmp_decoy = 0
            
            if(tmp_target > 0):
                fdr_2hit_unique = float((tmp_decoy / tmp_target) * prot_fdr_mayus)
            else:
                fdr_2hit_unique = 0
                
            two_hit_thresholds_unique.append(fdr_2hit_unique)
            
            print "Unique two peptide rule gives " + str(ntarget_2peprule_unique) + " target proteins and " + str(ndecoy_2peprule_unique) + " decoy proteins and " +str(fdr_2hit_unique) + " qvalue" 
            if(maxtarget < ntarget_2peprule_unique):
                maxtarget = ntarget_2peprule_unique
            if(maxdecoy < ndecoy_2peprule_unique):
                maxdecoy = ndecoy_2peprule_unique
                
            print "doing parsimony at " + str(t)
            parsimony_tmp = [prot for prot,pep in parsimony.getParsimony(peptides, t, qvalue)]
            ntarget_parsimony = len([prot for prot in parsimony_tmp if prot.find("random") == -1])
            ndecoy_parsimony = len([prot for prot in parsimony_tmp if prot.find("random") != -1])
            parsimony_proteins_target.append(ntarget_parsimony)
            parsimony_proteins_decoy.append(ndecoy_parsimony)
            
            if(num_targets_mayus > 0 ):
                tmp_target = float(ntarget_parsimony) / float(num_targets_mayus)
            else:
                tmp_target = 0
            
            if(num_decoys_mayus > 0 ):
                tmp_decoy = float(ndecoy_parsimony) / float(num_decoys_mayus)
            else:
                tmp_decoy = 0
            
            if(tmp_target > 0):
                fdr_parsimony = float((tmp_decoy / tmp_target) * prot_fdr_mayus)
            else:
                fdr_parsimony = 0
                
            parsimony_thresholds.append(fdr_parsimony)
            print "parsimony gives " + str(ntarget_parsimony) + " target proteins and " + str(ndecoy_parsimony) + " decoy proteins and " +str(fdr_parsimony) + " qvalue" 
            
            if(maxtarget < ntarget_parsimony):
                maxtarget = ntarget_parsimony
            if(maxdecoy < ndecoy_parsimony):
                maxdecoy = ndecoy_parsimony
                
        maxtarget = float(maxtarget)
        maxdecoy = float(maxdecoy)
        
        mayu_thresholds[len(mayu_thresholds)-1] =  thresholds[len(thresholds)-1]
        mayu_thresholds_2psm[len(mayu_thresholds_2psm)-1] =  thresholds[len(thresholds)-1]
        parsimony_thresholds[len(parsimony_thresholds)-1] =  thresholds[len(thresholds)-1]
        two_hit_thresholds_unique[len(two_hit_thresholds_unique)-1] =  thresholds[len(thresholds)-1]
        
        thresholds[0] = 0
        mayu_thresholds[0] =  thresholds[0]
        mayu_thresholds_2psm[0] =  thresholds[0]
        parsimony_thresholds[0] =  thresholds[0]
        two_hit_thresholds_unique[0] =  thresholds[0]  
         
        clf()   
        plot(thresholds,percolator_proteins_target, lw = '2', label = "Percolator", color = "blue")
        #plot(two_hit_thresholds,twohit_proteins_target, lw = '2', label = "two peptide rule", color = "red")
        plot(mayu_thresholds,mayu_proteins_target, lw = '2', label = "mayu all", color = "orange")
        #plot(mayu_thresholds_1psm,mayu_proteins_target_1psm, lw = '2', label = "mayu single PSM", color = "pink")
        plot(mayu_thresholds_2psm,mayu_proteins_target_2psm, lw = '2', label = "mayu two peptide rule", color = "brown")
        plot(parsimony_thresholds,parsimony_proteins_target, lw = '2', label = "mayu parsimony rule", color = "black")
        plot(two_hit_thresholds_unique,twohit_proteins_target_unique, lw = '2', label = "mayu two peptide rule(unique)", color = "yellow")
        
        v = [thresholds[0], thresholds[len(thresholds)-1], 0, maxtarget - (maxtarget * 0.25)] #AXES [x-min, x-max, y-min, y-max]
        axis(v)
        legend(loc = 'lower right')
        if(qvalue):
            xlabel("q values",fontsize=20) #X-label
        else:
            xlabel("pep values",fontsize=20) #X-label
        ylabel("number of target proteins",fontsize=20) #Y-label   
        savefig("target_" + outfile, format='png')  
 
        clf()   
        plot(thresholds,percolator_proteins_decoy, lw = '2', label = "Percolator", color = "blue")
        #plot(two_hit_thresholds,twohit_proteins_decoy, lw = '2', label = "two peptide rule", color = "red")
        plot(mayu_thresholds,mayu_proteins_decoy, lw = '2', label = "mayu all", color = "orange")
        #plot(mayu_thresholds_1psm,mayu_proteins_decoy_1psm, lw = '2', label = "mayu single PSM", color = "pink")
        plot(mayu_thresholds_2psm,mayu_proteins_decoy_2psm, lw = '2', label = "mayu two peptide rule", color = "brown")
        plot(two_hit_thresholds_unique,twohit_proteins_decoy_unique, lw = '2', label = "mayu two peptide rule(unique)", color = "yellow")
        plot(parsimony_thresholds,parsimony_proteins_decoy, lw = '2', label = "mayu parsimony rule", color = "black") 
        v = [thresholds[0], thresholds[len(thresholds)-1], 0, maxdecoy - (maxdecoy * 0.25)] #AXES [x-min, x-max, y-min, y-max]
        axis(v)
        legend(loc = 'lower right')
        if(qvalue):
            xlabel("q values",fontsize=20) #X-label
        else:
            xlabel("pep values",fontsize=20) #X-label
        ylabel("number of decoy proteins",fontsize=20) #Y-label   
        savefig("decoy_" + outfile, format='png')  
        
        clf()   
        plot(percolator_proteins_decoy,percolator_proteins_target, lw = '2', label = "Percolator", color = "blue")
        #plot(twohit_proteins_decoy,twohit_proteins_target, lw = '2', label = "two peptide rule hit", color = "red")
        plot(mayu_proteins_decoy,mayu_proteins_target, lw = '2', label = "mayu", color = "orange") 
        #plot(mayu_proteins_decoy_1psm,mayu_proteins_target_1psm, lw = '2', label = "mayu single PSM", color = "pink")
        plot(mayu_proteins_decoy_2psm,mayu_proteins_target_2psm, lw = '2', label = "mayu all but single PSM", color = "brown")
        plot(parsimony_proteins_decoy,parsimony_proteins_target, lw = '2', label = "mayu parsimony rule", color = "black") 
        plot(twohit_proteins_decoy_unique,twohit_proteins_target_unique, lw = '2', label = "mayu all but single PSM(unique)", color = "yellow")
        v = [0,maxdecoy - (maxdecoy * 0.5) , 0, maxtarget + (maxtarget * 0.25)] #AXES [x-min, x-max, y-min, y-max]
        axis(v)
        legend(loc = 'lower right')
        xlabel("number of decoy proteins",fontsize=20) #X-label
        ylabel("number of target proteins",fontsize=20) #Y-label   
        savefig("target_vs_decoy_" + outfile, format='png')  
        
#        ratio_percolator = list()
#        ratio_twopeptide = list()
#        ratio_twopeptide_unique = list()
#        ratio_parsimony = list()
#        for i in xrange(len(thresholds)):
#            if(float(percolator_proteins_target[i]) > 0):
#                ratio_percolator.append(float(percolator_proteins_decoy[i]) / float(percolator_proteins_target[i]))
#            else:
#                ratio_percolator.append(float(percolator_proteins_decoy[i]) / float(1))
#            if(float(twohit_proteins_target[i]) > 0):
#                ratio_twopeptide.append(float(twohit_proteins_decoy[i]) / float(twohit_proteins_target[i]))
#            else:
#                ratio_twopeptide.append(float(twohit_proteins_decoy[i]) / float(1))
#            if(float(twohit_proteins_target_unique[i]) > 0):
#                ratio_twopeptide_unique.append(float(twohit_proteins_decoy_unique[i]) / float(twohit_proteins_target_unique[i]))
#            else:
#                ratio_twopeptide_unique.append(float(twohit_proteins_decoy_unique[i]) / float(1))
#            if(float(parsimony_proteins_target[i]) > 0):
#                ratio_parsimony.append(float(parsimony_proteins_decoy[i]) / float(parsimony_proteins_target[i]))
#            else:
#                ratio_parsimony.append(float(parsimony_proteins_decoy[i]) / float(1))
#            
#        clf()   
#        plot(thresholds,ratio_percolator, lw = '2', label = "Percolator", color = "blue")
#        plot(thresholds,ratio_twopeptide, lw = '2', label = "two peptide rule", color = "red")
#        plot(thresholds,ratio_twopeptide_unique, lw = '2', label = "Unique two peptide rule", color = "yellow")
#        plot(thresholds,ratio_parsimony, lw = '2', label = "parsimony rule", color = "black") 
#         
#        v = [thresholds[0], thresholds[len(thresholds)-1], 0, 0.5] #AXES [x-min, x-max, y-min, y-max]
#        axis(v)
#        legend(loc = 'lower right')
#        if(qvalue):
#            xlabel("q values",fontsize=20) #X-label
#        else:
#            xlabel("pep values",fontsize=20) #X-label
#        ylabel("ratio decoy/target proteins",fontsize=20) #Y-label   
#        savefig("ratio_" + outfile, format='png')  
        
if __name__ == "__main__":
    main(sys.argv[1:]) 