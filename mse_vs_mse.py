#! /usr/bin/env python
# @Created by Jose Fernandez

from lxml import etree
import os
import getopt
import sys
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
    
def usage():
    print "this script generates reads a tab delimited file containing the values of the grid search of Percolator for fido and generates useful plots"
    print "Usage : mse_vs_mse.py <input.txt>  [-h, --help] [-v, --verbose]"
    
def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        verbose = False
        try:
            opts, args = getopt.getopt(sys.argv[2:], "hv", ["help", "verbose"])
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
            else:
                assert False, "unhandled option"
        
        if(os.path.isfile(argv[0])):
            infile = argv[0]
        else:
            sys.stderr.write("Error: file not found\n")
            sys.exit()
                
        if(verbose):
            print "Reading " + str(argv[0])   
        
        alpha_list = list()
        beta_list = list()
        gamma_list = list()
        
        mse_list1 = list() #MSE standard formula
        mse_list2 = list() # sum of area of absolute trapezoids
        mse_list3 = list() # sum of area of trapezoids  
        mse_list4 = list() # sum of area of trapezoids squared
        mse_list5 = list() # MSE1 + MSE2
        
        roc_list = list()
        
        score_list = list() 
        
        best_gamma_mse1 = 0.0
        best_gamma_mse2 = 0.0
        best_gamma_mse3 = 0.0
        best_gamma_mse4 = 0.0
        best_gamma_mse5 = 0.0
        
        best_alpha_mse1 = 0.0
        best_alpha_mse2 = 0.0
        best_alpha_mse3 = 0.0
        best_alpha_mse4 = 0.0
        best_alpha_mse5 = 0.0
        
        best_beta_mse1 = 0.0
        best_beta_mse2 = 0.0
        best_beta_mse3 = 0.0
        best_beta_mse4 = 0.0
        best_beta_mse5 = 0.0
        
        best_score_mse1 = -1000
        best_score_mse2 = -1000
        best_score_mse3 = -1000
        best_score_mse4 = -1000
        best_score_mse5 = -1000
        
        fn = open(infile)
        for line in fn.readlines():
            words = line.split()
            if(line.find("Grid searching") != -1):
                alpha = float(words[3])
                alpha_list.append(alpha)
                beta = float(words[5])
                beta_list.append(beta)
                gamma = float(words[7])
                gamma_list.append(gamma)
            
            if(line.find("ROC") != -1):
                #roc = float(line.split()[7].split(",")[1])
                roc = float(words[7])
                roc_list.append(roc)
            
            if(line.find("MSE") != -1):
                mse1 = abs(float(line.split()[7].split(",")[0]))
                mse2 = abs(float(line.split()[7].split(",")[1]))
                mse3 = abs(float(line.split()[7].split(",")[2]))
                mse4 = abs(float(line.split()[7].split(",")[3]))
                
                mse_list1.append(mse1)
                mse_list2.append(mse2)
                mse_list3.append(mse3)
                mse_list4.append(mse4)
                mse_list5.append(mse1 + mse2)
            
            if(line.find("Objective function") != -1):
                score = float(words[9])
                score_list.append(score)
                
        lamb = 0.15
        for i in xrange(len(gamma_list)):
            f1 = (lamb * roc_list[i]) - ((1 - lamb) * mse_list1[i])  
            f2 = (lamb * roc_list[i]) - ((1 - lamb) * mse_list2[i])  
            f3 = (lamb * roc_list[i]) - ((1 - lamb) * mse_list3[i])  
            f4 = (lamb * roc_list[i]) - ((1 - lamb) * mse_list4[i])
            f5 = (lamb * roc_list[i]) - ((1 - lamb) * mse_list5[i])
            
            if(f1 > best_score_mse1):
                best_score_mse1 = f1
                best_alpha_mse1 = alpha_list[i]
                best_beta_mse1 = beta_list[i]
                best_gamma_mse1 = gamma_list[i]
            if(f2 > best_score_mse2):
                best_score_mse2 = f2
                best_alpha_mse2 = alpha_list[i]
                best_beta_mse2 = beta_list[i]
                best_gamma_mse2 = gamma_list[i]
            if(f3 > best_score_mse3):
                best_score_mse3 = f3
                best_alpha_mse3 = alpha_list[i]
                best_beta_mse3 = beta_list[i]
                best_gamma_mse3 = gamma_list[i]
            if(f4 > best_score_mse4):
                best_score_mse4 = f4
                best_alpha_mse4 = alpha_list[i]
                best_beta_mse4 = beta_list[i]
                best_gamma_mse4 = gamma_list[i]
            if(f5 > best_score_mse5):
                best_score_mse5 = f5
                best_alpha_mse5 = alpha_list[i]
                best_beta_mse5 = beta_list[i]
                best_gamma_mse5 = gamma_list[i]
                
        print "MSE1 = (sum_i |x_i-y_i|^2 ) * 1/N"
        print "MAE2 = sum of the area of the trapezoids (integral of the abs)"
        print "MAE3 = sum of the area of the trapezoids (abs of the integral)"
        print "MSE4 = sum of the squared area of the trapezoids"
        #print "MSE5 = MSE1 + MSE2"
        
        print "The best alpha,beta,gamma and score for MSE1 are : " + str(best_alpha_mse1) + " " + str(best_beta_mse1) + " " + str(best_gamma_mse1) + " " + str(best_score_mse1)
        print "The best alpha,beta,gamma and score for MSE2 are : " + str(best_alpha_mse2) + " " + str(best_beta_mse2) + " " + str(best_gamma_mse2) + " " + str(best_score_mse2)
        print "The best alpha,beta,gamma and score for MSE3 are : " + str(best_alpha_mse3) + " " + str(best_beta_mse3) + " " + str(best_gamma_mse3) + " " + str(best_score_mse3)       
        print "The best alpha,beta,gamma and score for MSE4 are : " + str(best_alpha_mse4) + " " + str(best_beta_mse4) + " " + str(best_gamma_mse4) + " " + str(best_score_mse4)       
        #print "The best alpha,beta,gamma and score for MSE5 are : " + str(best_alpha_mse5) + " " + str(best_beta_mse5) + " " + str(best_gamma_mse5) + " " + str(best_score_mse5)       
        
        print "The standard deviation of MSE1 is " + str(np.std([x for x in mse_list1 if not np.isnan(x)])) 
        print "The standard deviation of MSE2 is " + str(np.std([x for x in mse_list2 if not np.isnan(x)])) 
        print "The standard deviation of MSE3 is " + str(np.std([x for x in mse_list3 if not np.isnan(x)])) 
        print "The standard deviation of MSE4 is " + str(np.std([x for x in mse_list4 if not np.isnan(x)]))
        #print "The standard deviation of MSE5 is " + str(np.std([x for x in mse_list5 if not np.isnan(x)]))
        
        clf()   
        scatter(mse_list2,mse_list1, lw = '2', label = "Mse2_vs_Mse1", color = "blue")
        v = [min(mse_list2), max(mse_list2), min(mse_list1), max(mse_list1)] #AXES [x-min, x-max, y-min, y-max]
        axis(v)
        xlabel("MSE2 = sum of the area of the trapezoids(abs)",fontsize=20) #X-label
        ylabel("MSE1 = (sum_i |x_i-y_i|^2 ) * 1/N",fontsize=20) #Y-label   
        legend(loc = 'lower right')
        savefig("mse2_vs_mse1.png", format='png')  
        
        clf()   
        scatter(mse_list2,mse_list3, lw = '2', label = "Mse2_vs_Mse3", color = "blue")
        v = [min(mse_list2), max(mse_list2), min(mse_list3), max(mse_list3)] #AXES [x-min, x-max, y-min, y-max]
        axis(v)
        xlabel("MSE2 = sum of the area of the trapezoids(abs)",fontsize=20) #X-label
        ylabel("MSE3 = sum of the area of the trapezoids",fontsize=20) #Y-label   
        legend(loc = 'lower right')
        savefig("mse2_vs_mse3.png", format='png')  
        
        clf()   
        scatter(mse_list2,mse_list4, lw = '2', label = "Mse2_vs_Mse4", color = "blue")
        v = [min(mse_list2), max(mse_list2), min(mse_list4), max(mse_list4)] #AXES [x-min, x-max, y-min, y-max]
        axis(v)
        xlabel("MSE2 = sum of the area of the trapezoids(abs)",fontsize=20) #X-label
        ylabel("MSE4 = sum of the sqarea of the trapezoids",fontsize=20) #Y-label   
        legend(loc = 'lower right')
        savefig("mse2_vs_mse4.png", format='png')  
        
#        clf()   
#        scatter(mse_list2,mse_list5, lw = '2', label = "Mse2_vs_Mse5", color = "blue")
#        v = [min(mse_list2), max(mse_list2), min(mse_list5), max(mse_list5)] #AXES [x-min, x-max, y-min, y-max]
#        axis(v)
#        xlabel("MSE2 = sum of the area of the trapezoids(abs)",fontsize=20) #X-label
#        ylabel("MSE5 = MSE1 + MSE2",fontsize=20) #Y-label   
#        legend(loc = 'lower right')
#        savefig("mse2_vs_mse5.png", format='png') 
        
        clf()   
        scatter(mse_list1,mse_list4, lw = '2', label = "Mse1_vs_Mse4", color = "blue")
        v = [min(mse_list1), max(mse_list1), min(mse_list4), max(mse_list4)] #AXES [x-min, x-max, y-min, y-max]
        axis(v)
        xlabel("MSE1 = (sum_i |x_i-y_i|^2 ) * 1/N",fontsize=20) #X-label
        ylabel("MSE4 = sum of the sq area of the trapezoids",fontsize=20) #Y-label   
        legend(loc = 'lower right')
        savefig("mse1_vs_mse4.png", format='png') 
        
        if(verbose):
            print "done"   
        
if __name__ == "__main__":
    main(sys.argv[1:]) 