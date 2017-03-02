'''
Created on Nov 15, 2012

@author: jfn
'''
#!/usr/bin/env python

import os
import random

def submit_range(v, f1, f2, v_vals, f1_val, f2_val, output_dir, job_file_path, 
                 percolator_path, pin_file, time, project): 
    '''Given variable names (v, f1, f2) and values (v_vals, fx_val),
    submit percolator with those values'''
    # Set some parameters
    pout_file = '/tmp/' + ''.join([random.choice(alphabet) for i in range(9)]) + '.xml'
    tab_file = '/dev/null'
    protein_opts = '-A -%s $v -%s %s -%s %s' % (v, f1, f1_val, f2, f2_val)
    v_vals_string = ''.join([str(val) + ' ' for val in v_vals]).strip()
    
    f = open(job_file_path, 'w')
    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write("#SBATCH -A %s\n" % (project))
    f.write("#SBATCH -p core\n")
    f.write("#SBATCH -n 1\n")
    f.write("#SBATCH -t %s\n" % (time))
    f.write("#SBATCH -J fido_%s\n" % (v))
    f.write("\n")
    f.write("date\n")
    f.write("percolator=\"%s\"\n" % (percolator_path))
    f.write("poutXml=\"%s\"\n" % (pout_file))
    f.write("poutTxt=\"%s\"\n" % (tab_file))
    f.write("pin=\"%s\"\n" % (pin_file))
    f.write("mkdir %s\n" % (output_dir))
    f.write("for v in %s; do\n" % (v_vals_string))
    f.write("stderrOut=\"%s/stderr_$v.txt\"\n" % (output_dir))
    f.write(" $percolator -v 3 %s -X $poutXml $pin 2> $stderrOut > $poutTxt\n" % (protein_opts))
    f.write("done\n")
    f.write("date\n")

    f.close()
    os.system("chmod +x %s" % (job_file_path))
    os.system("sbatch %s" % (job_file_path))

    print v


def main():
    '''
    Run Percolator-Fido with varied parameters, to generate 3 plots
    the distribution of MSE for each set of parameters.
    Parameters after grid_search: alpha=0.36, beta=0.01, gamma=0.5
    '''
    global pin_file, percolator_path, time, project, alphabet
    # Some necessary info
    pin_file = '/bubo/home/h14/viktorg/projects/fidointegration/dat/human_hoopmann/pin/112111-Human-2h-01.pin'
    percolator_path = '/bubo/home/h14/viktorg/bin/bin/percolator'
    dat_dir = '/bubo/home/h14/viktorg/s00111-213/dat/fidointegration/fido_parameters/human_hoopmann'
    time = 240  # TODO: Increase with ~2 hours
    project = 's00111-298'
    alphabet = map(chr, range(97, 123))
    # The parameter names, G must be written in capital, because -G is the option for percolator
    parameters = 'abG'
    # These different values are tested for the "fixed" parameters
    small_grid = {}
    small_grid['a'] = [0.1, 0.4, 0.7]
    small_grid['b'] = [0.0, 0.05, 0.5]
    small_grid['G'] = [0.1, 0.5, 0.9]
    # These values are tested for the variable parameter
    variable_values = [i*0.02 for i in range(37, 51)]  # TODO: range(51) is the right way

    # Alternate the variable parameter
    for variable in parameters:
        # Define the "fixed" parameters, but they're not really fixed
        fixed_params = parameters.replace(variable, '')
        fix_param1 = fixed_params[0]
        fix_param2 = fixed_params[1]
        # Try the values in small_grid for the "fixed" parameters
        for fix_param1_value in small_grid[fixed_params[0]]:
            for fix_param2_value in small_grid[fixed_params[1]]:
                # For these parameters, submit percolator job
                output_dir = '%s/%s/%s%s%s%s' % (dat_dir, variable, fix_param1, fix_param1_value, fix_param2, fix_param2_value)
                job_file_path = 'job_scripts/%s_%s%s_%s%s.sh' % (variable, fix_param1, fix_param1_value, fix_param2, fix_param2_value)
                submit_range(variable, fix_param1, fix_param2, variable_values, fix_param1_value, fix_param2_value, output_dir, job_file_path)
                

if __name__ == '__main__':
    main()
    
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