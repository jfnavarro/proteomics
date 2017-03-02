#!/usr/bin/env python

import sys
import os
import getopt
import fnmatch

global mass_tolerance 
mass_tolerance = "50.0ppm"
global instrument 
instrument = 1
global tda 
tda = 0
global enzyme 
enzyme = 1
global non_enzyme_terminy 
non_enzyme_terminy = 1
global min_length
min_length = 6
global max_length 
max_length = 40
global min_charge 
min_charge = 1
global max_charge 
max_charge = 6
global num_hits_psm 
num_hits_psm = 1
global aa_probs 
aa_probs = 1 
global msgfdb 
msgfdb = "/bubo/home/h22/navarro/bin/bin/MSGFDB.jar"
global local_path 
local_path = os.getcwd()
global project 
project = "s00111-298"
global time 
time = "20:00:00"

def getFilesinFolder(pathToFolder,extension):
    files = list()
    for file_ in os.listdir(pathToFolder):
        if fnmatch.fnmatch(file_, '*.' + extension):
            files.append(file_)
    return files

def make_dir(pathname):
    if not os.path.exists('%s' % pathname):
        cmd = 'mkdir %s' % pathname
        os.system(cmd)
    return pathname

def submit_job(spectra_file,partNum,type_,db,out_dir):
    
    spectra_output = local_path + "/" + spectra_file[:-5] + "_" + str(partNum) + "_output.txt"

    if type_ == "target":
        scriptName = "copyFileTarget%s" % (partNum)
    elif type_ == "decoy":
        scriptName = "copyFileDecoy%s" % (partNum)
        
    job_file = open(("%s/%s" % (out_dir + "/jobs",scriptName)), "w")
    job_file.write("#!/bin/bash\n")
    job_file.write("\n")
    job_file.write("#SBATCH -A %s\n" % (project))
    job_file.write("#SBATCH -p core\n")
    #job_file.write("#SBATCH -N 2\n")
    job_file.write("#SBATCH -n 1\n")
    job_file.write("#SBATCH -t %s\n" % (time))
    job_file.write("#SBATCH -J msgfdb_search\n")
    job_file.write("\n")
    job_file.write("msgfdb=\"%s\"\n" % ("java -Xmx2000M -jar " + msgfdb))
    job_file.write("$msgfdb %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % \
                   ("-s " + spectra_file, " -d " + db, " -t " + mass_tolerance,
                    " -o " + spectra_output, " -tda " + str(tda), " -e " + str(enzyme),
                    " -nnet " + str(non_enzyme_terminy), " -minLength " + str(min_length),
                    " -maxLength " + str(max_length), " -minCharge " + str(min_charge), 
                    " -maxCharge " + str(max_charge), " -n " + str(num_hits_psm), 
                    " -uniformAAProb " + str(aa_probs), " -inst " + str(instrument) ))
    if(type_ == "target"):
        job_file.write("mv %s %s" %(spectra_output ,out_dir + "/target"))
    if(type_ == "decoy"):
        job_file.write("mv %s %s" %(spectra_output ,out_dir + "/decoy"))
    job_file.close()
    os.chmod("%s/%s" % (out_dir + "/jobs", scriptName), 0755)
    command = "sbatch %s/%s" % (out_dir + "/jobs", scriptName)
    os.system(command)

    
def usage():
    print "this script runs msgfdb on a set of spectra files with the dabatabes and parameters given in UPPMAX"
    print "Usage : uppmax_msgfdb.py [-o, --output <pathfolder>] [-i, --input <pathfolder>] [-t, --target <pathname>] [d, --decoy <pathname>] [e, --extension <string>] [-h, --help] [-v, --verbose]"
    print "--output : the path of the directory where the results will be placed"
    print "--input : the path of the directory containing the spectra files"
    print "--target : the path of the target database in fasta format"
    print "--decoy : the path of the decoy database in fasta format"
    print "--extension : the path of the decoy database in fasta format"
    
def main(argv):
    if( len(argv) < 6):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        output_dir = ""
        input_dir = ""
        target_file = ""
        decoy_file = ""
        extension = ""
        verbose = False
        try:
            opts, args = getopt.getopt(sys.argv[1:], "o:i:t:d:e:hv", ["output=","input=","target=","decoy=","extension=" "help", "verbose"])
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
                output_dir = a
            elif o in ("-i", "--input"):
                input_dir = a
            elif o in ("-t", "--target"):
                target_file = a
            elif o in ("-d", "--decoy"):
                decoy_file = a
            elif o in ("-e", "--extension"):
                extension = a
            else:
                assert False, "unhandled option"
        
        if(os.path.isfile(target_file) and os.path.isfile(decoy_file) and os.path.isdir(output_dir) and os.path.isdir(input_dir)):
            if(verbose):
                print "Input parameters processed succesfully"
        else:
            sys.stderr.write("Error: file/s or directoriy/es not found\n")
            sys.exit()
            '''
    -s SpectrumFile (*.mzXML, *.mzML, *.mgf, *.ms2, *.pkl or *_dta.txt)
    -d DatabaseFile (*.fasta or .fa)
    -t ParentMassTolerance (e.g. 2.5Da, 30ppm, or 0.5Da,2.5Da)
       Use comma to set asymmetric values. E.g. "-t 0.5Da,2.5Da" will set 0.5Da to the left (expMass<theoMass) and 2.5Da to the right (expMass>theoMass).
    [-o outputFileName] (Default: stdout)
    [-thread NumOfThreads] (Number of concurrent threads to be executed, Default: Number of available cores)
    [-tda 0/1] (0: don't search decoy database (default), 1: search decoy database to compute FDR)
    [-m FragmentationMethodID] (0: as written in the spectrum or CID if no info (Default), 1: CID, 2: ETD, 3: HCD, 4: Merge spectra from the same precursor)
    [-inst InstrumentID] (0: Low-res LCQ/LTQ (Default for CID and ETD), 1: High-res LTQ (Default for HCD), 2: TOF)
    [-e EnzymeID] (0: No enzyme, 1: Trypsin (Default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 5: Glu-C, 6: Arg-C, 7: Asp-N, 8: aLP, 9: Endogenous peptides)
    [-c13 0/1/2] (Number of allowed C13, Default: 1)
    [-nnet 0/1/2] (Number of allowed non-enzymatic termini, Default: 1)
    [-mod ModificationFileName] (Modification file, Default: standard amino acids with fixed C+57)
    [-minLength MinPepLength] (Minimum peptide length to consider, Default: 6)
    [-maxLength MaxPepLength] (Maximum peptide length to consider, Default: 40)
    [-minCharge MinPrecursorCharge] (Minimum precursor charge to consider if not specified in the spectrum file, Default: 2)
    [-maxCharge MaxPrecursorCharge] (Maximum precursor charge to consider if not specified in the spectrum file, Default: 3)
    [-n NumMatchesPerSpec] (Number of matches per spectrum to be reported, Default: 1)
    [-uniformAAProb 0/1] (0: use amino acid probabilities computed from the input database (default), 1: use probability 0.05 for all amino acids)
    '''

        spectra_files = getFilesinFolder(input_dir,extension)
        ##create output target and decoy
        make_dir(output_dir + "/target")
        make_dir(output_dir + "/decoy")
        ##create jobs directory
        make_dir(output_dir + "/jobs")
        i = 0
        if(len(spectra_files) == 0):
            sys.stderr.write("Error: the folder does not contain any spectra file\n")
            sys.exit()
        for spectra in spectra_files:
            print "Reading spectra : " + input_dir + spectra 
            submit_job(input_dir + spectra,i,"target",target_file,output_dir)
            submit_job(input_dir + spectra,i,"decoy",decoy_file,output_dir)
            i+=1

if __name__ == "__main__":
    main(sys.argv[1:]) 
