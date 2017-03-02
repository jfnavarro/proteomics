#!/usr/bin/env python

import sys
import os
import getopt
import fnmatch

def usage():
    print "this script runs crux sequest-search on a set of spectra files with the dabatabes and parameters given in UPPMAX"
    print "YOU must use absolute paths,,,easy to fix this but lazy to do it now"
    print "Usage : uppmax_crux.py [-o, --output <pathfolder>] [-i, --input <pathfolder>] [-p, --params <pathname>] [-t, --target <pathname>] [d, --decoy <pathname>] [e, --extension <string>] [-h, --help] [-v, --verbose]"
    print "--output : the path of the directory where the results will be placed"
    print "--input : the path of the directory containing the spectra files"
    print "--params : the path of the directory containing the parameter file"
    print "--target : the path of the target database in fasta format"
    print "--decoy : the path of the decoy database in fasta format"
    print "--extension : the path of the decoy database in fasta format"

def make_dir(pathname):
    if not os.path.exists('%s' % pathname):
        cmd = 'mkdir %s' % pathname
        os.system(cmd)
    return pathname

def getFilesinFolder(pathToFolder,extension):
    files = list()
    for file_ in os.listdir(pathToFolder):
        if fnmatch.fnmatch(file_, '*.' + extension):
            files.append(file_)
    return files
    

def submit_job(spectraName,partNum,param_file,db_file,project,type_search,time,cruxPath,crux_mode,output_dir,input_dir,local_path):
    
    sqtName =     spectraName[:-5] + str(partNum) + ".sqt"
    pepXMLname =  spectraName[:-5] + str(partNum) + ".pep.xml"

    targetDir = output_dir + "target"
    decoyDir = output_dir + "decoy"
    jobFileDir = output_dir + "jobs"
    scripName = ""
    if type_search == "target":
        scriptName = "copyFileTarget%s" % (partNum)
    elif type_search == "decoy":
        scriptName = "copyFileDecoy%s" % (partNum)
    
    job_file = open(("%s/%s" % (jobFileDir,scriptName)), "w")
    job_file.write("#!/bin/bash\n")
    job_file.write("\n")
    job_file.write("#SBATCH -A %s\n" % (project))
    job_file.write("#SBATCH -p core\n")
    #job_file.write("#SBATCH -N 2\n")
    job_file.write("#SBATCH -n 1\n")
    job_file.write("#SBATCH -t %s\n" % (time))
    job_file.write("#SBATCH -J crux_search\n")
    job_file.write("\n")
    job_file.write("crux=\"%s\"\n" % (cruxPath))
    job_file.write("cd $SNIC_TMP/\n")
    job_file.write("cp %s%s .\n" % (input_dir,spectraName))
    job_file.write("cp %s database.fasta\n" % (db_file))
    job_file.write("cp %s params.txt\n" % (param_file))
    job_file.write("$crux create-index --overwrite T --parameter-file params.txt database.fasta database.index\n")
    job_file.write("$crux %s --overwrite T --parameter-file params.txt %s database.index\n" % (crux_mode, spectraName))
    if type_search == "target":
        if crux_mode == 'sequest-search':
            job_file.write("mv crux-output/sequest.target.sqt %s/%s\n" % (targetDir, sqtName))	
        elif crux_mode == 'search-for-matches':
            job_file.write("mv crux-output/search.target.txt %s/%s\n" % (targetDir, sqtName))
        job_file.write("mv crux-output/sequest.target.pep.xml %s/%s\n" % (targetDir, pepXMLname))
    elif type_search == "decoy":
        if crux_mode == 'sequest-search':
            job_file.write("mv crux-output/sequest.target.sqt %s/%s\n" % (decoyDir, sqtName))
        elif crux_mode == 'search-for-matches':
            job_file.write("mv crux-output/search.target.txt %s/%s\n" % (decoyDir, sqtName))
        job_file.write("mv crux-output/sequest.target.pep.xml %s/%s\n" % (decoyDir, pepXMLname))

    job_file.close()
    os.chmod("%s/%s" % (jobFileDir, scriptName), 0755)
    command = "sbatch %s/%s" % (jobFileDir, scriptName)
    os.system(command)

def cleanUp():
    os.system("rm *.mzXML")
    
########## MAIN ###############

def main(argv):
    if( len(argv) < 6):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        output_dir = ""
        input_dir = ""
        param_file = ""
        target_file = ""
        decoy_file = ""
        extension = ""
        verbose = False

        try:
            opts, args = getopt.getopt(sys.argv[1:], "o:i:p:t:d:e:hv", ["output=","input=","params=","target=","decoy=","extension=" "help", "verbose"])
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
            elif o in ("-p", "--params"):
                param_file = a
            elif o in ("-t", "--target"):
                target_file = a
            elif o in ("-d", "--decoy"):
                decoy_file = a
            elif o in ("-e", "--extension"):
                extension = a
            else:
                assert False, "unhandled option"
        
        if(os.path.isfile(param_file) and os.path.isfile(target_file) and os.path.isfile(decoy_file) and os.path.isdir(output_dir) and os.path.isdir(input_dir)):
            if(verbose):
                print "Input parameters processed succesfully"
        else:
            sys.stderr.write("Error: file/s or directoriy/es not found\n")
            sys.exit()
        
        spectra_files = getFilesinFolder(input_dir,extension)
        local_path = os.getcwd()
        project="s00111-298"
        time="20:00:00"
        cruxPath="/bubo/home/h22/navarro/bin/bin/crux"
        searchMode = 'sequest-search'
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
            print "Reading spectra : " + spectra 
            submit_job(spectra,i,param_file,target_file,project,"target",time,cruxPath,searchMode,output_dir,input_dir,local_path)
            submit_job(spectra,i,param_file,decoy_file,project,"decoy",time,cruxPath,searchMode,output_dir,input_dir,local_path)
            i+=1
        #cleanUp()
        
if __name__ == "__main__":
    main(sys.argv[1:]) 