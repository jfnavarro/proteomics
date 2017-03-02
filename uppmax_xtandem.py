#!/usr/bin/env python

import sys
import os
import getopt
import fnmatch

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

def createXTandemInput(spectra_file,output_file,database,parameters,type_,input_file_name):
    input_file= "<?xml version=\"1.0\"?>\n" \
                "	<bioml>\n" \
                "		<note>\n" \
                "			Each one of the parameters for x!tandem is entered as a labeled note node.\n" \
                "			Any of the entries in the default_input.xml file can be over-ridden by\n" \
                "			adding a corresponding entry to this file. This file represents a minimum\n" \
                "			input file, with only entries for the default settings, the output file\n" \
                "			and the input spectra file name.\n" \
                "			See the taxonomy.xml file for a description of how FASTA sequence list\n" \
                "			files are linked to a taxon name.\n" \
                "		</note>\n" \
                "		<note type=\"input\" label=\"list path, default parameters\">" + parameters + "</note>\n" \
                "		<note type=\"input\" label=\"list path, taxonomy information\">" + database + "</note>\n" \
                "		<note type=\"input\" label=\"protein, taxon\">" + type_ + "</note>\n" \
                "		<note type=\"input\" label=\"spectrum, path\">" + spectra_file + "</note>\n" \
                "		<note type=\"input\" label=\"output, path\">" + output_file + "</note>\n" \
            	"	</bioml>"
            
    f = open(input_file_name, "w")
    f.write(input_file)
    f.close()

def createXTandemTaxonomy(output_file,db,label):
    taxon_file = "<?xml version=\"1.0\"?>\n" \
                 "	<bioml label=\"x! taxon-to-file matching list\">\n" \
                 "		<taxon label=" + "\"" + label + "\"" + ">\n" \
                 "			<file format=\"peptide\" URL=\"" + db + "\"" + "/>\n" \
                 "		</taxon>\n" \
                 "	</bioml>"

    f = open(output_file, "w")
    f.write(taxon_file)
    f.close()

def submit_job(spectra_file,partNum,taxonomy,parameters,type_,project,time,xtandem,out_dir,input_dir,local_path):
    
    input_file = local_path + "/" + spectra_file[:-5] + "_" + str(partNum) + "_in.xml"
    spectra_output = local_path + "/" + spectra_file[:-5] + "_" + str(partNum) + "_output.xml"
    createXTandemInput(input_dir + spectra_file,spectra_output,taxonomy,parameters,type_,input_file)
    
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
    job_file.write("#SBATCH -J xtandem_search\n")
    job_file.write("\n")
    job_file.write("xtandem=\"%s\"\n" % (xtandem))
    job_file.write("$xtandem %s\n" % (input_file))
    if(type_ == "target"):
        job_file.write("mv %s* %s" %(spectra_output[:-4] ,out_dir + "/target"))
    if(type_ == "decoy"):
        job_file.write("mv %s* %s" %(spectra_output[:-4] ,out_dir + "/decoy"))
    job_file.close()
    os.chmod("%s/%s" % (out_dir + "/jobs", scriptName), 0755)
    command = "sbatch %s/%s" % (out_dir + "/jobs", scriptName)
    os.system(command)

def cleanUp():
    os.system("rm *in.xml")
    os.system("rm *taxonomy_target.xml")
    os.system("rm *taxonomy_decoy.xml")
    
def usage():
    print "this script runs x!tandem on a set of spectra files with the dabatabes and parameters given in UPPMAX"
    print "Usage : uppmax_xtandem.py [-o, --output <pathfolder>] [-i, --input <pathfolder>] [-p, --params <pathname>] [-t, --target <pathname>] [d, --decoy <pathname>] [e, --extension <string>] [-h, --help] [-v, --verbose]"
    print "--output : the path of the directory where the results will be placed"
    print "--input : the path of the directory containing the spectra files"
    print "--params : the path of the directory containing the parameter file"
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
        xtandem = "/bubo/home/h22/navarro/bin/bin/tandem"
        local_path = os.getcwd()
        taxonomy_target = local_path + "/taxonomy_target.xml"
        taxonomy_decoy = local_path + "/taxonomy_decoy.xml"
        project="s00111-298"
        time="20:00:00"
        ##create output target and decoy
        make_dir(output_dir + "/target")
        make_dir(output_dir + "/decoy")
        ##create jobs directory
        make_dir(output_dir + "/jobs")
        ##create taxonomy file for target and decoy databases
        createXTandemTaxonomy(taxonomy_target,target_file,"target")
        createXTandemTaxonomy(taxonomy_decoy,decoy_file,"decoy")
        i = 0
        if(len(spectra_files) == 0):
            sys.stderr.write("Error: the folder does not contain any spectra file\n")
            sys.exit()
        for spectra in spectra_files:
            print "Reading spectra : " + spectra 
            submit_job(spectra,i,taxonomy_target,param_file,"target",project,time,xtandem,output_dir,input_dir,local_path)
            submit_job(spectra,i,taxonomy_decoy,param_file,"decoy",project,time,xtandem,output_dir,input_dir,local_path)
            i+=1
        cleanUp()
        
if __name__ == "__main__":
    main(sys.argv[1:]) 
