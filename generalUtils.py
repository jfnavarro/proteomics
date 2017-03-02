'''
Created on Jul 18, 2012

@author: jfn
'''
from pylab import *
params = {'legend.fontsize': 10,'legend.linewidth': 1}
rcParams.update(params)
from operator import itemgetter, attrgetter
from curses.ascii import isprint
import hit as hit
import matplotlib.pyplot as plt
from numpy import array
from pylab import prctile

def writePsms(psms,filename):
    f = open(filename, "w")
    psms = sorted(psms,key=attrgetter('pep'),reverse=False)  
    f.write("PSM Peptide Probability Proteins\n")
    for psm in psms:
        prob = 1 - psm.pep
        sequence = psm.peptide
        proteins = psm.proteins
        f.write(str(psm.name) + "   " + str(sequence) + "   " + str(prob) + "   " )
        for protein in proteins:
            f.write(str(protein) + "   ")
        f.write("\n")
    f.close()

def writePeptides(peptides,outputfile):
    f = open(outputfile, "w")
    peptides = sorted(peptides,key=attrgetter('pep'),reverse=False)  
    f.write("Peptide Probability Proteins\n")
    for peptide in peptides:
        prob = 1 - peptide.pep
        sequence = peptide.peptide
        proteins = peptide.proteins
        f.write(str(sequence) + "   " + str(prob) + "   ")
        for protein in proteins:
            f.write(str(protein) + "   ")
        f.write("\n")
    f.close()
    
def writeProteins(proteins,outputfile):
    f = open(outputfile, "w")
    proteins = sorted(proteins,key=attrgetter('pep'),reverse=False)
    f.write("Protein Probability Peptides\n")
    for protein in proteins:
        prob = 1 - protein.pep
        sequence = protein.name
        peptides = protein.peptides
        f.write(str(sequence) + "   " + str(prob) + "   ")
        for peptide in peptides:
            f.write(str(peptide) + "   ")
        f.write("\n")
    f.close()  

def printable(inputfile):
    return ''.join([char for char in inputfile if isprint(char)])

def writeFidoInput(peptides,fidoOutputFile,fidoOutputFile2,decoy_prefix="random"):
#    e EEEMPEPK
#    r SW:TRP6_HUMAN
#    r GP:AJ271067_1
#    r GP:AJ271068_1
#    p 0.9849
#    e LLEIIQVR
#    r SW:TRP6_HUMAN
#    r GP:AJ271067_1
#    r GP:AJ271068_1
#    p 0.0
#    e FAFNNKPNLEWNWK
#    r gi|1574458|gb|AAC23247.1|
#    p 0.9750

#{ target1 , target2 , target3 , ... }
#{ decoy1 , decoy2 , decoy3 , ... }
    proteinDecoys = set()
    proteinTargets = set()
    
    f = open(fidoOutputFile, "w")
    for peptide in peptides:
        prots = peptide.proteins
        pepname = peptide.peptide
        pepprob = 1 - peptide.pep
        f.write("e " + pepname + "\n")
        for prot in prots:
            if(prot.find(decoy_prefix) == -1):
                proteinTargets.add(prot)
            else:
                proteinDecoys.add(prot)
            f.write("r " + prot + "\n")
        f.write("p " + str(pepprob) + "\n")     
    f.close()
    
    proteinDecoys = list(proteinDecoys)
    proteinTargets = list(proteinTargets)

    f = open(fidoOutputFile2, "w")
    f.write("{ ")
    for protein in proteinTargets[:-1]:
        f.write(printable(protein) + " , ")
    f.write(proteinTargets[-1])
    f.write(" }\n")
    f.write("{ ")
    
    for protein in proteinDecoys[:-1]:
        f.write(printable(protein) + " , ")
    f.write(proteinDecoys[-1])
    f.write(" }\n")
    
    f.close()
    
def writeMayuInput(psms,filename):
    #20060116_pep_IPG_5.2517.2517.1,LDFMGPK,rev_C44E4.6,4=147.192000,0.0609
    #20060116_pep_IPG_5.2607.2607.3,CQFGTSTSSIK,K04B12.1,1=160.190100,0.0584
    #20060116_pep_IPG_5.2609.2609.3,SEGAGGSMSLKHHLARK,Y59A8B.1a,8=147.192000,0.1133
    #run.scannr.scannr.charge , raw peptide sequence, protein name,  probability
    f = open(filename, "w")
    for psm in psms:
        proteins = psm.proteins
        pepname = str(psm.peptide[:-2][2:])
        pepprob = 1 - psm.pep
        scanitems = psm.name.split("_")
        run = ''.join(scanitems[:-3]) 
        scan = str(run + "." + scanitems[-3] + "." + scanitems[-3] + "." + scanitems[-2])
        for protein in proteins:
            f.write(str(scan) + "," + str(pepname) + "," + str(protein) + ",," + str(pepprob) + "\n" )
    f.close()

def writeMayuInputFromPeptides(peptides,filename):
    f = open(filename, "w")
    for peptide in peptides:
        proteins = peptide.proteins
        psms = peptide.psms
        pepname = peptide.peptide
        pepprob = 1 - peptide.pep
        scanitems = psms[0].split("_")
        run = ''.join(scanitems[:-3]) 
        scan = str(run + "." + scanitems[-3] + "." + scanitems[-3] + "." + scanitems[-2])
        for protein in proteins:
            f.write(str(scan) + "," + str(pepname) + "," + str(protein) + ",," + str(pepprob) + "\n" )
    f.close()
    
def writeMayuInputUniqueProteins(psms,filename):
    #20060116_pep_IPG_5.2517.2517.1,LDFMGPK,rev_C44E4.6,4=147.192000,0.0609
    #20060116_pep_IPG_5.2607.2607.3,CQFGTSTSSIK,K04B12.1,1=160.190100,0.0584
    #20060116_pep_IPG_5.2609.2609.3,SEGAGGSMSLKHHLARK,Y59A8B.1a,8=147.192000,0.1133
    #run.scannr.scannr.charge , raw peptide sequence, protein name,  probability
    f = open(filename, "w")
    for psm in psms:
        proteins = psm.proteins
        pepname = str(psm.peptide[:-2][2:])
        pepprob = 1 - psm.pep
        scanitems = psm.name.split("_")
        run = ''.join(scanitems[:-3]) 
        scan = str(run + "." + scanitems[-3] + "." + scanitems[-3] + "." + scanitems[-2])
        protein = proteins[0]
        f.write(str(scan) + "," + str(pepname) + "," + str(protein) + ",," + str(pepprob) + "\n" )
    f.close()

def writeMayuInputFromPeptidesUniqueProteins(peptides,filename):
    f = open(filename, "w")
    for peptide in peptides:
        proteins = peptide.proteins
        psms = peptide.psms
        pepname = peptide.peptide
        pepprob = 1 - peptide.pep
        scanitems = psms[0].split("_")
        run = ''.join(scanitems[:-3]) 
        scan = str(run + "." + scanitems[-3] + "." + scanitems[-3] + "." + scanitems[-2])
        protein = proteins[0]
        f.write(str(scan) + "," + str(pepname) + "," + str(protein) + ",," + str(pepprob) + "\n" )
    f.close()
       
def readFidoOuput(fidoInputFile,pattern="random"):
#    0.9988 { gi|1574458|gb|AAC23247 }
#    0.6788 { SW:TRP6_HUMAN , GP:AJ271067_1 , GP:AJ271068_1 }  
    f = open(fidoInputFile, "r")
    fidoproteins = dict()
    for line in f.readlines():
        words = line.split()
        prob = float(words[0])
        startpos = line.find("{") + len("{")
        endpos = line.find("}", startpos)
        proteins = line[startpos:endpos].split(",")
        for protein in proteins:
            decoy = True
            if(str(protein).find(pattern) == -1):  
                decoy = False
            if(fidoproteins.has_key(str(protein))):
                print "ERROR: protein " + str(protein) + " found more than once in the input file\n"
                if(fidoproteins[str(protein)].pep > (1 - prob)):
                    fidoproteins[str(protein)].pep = (1 - prob)
            else:
                fidoproteins[str(protein)] = hit.Protein(decoy,1 - prob,0.0,0.0,0.0,str(protein),[])
    return fidoproteins

def readMagnusrProteins(filein,pattern="random"):
#Q0050 0.101304
#Q0140 0.129315
    f = open(filein, "r")
    proteins = dict()
    for line in f.readlines():
        words = line.split()
        pep = 1 - float(words[1])
        protein = str(words[0])
        decoy = True
        if(str(protein).find(pattern) == -1):  
            decoy = False
        if(proteins.has_key(protein)):
            if(proteins[protein].pep > pep):
                proteins[protein].pep = pep
        else:
            proteins[protein] = hit.Protein(decoy,pep,0.0,0.0,0.0,str(protein),[])
    return proteins

def writeMagnusInput(psms,magnusOutputFile,hidden=False,limit = 1000):
#   PSM qvalue score PEP peptide proteins
    psms = sorted(psms,key=attrgetter('pep'),reverse=False)
    psmsdic_targets = dict()
    psmsdic_decoys = dict()
    f = open(magnusOutputFile, "w")
    f.write("PSM\tqvalue\tscore\tPEP\tPeptide\tProteins\n")
    if(hidden):
        psms = [ x for x in psms if not x.isdecoy]
    for psm in psms:
        cleanpeptide = str(psm.peptide[:-2][2:])
        if(psm.isdecoy):
            if(not psmsdic_decoys.has_key(cleanpeptide)):
                psmsdic_decoys[cleanpeptide] = 1
            else:
                psmsdic_decoys[cleanpeptide] += 1
        else:
            if(not psmsdic_targets.has_key(cleanpeptide)):
                psmsdic_targets[cleanpeptide] = 1
            else:
                psmsdic_targets[cleanpeptide] += 1     
        prots = psm.proteins
        name = printable(psm.peptide)
        prob = psm.pep
        score = psm.score
        qvalue = psm.qvalue
        scan = psm.name
        if( (psm.isdecoy and int(psmsdic_decoys[cleanpeptide]) <= int(limit)) or
            (not psm.isdecoy and int(psmsdic_targets[cleanpeptide]) <= int(limit))):
            f.write(str(scan) + "\t" + str(qvalue) + "\t" + str(score) + "\t" + str(prob) + "\t" + str(name) + "\t")
            for prot in prots:
                f.write(str(prot) + " ")
            f.write("\n")
  
    f.close()

def writeMagnusPeptides(peptides,psms,magnusOutputFile,hidden=False):
#   PSM qvalue score PEP peptide proteins
    peptides = sorted(peptides,key=attrgetter('pep'),reverse=False)
    f = open(magnusOutputFile, "w")
    #f.write("Peptide\tqvalue\tscore\tPEP\tPSM\tProteins\n")
    f.write("PSM\tqvalue\tscore\tPEP\tPeptide\tProteins\n")
    if(hidden):
        peptides = [ x for x in peptides if not x.isdecoy]
    for peptide in peptides:
        prots = peptide.proteins
        prob = peptide.pep
        score = peptide.score
        qvalue = peptide.qvalue
        psm = psms[peptide.peptide]
        scan = psm.name
        name = psm.peptide
        #f.write(str(name) + "\t" + str(qvalue) + "\t" + str(score) + "\t" + str(prob) + "\t" + str(scan) + "\t")
        f.write(str(scan) + "\t" + str(qvalue) + "\t" + str(score) + "\t" + str(prob) + "\t" + str(name) + "\t")
        for prot in prots:
            f.write(str(prot) + " ")
        f.write("\n")
  
    f.close()
    
    
def plotHist(elements,names,colors,xlabel_text,ylabel_text,filename,limit=0.1):
    clf()
    bound = 1
    for x in xrange(len(elements)):
        X = cumulate(elements[x])
        plot(X[1][0:1000], X[0], '-', lw = '2', label = names[x], color = colors[x])           
        newbound = len([y for y in elements[x] if y <= limit])
        bound = max(bound,newbound)
    v = [0, limit, 0, bound + bound/10] #AXES [x-min, x-max, y-min, y-max]
    axis(v)
    legend(loc = 'lower right')
    xlabel(xlabel_text,fontsize=20) #X-label
    ylabel(ylabel_text,fontsize=20) #Y-label   
    savefig(filename, format='png')

def plotCorrelation(elements,elements2,names,colors,xlabel_text,ylabel_text,filename,limit=0.1):
    clf()
    for x in xrange(len(elements)):
        plot(elements[x],elements2[x], '-', lw = '2', label = names[x], color = colors[x])
    x = y = [0,limit]
    plot(x,y)
    v = [0, limit, 0, limit] #AXES [x-min, x-max, y-min, y-max]
    axis(v)
    legend(loc = 'lower right')
    xlabel(xlabel_text,fontsize=20) #X-label
    ylabel(ylabel_text,fontsize=20) #Y-label   
    savefig(filename, format='png')   
    
def cumulate(x):
    figure()
    x = hist(x, bins = 1000, cumulative = True)
    cla()
    close()
    return x

def parseSqt(filename, hitsPerPSM = 1, isDecoy = False):

#H    SQTGenerator Crux
#H    SQTGeneratorVersion 1.0
#H    Comment Crux was written by...
#H    Comment ref...
#H    StartTime    Tue Dec  2 15:22:50 2008
#H    EndTime                               
#H    Database    test_crux_index/test-binary-fasta
#H    DBSeqLength    ?
#H    DBLocusCount    4
#H    PrecursorMasses    average
#H    FragmentMasses    mono
#H    Alg-PreMasTol    3.0
#H    Alg-FragMassTol    0.50
#H    Alg-XCorrMode    0
#H    Comment    preliminary algorithm sp
#H    Comment    final algorithm xcorr
#H    StaticMod    C=160.139
#H    Alg-DisplayTop    5
#H    EnzymeSpec    tryptic
#S    45894    45894    2    1    maccoss007    2038.59    9199.5    147.1    153628
#M      1     27    2040.244    0.0000     1.5881     245.6     11    34        V.YKCAADKQDATVVELTNL.T    U
#L    YCR102C
#M      2     68    2038.265    0.0116     1.5698     208.4     11    36        S.TQSGIVAEQALLHSLNENL.S    U
#L    YGR080W
#M      3     34    2039.247    0.1582     1.3369     239.3     11    36        I.NEKTSPALVIPTPDAENEI.S    U
#L    YLR035C
#M      4    322    2040.365    0.1699     1.3183     160.0      9    36        I.LKESKSVQPGKAIPDIIES.P    U
#L    YJL126W
#M      5     74    2039.453    0.2288     1.2248     203.6     10    32        D.MISVDLKTPLVIFKCHH.G    U
#L    YAL002W
#M     65      1    2041.246    0.4179     0.9245     370.2     13    32        S.CCGLSLPGLHDLLRHYE.E    U
#L    YLR403W
#S    45904    45904    3    1    maccoss007    2834.54    10103.7    246.4    152668
#M      1    237    2833.059    0.0000     1.9175     273.1     20    108            N.NSGSDTVDPLAGLNNLRNSIKSAGNGME.N    U
#L    YDR505C
#M      2    223    2834.278    0.1390     1.6510     274.8     18    96            G.HLSRISNIDDILISMRMDAFDSLIG.Y    U
#L    YLR247C
#M      3     52    2835.100    0.1503     1.6292     324.1     20    96      S.KSTTEPIQLNNKHDLHLGQELTEST.V    U
#L    YDR098C-B
#L    YDR365W-B
#L    YER138C
#L    YGR027W-B
#S     Spectrum     S [low scan] [high scan] [charge] [process time] [server] [experimental mass] [total ion intensity] [lowest Sp] [# of seq. matched]     
#M     Match     M [rank by Xcorr] [rank by Sp] [calculated mass] [DeltaCN] [Xcorr] [Sp] [matched ions] [expected ions] [sequence matched] [validation status U = unknown, Y = yes, N = no, M = Maybe]     yes
#L     Locus     L [locus name] [description if available]

    f = open(filename).readlines()
    psms = list()
    savePSM = False
    saveM = False
    tmp_name = str(filename).split("/")
    psm_name_file = tmp_name[len(tmp_name)-1][:-4]
    psm_name = ""
    hits = 1
    m = dict()
    l = dict()
    for line in f:
        if(line[0] == 'S'):
            if(savePSM):
                savePSM = False
                rank = 1
                for key,element in m.iteritems():
                    xcorr = element
                    peptide = key
                    proteins = l[key]
                    psm_name += "_" + str(rank)
                    psms.append(hit.PSM(isDecoy,0.0,0.0,0.0,0.0,psm_name,xcorr,peptide,proteins))
                    rank += 1
                m = dict()
                l = dict()
            words = line.split()
            scan_low = words[1]
            scan_up = words[2]
            charge = words[3]
            psm_name = psm_name_file + "_" + str(scan_low) + "_" + str(scan_up) + "_" + str(charge) 
            saveM = False
            hits = 1
        if(line[0] == 'M' and hits <= hitsPerPSM):
            words = line.split()
            xcorr = float(words[5])
            peptide = str(words[9])
            m[peptide] = xcorr
            hits +=1
            saveM = True
        elif(line[0] == 'M'):
            saveM = False
        if(line[0] == 'L' and saveM):
            words = line.split()
            if(l.has_key(peptide)):
                l[peptide].append(str(words[1]))
            else:
                l[peptide] = list()
                l[peptide].append(str(words[1]))
            savePSM = True
            
    if(savePSM):
        rank = 1
        for key,element in m.iteritems():
            xcorr = element
            peptide = key
            proteins = l[key]
            psm_name += "_" + str(rank)
            psms.append(hit.PSM(isDecoy,0.0,0.0,0.0,0.0,psm_name,xcorr,peptide,proteins))
            rank += 1

    return psms
        
def importer(path):
    hits = []
    f = open(path).readlines()
    firstLine = f.pop(0) #removes the first line
    firstline = firstLine.split()
    type = 0
    if( firstline[0] == "PSM" ):
        type = 1
    elif( firstline[0] == "Peptide" ):
        type = 2
    elif( firstline[0] == "Protein" ):
        type = 3
    else:
        sys.stderr.write("Error: could not recognize file format\n")
        sys.exit()
        
    for line in f:
        words = line.split()
        if type == 1 :
            prob = float(words[2])
            proteinseq = str(words[3])    
        elif type == 2:
            prob = float(words[1])
            proteinseq = str(words[2])  
        elif type == 3:
            if(len(words) > 1):
                prob = float(words[1])
            else:
                prob = 0.0
            proteinseq = str(words[0])
        hits.append(hit.Single_hit(proteinseq,1 - prob,prob))
    return hits

def plot_pvalues_calibration(pvalues,names,colors,outfile):
    #plot a list of (pvalues)
    #Quantiles
    clf()
    lengths = [len(p) for p in pvalues]
    binNumber = min(lengths)
    rank = array(range(1,binNumber+1)) / float(binNumber)
    quantiles = 100*(rank)

    #Plot
    for i in xrange(len(pvalues)):
        pvalue_list = pvalues[i]
        name = names[i]
        q = prctile(pvalue_list, p=quantiles)
        plt.scatter(rank, q, s=30, c=colors[i], edgecolor=colors[i], marker='o', label=name)
        
    #Diagonal
    linea = lineb = [0.000000000001,10]
    plt.plot(linea, lineb, c="black")
    linec = array(lineb)/2
    lined = array(lineb)*2
    plt.plot(linea, linec, '--', c="grey")
    plt.plot(linea, lined, '--', c="grey")
    
    #Axes
    ax = plt.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(min([1.0/l for l in lengths])*1, 1)
    ax.set_ylim(min([1.0/l for l in lengths])*1, 1)

    plt.legend(loc = 'upper left', scatterpoints=1)
    plt.xlabel(("Normalized rank / ideal $p$ values"),fontsize='x-large')
    plt.ylabel(('Reported $p$ values'),fontsize='x-large') 
    
    plt.savefig(outfile)
    