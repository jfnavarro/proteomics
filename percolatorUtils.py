#! /usr/bin/env python
# @Created by Jose Fernandez
'''
Created on Jul 18, 2012

@author: jfn

Some useful functions to process percolator files

'''
import hit as hit
        
def getPSMs(elems):     
    percolatorPSMs = []
    for elem in elems.iter("{http://per-colator.com/percolator_out/14}psm"):
        decoy = elem.get("{http://per-colator.com/percolator_out/14}decoy") == "true"    
        score = elem.findtext("{http://per-colator.com/percolator_out/14}svm_score")
        pep = elem.findtext("{http://per-colator.com/percolator_out/14}pep")
        q = elem.findtext("{http://per-colator.com/percolator_out/14}q_value")
        p = elem.findtext("{http://per-colator.com/percolator_out/14}p_value")
        proteins = elem.findall("{http://per-colator.com/percolator_out/14}protein_id")
        scan = elem.get("{http://per-colator.com/percolator_out/14}psm_id")
        name = elem.findall("{http://per-colator.com/percolator_out/14}peptide_seq")[0].get("seq")
        nterm = elem.findall("{http://per-colator.com/percolator_out/14}peptide_seq")[0].get("n")
        cterm = elem.findall("{http://per-colator.com/percolator_out/14}peptide_seq")[0].get("c")
        name = str(nterm) + "." + str(name) + "." + str(cterm)
        psm = hit.PSM(decoy,float(pep), float(q), float(p), 0.0, str(scan), float(score), str(name),[x.text for x in proteins])
        percolatorPSMs.append(psm)
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    return percolatorPSMs

def getPeptides(elems):            
    percolatorPeptides = []
    peptideUniqueTarget = dict()
    peptideUniqueDecoy = dict()
    for elem in elems.iter("{http://per-colator.com/percolator_out/14}peptide"):
        decoy = elem.get("{http://per-colator.com/percolator_out/14}decoy") == "true" 
        score = elem.findtext("{http://per-colator.com/percolator_out/14}svm_score")
        pep = elem.findtext("{http://per-colator.com/percolator_out/14}pep")
        q = elem.findtext("{http://per-colator.com/percolator_out/14}q_value")
        p = elem.findtext("{http://per-colator.com/percolator_out/14}p_value")
        name = elem.get("{http://per-colator.com/percolator_out/14}peptide_id")
        proteins = elem.findall("{http://per-colator.com/percolator_out/14}protein_id")
        psms = elem.find("{http://per-colator.com/percolator_out/14}psm_ids")
        psms = [str(id.text) for id in psms.findall("{http://per-colator.com/percolator_out/14}psm_id")]
        peptide = hit.Peptide(decoy,float(pep),float(q),float(p),0.0,str(name),float(score),psms,[x.text for x in proteins])
        percolatorPeptides.append(peptide)
        
        ##uff this check of the uniqueness of the peptides is really bad programming wise
        if(decoy):
            if(not peptideUniqueDecoy.has_key(str(name)) ):
                peptideUniqueDecoy[str(name)] = peptide
            else:
                print "Repeated decoy peptide : " + name + " Found in Percolator files"
        else:  
            if(not peptideUniqueTarget.has_key(str(name)) ):
                peptideUniqueTarget[str(name)] = peptide
            else:
                print "Repeated target peptide : " + name + " Found in Percolator files"
            
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    return percolatorPeptides
 
def getProteins(elems):   
    percolatorProteins = []
    for elem in elems.iter("{http://per-colator.com/percolator_out/14}protein"):
        decoy = elem.get("{http://per-colator.com/percolator_out/14}decoy") == "true" 
        pep = elem.findtext("{http://per-colator.com/percolator_out/14}pep")
        q = elem.findtext("{http://per-colator.com/percolator_out/14}q_value")
        p = elem.findtext("{http://per-colator.com/percolator_out/14}p_value")
        name = elem.get("{http://per-colator.com/percolator_out/14}protein_id")
        qmp = elem.findtext("{http://per-colator.com/percolator_out/14}q_value_emp")
        peptides = [peptide.get("seq") for peptide in elem.findall("{http://per-colator.com/percolator_out/14}peptide_seq")]
        if(qmp is None):
            qmp = 1.0
        if(p is None):
            p = 1.0
        protein = hit.Protein(decoy,float(pep),float(q),float(p),float(qmp),str(name),peptides)
        percolatorProteins.append(protein)
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    return percolatorProteins

