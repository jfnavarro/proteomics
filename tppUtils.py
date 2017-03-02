#! /usr/bin/env python
# @Created by Jose Fernandez
'''
Created on Jul 18, 2012

@author: jfn

Some useful functions for processing iProphet and Protein Prophet files

'''
import hit as hit

def readIprophetProteins(elements,decoy_prefix="random"):
    tpp = dict()
    for elem in elements.iter("{http://regis-web.systemsbiology.net/protXML}protein_group"):
        proteins = elem.findall("{http://regis-web.systemsbiology.net/protXML}protein")
        for protein in proteins:
            pname = protein.get("protein_name")
            decoy = pname.find(decoy_prefix) != -1 
            prob = protein.get("probability")
            if(prob is not None):
                pep = 1 - float(prob)
                if(tpp.has_key(pname)):
                    if(tpp[pname].pep > pep):
                        tpp[pname].pep = pep
                    print "ERROR found in Protein Prophet repeated protein " + str(pname)    
                else:
                    tpp[pname] = hit.Protein(decoy,pep,0.0,0.0,0.0,pname,[])

        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0] 
            
    return [v for k,v in tpp.items()]

def readIprophetPSMsPeptides(elemets,decoy_prefix="random"):
    tppPSMs = []
    tppPeptides = dict()
    tppPeptidesDecoy = dict()
    
    for elem in elemets.iter("{http://regis-web.systemsbiology.net/pepXML}spectrum_query"):
        scan = elem.get("spectrum")
        result = elem.find("{http://regis-web.systemsbiology.net/pepXML}search_result")
        #I take only the first search hit
        match = result.find("{http://regis-web.systemsbiology.net/pepXML}search_hit")
        if(match is not None):
            pname = match.get("protein")
            peptidename = match.get("peptide")
            decoy = pname.find(decoy_prefix) != -1
            proteins = match.findall("{http://regis-web.systemsbiology.net/pepXML}alternative_protein")
            proteins = [str(p.get("protein")) for p in proteins]
            proteins.insert(0,pname)
            score = match.findall("{http://regis-web.systemsbiology.net/pepXML}search_score")[1].get("value")
            peps = match.findall("{http://regis-web.systemsbiology.net/pepXML}analysis_result")
        
            if(len(peps)>0):
                if(len(peps) == 1): ##peptide prophet
                    pepww = peps[0].find("{http://regis-web.systemsbiology.net/pepXML}peptideprophet_result").get("probability")
                    pep = 1 - float(pepww)
                    prob = float(pepww)
                    log_prob = 0.0
                else:
                    pepIprophet = peps[1].find("{http://regis-web.systemsbiology.net/pepXML}interprophet_result").get("probability")
                    pep = 1 - float(pepIprophet)  
                    prob = float(pepIprophet)
                    pepIprophet_log = peps[1].find("{http://regis-web.systemsbiology.net/pepXML}interprophet_result").get("probability_log")
                    if(pepIprophet_log is not None):
                        log_prob = float(pepIprophet_log)
                    else:
                        log_prob = 0.0
                
                pep_prophet = peps[0].find("{http://regis-web.systemsbiology.net/pepXML}peptideprophet_result")
                score_summary = pep_prophet.find("{http://regis-web.systemsbiology.net/pepXML}search_score_summary")
                if(score_summary is not None):
                    fval = score_summary[0].get("value")
                else:
                    fval = 0.0

                psm = hit.PSM(decoy,float(pep),0.0,0.0,0.0,str(scan),float(score),str(peptidename),proteins)
                peptide = hit.Peptide(decoy,float(pep),0.0,0.0,0.0,str(peptidename),float(score),[str(scan)],proteins)
                tppPSMs.append(psm)
                
                if(not decoy):
                    if(tppPeptides.has_key(peptidename)):
                        if(tppPeptides[peptidename].pep > pep):
                            tppPeptides[peptidename].score = float(score)
                            tppPeptides[peptidename].pep = pep  
                            tppPeptides[peptidename].psms.append(str(scan))
                    else:
                        tppPeptides[str(peptidename)] = peptide
                else:
                    if(tppPeptidesDecoy.has_key(peptidename)):
                        if(tppPeptidesDecoy[peptidename].pep > pep):
                            tppPeptidesDecoy[peptidename].score = float(score)
                            tppPeptidesDecoy[peptidename].pep = pep 
                            tppPeptidesDecoy[peptidename].psms.append(str(scan))
                    else:
                        tppPeptidesDecoy[str(peptidename)] = peptide
            
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    
    tppPeptides = [v for k,v in tppPeptides.items()] + [v for k,v in tppPeptidesDecoy.items()]
    return tppPSMs,tppPeptides
