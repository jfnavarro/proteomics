'''
Created on Jul 18, 2012

@author: jfn

'''

class hit(object):
    def __init__(self, isdecoy=False,pep=0.0,q=0.0,empq=0.0,p=0.0):
        self.isdecoy = isdecoy
        self.pep = pep
        self.qvalue = q
        self.pvalue = p
        self.empqvalue = empq
        
class PSM(hit):
    def __init__(self,isdecoy=False,pep=0.0,q=0.0,empq=0.0,p=0.0,pms_id="",score=0.0,sequence="",proteins=[]):
        super(PSM, self).__init__(isdecoy,pep,q,empq,p)
        self.name = pms_id
        self.score = score 
        self.peptide = sequence 
        self.proteins = proteins 
    
class Peptide(hit):
    def __init__(self,isdecoy=False,pep=0.0,q=0.0,empq=0.0,p=0.0,sequence="",score=0.0,psms=[],proteins=[]):
        super(Peptide, self).__init__(isdecoy,pep,q,empq,p)
        self.peptide = sequence
        self.score = score 
        self.psms = psms 
        self.proteins = proteins 
        
class Protein(hit):
    def __init__(self,isdecoy=False,pep=0.0,q=0.0,empq=0.0,p=0.0,name="",peptides=[]):
        super(Protein, self).__init__(isdecoy,pep,q,empq,p)
        self.name = name
        self.peptides = peptides 
        
class Single_hit:
    def __init__(self, protein="",pep=0.0,score=0.0):
        self.protein = protein
        self.pep = pep
        self.score = score