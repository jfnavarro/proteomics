#!/usr/bin/python
#$Id: PepHype.py,v 1.20 2006/02/24 08:05:41 gfinney Exp $

#the string below is available as the __doc__ variable
'''
PepHype.py :
 
Motivation: Create different decoy database types for false-discovery estimates.

The general practice in proteomics MS/MS spectral database
(Peng 2003) searching
has been to append a reversed copy of the original database to the
original. We are investigating whether other methods for creating decoy
peptides will provide a better false positive estimate rate. Reversed
databases may overestimate the false positive rate because:
*Palindromic portions of the sequence could be repeated
*Homologous portions may be retained across reversing the database
*The deltas in masses between fragment series between reversed peptides
 are the same - looking into whether this affects correlations

Shuffling will essentially produce proteins with the same amino acid
distribution as the original proteins, but with a very low chance of
the original patterns being present.

A Markov-chain generated seqeuence (Allet N. Proteomics 2004) may generate
a random sequence that still retains ... biological similarity as preserved
by adjacent strings of residues compared to a shuffled set of sequences ...
(need to rewrite)

The program reads in a fasta file (reading from stdin is under progress),
and applies one of: 
*Shuffling each sequence,
*Reversing each sequence,
*Creating a sequence of the same length drawn from the amino acid composition of
the original sequence,
*Generating protein sequences drawn from a Markov chain of specified order (tested
up to 5) generated from the observed distribution of n-mers and the n-mer towards
which they transition in the protein database
'''

__license = '''
Copyright (c) 2004-2005, University of Washington
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.
* Neither the name of the University of Washington nor the names of its
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''




#general utility functions we need
import sys,os,random,getopt,tempfile, re

#IF YOU WOULD LIKE THE SCRIPT TO BE ~5X FASTER TRY INSTALLING
#THE PSYCO LIBRARY ( http://psyco.sourceforge.net )
#--------------------------------------------------------------------------------
try :
    import psyco
    psyco.full()
except :  
     pass
#    sys.stderr.write("Wasn't able to import the psyco optimizer, you should install this from\n")
#    sys.stderr.write("psyco.sourceforge.net. It significantly decreases execution time\n")

#for debugging to make sure seed is the same
#counter for random sequences produced
randseqcount = 0
#flag for the name of files, blast sequences
SIMPLE_NAMES = 1
#flag to check shuffled results against the sequence database
CHECKBACK = 0
#blast e-value threshold
BLAST_EVALUE = 0.001
#check prematurely terminating transitions
STOP_EMPTY_CHAINS = 1
#maximum Markov chain level
MAX_MARKOV_ORDER = 4
#output a transition histogram
TRANSITION_HIST = 0
#chain limit ... number of uniquely transferred
UNIQ_CHAIN_LIMIT = 3


FASTA_LINE_LEN = 80
RAND_PREFIX = "random_seq_"
REVERSE_PREFIX = "random_seq_"

USAGE = """PepHype.py [-t -p (prefix) -s -r -m (order) -d -h -n -c -l (linelen) -m (order) ] inputfile.faa [output.faa]

[-t | --showname] : include the name of the original target protein id in the in the name of the decoy

[-p | --prefix] [string]: prefix to use to name the decoy entries

[-z | --combined] : generates a combined target-decoy file as output

[-s | --shuffle] : shuffle the sequences in inputfile.faa themselves
     [default]

[-r | --reverse ] : decoy sequences are created as the reverse of the input sequence

[-m 1-5 | --markov 1-5 ] : uses an Nth order Markov model constructed from the observed n-mer
                   residue frequencies to output the decoy sequences. 

[-n | --nothing ] : output the original database, although with the formatting
                    standardized (i.e. identical line lengths)

[-d | --distribution] : create sequence database with lengths drawn
from the distribution of sequence lengths and residue frequencies
from the whole database

[ --per_seq_distribution ] : calculate the amino acid distribution in each sequence,
     recreate a sequence of that length drawn from that
     distribution (slower, not tested much)

[-l | --line-length ] [integer] : sets the output FASTA line length to the supplied value

[-c | --concat ] : output consists of the original database with a decoy version
                   appended to it 

[-h | --help ] : print this message, then exit.

[--tryp_pep_shuf]  : maintain tryptic-determinant residues, shuffle each tryptic peptide individually

[--tryp_prot_shuf] : maintain tryptic-determinant residues, shuffle all other residues

[--tryp_rev]       : maintain tryptic-determinant residues, reverse the sequence of each tryptic peptide

inputfile.faa : fasta database file which will be outputted as
                a permutation

[output.faa]  : name of output file to which the permuted sequences will be
                written. if left empty, output goes to standard out
"""

def usage(errmsg=None) :
    if errmsg is not None :
        sys.stderr.write(errmsg + os.linesep)
    sys.stderr.write(__doc__)
    sys.stderr.write(USAGE)
    sys.exit(-1)

def main() :
    #Process options
    opts = None
    argv = sys.argv[1:]
    global DEFAULT_FASTA_LINELEN

    try:
        opts,args = getopt.getopt(argv, "m:tp:zrsdhncl:",
                                  ['markov','showname','prefix','combined','reverse','shuffle','distribution',
                                   'help','nothing','concat','per_seq_distribution','seed',
                                   'transition_hist','tryp_rev','tryp_pep_shuf','tryp_prot_shuf'])
    except getopt.GetoptError, ge :
        usage(str(ge))

    #Setup our input and output filehandles
    #open the first argument as input, optional second argument as output
    #or send output to stderr
    #A limitation in the fasta reading function makes it difficult
    #to read from standard input, but that should be fixed soon.
    try :
	fsa_fname = args[0]
        fsafh = open(fsa_fname,'r')
        if ( len(args) == 2 ) :
            outfh = open(args[1],'w')
        else :
            outfh = sys.stdout
    except IOError,ie :
        sys.stderr.write("Error opening file: %s" % str(ie))
        sys.exit(-1)
    except :
        usage()

    #shuffle by default
    concat = 0
    markov = 0 #order of markov model
    #flags for determining decoy technique
    RAND_PREFIX = "random_seq_"
    REVERSE_PREFIX = "random_seq_"
    SIMPLE_NAMES = 1
    COMBINED = 0
    (shuf,whole_seq_dist,per_seq_dist,reverse,markov_on,nothing,tryp_rev,tryp_pep_shuf,tryp_prot_shuf) = (1,0,0,0,0,0,0,0,0)
    for ( opt, val ) in opts :
        if opt in ('-l','--line-length') :
            try :
                FASTA_LINE_LEN = int(val)
                assert(FASTA_LINE_LEN > 0)
                sys.stderr.write("fasta line length is %d" % FASTA_LINE_LEN)
                #set the default fasta line length of the Sequence class to the argument
                Sequence.default_fastalinelen = FASTA_LINE_LEN
            except ValueError :
                usage("Couldn't understand line-length of %s" % val)
        
        if opt in ('-t','--showname'):
            SIMPLE_NAMES = 0
        if opt in ('-p','--prefix'):
            if(val != ""):
                RAND_PREFIX = val
                REVERSE_PREFIX = val
        if opt in ('-z','--combined'):
            COMBINED = 1
            print "Generating Combined file..\n"
        if opt in ('-c','--concat') :
            concat = 1
        if opt in ('-s','--shuffle') :
            shuf = 1
        if opt in ('-d','--distribution') :
            shuf = 0
            whole_seq_dist = 1
        if opt in ('--tryp_rev',) :
            tryp_rev = 1
            shuf = 0
        if opt in ('--tryp_pep_shuf',) :
            tryp_pep_shuf = 1
            shuf = 0
        if opt in ('--tryp_prot_shuf',) :
            tryp_prot_shuf = 1
            shuf = 0
        if opt in ('--per_seq_distribution',) :
            shuf = 0
            per_seq_dist = 1
        if opt in ('--markov','-m') :
            shuf = 0
            markov_on = 1
            try :
                markov = int(val)
            except :
                usage("failed to parse markov level %s" % val)
            if markov > MAX_MARKOV_ORDER :
                usage("Maximum Markov Chain Level is %d" % MAX_MARKOV_ORDER)
        if opt in ('-r','--reverse') :
            shuf = 0
            reverse = 1
        if opt in ('-n','--nothing') :
            shuf = 0
            nothing = 1
        if opt in ('-h','--help') :
            usage()
            sys.exit(0)
        if opt == '--transition_hist' :
            global TRANSITION_HIST
            TRANSITION_HIST = 1

    if ( (shuf + whole_seq_dist + per_seq_dist + reverse + markov_on + nothing) > 1 ) :
        usage("Cannot invoke a combination of -m, ,-n ,-s, -d, and -r options" + os.linesep)

    #counter for random sequences produced
    randseqcount = 0

    #flag for the name of files, blast sequences


    #spit out our initial input database 1st, if that is the case
    if ( concat ) :
        while 1 :
            seq = fasta_to_seq(fsafh)
            if seq is None :
                break
            outfh.write(seq.to_fasta())
        fsafh.seek(0)

    #whole_seq_dist, markov both read in all the data into a set of seq objects to
    #build their decoy databases
    if ( whole_seq_dist ) :
      sys.stderr.write("Generating frequency distribution from the whole database")
      oldseqs = fasta_to_seqs(fsafh)
      newseqs = generate_seqs_from_db_dist(oldseqs)
      for seq in newseqs :
        outfh.write(seq.to_fasta())

    elif markov :
        #TODO -- rewrite this to be less memory hungry -- i.e.
        #generate the model incrementally from each sequence, then generate
        #the output Markov sequences
        sys.stderr.write("Generating with Markov Model of Order %d\n" % markov)
        oldseqs = fasta_to_seqs(fsafh)
        newseqs = generate_new_db( oldseqs, mm_len = markov )
        for seq in newseqs :
            outfh.write(seq.to_fasta())

    else :

      if shuf :
          sys.stderr.write("Shuffling Sequences (drawing from seq dist without replacement)\n")
      elif reverse :
          sys.stderr.write("Reversing Sequences\n")
      elif per_seq_dist:
          sys.stderr.write("Drawing from seq distribution with replacement\n")
      

    #loop over input file
      seqcnt = 0
      #loop while we can create sequences from the fasta filehandle
      while 1 :
        seq = fasta_to_seq(fsafh)
        if seq is None :
            break
        seqcnt += 1
        seqstr = ''.join(seq.seq)
        newseqlist = []
        if(COMBINED):
            outfh.write(seq.to_fasta())
        tmp_seq = seq.seq
        if shuf :
            tmp_seq = shuffle_list(seq.seq)
        elif per_seq_dist :
            tmp_seq = randomize_single_seqstr(seqstr)
        elif reverse :
            tmp_seq = seq.seq
            tmp_seq.reverse()
        elif tryp_pep_shuf :
            tmp_seq = tryptic_pep_shuffle(seq.seq)
        elif tryp_prot_shuf :
            tmp_seq = tryptic_prot_shuffle(seq.seq)
        elif tryp_rev :
            tmp_seq = tryptic_reverse(seq.seq)
        elif nothing :
            pass
        else :
            raise Exception("Error in parsing arguments")
        seq.seq = tmp_seq

        randseqcount += 1

        #determine whether our output names contain the originating protein ID
        #or just numbers (only relevant for reverse, shuffled sequences)
	if SIMPLE_NAMES :
            if reverse :
                seq.name = REVERSE_PREFIX + str(randseqcount)
                seq.anno = ""
            elif nothing :
                pass
            else :
                seq.name = RAND_PREFIX + str(randseqcount)
                seq.anno = ""
        else :
            if reverse :
                seq.name = REVERSE_PREFIX + seq.name
                seq.anno = ''
            elif nothing :
                pass
            else :
                seq.name = RAND_PREFIX + seq.name
                seq.anno = ''
        outfh.write(seq.to_fasta())

def sort_by_0 (x,y) :
    '''Comparison function for comparing a list'''
    return cmp(y[0],x[0])

#TODO -- recheck shuffling code....
def shuffle_list(l) :
    '''Does an in-place Fisher-Yates shuffle of a sequence
    Described in some detail at
    http://www.nist.gov/dads/HTML/fisherYatesShuffle.html

    exchange each element in an array l of size n,
    starting at l[n-1], going down to 0, let the current index be i
    pick a random index from 0 to i to swap with i.

    Notes on randomization:
    This uses pythons builtin random library, which uses a Wichmann-Hill generator
    According to the documentation:
    "While of much higher quality than the rand() function supplied by most C libraries,
    the theoretical properties are much the same as for a single linear congruential generator
    of large modulus. It is not suitable for all purposes, and is completely unsuitable for
    cryptographic purposes."
    So, it should be fine for shuffling sequence databases
    The seed is automatically created from the system time,
    and the generator has a period of ~ 7e12.
    '''
    #i ranges from the last index to zero, decrementing by 1
    for i in range(len(l)-1,0,-1) :
        #j (swap index) goes from 0 to i

        j = int(random.random() * (i+1))
        #swap indices i,j
        (l[i],l[j]) = (l[j],l[i])
    return l


def do_count_map(sequence, seq_freq_map=None ) :
    '''
    Produces a map containing a count of amino acids in a sequence, ignoring 'X' and '*'
    Can also add onto an existing map.
    '''
    count_map = {}
    #optional prexisting map we could add to
    if seq_freq_map :
        count_map = seq_freq_map
    t = 0
    for s in sequence :
        if s in ('X','*') :
            continue
        if not count_map.has_key(s) :
            count_map[s] = 1
        else :
            count_map[s] += 1
        t += 1
    return (count_map,t)

def do_freq_list(count_map,total) :
    '''
    Creates a list of two-element tuples containing the cumulative proportion
    of a residues frequency and the residue.
    i.e for k residues, element n contains the cumulative frequency of
    the residues from 1 to n - 1 plus the frequency of residue n.
    This allows a residue to be selected proportionally from a
    random number. By searching through the ordered cumulative frequencies
    from 0 to 1, when a random number r is less than the current cumulative
    frequency, the corresponding residue is selected.

    [ ( 0.35 , 'A'),(0.65,'K'),(0.85,'Q'),(1,'Y') ]
    will select Alanine with 35% probability, Lysine with 30%, Glutamate with 20%,
    and Tyrosine with 15%
    '''
    #create an ordered array of two elements, one being the proportion, the other the letter
    freq_elements = []
    for (letter,freq) in count_map.items() :
        fe = [ freq * 1.0 / total , letter ]
        freq_elements.append(fe)
    freq_elements.sort(sort_by_0)
    #convert frequency counts to a 0-to-1 scale
    lh = 0.0
    for fe in freq_elements :
        lh = lh + fe[0]
        fe[0] = lh
    #the line below needs to be worked out a bit :
    #frequently summing up the proportions does not equal one,
    #but I was looking for the tolerance to be within 10**-7
    try :
        assert(abs(freq_elements[-1][0] - 1.0) < 0.0000001 )
    except :
        sys.stderr.write( "last element freq = %f\n" % freq_elements[-1][0])
    return freq_elements

#We pick an amino acid by selecting a random number from 0 to 1,
#checking each elements 'proportion boundary' as we go along.
#this could be done more efficiently using a data structure that provides you
#with an element lt/gt the index, like a red-black tree. i.e. search through
#the tree using the random number, where the indices in the tree are
#the frequency bounds of the amino acids

def randomize_single_seqstr (seqstr) :
    '''helper function for generate sequence above -
    creates a frequency distribution and passes it along'''
    (cmap,t) = do_count_map(seqstr)
    seq_freq_list = do_freq_list(cmap,t)
    l = len(seqstr)
    return (generate_sequence(seq_freq_list,l))

def all_seqs_freqmap ( seqs ) :
    count_map = {}
    lengths = []
    tl = 0
    for s in seqs :
	seqstr = ''.join(s.seq)
	(count_map,t) = do_count_map(seqstr,count_map)
	tl += t
        lengths.append(t)
    return (count_map,lengths,tl)


def draw_seq_length(lengths) :
    i = random.randint(0,len(lengths)-1)
    sys.stderr.write("%d\n" % i)
    return lengths[i]


def generate_seqs_from_db_dist ( seqs, length=0 ) :
    '''generates a set of sequences from the distribution of aa composition and
    lengths in seqs
    '''
    (count_map,lengths,tl) = all_seqs_freqmap(seqs)
    if length <= 0 :
      genlen = 0
      for l in lengths :
        genlen += l
    else :
      genlen = tl

    freq_list = do_freq_list(count_map,tl)  
    newseqs = []
    tn = 0
    seq_cnt = 1
    while (tn < genlen) :
        l = draw_seq_length(lengths)
        tn += l
        res_list = generate_sequence(freq_list,l)
        try :
            new_seq = Sequence("random_seq_" + str(seq_cnt),''.join(res_list))
        except :
            print res_list
        newseqs.append(new_seq)
        seq_cnt += 1
    return newseqs
        

def transition_count_per_seq (seq,n,count_map,starting_map) :
    pass

def transition_count_n (seqs,n) :
    count_map = {}   #contains a map of nmers to maps containing following nmers and their counts
    starting_map = {}  #map and counts for starting nmers
    #s_nmer -- starting end mer
    #t_nmer -- transition end mer
    trans = 0

    for seq in seqs :
        seqlist = seq.seq
        #grab starting nmer (s_nmer)
        s_nmer = tuple(seqlist[0:n])
  	#generate a count of triplets at the start
        if starting_map.has_key(s_nmer) :
            starting_map[s_nmer] += 1
        else :
            starting_map[s_nmer] = 1
        _assumption = '''
        If I transit to the last set of residues....
        OK - I was worried that if I set my transitions to jump to the last set of residues
        then, they might be left with nowhere to go -- but that could happen with any kind of
        transition of this nature.

        What kind of things are reasonable to do if I reach a dead end?

        '''

        #start out with the starting nmer

        #probably memory-inefficient to be doing this with nmers
        
        for i in range(0,len(seqlist)-(n)) :

            trans += 1
            s_nmer = tuple(seqlist[i:i+n])   #starting nmer 
            t_nmer = tuple(seqlist[i+1:i+(n+1)]) #the one which follows by one residue

            #old bug had None in s_nmer, t_nmer
            assert(None not in s_nmer)
            assert(None not in t_nmer)
                        
            if count_map.has_key(s_nmer) :
                trans_map = count_map[s_nmer]
                if trans_map.has_key(t_nmer) :
                    try :
                        trans_map[t_nmer] += 1
                    #what is going on with this KeyError for debugging
                    except KeyError :
                        #to facilitate debugging
                        sys.stderr.write("Reached KeyError")
                        pass
                else :
                    try :
                        trans_map[t_nmer] = 1
                    except KeyError :
                        #to facilitate debugging
                        sys.stderr.write("Reached KeyError")
                        pass
            else :
                count_map[s_nmer] = { t_nmer : 1 }

    #an option here would be to terminate the protein instead if a breaking point was found
    #do empty count is a flag FOR WHAT GF

    #make sure the chain doesn't terminate early by reaching a

    deleted_trans_tuples = 0
    deleted_start_tuples = 0
    

    if STOP_EMPTY_CHAINS :
        first_time = 1
        do_empty_cnt = 1
        while do_empty_cnt:
            if not first_time :
                sys.stderr.write("Repeating search for non-terminal maps")
            do_empty_cnt = 0
            first_time = 0

            for (s_nmer,t_nmer_map) in count_map.items() :
                #for each of the following nmers
                t_nmers = t_nmer_map.keys()
                #first remove anything which is self-referential (exclusively)
                #for each starting nmer, and the map of the following nmers and their counts:
                if len(t_nmers) == 1 :
                    if t_nmers[0] == s_nmer :
                        del(count_map[t_nmer[0]])
                    else :
                        uniq_chain_count = 1
                        s_mer = s_nmer
                        t_mers = t_nmers
                        t_mers_to_del = [ s_mer ]
                        while len(t_mers) == 1 :
                            s_mer = t_mers[0]
                            try :
                                t_mers = count_map[s_mer]
                                if len(t_mers) == 1 :
                                    tmers_to_del.append(s_mer)
                                    uniq_chain_count += 1
                                    continue
                                else :
                                    break
                            except :
                                break
                        if uniq_chain_count >= UNIQ_CHAIN_LIMIT :
                            sys.stderr.write("Unique Chain of >= %d links: %s" %
                                             (UNIQ_CHAIN_LIMIT, str(tmers_to_del)))
                            for s_mer in tmers_to_del :
                                del(count_map[s_mer])
                            
                            
                    
                
                
                for t_nmer in t_nmer_map.keys() :
                    #if no entry for a starting nmer (i.e. not the last one)
                    #then delete it as a transition option for the map containing the nmers we transition to
                    if not (count_map.has_key(t_nmer)) :  
                        del(t_nmer_map[t_nmer])
                        sys.stderr.write("Deleted tuple %s -> %s" % (str(s_nmer),str(t_nmer)))
                        deleted_trans_tuples += 1
                        
                        if(len(t_nmer_map.keys()) == 0) :
                            sys.stderr.write("Starting nmer %s no longer has any transitions leaving it" % str(s_nmer))
                            del(count_map[s_nmer])
                            deleted_start_tuples += 1
                            do_empty_cnt = 1
    
                            #    for (nmer_from, nmer_to_map) in count_map.items() :
                            #        #ok, this is a frequency count
                            #        
                            #        if len(nmer_to_map.items()) == 0 :
                            #            empty_cnt += 1
                            #        filter_out_self_transitions_n ( count_map, nmer_from, n )
    return (count_map,starting_map,trans)


def filter_out_self_transitions_n ( count_map, nmer_from, n ) :
    cntdown = n
    nmer = nmer_from

    nmers_to_delete = []

    #tmap is a map of nmers and counts
    while 1 :
        try :
            tmap = count_map[nmer]
        except :
            #simulate that there are no exits from this state
            tmap == {}
        nmers_to_delete.append(nmer)
        if len(tmap.keys()) == 1 :
            cntdown-=1
            if cntdown == 0 :
                break
            nmer = tmap.keys()[0]
        elif len(tmap.keys()) == 0 :
            if cntdown :
                sys.stderr.write("deleting nmers:\n")
                for nmer in nmers_to_delete :
                    sys.stderr.write(count_map[nmer] + "\n")
                    del(count_map[nmer])
        else :
            break

def calculate_transitions (count_map) :
    transition_map = {}
    for res in count_map.keys() :
        trans_map = count_map[res]
        t = 0.0
        for v in trans_map.values() :
            t += v
        freq_list = do_freq_list(trans_map,t)
        transition_map[res] = freq_list
    return transition_map

def generate_sequence_mm (transition_map,start_dist,length) :
    seqlist = [None] * length
    r = random.random()
    seqlist[0] = 'M'
    for sres in start_dist :
        if r < sres[0] :
            seqlist[1] = sres[1]
            break
    l = 2
    while l < length :
        r = random.random()
        last_res = seqlist[l-1]
        transition_dist = transition_map[last_res]
        for tres in transition_dist :
            if r < tres[0] :
                seqlist[l] = tres[1]
                break
        l += 1
    return seqlist





def generate_sequence_mm_n( transition_map, start_dist, length, n ) :
    seqlist = ['Z'] * length
    r = random.random()
    #draw a starting pair from the start distribution
    last_pair = None
    for s_pair_set in start_dist :
        if r < s_pair_set[0] :
            s_pair = s_pair_set[1]
            try :
              seqlist[0:n] = s_pair
            except Exception,e :
              raise(e)

            last_pair = s_pair
            break
    l = n
    while l < length :
        found = 0
        r = random.random()
        try :
            transition_dist = transition_map[last_pair]
        except :
            #problem is what to do when I reach something which only exists
            #as a terminal group of sequences
            #I would want to make sure this one is not chosen!!
            sys.stderr.write("failed to find transitions for %s\n" % str(last_pair))
            break
#        print transition_map
        for t_pair_set in transition_dist :
            if r < t_pair_set[0] :
 #               print "found t_pair_set!!"
  #              print t_pair_set
                t_pair = t_pair_set[1]
#                print t_pair_set
                last_pair = t_pair
                seqlist[l] = t_pair[-1]
                found = 1
                break
        if found == 0 :
            sys.stderr.write("not found!\n" + str(t_pair_set))
            sys.exit(-1)
            
        l += 1
    return seqlist

def calculate_unique_transition_hist ( transition_map ) :
    leaving_counts = {}
    for ( s_nmer, t_list ) in transition_map.items() :
        t_len = len(t_list)

        if t_len == 1 :
            sys.stderr.write("starting group %s had unique exit %s\n" % (str(s_nmer),str(t_list[0])))
        
        if leaving_counts.has_key(t_len) :
            leaving_counts[t_len] += 1
        else :
            leaving_counts[t_len] = 1
    #convert to a list and return
    hval = max(leaving_counts.keys())
    retl = [0] * hval
    for ( lcount, freq ) in leaving_counts.items() :
        assert(lcount != 0)
        retl[lcount-1] = freq
    return retl




def get_len_dist (seqs) :
    lens = []
    for s in seqs :
        lens.append(s.length)
    return lens

def generate_new_db ( seqs, length = 0, mm_len=1 ) :
    lengths = get_len_dist(seqs)
    if length :
        pass
    else :
        for l in lengths :
            length += l
#    (count_map,starting_map,trans) = trans_func(seqs)

    (count_map, starting_map, trans) = transition_count_n(seqs,mm_len)

    #output transition count map
    
    for k in count_map.keys() :
        assert(None not in k)

    trans_dist_map = calculate_transitions(count_map)

    if TRANSITION_HIST :
        trans_hist = calculate_unique_transition_hist ( trans_dist_map )
        print "transition count histogram"
        for c in trans_hist :
            print c
        print "---"
    
    start_t = 0.0
    for sval in starting_map.values() :
        start_t += sval
    start_dist_map = do_freq_list(starting_map,start_t)

    tl = 0
    seqs = []
    seq_cnt = 1
    while tl < length :
        l = draw_seq_length(lengths)
        tl += l
        res_list = generate_sequence_mm_n(trans_dist_map,start_dist_map,l,mm_len)
        #        print res_list
        try :
            new_seq = Sequence("random_seq_" + str(seq_cnt),''.join(res_list))
        except :
            sys.stderr.write("res_list_error:" + str(res_list) + "\n")
        seqs.append(new_seq)
        seq_cnt += 1
    return seqs
            
      
    

    
	

        

  


#------------------------------------------------------------
# Generic code for handling sequences
#------------------------------------------------------------

class Sequence :
  DNA = 1
  PROTEIN = 2
  types = ( DNA, PROTEIN )
  #CHANGE fastalinelen below to use a different sequence length
  default_fastalinelen = 80
   
  def __init__ (self,name,seq="",anno="",type=DNA,
                stash=None,line_length=None) :
    self.name = name
    self.fastalinelen = line_length or self.default_fastalinelen
    #remove '*' from the end of the seq
    seqlist = list(seq)
    #    if seqlist[-1] == '*' :
    #        seqlist.pop()
    #    if filter(lambda(x) : x == '*',seqlist ) :
    #        raise Exception("Found a '*' character before the end of sequence ID %s" % self.name)
    seqlist = filter(lambda(x) : x not in '*' , seqlist)
    self.seq = seqlist
    self.type = type
    self.anno = anno
    self.length = len(self.seq)
    self.stash = stash
    
  def __getitem__(self,pos) :
    return self.seq[pos]
    
  def _output_seq (self,linelen) :
    outputted = 0
    outstr = ""
    while ( outputted < self.length ) :
      outstr += ''.join(self.seq[outputted:outputted+linelen])
      outstr += os.linesep
      outputted += linelen
    return outstr

  def frac_diff(self,othseq) :
    '''returs the fractional difference between the two sequences'''
    l = self.length
    dsum = 0
    if ( self.length != othseq.length ) :
      raise Exception("Expect the sequence lengths of self and other to be the same")
    for i in range(self.length) :
      if self.seq[i] != othseq.seq[i] :
        dsum += 1
    return dsum * 1.0 / l
  
    
  
  def set_seq(self,sequence) :
    #copies a list or converts a string
    self.seq = list(sequence)
    self.length = len(self.seq)
    
  def to_fasta (self,linelen=None) :
    if linelen is None :
      linelen = self.fastalinelen
    outstr = ""
    outstr += ">" + self.name
    if self.anno :
      outstr += "\t" + self.anno
    outstr += os.linesep
    outstr += self._output_seq(linelen)
    return outstr


def fasta_to_seq(fhndl,type='dna' ) :
  '''Expects the next sequence to begin with ">"
  '''
  line = fhndl.readline()
  line = line.replace("\n","")
  line = line.replace("\r","")
  if len(line) == 0 :
    return None
  try :
    if line[0] != '>' :
      raise Exception("Incorrect line beginning for FASTA seq: %s",line)
  except Exception, e :
    raise(e)
    
  lastpos = fhndl.tell()
  sequence = ''
  annotation = ''

  #it would be better to split the ID by whitespace rather than tab
  #use a regex here to split on whitespace

  try :
      assert(line[0] == '>')
  except Exception, e :
      print >> sys.stderr, "expected > at line beginning", line
      raise(e)
  fields = re.split("\s+",line)
  id = fields[0][1:]
  annotation = line[len(fields[0]):]
  if annotation is not None :
      annotation = annotation.strip()

  if type == 'dna' :
      seqtype = Sequence.DNA
  elif type == 'protein' :
      seqtype = Sequence.PROTEIN
  
  while 1 :
    line = fhndl.readline()
    if line is None or line == '' :
      #return finished sequence
      return Sequence(id,sequence,annotation)
    elif line[0] == '>' :
      #backtrack, return the sequence
      fhndl.seek(lastpos)
      #return finished sequence
      return Sequence(id,sequence,annotation,type=seqtype)
    
    lastpos = fhndl.tell()
    line = line.replace("\n","")
    line = line.replace("\r","")

    sequence += line
    

def fasta_to_seqs (fhndl) :
    seqs = []
    while 1 :
        seq = fasta_to_seq(fhndl)
        if seq is None :
            break
        seqs.append(seq)
    #sys.stderr.write("Num of seqs: %d") % len(seqs)
    return seqs

def insilico_trypsinized(seq) :
    segments = []
    seg = []
    for i in range(len(seq)) :
        if seq[i] in ('K','R') :
            if i == len(seq)-1 :
                seg.append(seq[i])
            elif seq[i+1] == 'P' :
                seg.append(seq[i])
            else :
            #found first tryptic site
                if len(seg) :
                    segments.append(seg)
                segments.append( [seq[i]] )
                seg = []
        else :
            seg.append(seq[i])
    if len(seg) :
        segments.append(seg)
    segs_len = reduce(lambda x,y : x+len(y) , segments, 0 )
    try :
        assert(segs_len == len(seq))
    except Exception, e :
        segged_seq = []
        for s in segments :
            segged_seq.extend(s)
        print >> sys.stderr , "lens:" , len(seq), len(segged_seq)
        print >> sys.stderr , "original_seq:"
        print >> sys.stderr , "".join(seq)
        print >> sys.stderr , "new_seq:"
        print >> sys.stderr , "".join(segged_seq)
        raise(e)
        
    return segments

def tryptic_prot_shuffle(seq) :
    segments = insilico_trypsinized(seq)
    pool = []
    for s in segments :
        if len(s) > 1 or s[0] not in ('K','R') :
            pool.extend(s)
    pool = shuffle_list(pool)
    pool_idx = 0
    new_seq = []
    for s in segments :
        l = len(s)
        assert(l)
        if l > 1 or s[0] not in ('K','R') :
            #new shuffled sets are the first l items of shuffled pool
            new_seq.extend(pool[pool_idx:pool_idx+l])
            pool_idx += l
        else :
            #just reinsert K,R
            new_seq.append(s[0])
    assert(len(new_seq) == len(seq))
    return new_seq

def tryptic_pep_shuffle(seq) :
    segments = insilico_trypsinized(seq)
    final_seq = []
    for s in segments :
        if len(s) > 1 :
            new_s = shuffle_list(s)
        else :
            new_s = s
        final_seq.extend(new_s)
    assert(len(final_seq) == len(seq))
    return final_seq

def tryptic_reverse(seq) :
    segments = insilico_trypsinized(seq)
    final_seq = []
    for s in segments :
        if len(s) > 1 :
            new_s = shuffle_list(s)
        else :
            new_s = s
        final_seq.extend(new_s)
    return final_seq

#convention for calling the main function
if __name__ == '__main__' :
  main()
