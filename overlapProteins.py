#! /usr/bin/env python
#@ Created by L. Moruz 
#June 29th, 2012

import sys 
import os
import getopt
import generalUtils as util
import pylab
from matplotlib.patches import Circle, Ellipse
from itertools import chain
from collections import Iterable

#--------------------------------------------------------------------
alignment = {'horizontalalignment':'center', 'verticalalignment':'baseline'}


def Get2Venn(sA, lA, sB, lB):
    print "--------------"
    a1 = len(sA.difference(sB))
    print "1 = (%s - %s) = %d" %(lA, lB, a1)
    a2 = len(sA.intersection(sB))
    print "2 = (%s int %s) = %d" %(lA, lB, a2)
    a3 = len(sB.difference(sA))
    print "3 = (%s - %s) = %d" %(lB, lA, a3)
    print "--------------"
    print "Verification: %s = %d = %d" % (lA, len(sA), a1 + a2)
    print "              %s = %d = %d" % (lB, len(sB), a2 + a3)  
    print "--------------"
    print "The common represent: %.0f%% from %s, and %.0f%% from %s" % \
    (float(a2) / len(sA) * 100, lA, float(a2) / len(sB) * 100, lB)



def Get3Venn(sA, lA, sB, lB, sC, lC):
    print "--------------"
    a1 = len(sA.difference(sB.union(sC)))
    print "1 = %s - (%s U %s) = %d" %(lA, lB, lC, a1)
    a2 = len((sA.intersection(sB)).difference(sC))
    print "2 = (%s int %s) -%s = %d" %(lA, lB, lC, a2)
    a3 = len(sA.intersection(sB, sC))
    print "3 = %s int %s int %s = %d " % (lA, lB, lC, a3)
    a4 = len((sA.intersection(sC)).difference(sB)) 
    print "4 = (%s int %s) - %s = %d" % (lA, lC, lB, a4)
    a5 = len(sB.difference(sA.union(sC))) 
    print "5 = %s - (%s U %s) = %d" % (lB, lA, lC, a5)
    a6 = len((sB.intersection(sC)).difference(sA)) 
    print "6 = (%s int %s) - %s = %d" % (lB, lC, lA, a6)
    a7 = len(sC.difference(sA.union(sB))) 
    print "7 = %s - (%s U %s) = %d" % (lC, lA, lB, a7)

    print "--------------"
    print "Verification: %s = %d = %d" % (lA, len(sA), a1 + a2 + a3 + a4)
    print "              %s = %d = %d" % (lB, len(sB), a2 + a3 + a5 + a6)
    print "              %s = %d = %d" % (lC, len(sC), a3 + a4 + a6 + a7)    


def Get4Venn(sA, lA, sB, lB, sC, lC, sD, lD):
    print "--------------"
    a1 = len( sA.difference(sB.union(sC, sD)) )
    print "1 = %s - (%s U %s U %s) = %d" %(lA, lB, lC, lD, a1) 
    a2 = len( (sA.intersection(sB)).difference(sC.union(sD)) )
    print "2 = (%s int %s) - (%s U %s) = %d" %(lA, lB, lC, lD, a2)
    a3 = len( sB.difference(sA.union(sC, sD)) )
    print "3 = %s - (%s U %s U %s) = %d" %(lB, lA, lC, lD, a3)
    a4 = len( (sA.intersection(sC)).difference(sB.union(sD)) )
    print "4 = (%s int %s) - (%s U %s) = %d" %(lA, lC, lB, lD, a4)
    a5 = len( (sA.intersection(sB, sC)).difference(sD) )
    print "5 = (%s int %s int %s) - %s = %d" %(lA, lB, lC, lD, a5)
    a6 = len( (sB.intersection(sC)).difference(sA.union(sD)) )
    print "6 = (%s int %s) - (%s U %s) = %d" %(lB, lC, lA, lD, a6)
    a7 = len( sC.difference(sA.union(sB, sD)) )
    print "7 = %s - (%s U %s U %s) = %d" %(lC, lA, lB, lD, a7)
    a8 = len( (sA.intersection(sC, sD)).difference(sB) )
    print "8 = (%s int %s int %s) - %s = %d" %(lA, lC, lD, lB, a8)
    a9 = len( sA.intersection(sB, sC, sD) )
    print "9 = %s int %s int %s int %s = %d" %(lA, lB, lC, lD, a9)
    a10 = len( (sB.intersection(sC, sD)).difference(sA) )
    print "10 = (%s int %s int %s) - %s = %d" %(lB, lC, lD, lA, a10)
    a11 = len( (sC.intersection(sD)).difference(sA.union(sB)) )
    print "11 = (%s int %s) - (%s U %s) = %d" %(lC, lD, lA, lB, a11)
    a12 = len( (sA.intersection(sD)).difference(sB.union(sC)) )
    print "12 = (%s int %s) - (%s U %s) = %d" %(lA, lD, lB, lC, a12)
    a13 = len( (sA.intersection(sB, sD)).difference(sC) )
    print "13 = (%s int %s int %s) - %s = %d" %(lA, lB, lD, lC, a13)
    a14 = len( (sB.intersection(sD)).difference(sA.union(sC)) )
    print "14 = (%s int %s) - (%s U %s) = %d" %(lB, lD, lA, lC, a14)
    a15 = len( sD.difference(sA.union(sB, sC)) )
    print "15 = %s - (%s U %s U %s) = %d" %(lD, lA, lB, lC, a15)

    print "--------------"
    print "Verification: %s = %d = %d" % (lA, len(sA), a1 + a2 + a4 + a5 + a8 + a9 + a12 + a13)
    print "              %s = %d = %d" % (lB, len(sB), a2 + a3 + a5 + a6 + a9 + a10 + a13 + a14)
    print "              %s = %d = %d" % (lC, len(sC), a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11)
    print "              %s = %d = %d" % (lD, len(sD), a8 + a9 + a10 + a11 + a12 + a13 + a14 + a15)        

#--------------------------------------------------------------------
def get_labels(data, fill="number"):
    """
to get a dict of labels for groups in data

input
data: data to get label for
fill = ["number"|"logic"|"both"], fill with number, logic label, or both

return
labels: a dict of labels for different sets

example:
In [12]: get_labels([range(10), range(5,15), range(3,8)], fill="both")
Out[12]:
{'001': '001: 0',
'010': '010: 5',
'011': '011: 0',
'100': '100: 3',
'101': '101: 2',
'110': '110: 2',
'111': '111: 3'}
"""

    N = len(data)

    sets_data = [set(data[i]) for i in range(N)] # sets for separate groups
    s_all = set(chain(*data)) # union of all sets

    # bin(3) --> '0b11', so bin(3).split('0b')[-1] will remove "0b"
    set_collections = {}
    for n in range(1, 2**N):
        key = bin(n).split('0b')[-1].zfill(N)
        value = s_all
        sets_for_intersection = [sets_data[i] for i in range(N) if key[i] == '1']
        sets_for_difference = [sets_data[i] for i in range(N) if key[i] == '0']
        for s in sets_for_intersection:
            value = value & s
        for s in sets_for_difference:
            value = value - s
        set_collections[key] = value

    if fill == "number":
        labels = {k: len(set_collections[k]) for k in set_collections}
    elif fill == "logic":
        labels = {k: k for k in set_collections}
    elif fill == "both":
        labels = {k: ("%s: %d" % (k, len(set_collections[k]))) for k in set_collections}
    else: # invalid value
        raise Exception("invalid value for fill")

    return labels

def Plotvenn2(data=None, names=None, fill="number", show_names=True, show_plot=True, **kwds):

    if (data is None) or len(data) != 2:
        raise Exception("length of data should be 2!")
    if (names is None) or (len(names) != 2):
        names = ("set 1", "set 2")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (8, 8)

    fig = pylab.figure(figsize=figsize)
    ax = fig.gca(); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_xlim(0, 8); ax.set_ylim(0, 8)

    # r: radius of the circles
    # (x1, y1), (x2, y2): center of circles
    r, x1, y1, x2, y2 = 2.0, 3.0, 4.0, 5.0, 4.0

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 2:
        colors = kwds['colors']
    else:
        colors = ['red', 'green']

    c1 = Circle((x1,y1), radius=r, alpha=0.5, color=colors[0])
    c2 = Circle((x2,y2), radius=r, alpha=0.5, color=colors[1])

    ax.add_patch(c1)
    ax.add_patch(c2)

    ## draw text
    #1
    pylab.text(x1-r/2, y1, labels['10'], **alignment)
    pylab.text(x2+r/2, y2, labels['01'], **alignment)
    # 2
    pylab.text((x1+x2)/2, y1, labels['11'], **alignment)
    # names of different groups
    if show_names:
        pylab.text(x1, y1-1.2*r, names[0], fontsize=16, **alignment)
        pylab.text(x2, y2-1.2*r, names[1], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    if show_plot:
        pylab.show()

#--------------------------------------------------------------------
def Plotvenn3(data=None, names=None, fill="number", show_names=True, show_plot=True, **kwds):

    if (data is None) or len(data) != 3:
        raise Exception("length of data should be 3!")
    if (names is None) or (len(names) != 3):
        names = ("set 1", "set 2", "set 3")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (10, 10)

    fig = pylab.figure(figsize=figsize) # set figure size
    ax = fig.gca()
    ax.set_aspect("equal") # set aspect ratio to 1
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_xlim(0, 8); ax.set_ylim(0, 8)

    # r: radius of the circles
    # (x1, y1), (x2, y2), (x3, y3): center of circles
    r, x1, y1, x2, y2 = 2.0, 3.0, 3.0, 5.0, 3.0
    x3, y3 = (x1+x2)/2.0, y1 + 3**0.5/2*r

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 3:
        colors = kwds['colors']
    else:
        colors = ['red', 'green', 'blue']

    c1 = Circle((x1,y1), radius=r, alpha=0.5, color=colors[0])
    c2 = Circle((x2,y2), radius=r, alpha=0.5, color=colors[1])
    c3 = Circle((x3,y3), radius=r, alpha=0.5, color=colors[2])
    for c in (c1, c2, c3):
        ax.add_patch(c)

    ## draw text
    # 1
    pylab.text(x1-r/2, y1-r/2, labels['100'], **alignment)
    pylab.text(x2+r/2, y2-r/2, labels['010'], **alignment)
    pylab.text((x1+x2)/2, y3+r/2, labels['001'], **alignment)
    # 2
    pylab.text((x1+x2)/2, y1-r/2, labels['110'], **alignment)
    pylab.text(x1, y1+2*r/3, labels['101'], **alignment)
    pylab.text(x2, y2+2*r/3, labels['011'], **alignment)
    # 3
    pylab.text((x1+x2)/2, y1+r/3, labels['111'], **alignment)
    # names of different groups
    if show_names:
        pylab.text(x1-r, y1-r, names[0], fontsize=16, **alignment)
        pylab.text(x2+r, y2-r, names[1], fontsize=16, **alignment)
        pylab.text(x3, y3+1.2*r, names[2], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    if show_plot:
        pylab.show()

#--------------------------------------------------------------------
def Plotvenn4(data=None, names=None, fill="number", show_names=True, show_plot=True, **kwds):

    if (data is None) or len(data) != 4:
        raise Exception("length of data should be 4!")
    if (names is None) or (len(names) != 4):
        names = ("set 1", "set 2", "set 3", "set 4")

    labels = get_labels(data, fill=fill)

    # set figure size
    if 'figsize' in kwds and len(kwds['figsize']) == 2:
        # if 'figsize' is in kwds, and it is a list or tuple with length of 2
        figsize = kwds['figsize']
    else: # default figure size
        figsize = (10, 10)

    # set colors for different Circles or ellipses
    if 'colors' in kwds and isinstance(kwds['colors'], Iterable) and len(kwds['colors']) >= 4:
        colors = kwds['colors']
    else:
        colors = ['r', 'g', 'b', 'c']

    # draw ellipse, the coordinates are hard coded in the rest of the function
    fig = pylab.figure(figsize=figsize) # set figure size
    ax = fig.gca()
    patches = []
    width, height = 170, 110 # width and height of the ellipses
    patches.append(Ellipse((170, 170), width, height, -45, color=colors[0], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -45, color=colors[1], alpha=0.5))
    patches.append(Ellipse((200, 200), width, height, -135, color=colors[2], alpha=0.5))
    patches.append(Ellipse((230, 170), width, height, -135, color=colors[3], alpha=0.5))
    for e in patches:
        ax.add_patch(e)
    ax.set_xlim(80, 320); ax.set_ylim(80, 320)
    ax.set_xticks([]); ax.set_yticks([]);
    ax.set_aspect("equal")

    ### draw text
    # 1
    pylab.text(120, 200, labels['1000'], **alignment)
    pylab.text(280, 200, labels['0100'], **alignment)
    pylab.text(155, 250, labels['0010'], **alignment)
    pylab.text(245, 250, labels['0001'], **alignment)
    # 2
    pylab.text(200, 115, labels['1100'], **alignment)
    pylab.text(140, 225, labels['1010'], **alignment)
    pylab.text(145, 155, labels['1001'], **alignment)
    pylab.text(255, 155, labels['0110'], **alignment)
    pylab.text(260, 225, labels['0101'], **alignment)
    pylab.text(200, 240, labels['0011'], **alignment)
    # 3
    pylab.text(235, 205, labels['0111'], **alignment)
    pylab.text(165, 205, labels['1011'], **alignment)
    pylab.text(225, 135, labels['1101'], **alignment)
    pylab.text(175, 135, labels['1110'], **alignment)
    # 4
    pylab.text(200, 175, labels['1111'], **alignment)
    # names of different groups
    if show_names:
        pylab.text(110, 110, names[0], fontsize=16, **alignment)
        pylab.text(290, 110, names[1], fontsize=16, **alignment)
        pylab.text(130, 275, names[2], fontsize=16, **alignment)
        pylab.text(270, 275, names[3], fontsize=16, **alignment)

    leg = ax.legend(names, loc='best', fancybox=True)
    leg.get_frame().set_alpha(0.5)

    if show_plot:
        pylab.show()
        
        
def usage():
    print "Given a meta file including the names of 2, 3 or 4 tab_files, give the information necessary to build a Venn diagram"
    print "Usage : overlapProteins.py <metaFile.txt> [-d, --decoy] [-h, --help] [-v, --verbose]"
    print "decoy : prefix to identify decoys"
    
def main(argv):
    if( len(argv) < 1):
        sys.stderr.write("Error: Number of arguments incorrect\n")
        usage()
        sys.exit()
    else:
        verbose = False
        decoy_prefix = "random"
        try:
            opts, args = getopt.getopt(sys.argv[2:], "d:hv", ["decoy=", "help", "verbose"])
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
            elif o in ("-d", "--decoy"):
                decoy_prefix = a
            else:
                assert False, "unhandled option"
        
        if(os.path.isfile(argv[0])):
            infile = argv[0]
        else:
            sys.stderr.write("Error: XML file not found\n")
            sys.exit()
        
        sets = []
        labels = []
        lines = [l for l in open(infile).readlines() if not l.startswith("%")]
        for line in lines:
            words = line.split() 
            if(line != ""):
                prot_file = words[0]
                labels.append(words[1])
                if(os.path.isfile(prot_file)):
                    
                    proteins = util.importer(prot_file)
                    sets.append(set([str(prot.protein) for prot in proteins if prot.protein.find(decoy_prefix) == -1]))
                    if(verbose):
                        print "reading file " +str(prot_file)  

        if len(sets) < 2 or len(sets) > 4:
            print "Error! Cant calculate Venn diagrams for this number of datasets = " + str(len(sets))
            sys.exit(1)

        if len(sets) == 2:
            Get2Venn(sets[0], labels[0], sets[1], labels[1])
            Plotvenn2(sets, labels, fill="both", show_names=False)
        elif len(sets) == 3:
            Get3Venn(sets[0], labels[0], sets[1], labels[1], sets[2], labels[2])
            Plotvenn3(sets, labels, fill="both", show_names=False)
        else: 
            Get4Venn(sets[0], labels[0], sets[1], labels[1], sets[2], labels[2], sets[3], labels[3])
            Plotvenn4(sets, labels, fill="both", show_names=False)

if __name__ == "__main__":
    main(sys.argv[1:]) 
