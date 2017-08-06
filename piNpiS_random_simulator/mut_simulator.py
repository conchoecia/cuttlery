#!/usr/bin/env python

import numpy as np
import scipy.stats as stats
import pandas as pd


#plotting stuff
import argparse
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as mplpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path
import urllib.request
import sys, os
import copy
from Bio import SeqIO

import progressbar
#pd.options.mode.chained_assignment = None

# This program needs to do a few things
# 0. Select a codon usage table
# 1. Read in a fasta alignment
# 2. Determine all of the polymorphic sites in the alignment
# 3. Calculate pi, piN, and piS for the alignment

def printer():
    B = mutate(np.copy(A), 35)
    hets = (A==B)
    A_gap= (A=='-')
    B_gap= (B=='-')
    hets_A_and = np.logical_or(hets, A_gap)
    hets_all   = np.logical_or(hets_A_and, B_gap)
    # convert to show hets instead of matches
    hets_final = np.logical_not(hets_all).astype(int)
    myCodons = CodonFetch("gc.prt.txt")
    thisCodonObj = myCodons.codonObjects[options.tt_code]
    print(thisCodonObj.description)
    #print(thisCodonObj.table)
    #print(thisCodonObj.AAs)

    Astr = "".join(A)
    Bstr = "".join(B)
    print(len(hets_final))
    Hetstr = ''.join(hets_final.astype(int).astype(str))
    print(len(Hetstr))
    Acodons = [Astr[i:i+3] for i in range(0,len(Astr),3)]
    Bcodons = [Bstr[i:i+3] for i in range(0,len(Bstr),3)]
    Hcodons = [Hetstr[i:i+3] for i in range(0,len(Hetstr),3)]
    Aprint = ""
    Aaa    = ""
    Bprint = ""
    Baa    = ""
    Hprint = ""
    synonymous = 0
    nonsynonymous = 0
    for i in range(len(Acodons)):
        if int(Hcodons[i]) > 0:
            thisAaa = thisCodonObj.table[Acodons[i]]
            thisBaa = thisCodonObj.table[Bcodons[i]]
            Aprint = Aprint + Acodons[i] + " "
            Aaa    = Aaa + " {}  ".format(thisAaa)
            Bprint = Bprint + Bcodons[i] + " "
            Baa    = Baa + " {}  ".format(thisBaa)
            Hprint = Hprint + Hcodons[i] + " "
            if thisAaa == thisBaa:
                synonymous += 1
            else:
                if Hcodons[i] in ["001", "010", "100"]:
                    nonsynonymous +=1
                elif Hcodons[i] in ["011", "101"]:
                    synonymous +=1
                    nonsynonymous +=1
                elif Hcodons[i] == "110":
                    nonsynonymous +=2
                elif Hcodons[i] == "111":
                    nonsynonymous += 3
    print(Aprint)
    print(Aaa)
    print(Baa)
    print(Bprint)
    print(Hprint)
    print("there are {} het sites out of {}, or {:.4f}".format(sum(hets_final), len(hets_final), sum(hets_final)/len(hets_final)))
    print("there are {} S sites out of {}".format(synonymous, synonymous+nonsynonymous))
    print("there are {} N sites out of {}".format(nonsynonymous, synonymous+nonsynonymous))
    nSites = len(Astr)
    piN = nonsynonymous/nSites
    piS = synonymous/nSites
    piNpiS = piN/piS
    print("piN/piS is {:.4f}".format(piNpiS))

def parse_arguments():
    """This method parses the arguments provided during command-line input
    to control how the kmers are counted. Please see the documentation for help.
    """

    #using Kevin's idea here to make te docstring the help file
    parser=argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta",
                        type=str,
                        help="""the input fasta file to use. Don't use an alignment""")
    parser.add_argument("--translation_table", type=str,
                        default="ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt",
                        help="""The genetic code table to read from. The default
                        is from ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt""")
    parser.add_argument("--outprefix", type=str,
                        default=sys.stdout,
                        help="""the output will be saved here""")
    parser.add_argument("--tt_code", type = int,
                        default = 4,
                        help="""Select which gene code to use. Default is
                        Standard""")
    parser.add_argument("--tt_options", action = 'store_true',
                        default = False,
                        help="""Display the optional gene code names that you
                        may pass to the <--tt_code> argument in a subsequent
                        run""")

    args = parser.parse_args()
    return args

def forward_sum(array, window):
    forward = np.zeros(len(array), dtype=int)
    for i in range(len(array)-window):
        #print(i, array[i:i+window])
        forward[i] = sum(array[i:i+window])
    return forward

def reverse_sum(array, window):
    reverse = np.zeros(len(array), dtype=int)
    for i in range(window, len(array)):
        #print(i, array[i:i+window])
        reverse[i] = sum(array[i-window+1:i+1])
    return reverse

def middle_sum(array, window):
    middle = np.zeros(len(array), dtype=int)
    windivtwo = int((window-1)/2)
    for i in range(windivtwo+1, len(array)-windivtwo):
        #print(i, array[i:i+window])
        middle[i] = sum(array[i-windivtwo:i+windivtwo+1])
    return middle

class CodonObject:
    """This object stores a codon table and performs calculations of
    degenerate codon imbalance, as well as calculates an optimized degenerate
    codon string from a set of amino acids"""

    def __init__(self, name, tt_id, description, ttDict):
        # the name of the codon usage table. often sam as description
        self.name = name
        #this is the id of the codon usage table taken from the file or NCBI website
        self.tt_id = tt_id
        #the long description of the codon usage table
        self.description = description

        #This is the codon usage table
        #  key = codon as <str>
        #  value = frequency per thousand as <str>
        self.codonTable = ttDict

        #self.table
        #  key = codon s.t. each x in codon in ["G","A","T","C"]
        #  value = resultant amino acid from codon
        self.table = ttDict
        self.AAs = "".join(sorted({x for x in self.table.values()}))

        #one codon usage table per sequence
        self.codonUsageTables = []

    def add_CUT(self, seq, name):
        # count the info from a sequence and put it into a CodonUsageTable object
        thisUsageTable = CodonUsageTable(self.table)
        thisUsageTable.fill_table(seq, name)
        self.codonUsageTables.append(thisUsageTable)

class CodonFetch:
    """This class is a parser that fetches a codon useage table using the
    NCBI format.

    Use the parse_html method to return a dictionary with the codon
    names as keys, and CodonObjects as the values.

    The initialization parameters for the codon usage table can be a
    file, or a URL that contains the information.

    The input format is text, using the string '--' as comments to ignore.
    An example input is:
    ```example start
     {
      name "Standard" ,
      name "SGC0" ,
      id 1 ,
      ncbieaa  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
      sncbieaa "---M---------------M---------------M----------------------------"
      -- Base1  TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
      -- Base2  TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
      -- Base3  TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
     },
    ```example end

    - The first `name` field that occurs may be on one line or split to two. This
       string shall be used to select which codon usage with the --tt_code option.
       The default for --tt_code is 'Standard'.
    - The second `name` field is not important, but a necessary placeholder.
       Please leave the field blank after the keyword if there is no
       associated value.
    - The `id` field is also not important. Please leave the field blank after
       the keyword if there is no associated value.
    - The `ncbieaa` field is a positional string of the amino acids resulting
        from a specific codon
    - The `sncbieaa` field is not important for this analysis. Please leave the
       field blank after the keyword if there is no associated value.
    - The `Base1`, `Base2`, and `Base3` fields are the positional strings
       representing the nucleotides in positions 1, 2, and 3 in the codon.
       The indices of these nucleotides correspond to the indices of the
       `ncbieaa` field.
    """

    def __init__(self, codonURL = None):
        ncbiCodonUsage=''
        if os.path.exists(os.path.abspath(codonURL)):
             # print("""You have selected to parse the codon usage information
             # from the following file: {}""".format(codonURL), file=sys.stderr)
             handle = open(codonURL, "r")
             ncbiCodonUsage = "".join(list(handle))
        else:
            # print("""You have selected to parse the codon usage information
            # from the following url: {}""".format(codonURL), file=sys.stderr)
            with urllib.request.urlopen(codonURL) as response:
               ncbiCodonUsage = response.read().decode("utf-8")
        #self.usageTable = usageTable
        self.codonObjects = self.parse_ncbiCodonUsage(ncbiCodonUsage)

    def parse_ncbiCodonUsage(self, ncbiCodonUsage):
        """Parse the ncbiCodonUsage file and return a dictionary with codon names as keys
        and CodonObjects as values"""
        split = ncbiCodonUsage.split("{")
        #only keep the entries with a codon table
        codes = [x for x in split if x.strip().startswith("name ")]
        codonObjects = {}
        for each in codes:
            if each:
                thisCodonObject = self.parse_single(each)
                if thisCodonObject:
                    codonObjects[thisCodonObject.tt_id] = thisCodonObject
        return codonObjects

    def parse_single(self, codonRawString):
        """Parse a single block of text that will become a single CodonObject.
        Return a single CodonObject"""
        lines = codonRawString.strip().split("\n")
        #this then splits on commas, removes blank lines, and flattens list
        lines = [item.strip()
                 for sublist in lines
                     for item in sublist.strip().split(",")
                 if item.strip()]
        if lines[1].split()[0] != "name":
            lines[0] = lines[0] + " " + lines[1]
            del lines[1]
        #there shouldn't be anything after the " character on the first line,
        #  so split those up
        if not lines[0].endswith("\""):
            idIndex = lines[0].rfind("\"")
            lines.insert(1, lines[0][idIndex+1:])
            lines[0] = lines[0][0:idIndex+1]
        #remove all the entries that are just '}' characters
        lines = [x.replace("\"", "") for x in lines if x.replace("}", "")]
        description = " ".join(lines[0].split()[1:])
        i=0
        if len(lines) < 8:
            name = description
            i = 1
        else:
            name = lines[1].split()[1].replace("\"", "")
        tableid = int(lines[2-i].split()[1].replace("\"", ""))
        AAs = lines[3-i].split()[1].replace("\"", "")
        base1 = lines[5-i].split()[2]
        base2 = lines[6-i].split()[2]
        base3 = lines[7-i].split()[2]
        ttDict = {"{}{}{}".format(base1[i], base2[i], base3[i]):AAs[i]
                     for i in range(len(AAs))}
        #here, I verified that the file parsers all return 21 AAs
        #return CodonObject(name, description, codonDict, self.usageTable)
        return CodonObject(name, tableid, description, ttDict)

def mutate(sequence_array, n_muts):
    bases = ['G', 'A', 'T', 'C']
    mut_sites = np.random.choice(range(3, len(sequence_array)), n_muts, replace=False)
    #print(mut_sites)
    for i in mut_sites:
        options = [x for x in bases if x != sequence_array[i]]
        prev = sequence_array[i]
        sequence_array[i] = np.random.choice(options, 1)[0]
        #print(prev, sequence_array[i])
    return sequence_array

def calcpiNpiS(thisCodonObj, sequenceA, sequenceB):
    A = sequenceA
    B = sequenceB
    hets = (A==B)
    A_gap= (A=='-')
    B_gap= (B=='-')
    hets_A_and = np.logical_or(hets, A_gap)
    hets_all   = np.logical_or(hets_A_and, B_gap)
    # convert to show hets instead of matches
    hets_final = np.logical_not(hets_all).astype(int)

    Astr = "".join(A)
    Bstr = "".join(B)
    Hetstr = ''.join(hets_final.astype(int).astype(str))
    Acodons = [Astr[i:i+3] for i in range(0,len(Astr),3)]
    Bcodons = [Bstr[i:i+3] for i in range(0,len(Bstr),3)]
    Hcodons = [Hetstr[i:i+3] for i in range(0,len(Hetstr),3)]
    synonymous = 0
    nonsynonymous = 0
    for i in range(len(Acodons)):
        if int(Hcodons[i]) > 0:
            thisAaa = thisCodonObj.table[Acodons[i]]
            thisBaa = thisCodonObj.table[Bcodons[i]]
            if thisAaa == thisBaa:
                if Hcodons[i] in ["001", "010", "100"]:
                    synonymous +=1
                elif Hcodons[i] in ["011", "101", "110"]:
                    synonymous +=2
                elif Hcodons[i] == "111":
                    synonymous += 3
            else:
                if Hcodons[i] in ["001", "010", "100"]:
                    nonsynonymous +=1
                elif Hcodons[i] in ["011", "101"]:
                    synonymous +=1
                    nonsynonymous +=1
                elif Hcodons[i] == "110":
                    nonsynonymous +=2
                elif Hcodons[i] == "111":
                    nonsynonymous += 3
    nSites = len(Astr)
    piN = nonsynonymous/nSites
    piS = synonymous/nSites
    try:
        piNpiS = piN/piS
        return piNpiS
    except:
        return -1

def pi(seqs):
    """
    there is some info here http://www.columbia.edu/cu/biology/courses/c3020/solutions-2.html
    """
    df_list = []
    names_list = []
    for index in range(len(seqs)):
        names_list.append("seq{}".format(index))
    df = pd.DataFrame.from_items(zip(names_list, seqs))
    #df.apply(rowwise_unique, axis=1)

    x = []
    #this block calculates the frequencies of each sequence
    for i in range(len(seqs)):
        this_seq = seqs[i]
        this_x = 0
        for j in range(len(seqs)):
            if np.array_equal(seqs[i], seqs[j]):
                this_x += 1
        x.append(this_x/len(seqs))

    print(x)
    running_sum = 0
    for i in range(len(seqs)):
        for j in range(i+1, len(seqs)):
            comp = df.copy().iloc[:,[i,j]]
            comp.replace('-', np.nan, inplace=True)
            comp.dropna(axis=0,how='any', inplace=True)
            #print(comp)
            num_difs = sum(comp.iloc[:,0] != comp.iloc[:,1])
            len_seqs = len(comp.iloc[:,0])
            this_sum = x[i] * x[j] * (num_difs/len_seqs)
            running_sum += this_sum
            #print("x[i] * x[j] * pi = {} * {} * {}/{}".format(x[i], x[j], num_difs, len_seqs))
    print("running sum: {}".format(running_sum))
    return running_sum

def main():
    #First, read in the options
    options = parse_arguments()
    print(options)
    #Second, get the translation table
    myCodons = CodonFetch(options.translation_table)
    # If the user requests to see the options, print them and exit
    if options.tt_options:
        print("Options for codon usage tables in --tt_code are:")
        for key in sorted(myCodons.codonObjects):
            print("{}: {}".format(key, myCodons.codonObjects[key].description))
        sys.exit()

    # Third, select the codon table that was passed as an option
    thisCodonObj = myCodons.codonObjects[options.tt_code]
    print("Selected", thisCodonObj.description)
    #print(thisCodonObj.table)
    #print(thisCodonObj.AAs)

    seqA=[]
    #Fourth,read in the sequences to mutate
    for record in SeqIO.parse(options.fasta, "fasta"):
        seqA.append(np.array(record.seq))

    pi(seqA)

    #piNpiS = []
    #numSimulations = 100001
    #bar = progressbar.ProgressBar()
    ##Third, mutate this sequence 35 times
    #for i in bar(range(numSimulations)):
    #    for j in range(len(seqA)):
    #        A = seqA[j]
    #        B = mutate(copy.deepcopy(A), 21)
    #        thisVal =calcpiNpiS(thisCodonObj, copy.deepcopy(A), B)
    #        if thisVal != -1:
    #            piNpiS.append(thisVal)
    #print(piNpiS)

    ## plot normed histogram
    #plt.hist(piNpiS, bins=100, normed=True)
    ## find minimum and maximum of xticks, so we know
    ## where we should compute theoretical distribution
    #xt = plt.xticks()[0]
    #xmin, xmax = min(xt), max(xt)
    #lnspc = np.linspace(xmin, xmax, len(piNpiS))
    ## exactly same as above
    #ag,bg,cg = stats.gamma.fit(piNpiS, loc=0)
    #print(ag, bg, cg)
    #pdf_gamma = stats.gamma.pdf(lnspc, ag, bg,cg)
    #print(pdf_gamma)
    #plt.plot(lnspc, pdf_gamma, label="Gamma")
    #plt.xlim([0,10])
    #plt.show()

if __name__ == "__main__":
    sys.exit(main())
