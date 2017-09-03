#!/usr/bin/env python3

"""file: degenerate(.py)
author: Darrin Schultz
date: 20161111
SOE location: /soe/dschultz/BME205_hw9/degenerate

NOTE! Please set your python3 variable, and execute (export PATH=$PATH:.) in
      your shell prior to executing this script. This will avoid having to call
      the script with the ./ convention.

Program Abstract:
  This program generates the possible sets of amino acids that can be generated
  with all degenerate codons given a codon usage table. The program outputs
  calculated codon imbalance values and codon frequency values, as well as the
  optimal degenerate codon to use for that specific species and translated
  table. Lastly, the program can calculate the optimum degenerate nucleotide
  code to use for a given amino acid sequence, coding table, and codon
  usage frequency.

Progam Deatiled Description:
  Args:
    <--translation_table>
      - (type = str)
      - (default = "gc.prt")
      - This is the file or website from which the codon tables are
        drawn. Every codon table on the page is parsed, then the user
        is given the option to choose one of the codon tables via a
        command prompt. The default is the information found at the url,
        "ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt".
      - You may select which codon table you would like by passing in the string
        under the first `name` field in each curly braces-delimited codon table.
      - If you would like to see a list of options for the codon table URL or
        file, you may do so with the <--gene_code_options> argument.

    <--codon>
      - (type = str)
      - (default = "species199310.html")
      - This is the website or html file from which the
        species-specific codon usage table will be parsed. This parser
        is only designed to work with http://www.kazusa.or.jp
        associated websites.
      - The default is the data found at:
        "http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=199310"
      - You may input any similarly-formatted species-specific codon
        usage frequency table from http://www.kazusa.or.jp

    <--gene_code>
      - (type = str)
      - (default = "Standard")
      - This selects the codon usage type from the codon tables parsed from
        <--translation_table>'s input.
      - The default is 'Standard', this assumes that the user is using the
        default NCBI codon usage ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
      - You may discover other options for <--gene_code> by using the
        <--gene_code_options> command. Copy the entire string (not including
        integer index) of the output of the desired table to pass as an argument
        for <--gene_code>.

    <--gene_code_options>
      - This prints all of the possible options for <--gene_code> that you may
        use in subsequent runs. Prints to STDERR.
      - If you use the default "gc.prt" file or the corresponding url at
        ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt, these are the options
        that will appear:
        0. "Alternative Flatworm Mitochondrial"
        1. "Alternative Yeast Nuclear"
        2. "Ascidian Mitochondrial"
        3. "Bacterial Archaeal and Plant Plastid"
        4. "Blepharisma Macronuclear"
        5. "Candidate Division SR1 and Gracilibacteria"
        6. "Chlorophycean Mitochondrial"
        7. "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"
        8. "Echinoderm Mitochondrial; Flatworm Mitochondrial"
        9. "Euplotid Nuclear"
        10. "Invertebrate Mitochondrial"
        11. "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma"
        12. "Pachysolen tannophilus Nuclear"
        13. "Pterobranchia Mitochondrial"
        14. "Scenedesmus obliquus Mitochondrial"
        15. "Standard"
        16. "Thraustochytrium Mitochondrial"
        17. "Trematode Mitochondrial"
        18. "Vertebrate Mitochondrial"
        19. "Yeast Mitochondrial"

    <--output>
      - (type = str)
      - This is the filename where the output from the program will be saved.
        Both the output from the print option and --optimize_this_protein go
        here. Do not use unix pipes with this program, as you will not be
        prompted to choose a codon table.

  Example usage:
    - All examples first prompt the user for the codon table to use.

    - `degenerate`
      - uses the "minimal" printing style described above
    - `degenerate --optimize_this_protein A[AEV]G[GRWY*] --output_format=full`
      - prints out the "full" printing style described above.
      - outputs the optimized protein sequence
    - `degenerate --optimize_this_protein A[AEV]G[GRWY*]`
      - prints out the "minimal" printing style described above
      - outputs the optimized protein seq
    - `degenerate --output foo.txt`
      - prints the "minimal" printing style to the file 'foo.txt'

"""
#this is here in case this script is run in python 2 for some reason
from __future__ import print_function, division

import argparse
from collections import Counter
import itertools
import operator
import os
import sys
import pandas as pd
pd.options.display.float_format = '{:,.3f}'.format
import urllib.request
from Bio import SeqIO

def parse_arguments():
    """This method parses the arguments provided during command-line input
    to control how the kmers are counted. Please see the documentation for help.
    """

    #using Kevin's idea here to make te docstring the help file
    parser=argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta",
                        type=str,
                        help='the input fasta file to use')
    parser.add_argument("--translation_table", type=str,
                        default="ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt",
                        help="""The genetic code table to read from. The default
                        is from ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt""")
    parser.add_argument("--codon", type=str,
                        default="http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=199310",
                        help="""The codon usage table to read from. The default
                        table is from E. coli.
                        http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=199310
                        You may use any species from the www.kazusa.or.jp
                        website with this parser.""")
    parser.add_argument("--output", type=str,
                        default=sys.stdout,
                        help="""the output will be saved here""")
    parser.add_argument("--tt_code", type = int,
                        default = 1,
                        help="""Select which gene code to use. Default is
                        Standard""")
    parser.add_argument("--tt_options", action = 'store_true',
                        default = False,
                        help="""Display the optional gene code names that you
                        may pass to the <--tt_code> argument in a subsequent
                        run""")

    args = parser.parse_args()
    return args

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

class CodonUsageTable:
    def __init__(self, transTable):
        self.codondf = pd.DataFrame.from_dict(transTable, orient='index')
        self.codondf.index.name = 'codon'
        self.codondf.reset_index(inplace=True)
        self.codondf.columns = ['codon', 'AA']
        self.codondf['Fraction'] = 0
        self.codondf['Frequency'] = 0
        self.codondf['Number'] = 0
        self.codondf.sort_values(['AA', 'codon'], ascending=[True, True], inplace=True)
        self.codondf.reset_index(drop=True,inplace=True)
        self.seq = ''
        self.name = ''
        self.codonCount = 0
        self.codingGC = 0.0
        self.firstGC  = 0.0
        self.secondGC = 0.0
        self.thirdGC  = 0.0

    def fill_table(self, seq, name):
        self.name = name
        self.seq = seq.upper()
        codons = [seq[i:i+3] for i in range(0,len(seq),3)]
        self.codonCount = len(codons)
        self.codingGC = (seq.count('G') + seq.count('C'))/len(seq)
        firstSeq = "".join([codon[0] for codon in codons])
        self.firstGC = (firstSeq.count('G') + firstSeq.count('C'))/len(firstSeq)
        secondSeq = "".join([codon[1] for codon in codons])
        self.secondGC = (secondSeq.count('G') + secondSeq.count('C'))/len(secondSeq)
        thirdSeq = "".join([codon[2] for codon in codons])
        self.thirdGC = (thirdSeq.count('G') + thirdSeq.count('C'))/len(thirdSeq)
        for each in codons:
            self.codondf.loc[self.codondf['codon'] == each, 'Number'] += 1
        self.codondf['Frequency'] = (self.codondf['Number']/self.codonCount) * 1000
        uniqueAA = self.codondf['AA'].unique()
        self.AA_counts = {AA: sum(self.codondf.loc[self.codondf['AA'] == AA, 'Number'])
                        for AA in uniqueAA}
        print(self.AA_counts)
        self.codondf['Fraction'] = self.codondf.apply(self.calc_fraction, axis=1)

    def calc_fraction(self, row):
        if not row['Number'] == 0:
            return row['Number']/self.AA_counts[row['AA']]
        else:
            return 0.0

    def __repr__(self):
        joinThese = []
        joinThese.append("#Name: {}\n\n".format(self.name))
        joinThese.append("#NumberOfCodons: {}\n\n".format(self.codonCount))
        joinThese.append("#Coding GC {0:0.2f}%\n".format(self.codingGC * 100))
        joinThese.append("#1st letter GC {0:0.2f}%\n".format(self.firstGC * 100))
        joinThese.append("#2nd letter GC {0:0.2f}%\n".format(self.secondGC * 100))
        joinThese.append("#3rd letter GC {0:0.2f}%\n\n".format(self.thirdGC * 100))
        #joinThese.append(str(self.codondf))
        return "".join(joinThese)

class PrepCUTforPlot:
    def __init__(self, CUTS):
        #initialize the first df based on the 0th element,
        # then just add everything else.
        self.fracdf = CUTS[0].codondf.copy(deep=True)
        self.fracdf.rename(columns={'Fraction':CUTS[0].name},inplace=True)
        self.fracdf.drop(['Frequency', 'Number'], axis = 1, inplace=True)
        for i in range(1,len(CUTS)):
            self.fracdf[CUTS[i].name] = CUTS[i].codondf['Fraction']

def plot_codon_usage(cleandf, title):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mplpatches
    import matplotlib.colors as colors
    from scipy.stats import mannwhitneyu as mu
    from collections import Counter
    from matplotlib import rcParams
    import numpy as np
    import sys
    import os
    import time
    # This mpl style is from the UCSC BME163 class.
    rcParams.update({
    'font.size'           : 8.0      ,
    'font.sans-serif'     : 'Helvetica'    ,
    'xtick.major.size'    : 2        ,
    'xtick.major.width'   : 0.75     ,
    'xtick.labelsize'     : 8.0      ,
    'xtick.direction'     : 'out'      ,
    'ytick.major.size'    : 2        ,
    'ytick.major.width'   : 0.75     ,
    'ytick.labelsize'     : 8.0      ,
    'ytick.direction'     : 'out'      ,
    'xtick.major.pad'     : 2        ,
    'xtick.minor.pad'     : 2        ,
    'ytick.major.pad'     : 2        ,
    'ytick.minor.pad'     : 2        ,
    'savefig.dpi'         : 600      ,
    'axes.linewidth'      : 0.75     ,
    'text.usetex'         : True     ,
    'text.latex.unicode'  : True     ,
    'mathtext.fontset'    : 'custom' ,
    'mathtext.rm'         :'Helvetica'})

    df = cleandf.transpose()

    #set the figure dimensions
    figWidth = 5
    figHeight = 3
    plt.figure(figsize=(figWidth,figHeight))

    #set the panel dimensions
    panelWidth = 4
    panelHeight = 2

    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2)

    panel0=plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='off', labelbottom='off',\
                       left='on', labelleft='on', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')


    #dfnoNAD2L = df.drop(['NAD2L', 'RCNAD2LNAD21', 'RCNAD2LNAD22','RCNAD2LNAD23'], axis = 1)
    dfnoNAD2L = df
    print(dfnoNAD2L)
    observations = []
    for n, col in enumerate(dfnoNAD2L.columns):
        observations.append(list(dfnoNAD2L.iloc[2:,col]))

    bp=panel0.boxplot(observations, \
                      positions=np.arange(1, len(observations)+1, 1), \
                      patch_artist=True,widths=0.5)

    for box in bp['boxes']:
        box.set(edgecolor='black',facecolor=(1, 1, 1, 0.0),linewidth=0.5)
    for whisker in bp['whiskers']:
        whisker.set(color='black', linestyle='-',linewidth=0.25)
    for median in bp['medians']:
        median.set(color='black', linestyle='-',linewidth=0.5)
    for flier in bp['fliers']:
        flier.set(markersize=0)

    #panel0.set_xticks(np.arange(1, len(observations)+1, 1), minor=False)
    #panel0.set_xticklabels(rotlabels, fontdict=None, minor=False, fontsize=4)
    panel0.set_xlim([0,len(observations)+2])
    #panel0.set_ylim([-0.25,1.1])
    panel0.set_yscale("log")
    panel0.set_yticks(np.arange(0, 1.25, 0.25))
    #panel0.set_yticks(np.arange(0,16, 2))
    #panel0.set_xticks(np.arange(1, 5, 1))
    panel0.set_ylabel("Fraction", \
                       fontsize = 8)
    #panel0.tick_params(axis='both', which='major', labelsize=8)
    for each in panel0.spines:
        panel0.spines[each].set_visible(False)

    # colors mostly taken from https://biologywarakwarak.files.wordpress.com/2012/01/amino-acid-table-singlet-code.png
    colormap = {"F": '#c0d9ef',
                "L": '#e8f1df',
                "I": '#f7dedc',
                "M": '#ff2e00',
                "V": '#ffc239',
                "S": '#ffff54',
                "P": "#7fce66",
                "T": "#00ae60",
                "A": "#00aeec",
                "Y": "#006fbb",
                "*": "#ffffff",
                "H": "#ded9c5",
                "Q": "#ffc294",
                "N": "#b5a2c4",
                "K": "#968b5a",
                "D": "#00fc65",
                "E": "#00dcf0",
                "C": "#ff994e",
                "W": "#dc31e6",
                "R": "#d8d8d8",
                "G": "#abdce7"}
    AAcoords={}
    uniqueAAs = df.loc['AA'].unique()
    AAlist = list(df.loc['AA'])

    for each in uniqueAAs:
        # add one to the end of these two because we start plotting at 1, not 0
        firstx= AAlist.index(each) + 1
        lastx = len(AAlist) - 1 - AAlist[::-1].index(each) + 1
        AAcoords[each] = (firstx, lastx)

    for key in AAcoords:
        xStart = AAcoords[key][0]-0.5
        xStop = AAcoords[key][1]+0.5
        print(key, xStart, xStop)
        yStart = -0.25
        yStop = -0.05
        width = xStop - xStart
        height = abs(yStop - yStart)
        facec = colors.to_rgba(colormap[key])
        if key != '*':
            rectangle1=mplpatches.Rectangle((xStart,yStart),width,height,\
                                linewidth=0.0,\
                                facecolor=facec,\
                                edgecolor=(0,0,0),\
                                alpha=0.8)
        else:
            rectangle1=mplpatches.Rectangle((xStart,yStart),width,height,\
                                linewidth=0.25,\
                                facecolor=facec,\
                                edgecolor=(0,0,0),\
                                alpha=0.8)

        panel0.add_patch(rectangle1)
    #now add the text labels on top of the boxes
    labels = list(dfnoNAD2L.loc['codon',])
    rotlabels = ["{0}\n{1}\n{2}".format(x[0],x[1],x[2]) for x in labels]
    for i in range(len(rotlabels)):
        label = rotlabels[i]
        #add one since we start plotting at 1
        xpos = i + 1
        ypos = -0.06
        panel0.text(xpos, ypos, label,
                    horizontalalignment='center',
                    verticalalignment='top',
                    size = 4)
                    #transform=panel0.transAxes)
    threeletcode = {'*': 'Sto',
                    'A': 'Ala',
                    'R': 'Arg',
                    'N': 'Asn',
                    'D': 'Asp',
                    'C': 'Cys',
                    'E': 'Glu',
                    'Q': 'Gln',
                    'G': 'Gly',
                    'H': 'His',
                    'I': 'Ile',
                    'L': 'Leu',
                    'K': 'Lys',
                    'M': 'Met',
                    'F': 'Phe',
                    'P': 'Pro',
                    'S': 'Ser',
                    'T': 'Thr',
                    'W': 'Trp',
                    'Y': 'Tyr',
                    'V': 'Val'}
    #now add the one-leter or three-letter codes
    for key in AAcoords:
        width = AAcoords[key][1]-AAcoords[key][0]
        xpos = AAcoords[key][1] - (width/2)
        ypos = -0.185
        if width > 1:
            label = threeletcode[key]
        else:
            label = key
        panel0.text(xpos, ypos, label,
                    horizontalalignment='center',
                    verticalalignment='top',
                    size = 6)

    masterPP = []
    start = time.time()
    for i in np.arange(0,len(observations),1):
        center=1 + i
        left_bound= i
        right_bound=2 + i
        step_size=0.04
        cutoff=0.06
        placed_points=[]
        counter=0
        for y_value in observations[i]:
            counter+=1
            if len(placed_points)==0:
                placed_points.append((center, y_value))
            else:
                potential_x_position=[]
                left_search = np.arange(center, left_bound, -1 * step_size)
                right_search = np.arange(center,right_bound, step_size)
                xBinarySearch = [val for pair in zip(left_search, right_search) for val in pair][1:]
                for x_position in xBinarySearch:
                    distances=[]
                    binaryPoints = [tup for tup in placed_points if abs(y_value - tup[1]) < cutoff*2]
                    if len(binaryPoints) == 0:
                        potential_x_position.append(x_position)
                    else:
                        for placed_point in binaryPoints:
                            distance=((x_position-placed_point[0])**2+(y_value-placed_point[1])**2)**0.5
                            distances.append(distance)
                        if min(distances)>cutoff:
                            potential_x_position.append(x_position)
                if len(potential_x_position)>0:
                     best_x_position=sorted(potential_x_position,key=lambda x: np.absolute(x-center))[0]
                     placed_points.append((best_x_position,y_value)) 
                else:
                     print('point not placed: ',y_value)
        masterPP += placed_points

    end = time.time()
    darrintime = end-start
    print("dot algorithm took {} seconds".format(end - start))
    facec = colors.to_rgba("#00aeec", 0.6)
    for point in masterPP:
        panel0.plot(point[0],point[1],marker='o',ms=2,mfc=facec, mew=0,linewidth=0)


    GOI = {'NAD2L': {'marker': '*', 'mfc': 'red'},
           'RCNAD2LNAD21'  : {'marker': '^', 'mfc': 'purple'},
           'RCNAD2LNAD22'  : {'marker': 'v', 'mfc': 'purple'},
           'RCNAD2LNAD23'  : {'marker': '>', 'mfc': 'purple'}}

    for key in GOI:
        #now plot the points we are interested in for comparison
        tempdf = df.loc[key,]
        pointlist = zip(np.arange(1, len(observations)+1, 1), tempdf)
        for point in pointlist:
            panel0.plot(point[0],point[1],
                        marker=GOI[key]['marker'],ms=4,
                        mfc=GOI[key]['mfc'],alpha =0.75,
                        mew=0,linewidth=0)

    #plt.show()
    #plt.savefig('Schultz_Darrin_BME163_Assignment_Week4.png')
    plt.savefig('Schultz_Darrin_BME163_Assignment_Week4.pdf')


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

def main():
    """ Controls argument parsing, creating markov model from stdin, calculating
    cost of seqs in positional <file> argument. Calculates encoding costs and
    prints to stdout."""
    options = parse_arguments()
    print(options)

    #prompt user for input
    if options.tt_options:
        print("""The options for codon table are below. Please use the full name
        when passing the argument to:""")
        myCodons = CodonFetch(options.translation_table)
        print("{0:>2}  {1:}".format('ID', 'DESCRIPTION'))
        for key in myCodons.codonObjects:
            codonObject = myCodons.codonObjects[key]
            print("{0:>2}) {1:}".format(codonObject.tt_id, codonObject.description))
    else:
        myCodons = CodonFetch(options.translation_table)
        thisCodonObj = myCodons.codonObjects[options.tt_code]
        print(thisCodonObj.description)
        print(thisCodonObj.table)
        print(thisCodonObj.AAs)
        for record in SeqIO.parse(options.fasta, "fasta"):
           thisCodonObj.add_CUT(record.seq, record.name)
        #for CUT in thisCodonObj.codonUsageTables:
        #    print(CUT)
        thisCleaned = PrepCUTforPlot(thisCodonObj.codonUsageTables)
        plot_codon_usage(thisCleaned.fracdf, 'temp')


    #analyze the optimal degenerate codon for each amino acid string
    #analyze_this = fetch.codonObjects[options.gene_code]
    #analyze_this.generate_degen_to_aa_counts()
    #analyze_this.printer(options.output_format, output=options.output)

if __name__ == "__main__":
    sys.exit(main())
