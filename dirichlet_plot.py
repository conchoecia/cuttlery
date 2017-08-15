#!/usr/bin/env python

#mathy stuff
import numpy as np
from numpy import cumsum
from numpy.random import rand
import scipy.stats as stats
import pandas as pd
import time

# functools because apply doesn't work with two arguments
from functools import partial
from itertools import permutations

#plotting stuff
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import matplotlib.patches as mplpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path
import urllib.request
import sys, os
import copy
from Bio import SeqIO
from Bio.codonalign import CodonSeq
import Bio.Data.CodonTable
from Bio.codonalign.codonalphabet import get_codon_alphabet
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import _get_pi

# Also import some stufffff to parallelize the function
#from multiprocessing import cpu_count
#from multiprocessing import Pool
#from multiprocessing.dummy import Pool as ThreadPool
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# don't print out options
pd.options.mode.chained_assignment = None

from matplotlib import rc
rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}',    # helvetica font
        r'\usepackage{sansmath}',   # math-font matching  helvetica
        r'\sansmath'                # actually tell tex to use it!
        r'\usepackage{siunitx}',    # micro symbols
        r'\sisetup{detect-all}',    # force siunitx to use the fonts
        ]

def main():
    results_file = "/Users/darrin/git/cutlery/dirichlet_results_2k.csv"
    results_df = pd.read_csv(results_file)
    plot_results(results_df)

def plot_results(results):
    print(results.head())
    plt.style.use('BME163')
    #set the figure dimensions
    figWidth = 5
    figHeight = 4
    plt.figure(figsize=(figWidth,figHeight))
    #set the panel dimensions
    panelWidth = 4
    panelHeight = 2.5
    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2 + 0.125
    bottomMargin = ((figHeight - panelHeight - 0.5)/2) + 0.75
    panel0 =plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panel0.spines['top'].set_visible(False)
    panel0.spines['right'].set_visible(False)
    panel0.spines['left'].set_visible(False)

    noncodings = results[results['real_val'] == 'noncoding']
    codings = results[results['real_val'] == 'coding']
    tests = results[results['analysis_type'] == 'test']

    noncoding_seqnames = sorted(noncodings['seqname'].unique())
    coding_seqnames = sorted(codings['seqname'].unique())
    tests_seqnames = sorted(tests['seqname'].unique())

    print(noncoding_seqnames)
    print(coding_seqnames)
    print(tests_seqnames)


    # autumn colors for noncoding
    cmap = cm.get_cmap('autumn')
    autumn = {noncoding_seqnames[i]: cmap(1 - (i/len(noncoding_seqnames)))
            for i in range(len(noncoding_seqnames)) }
    # winter colors for coding
    cmap = cm.get_cmap('winter')
    winter = {coding_seqnames[i]: cmap(1 - (i/len(coding_seqnames)))
            for i in range(len(coding_seqnames)) }
    # magma colors for tests
    cmap = cm.get_cmap('magma')
    magma = {tests_seqnames[i]: cmap(((i/len(coding_seqnames))) ) for i in range(len(tests_seqnames)) }

    xmin= int(min(results['ll_ratio']))
    xmax= np.ceil(max(results['ll_ratio']))
    bins = np.linspace(xmin, xmax, 500)
    print("autumn", autumn)
    # first get the noncoding data
    noncoding_colors = [autumn[noncoding_seqname] for noncoding_seqname in noncoding_seqnames]
    noncoding_data = [noncodings.query("seqname == '{}'".format(nc_seqname))['ll_ratio'] \
                      for nc_seqname in noncoding_seqnames]

    out = panel0.hist(noncoding_data, bins, stacked=True, color = noncoding_colors,
                     alpha = 0.66, label = noncoding_seqnames, normed = True)
    # then get the coding data
    coding_colors = [winter[coding_seqname] for coding_seqname in coding_seqnames]
    coding_data = [codings.query("seqname == '{}'".format(coding_seqname))['ll_ratio'] \
                      for coding_seqname in coding_seqnames]
    panel0.hist(coding_data, bins, stacked=True, color = coding_colors,
                alpha = 0.66, label = coding_seqnames, normed = True)
    # then get the test data
    tests_colors = [magma[tests_seqname] for tests_seqname in tests_seqnames]
    tests_data = [tests.query("seqname == '{}'".format(tests_seqname))['ll_ratio'] \
                      for tests_seqname in tests_seqnames]
    for i in range(len(tests_data)):
        panel0.hist(tests_data[i], bins, color = tests_colors[i],
                    alpha = 0.25, label = tests_seqnames[i],
                    normed = True)

    #data = noncoding_data + coding_data
    #colors = noncoding_colors + coding_colors
    #panel0.hist(data, bins, stacked=True, color = colors,
    #            alpha = 0.75, label = noncoding_seqnames + coding_seqnames)
    panel0.set_xlim([xmin * 1.1, xmax * 1.1])
    panel0.set_xlim([xmin * 1.1, xmax * 1.1])
    panel0.legend(loc = 'upper right', fontsize = 8,
                  ncol = 5, bbox_to_anchor=(1.05, -0.15))

    panel0.set_xlabel("Log-likelihood ratio")
    panel0.set_ylabel("Normalized Observation Proportions")
    panel0.set_title("Codon Usage Log-likelihood Ratios")
    plt.savefig("dirichlet_histogram_{}.png".format(timestamp()), dpi=600, transparent=False)

def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")

if __name__ == "__main__":
    sys.exit(main())
