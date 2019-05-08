#!/usr/bin/env python

# cuttlery - 'Codon Usage Table Tools-lery.'
# Copyright (c) 2016-2017 Darrin T. Schultz. All rights reserved.
#
# This file is part of cuttlery.
#
# cuttlery is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cuttlery is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cuttlery.  If not, see <http://www.gnu.org/licenses/>.
#
#Created on Thu Jun  1 14:54:48 2017
"""
title: cuttlery dirichlet
authors: - Darrin T Schultz (github@conchoecia)
         - Jordan M Eizenga (github@jeizenga)
conception: Jordan M Eizenga (and a little bit DTS)


In short, this test asks if the trinucleotide frequency of a test ORF
more closely matches the trinucleotides frequencies of coding or
noncoding DNA from the same species.
"""

#info about confusion matrix here: http://www.dataschool.io/simple-guide-to-confusion-matrix-terminology/

# - TODO
#   - Streamline how sequences are stored. Use a dictionary or something
#   - Add the ability to perform tests independently.
#     - Instead of randomly selecting sequences add the ability to permute
#       through all possible combinations.


import argparse
import collections
import math
import itertools
import sys
import random
from Bio import SeqIO
from builtins import range
from scipy import stats
import progressbar
import numpy as np
import matplotlib
matplotlib.use('agg')
#matplotlib.use('webagg')
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import cm
import os
import pandas as pd
from sklearn.metrics import confusion_matrix as cmcalc

# multiprocessing stuff
from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
import time
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# cuttlery.codonfunctions stuff
from cuttlery.codonfunctions import fasta_dir_to_gene_filelist, timestamp,\
    print_images


from matplotlib import rc
rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}',    # helvetica font
        r'\usepackage{sansmath}',   # math-font matching  helvetica
        r'\sansmath'                # actually tell tex to use it!
        r'\usepackage{siunitx}',    # micro symbols
        r'\sisetup{detect-all}',    # force siunitx to use the fonts
        ]


global nts
global starts
global stops
nts = ['A', 'C', 'G', 'T']
starts = ['TTA', 'TTG','CTG',
          'ATT', 'ATC','ATA',
          'ATG','GTG']
stops = ['TAA','TAG']
global analyses
analyses = set()

def parallel_process(array, function, n_jobs=2, use_kwargs=False, front_num=3):
    """
        A parallel version of the map function with a progress bar.

        Args:
            array (array-like): An array to iterate over.
            function (function): A python function to apply to the elements of array
            n_jobs (int, default=16): The number of cores to use
            use_kwargs (boolean, default=False): Whether to consider the elements of array as dictionaries of 
                keyword arguments to function
            front_num (int, default=3): The number of iterations to run serially before kicking off the parallel job. 
                Useful for catching bugs
        Returns:
            [function(array[0]), function(array[1]), ...]

        For jobs that have a quick compute time and for many threads, it is
         probably more efficient to split this process up so that it runs in
         chunks rather than take a lot of overhead creating new jobs. This also
         needs to be optimized for list unwrapping

        I optimized for this program and it works best in chunks between 20 and
         40 analyses per process. Using one analysis per pcoess is 7x slower.
    """
    #We run the first few iterations serially to catch bugs
    if front_num > 0:
        front = [function(**a) if use_kwargs else function(a) for a in array[:front_num]]
    #If we set n_jobs to 1, just run a list comprehension. This is useful for benchmarking and debugging.
    if n_jobs==1:
        return front + [function(**a) if use_kwargs else function(a) for a in tqdm(array[front_num:])]
    #Assemble the workers
    with ProcessPoolExecutor(max_workers=n_jobs) as pool:
        #Pass the elements of array into function
        if use_kwargs:
            futures = [pool.submit(function, **a) for a in array[front_num:]]
        else:
            futures = [pool.submit(function, a) for a in array[front_num:]]
        kwargs = {
            'total': len(futures),
            'unit': 'it',
            'unit_scale': True,
            'leave': True
        }
        #Print out the progress as tasks complete
        for f in tqdm(as_completed(futures), **kwargs):
            pass
    out = []
    #Get the results from the futures.
    for i, future in tqdm(enumerate(futures)):
        try:
            out.append(future.result())
        except Exception as e:
            out.append(e)
    return front + out

def interpret_hypothesis_test(log_likelihood_ratio):
    """
     Author:
    - Jordan M Eizenga (github@jeizenga)
    """
    if log_likelihood_ratio < 0:
        favored = 2
        unfavored = 1
        lookup_ratio = -log_likelihood_ratio
    else:
        favored = 1
        unfavored = 2
        lookup_ratio = log_likelihood_ratio

    lookup_ratio = lookup_ratio

    if lookup_ratio / math.log(10) < 10**0.5:
        jeffreys_support = "barely worth mentioning"
    elif lookup_ratio / math.log(10) < 10**1.0:
        jeffreys_support = "substantial"
    elif lookup_ratio / math.log(10) < 10**1.5:
        jeffreys_support = "strong"
    elif lookup_ratio / math.log(10) < 10**2.0:
        jeffreys_support = "very strong"
    else:
        jeffreys_support = "decisive"

    if lookup_ratio < 1.0:
        kass_support = "not worth more than a bare mention"
    elif lookup_ratio < 3.0:
        kass_support = "positive"
    elif lookup_ratio < 5.0:
        kass_support = "strong"
    else:
        kass_support = "very strong"

    #print("Bayes factor is " + str(math.exp(log_likelihood_ratio)))
    #print("Hypothesis test favors hypothesis " + str(favored) + " and disfavors hypothesis " + str(unfavored))
    #print("Level of support is '" + jeffreys_support + "' on Jeffreys' scale and '" + kass_support + "' on Kass and Raferty's scale")

    return favored

def dirichlet_log_marginal_likelihood(counts_test, counts_background,
                                      pseudocounts):
    """
    Author:
    - Jordan M Eizenga (github@jeizenga)
    """

    #print(counts_test)
    #print(counts_background)
    codons = ["".join(codon) for codon in itertools.product(nts, repeat = 3)]
    log_like = 0
    log_like += math.lgamma(sum(counts_test[codon] for codon in codons))
    log_like += math.lgamma(sum(pseudocounts + counts_background[codon] for codon in codons))
    log_like -= math.lgamma(sum(pseudocounts + counts_background[codon] + counts_test[codon] for codon in codons))
    for codon in codons:
        log_like -= math.lgamma(pseudocounts + counts_test[codon])
        log_like -= math.lgamma(pseudocounts + counts_background[codon])
        log_like += math.lgamma(pseudocounts + counts_background[codon] + counts_test[codon])
    return log_like

def gen_codons(seq):
    """this generates codons as a list of strings. Each element of the list is
    a string of length 3

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    return [seq[i:i+3] for i in np.arange(0, len(seq), 3)]

def codon_dirichlet_log_likelihood_ratio(test_seq, seq_list_1, seq_list_2,
                                         pseudocounts = 0.5, verbose = True):
    """
    This method takes one test sequence and two lists of sequences representing
    the two hypothesis groups. This method determines the log likelihood
    probability, then log likelihood ratio of whether the trinucleotide frequency
    of the test sequence more closely matches the trinucleotide frequency of
    hypothesis sequence list 1 or hypothesis sequence list 2.

    Ambiguities will be skipped over

    Parameters:
     <test_seq>
     <seq_list_1>
     <seq_list_2>
     <pseudocounts>
     <verbose> -  This flag toggles whether test will be interpreted for you on
       stdout or not.

    Author:
    - Jordan M Eizenga (github@jeizenga)
    """
    # check for data invariants
    #  nt invariants in codons are not added later
    assert len(test_seq) % 3 == 0
    for seq in seq_list_1:
        assert len(seq) % 3 == 0
    for seq in seq_list_2:
        assert len(seq) % 3 == 0

    # count up codons
    codon_counts_1 = collections.Counter()
    codon_counts_2 = collections.Counter()
    codon_counts_test = collections.Counter()

    for codon in gen_codons(test_seq):
        codon_counts_test[codon] += 1

    for seq in seq_list_1:
        for codon in gen_codons(seq):
            add = True
            for nt in codon:
                if nt not in nts:
                    add = False
            if add:
                codon_counts_1[codon] += 1

    for seq in seq_list_2:
        for codon in gen_codons(seq):
            add = True
            for nt in codon:
                if nt not in nts:
                    add = False
            if add:
                codon_counts_2[codon] += 1

    # compute bayes factor
    log_likelihood_ratio = dirichlet_log_marginal_likelihood(codon_counts_test, codon_counts_1, pseudocounts) \
                           - dirichlet_log_marginal_likelihood(codon_counts_test, codon_counts_2, pseudocounts)

    # interpret as indicated
    if verbose:
        interpret_hypothesis_test(log_likelihood_ratio)

    return log_likelihood_ratio

def gen_noncodingseq_dict(fasta_directory):
    """This function takes a filepath of the fasta file that contains the
    noncoding sequences to use in the analysis. The method then
    generates all three artificial reading frames and returns in dict format.
    The keys are the sequence name and the values are a list of sequences.

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    seqdict = {}
    filelist = fasta_dir_to_gene_filelist(fasta_directory)

    for genename in filelist:
        seqs = []
        for record in SeqIO.parse(filelist[genename], "fasta"):
            seq = str(record.seq)
            for i in range(3):
                #this block creates sequences in all three frames
                frameshifted_seq = seq[i:]
                mod = len(frameshifted_seq) % 3
                trimmed_seq = frameshifted_seq[0:len(frameshifted_seq)-mod]
                # first make sure that the sequence is divisible by three, It should
                if len(trimmed_seq) % 3 != 0:
                    raise IOError("""One of the sequences ({}) in ({}) is not
                    divisible by three. Are there extra bases? Please only save
                    in-frame fasta sequences.""".format(
                        genename, filelist[genename]))
                nostops = remove_all_stops(trimmed_seq)
                no5primestarts = remove_5p_starts(nostops)
                nodeletions = remove_deletions(no5primestarts)
                # just to be pedantic, let's reaffirm that the sequence is divisible
                #  by three after some processing
                if len(nodeletions) % 3 != 0:
                    raise IOError("""One of the sequences ({}) in ({}) is not
                    divisible by three after filtering for 3' stop codons, 5'
                    start codons, and alignment-based deletions""".format(
                        genename, filelist[genename]))
                seqs.append(nodeletions)
            seqdict[genename] = seqs
    return seqdict

def gen_codingseqs_dict(fasta_directory):
    """This opens a directory of fasta files and reads them into a dictionary.
    This method does not recirsively search for fasta files, but only
    The keys are the sequence name and the values are a list of sequences

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    seqdict = {}
    filelist = fasta_dir_to_gene_filelist(fasta_directory)
    for genename in filelist:
        seqs = []
        for record in SeqIO.parse(filelist[genename], "fasta"):
            seq = str(record.seq)
            # first make sure that the sequence is divisible by three, It should
            if len(seq) % 3 != 0:
                raise IOError("""One of the sequences ({}) in ({}) is not
                divisible by three. Are there extra bases? Please only save
                in-frame fasta sequences.""".format(
                    genename, filelist[genename]))
            noendstops = remove_3p_stops(seq)
            no5primestarts = remove_5p_starts(noendstops)
            nodeletions = remove_deletions(no5primestarts)
            # just to be pedantic, let's reaffirm that the sequence is divisible
            #  by three after some processing
            if len(nodeletions) % 3 != 0:
                raise IOError("""One of the sequences ({}) in ({}) is not
                divisible by three after filtering for 3' stop codons, 5'
                start codons, and alignment-based deletions""".format(
                    genename, filelist[genename]))
            seqs.append(nodeletions)
        seqdict[genename] = seqs
    return seqdict

def remove_3p_stops(seq):
    """This function removes the 3' terminal stop codons from a DNA sequence
    formatted as a python string. It just returns the original sequence
    otherwise.

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    codons = [seq[i:i+3] for i in np.arange(0, len(seq), 3)]
    if codons[-1] in stops:
        new_codons = codons[0:-1]
    else:
        new_codons = codons
    return "".join([item for sublist in new_codons for item in sublist])

def remove_deletions(seq):
    """This function removes codons that contain deletions. For example, the
    following codons would be removed '---', 'A--', '-G-', 'T-A', et cetera.

    This is important in the analysis because deletions/gaps in a codon are
    often ambiguous for frameshifts or unknown bases. In this case, we assume
    that the user has provided a completely in-frame sequence and that any '-'
    characters are unknown bases. Including such bases in the analysis would
    interfere with the dirichlet distribution.

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    codons = [seq[i:i+3] for i in np.arange(0, len(seq), 3) if '-' not in seq[i:i+3]]
    return "".join([item for sublist in codons for item in sublist])

def remove_5p_starts(seq):
    """This function removes 5' start codons from sequences if they are present.
    Returns the original sequence otherwise.

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    codons = [seq[i:i+3] for i in np.arange(0, len(seq), 3)]
    if codons[0] in starts:
        new_codons = codons[1:]
    else:
        new_codons = codons
    return "".join([item for sublist in new_codons for item in sublist])

def remove_all_stops(seq):
    """This function removes all stop codons from a sequence.
    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    codons = [seq[i:i+3] for i in np.arange(0, len(seq), 3) if seq[i:i+3] not in stops]
    return "".join([item for sublist in codons for item in sublist])

def estimate_power_and_false_pos(log_likelihood_score, positive_simulations,
                                 negative_simulations):
    """
    This method estimates the power and false positive rate for a single test
    based on a set of simulated positive and negative results.

    As inputs, the method takes the log likelihood value of a test, a list of
    log likelihood values from leave-one-out positive simulations, and a
    list of log likelihood values from leave-one-out negative simulations.

    Parameters:
    <log_likelihood_score> -
    <positive_simulations> -
    <negative_simulations -

    Author:
    - Jordan M Eizenga (github@jeizenga)
    """
    false_pos = np.mean(np.greater(negative_simulations, log_likelihood_score))
    power = np.mean(np.greater(positive_simulations, log_likelihood_score))
    return false_pos, power

def estimate_power_and_false_pos_from_distr(log_likelihood_scores,
                                            positive_simulations,
                                            negative_simulations):
    """
    Author:
    - Jordan M Eizenga (github@jeizenga)
    """
    np.sort(log_likelihood_scores)
    np.sort(positive_simulations)
    np.sort(negative_simulations)

    pos_num = float(len(positive_simulations))
    neg_num = float(len(negative_simulations))

    j = 0
    k = 0

    total_power = 0.0
    total_false_pos = 0.0

    for i in range(len(log_likelihood_scores)):
        if j < len(positive_simulations):
            while positive_simulations[j] < log_likelihood_scores[i]:
                j += 1
                if j == len(positive_simulations):
                    break
        if k < len(negative_simulations):
            while negative_simulations[k] < log_likelihood_scores[i]:
                k += 1
                if k == len(negative_simulations):
                    break

        total_power += 1.0 - j / pos_num
        total_false_pos += 1.0 - k / neg_num

    power = total_power / len(log_likelihood_scores)
    false_pos = total_false_pos / len(log_likelihood_scores)
    return false_pos, power

def plot_ratios(results_df):
    """
    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    coding_llratios = np.array([llratio_list[i]
                           for i in range(len(true_value_list))
                           if true_value_list[i] == 1])
    nc_llratios = np.array([llratio_list[i]
                           for i in range(len(true_value_list))
                           if true_value_list[i] == 2])
    test_seq_sets = []
    for name in test_seq_names:
        these_llratios = np.load("{}_llratio.npy".format(name))
        test_seq_sets.append(these_llratios)
        false_pos, power = estimate_power_and_false_pos_from_distr(these_llratios,
                                                coding_llratios,
                                                nc_llratios)
        print("{}: power = {}, false_pos = {}".format(name, power, false_pos))

    fig, ax = plt.subplots()
    ax.hist(coding_llratios, bins='auto', label = 'coding',
            color='lightblue', alpha=0.5, normed=True)
    ax.hist(nc_llratios, bins='auto', label = 'noncoding',
            color='salmon', alpha=0.5, normed=True)
    ax.hist(test_seq_sets[0], bins='auto', label = 'unk1',
            color='purple', alpha=0.5, normed=True)
    ax.hist(test_seq_sets[1], bins='auto', label = 'nd2l',
            color='gold', alpha=0.5, normed=True)

    #create legend
    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in ["lightblue", "salmon", "purple", "gold"]]
    labels= ["coding", "nc", "unk1", "nd2l"]
    ax.legend(handles, labels)
    ax.set(title='Histogram Comparison', ylabel='% of Dataset in Bin')
    ax.margins(0.05)
    ax.set_ylim(bottom=0)
    plt.savefig("dirichlet_results_{}.png".format(timestamp()), dpi=600, transparent=False)

def LOO_analysis_chunk(args):
    """randomly performs one leave-one-out analysis

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    codingseqs_dict = args['codingseqs_dict']
    ncseqs_dict = args['ncseqs_dict']
    n_iterations = args['n_iterations']

    results = []
    for iteration in range(n_iterations):
        real_val = np.random.choice(['coding', 'noncoding'])
        if real_val == 'coding':
            # remove one of the coding genes and test it
            dict_keys = [x for x in codingseqs_dict.keys()]
            test_gene_name = np.random.choice(dict_keys)
            other_gene_names = [x for x in dict_keys if x != test_gene_name]
            random_test_seq = np.random.choice(codingseqs_dict[test_gene_name])
            coding_seqlist = [np.random.choice(codingseqs_dict[genename])
                              for genename in other_gene_names]
            nc_seqlist = [np.random.choice(ncseqs_dict[genename])
                          for genename in ncseqs_dict]
        elif real_val == 'noncoding':
            dict_keys = [x for x in ncseqs_dict.keys()]
            test_gene_name = np.random.choice(dict_keys)
            other_gene_names = [x for x in dict_keys if x != test_gene_name]
            random_test_seq = np.random.choice(ncseqs_dict[test_gene_name])
            coding_seqlist = [np.random.choice(codingseqs_dict[genename])
                              for genename in codingseqs_dict]
            nc_seqlist = [np.random.choice(ncseqs_dict[genename])
                          for genename in other_gene_names]

        log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(random_test_seq,
                                    coding_seqlist, nc_seqlist, pseudocounts = 0.1,
                                    verbose = False)

        if log_likelihood_ratio < 0:
            favored = 'noncoding'
        else:
            favored = 'coding'

        result = {"seqname": test_gene_name, "ll_ratio": log_likelihood_ratio,
                "analysis_type": "LOO", "real_val": real_val,
                "observed": favored}
        results.append(result)
    return results


def unknown_seq_analysis_chunk(args):
    """
    Author:
    - Darrin T Schultz (github@conchoecia)
    """

    testseqs_dict = args['testseqs_dict']
    codingseqs_dict = args['codingseqs_dict']
    ncseqs_dict = args['ncseqs_dict']
    n_iterations = args['n_iterations']

    # args should be: testseqs_dict, codingseqs_dict, ncseqs_dict
    # for one simulation, we randomly choose:
    #   - one of the seqs in the test seqlist
    #   - one of the seqs in each of the genes in the noncoding seqlist
    #   - one of the seqs in each of the genes in the coding seqlist

    # each entry of llratio is one analysis
    # {seqname: <test_gene_name>, ll_ratio: <log_likelihood_ratio>,
    #  analysis_type: 'test'}

    # This test returns a list of dictionaries of individual values

    testseqs_dict_keys = [x for x in testseqs_dict.keys()]
    results = []
    for iteration in range(n_iterations):
        test_gene_name = np.random.choice(testseqs_dict_keys)
        random_test_seq = np.random.choice(testseqs_dict[test_gene_name])
        coding_seqlist = [np.random.choice(val)
                          for seqname_key, val in codingseqs_dict.items()]
        nc_seqlist = [np.random.choice(val)
                          for seqname_key, val in ncseqs_dict.items()]

        log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(random_test_seq,
                                    coding_seqlist, nc_seqlist, pseudocounts = 0.1,
                                    verbose = False)

        if log_likelihood_ratio < 0:
               favored = 'noncoding'
        else:
               favored = 'coding'

        result = {"seqname": test_gene_name, "ll_ratio": log_likelihood_ratio,
                "analysis_type": "test", "real_val": None,
                "observed": favored}
        results.append(result)
    return results

def _calculate_combination_space(noncoding_dict, coding_dict,
                                 test_dict):
    """Calculates how many possible combinations of genes there are given
    the data.

    Edit: The combinatorics of the below equations is not correct. Fix those.

    For the total number of possible combinations for the LOO analyses:
      - For loci that are independent (nuclear genomes), the number of possible
         combinations is: ( <num of individuals ^ num coding loci> *\
                          <(num of individuals * 3) ^ num noncoding loci>)
      - For loci that are not independent (mito genomes), the number of
         possible combos is: <num of individuals * (3 ^ num noncoding loci)>
        - The reasoning behind the combinatorics for independent loci is that
           there are no ways to arrange the coding loci between individuals
           since they are not independent. For the noncoding loci we are able
           to make rearrangements between the pseudocodons in each individual,
           but not to make pseudocodon rearrangements between individuals.

    For the total number of possible combinations for the unknown test analyses:
      - For independent loci:
         <num of individuals ^ (num coding loci + 1)> * \
           <(num of individuals * 3) ^ (num noncoding loci)>
        - We simply add one to the first power since we are testing each
           possible combination of genes <num individuals> more times for
           a single test gene
      - For non-independent loci: <num of individuals * (3 ^ num noncoding loci)>
        - This is the same as the LOO analysis number of combinations since
           we have not added any possible way to make rearrangements when
           adding test loci
    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    results_dict = {}
    results_dict["num_individuals"]     = len(coding_dict[random.choice(list(coding_dict.keys()))])
    results_dict["num_coding_loci"]     = len(coding_dict)
    results_dict["num_nc_loci"]         = len(noncoding_dict)
    results_dict["num_test_loci"]       = len(test_dict)

    ni   = results_dict["num_individuals"]
    ncl  = results_dict["num_coding_loci"]
    nncl = results_dict["num_nc_loci"]
    ntl  = results_dict["num_test_loci"]

    results_dict["LOO_independent"]     = ((ni ** ncl) * ( (ni * 3) ** nncl ))
    results_dict["LOO_nonindependent"]  = ni * ( 3 ** nncl)
    results_dict["test_independent"]    = ( ni ** (ncl + 1) ) * ( (ni * 3) ** nncl )
    results_dict["test_nonindependent"] = results_dict["LOO_nonindependent"]

    return results_dict

def dirichlet(args):
    """
    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    outfile = sys.stderr
    print("{} - cuttlery dirichlet.".format(timestamp()), file = outfile)
    print("{} - reading options.".format(timestamp()), file = outfile)
    #First, read in the options
    global options
    options = args
    #print(options)
    ## FOR INTERNAL CONSISTENCY, ANY SET_1 IS CODING and SET_2 is NONCODING!
    results_file = options.results_file
    if results_file == None:
        raise IOError("""You must specify a name for a results file. If you already
        have results and just want to plot your data, pass the existing file to --results_file""")
    chunksize = 60 # this is the optimum chunk size for parallelization
    # If we've already done an analysis and the file is there already
    #  don't bother to do it again, but instead just plot
    seqsize = {}
    print("""{} - Couldn't find the results file, {}""".format(
        timestamp(), results_file), file = outfile)
    print("""{} - Starting analysis""".format(
        timestamp()), file = outfile)
    print("""{} - Looking for coding, noncoding, test fasta directories.""".format(
        timestamp()), file = outfile)
    # check that the test_dir, coding dir, and noncoding dir exist and that
    # they have fasta files inside
    look_at_these = [["noncoding_dir", options.noncoding_dir],
                     ["coding_dir", options.coding_dir]]
    if options.test_dir != None:
        look_at_these.append(["test_dir", options.test_dir])
    for dir_type, check_dir in look_at_these:
        if not os.path.isdir(check_dir):
            raise IOError("""Your results file {} didn't exist so cuttlery dirichlet looked
            for the fasta files for analysis. We couldn't find the {}: {}.
            Please make sure that this directory exists.""".format(
                results_file, dir_type, check_dir))
    ##first generate a list of individuals of the sequence to test
    print("""{} - Processing test sequences.""".format(timestamp()), file = outfile)
    if options.test_dir != None:
        testseqs_dict = gen_codingseqs_dict(options.test_dir)
    else:
        testseqs_dict = {}
    #Then get a list of the known sequences
    print("""{} - Processing known coding sequences.""".format(timestamp()), file = outfile)
    codingseqs_dict = gen_codingseqs_dict(options.coding_dir)
    print("""{} - Processing noncoding sequences and generating pseudocodons.""".format(timestamp()), file = outfile)
    ncseqs_dict = gen_noncodingseq_dict(options.noncoding_dir)

    print("""{} - Calculating the combination space.""".format(timestamp()), file = outfile)
    cspace = _calculate_combination_space(ncseqs_dict, codingseqs_dict, testseqs_dict)
    print("""{} - There are {:,} individuals.""".format(
        len(timestamp())*" ", cspace["num_individuals"]), file = outfile)
    print("""{} - There are {:,} coding loci.""".format(
        len(timestamp())*" ", cspace["num_coding_loci"]), file = outfile)
    print("""{} - There are {:,} noncoding loci.""".format(
        len(timestamp())*" ", cspace["num_nc_loci"]), file = outfile)
    print("""{} - There are {:,} test loci.""".format(
        len(timestamp())*" ", cspace["num_test_loci"]), file = outfile)
    print("""{} - There are {:,} combinations for LOO independent analyses.""".format(
        len(timestamp())*" ", cspace["LOO_independent"]), file = outfile)
    print("""{} - There are {:,} combinations for LOO nonindependent analyses.""".format(
        len(timestamp())*" ", cspace["LOO_nonindependent"]), file = outfile)
    print("""{} - There are {:,} combinations for test independent analyses.""".format(
        len(timestamp())*" ", cspace["test_independent"]), file = outfile)
    print("""{} - There are {:,} combinations for test nonindependent analyses.""".format(
        len(timestamp())*" ", cspace["test_nonindependent"]), file = outfile)

    # then get the list of lists of noncoding sequences.
    # We first have to pass the argument of how many copies of each gene there
    #  are so that the program adds the sequences in groups, this determines how many
    #  sequences to add in each sublist.
    results_dict_list = []
    seqs_dicts_args = {'testseqs_dict': testseqs_dict,
                       'codingseqs_dict': codingseqs_dict,
                       'ncseqs_dict': ncseqs_dict,
                       'n_iterations': chunksize}

    #get the mean sequence size
    for thisdict in [testseqs_dict, codingseqs_dict, ncseqs_dict]:
        for k in thisdict:
            seqsize[k] = np.mean([len(seq) for seq in thisdict[k]])

    if not os.path.exists(results_file):

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        #  This block performs the test analyses of the unknown seqs.
        #
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        print("""{} - Performing test_seqs analyses.""".format(timestamp()), file = outfile)
        num_simulations = int(options.numsims/chunksize) * len(testseqs_dict.keys())
        results = parallel_process([seqs_dicts_args for x in range(num_simulations)],
                                   unknown_seq_analysis_chunk, n_jobs = options.threads,
                                   use_kwargs = False, front_num=3)
        flat_results = [item for sublist in results for item in sublist]
        results_dict_list += flat_results

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        #  now perform the LOO analyses of all the known coding and nc seqs
        #
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        print("""{} - Performing leave-one-out analyses.""".format(timestamp()), file = outfile)
        num_gene_simulations = int(options.numsims/chunksize) * (len(codingseqs_dict.keys()) + len(ncseqs_dict.keys()))
        results = parallel_process([seqs_dicts_args for x in range(num_gene_simulations)],
                                   LOO_analysis_chunk, n_jobs = options.threads,
                                   use_kwargs = False, front_num=3)
        flat_results = [item for sublist in results for item in sublist]

        results_dict_list += flat_results
        print("""{} - Converting results to dataframe.""".format(timestamp()), file = outfile)
        results_df = pd.DataFrame.from_dict(results_dict_list)
        print("""{} - Saving results to {}.""".format(
            timestamp(), results_file), file = outfile)
        results_df.to_csv(results_file, index=False)
    else:
        print("""{} - Found the results file: {}.""".format(
            timestamp(), results_file), file = outfile)
        print("""{} - Loading the results file.""".format(
            timestamp()), file = outfile)
        results_df = pd.read_csv(results_file, index_col = False)

    print("""{} - Plotting with all genes separated.""".format(timestamp()), file = outfile)
    plot_results_simple(results_df, **vars(options))
    print("""{} - Plotting with all genes overlaid.""".format(timestamp()), file = outfile)
    plot_results(results_df, **vars(options))
    print("""{} - Plotting violinplots with locus size.""".format(timestamp()), file = outfile)
    plot_results_size(results_df, seqsize, **vars(options))

    # now output statistics about the data
    print("""{} - Calculating the confusion matrix.""".format(timestamp()), file = outfile)
    real_val = results_df.loc[results_df['analysis_type'] == 'LOO', 'real_val']
    observed_val = results_df.loc[results_df['analysis_type'] == 'LOO', 'observed']
    matrix = cmcalc(real_val, observed_val)
    print(matrix)

    # now calculate the false positive rate and the power of the simulation for
    #  each test data type
    print("""{} - Calculating false positive and power rate.""".format(timestamp()), file = outfile)
    tests = results_df[results_df['analysis_type'] == 'test']
    tests_seqnames = sorted(tests['seqname'].unique())
    for seqname in tests_seqnames:
        test_ll = tests.loc[tests['seqname'] == seqname, 'll_ratio'].as_matrix()
        noncodings = results_df.loc[results_df['real_val'] == 'noncoding', 'll_ratio'].as_matrix()
        codings = results_df.loc[results_df['real_val'] == 'coding', 'll_ratio'].as_matrix()
        false_pos, power = estimate_power_and_false_pos_from_distr(test_ll,
                                            codings,
                                            noncodings)
        print("{} marginalized false positive rate: {}".format(seqname, false_pos))
        print("{} power: {}".format(seqname, power))

def plot_results(results, **kwargs):
    """
    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    if kwargs.get("debug"):
        print(results.head())
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

    if kwargs.get("debug"):
        print(noncoding_seqnames)
        print(coding_seqnames)
        print(tests_seqnames)

    # autumn colors for noncoding
    cmap = cm.get_cmap('autumn')
    autumn = {noncoding_seqnames[i]: cmap(1 - ((i+1)/(len(noncoding_seqnames) + 1)))
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
    if kwargs.get("debug"):
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
    # Print image(s)
    print_images(
        base_output_name= kwargs["output_basename"] + "_histogram",
        image_formats=kwargs["fileform"],
        dpi=kwargs["dpi"],
        no_timestamp = kwargs["no_timestamp"],
        transparent=kwargs["transparent"])

def plot_results_simple(results, **kwargs):
    """
    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    #set the figure dimensions
    figWidth = 5
    figHeight = 4
    plt.figure(figsize=(figWidth,figHeight))
    #set the panel dimensions
    panelWidth = 4
    panelHeight = 2.5
    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2)
    panel0 =plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='off', labelleft='off', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panel0.spines['top'].set_visible(False)
    panel0.spines['right'].set_visible(False)
    panel0.spines['left'].set_visible(False)

    # Subset the dataframes to separate the noncoding seqs, et cetera.
    #  This makes them easier to process.
    noncodings = results[results['real_val'] == 'noncoding']
    codings = results[results['real_val'] == 'coding']
    tests = results[results['analysis_type'] == 'test']

    # Get the seqnames and organize somehow
    noncoding_seqnames = sorted(noncodings['seqname'].unique())
    coding_seqnames = sorted(codings['seqname'].unique())
    tests_seqnames = sorted(tests['seqname'].unique())
    if kwargs.get("sort_type") in ["meaninc", "meandec", "medinc", "meddec"]:
        # This block sorts the sequences based on the arguments to cuttlery dirichlet
        decreasing = False
        if kwargs.get("sort_type") in ["meandec","meddec"]:
            decreasing = True
        if kwargs.get("sort_type") in ["meaninc","meandec"]:
            sorton = "mean"
        elif kwargs.get("sort_type") in ["medinc","meddec"]:
            sorton = "median"
        noncoding_seqnames = _sorted_by_mean_ll(noncoding_seqnames,
                                                noncodings, sorton,
                                                decreasing = decreasing)
        coding_seqnames = _sorted_by_mean_ll(coding_seqnames,
                                             codings, sorton,
                                             decreasing = decreasing)
        tests_seqnames = _sorted_by_mean_ll(tests_seqnames,
                                            tests, sorton,
                                            decreasing = decreasing)

    all_seqnames = noncoding_seqnames + tests_seqnames + coding_seqnames

    # make a list of every gene's ll_ratios given the sort order we have just
    #  selected
    data_lists = [list(results.loc[results['seqname'] == seqname, 'll_ratio']) for
                  seqname in all_seqnames]

    if kwargs.get("debug"):
        print(noncoding_seqnames)
        print(coding_seqnames)
        print(tests_seqnames)


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # This block handles selecting the colors of all of the violinplots
    #
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # color the noncodings warm by default or according to argparse
    if kwargs.get("color_nc"):
        nc_colordict = {x: kwargs.get("color_nc") for x in noncoding_seqnames }
    else:
        cmap = cm.get_cmap('autumn')
        nc_colordict = {noncoding_seqnames[i]: cmap(1 - ((i + 1)/(len(noncoding_seqnames) +1) ))
                  for i in range(len(noncoding_seqnames)) }
    # color the codings cool by default or according to argparse
    if kwargs.get("color_coding"):
        coding_colordict = {x: kwargs.get("color_coding") for x in coding_seqnames }
    else:
        cmap = cm.get_cmap('winter')
        coding_colordict = {coding_seqnames[i]: cmap(1 - (i/len(coding_seqnames)))
                  for i in range(len(coding_seqnames)) }
    # magma colors for tests
    # color the test sequences dark by default or according to argparse
    if kwargs.get("color_test"):
        test_colordict = {x: kwargs.get("color_test") for x in tests_seqnames }
    else:
        cmap = cm.get_cmap('magma')
        test_colordict = {tests_seqnames[i]: cmap(((i/len(coding_seqnames))) )
                          for i in range(len(tests_seqnames)) }
    ## We need to plot each of the sections individually to accommodate the gaps
    ## Put all the colors into a list with the same indices as all_seqnames
    noncoding_colors = [nc_colordict[noncoding_seqname]
                        for noncoding_seqname in noncoding_seqnames]
    coding_colors    = [coding_colordict[coding_seqname]
                        for coding_seqname in coding_seqnames]
    tests_colors     = [test_colordict[tests_seqname]
                        for tests_seqname in tests_seqnames]
    all_colors = noncoding_colors + tests_colors + coding_colors
    if kwargs.get("debug"):
        print(all_colors)

    plt.axvline(x=0, color=(0,0,0,0.5), linewidth=0.5)

    positions = []
    counter = -1
    for this_list in [noncoding_seqnames, tests_seqnames, coding_seqnames]:
        counter += 1
        positions = positions + list(np.arange(counter, counter+len(this_list), 1))
        counter += len(this_list)

    nc=panel0.violinplot(data_lists,
                  positions=positions,
                  widths=0.8, vert = False,
                  showmeans=False,showmedians=False,showextrema=False,
                  bw_method=0.05)
    for i in range(len(nc['bodies'])):
        nc['bodies'][i].set_facecolor(all_colors[i])
        nc['bodies'][i].set_alpha(0.75)
        nc['bodies'][i].set_linewidth(0)

    #panel0.set_ylim([len(data_lists), -1 ])
    plt.gca().invert_yaxis()
    xmin  = results['ll_ratio'].min() * 1.1
    xmax = results['ll_ratio'].max() * 1.1
    ax_width = abs(xmax - xmin)
    panel0.set_xlim(xmin, xmax)

    smart_ticks = True
    if smart_ticks:
        # first plot text for noncodings
        # find the max ll_ratio in noncodings and left align all labels there
        noncoding_max = noncodings['ll_ratio'].max() * 1.05
        pos_ix = 0
        noncodingtext_xmaxes = []
        noncodingtext_ymaxes = []
        noncodingtext_ymins = []
        for i in range(len(noncoding_seqnames)):
            tx = panel0.text(noncoding_max, positions[pos_ix], all_seqnames[pos_ix],
                        verticalalignment='center',
                        horizontalalignment='left',
                        color='black', fontsize=8)
            plt.draw()
            txcoords = tx.get_window_extent().transformed(panel0.transData.inverted())
            noncodingtext_xmaxes.append(txcoords.xmax)
            noncodingtext_ymaxes.append(txcoords.ymax)
            noncodingtext_ymins.append(txcoords.ymin)
            pos_ix += 1
        # Now draw a bounding box for the noncoding region
        plot_left_facing_bracket(panel0, noncodingtext_xmaxes, noncodingtext_ymaxes,
                                 noncodingtext_ymins, ax_width, "noncoding")

        # Put the test sequences on individually
        for i in range(len(tests_seqnames)):
            if kwargs.get("debug"):
                print(coding_seqnames[i])
            thismin = tests.loc[tests['seqname'] == tests_seqnames[i], 'll_ratio'].min()
            thismax = tests.loc[tests['seqname'] == tests_seqnames[i], 'll_ratio'].max()
            leftdif = abs(xmax - thismax)
            rightdif = abs(thismin - xmin)
            if leftdif <= rightdif:
                panel0.text(thismin * 0.95, positions[pos_ix], all_seqnames[pos_ix],
                        verticalalignment='center',
                        horizontalalignment='right',
                        color='black', fontsize=8)
            else:
                panel0.text(thismax * 1.05, positions[pos_ix], all_seqnames[pos_ix],
                        verticalalignment='center',
                        horizontalalignment='left',
                        color='black', fontsize=8)
            pos_ix += 1
        # and finish up with the coding seqs
        codingtext_xmins = []
        codingtext_ymaxes = []
        codingtext_ymins = []
        coding_min = codings['ll_ratio'].min()
        if coding_min < 0:
            coding_min *= 1.05
        else:
            coding_min *= 0.95
        for i in range(len(coding_seqnames)):
            tx = panel0.text(coding_min, positions[pos_ix], all_seqnames[pos_ix],
                        verticalalignment='center',
                        horizontalalignment='right',
                        color='black', fontsize=8)
            plt.draw()
            txcoords = tx.get_window_extent().transformed(panel0.transData.inverted())
            codingtext_xmins.append(txcoords.xmin)
            codingtext_ymaxes.append(txcoords.ymax)
            codingtext_ymins.append(txcoords.ymin)
            pos_ix += 1
        # Now draw a bounding box for the coding region
        plot_right_facing_bracket(panel0, codingtext_xmins, codingtext_ymaxes,
                                 codingtext_ymins, ax_width, "coding")

    else:
        # this is here because I need to align the tick labels
        # https://stackoverflow.com/questions/15882249/
        plt.draw()
        panel0.set_yticks(positions)
        panel0.set_yticklabels(all_seqnames)
        # this stuff is to align the tick labels but I think it looks OK now
        #yax = panel0.get_yaxis()
        ## find the maximum width of the label on the major ticks
        #pad = max(T.label.get_window_extent().width for T in yax.majorTicks)
        #yax.set_tick_params(pad=pad/3)

    panel0.set_xlabel("Log-likelihood ratio")
    panel0.set_title("Codon Usage Log-likelihood Ratios")

    #plt.show()
    # Print image(s)
    print_images(
        base_output_name= kwargs["output_basename"] + "_violins",
        image_formats=kwargs["fileform"],
        dpi=kwargs["dpi"],
        no_timestamp = kwargs["no_timestamp"],
        transparent=kwargs["transparent"])

def _sorted_by_mean_ll(genenames, df, sorton, **kwargs):
    """This method returns a list of gene names organized by the mean value
    of log-likelihoods of all simulations for that gene.

    As input, this method takes an iterable of gene names and a pandas
    dataframe containing the following information (for example).

      analysis_type   ll_ratio observed real_val seqname
    0          test  29.198217   coding      NaN     UNK
    1          test  25.050107   coding      NaN    ND2L
    2          test  31.772205   coding      NaN     UNK
    3          test  21.590106   coding      NaN    ND2L
    4          test  33.697691   coding      NaN     UNK

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    try:
        kwargs["decreasing"]
    except:
        raise Exception("""You need to specify if the sort is
                        increasing or decreasing""")
    gene_mean_dict = {x:getattr(df.loc[df['seqname'] == x, 'll_ratio'],sorton)()
                     for x in genenames}
    sorted_values = sorted(gene_mean_dict, key = gene_mean_dict.get,
                           reverse = kwargs.get("decreasing"))
    return sorted_values

def plot_results_size(results, seqsize, **kwargs):
    """
    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    #set the figure dimensions
    figWidth = 5
    figHeight = 4
    plt.figure(figsize=(figWidth,figHeight))
    #set the panel dimensions
    panelWidth = 4
    panelHeight = 2.5
    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2)
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
    #panel0.spines['left'].set_visible(False)

    # Get the seqnames and organize somehow
    seqnames = sorted(results['seqname'].unique())

    # Subset the dataframes to separate the noncoding seqs, et cetera.
    #  This makes them easier to process.
    noncodings = results[results['real_val'] == 'noncoding']
    codings = results[results['real_val'] == 'coding']
    tests = results[results['analysis_type'] == 'test']

    # Get the seqnames and organize somehow
    noncoding_seqnames = sorted(noncodings['seqname'].unique())
    coding_seqnames = sorted(codings['seqname'].unique())
    tests_seqnames = sorted(tests['seqname'].unique())
    colortypes = {"noncoding": "red",
                  "coding": "blue",
                  "test": "black"}

    seqname_to_type = {}
    for thisname in noncoding_seqnames:
        seqname_to_type[thisname] = "noncoding"
    for thisname in coding_seqnames:
        seqname_to_type[thisname] = "coding"
    for thisname in tests_seqnames:
        seqname_to_type[thisname] = "test"

    # make a list of every gene's ll_ratios given the sort order we have just
    #  selected
    data_lists = [list(results.loc[results['seqname'] == seqname, 'll_ratio']) for
                  seqname in seqnames]
    data_lists_to_color_type = [colortypes[seqname_to_type[seqname]] for
                  seqname in seqnames]

    if kwargs.get("debug"):
        print(seqnames)

    positions = []
    print(seqnames)
    print(seqsize)
    for thisname in seqnames:
        positions.append(seqsize[thisname])
    xmax = np.max([seqsize[thisseq] for thisseq in seqnames])

    nc=panel0.violinplot(data_lists,
                  positions=positions,
                  widths=xmax/80, vert = True,
                  showmeans=False,showmedians=False,showextrema=False,
                  bw_method=0.05)
    for i in range(len(nc['bodies'])):
        nc['bodies'][i].set_facecolor(data_lists_to_color_type[i])
        nc['bodies'][i].set_alpha(0.75)
        nc['bodies'][i].set_linewidth(0)

    # this plots the seqnames above the text
    for i in range(len(seqnames)):
        text =  seqnames[i]
        x = seqsize[text]
        y = np.mean(data_lists[i])
        text = "     " + text
        panel0.text(x, y, text, fontsize=4,
                    withdash=False,
                    ha='left', va='center',
                    color=(0,0,0,0.75))

    # LINEAR FIT FOR CODING POINTS
    xs = []
    ys = []
    for i in range(len(seqnames)):
        if seqnames[i] in coding_seqnames:
            xs = xs + [seqsize[seqnames[i]]]*len(data_lists[i])
            ys = ys + data_lists[i]
    c_slope, c_intercept, c_r_value, c_p_value, c_std_err = stats.linregress(xs,ys)
    print("coding slope: {}".format(c_slope))
    print("coding intercept: {}".format(c_intercept))
    print("coding r_value: {}".format(c_r_value))
    print("coding p_value: {}".format(c_p_value))
    print("coding std_err: {}".format(c_std_err))
    print()
    xvals = np.array(panel0.get_xlim())
    xvals[0] = 0
    yvals = c_intercept + c_slope * xvals
    panel0.plot(xvals, yvals, '--', color = "blue", lw = 1)

    # LINEAR FIT FOR NONCODING POINTS
    xs = []
    ys = []
    for i in range(len(seqnames)):
        if seqnames[i] in noncoding_seqnames:
            xs = xs + [seqsize[seqnames[i]]]*len(data_lists[i])
            ys = ys + data_lists[i]
    nc_slope, nc_intercept, nc_r_value, nc_p_value, nc_std_err = stats.linregress(xs,ys)
    print("noncoding slope: {}".format(nc_slope))
    print("noncoding intercept: {}".format(nc_intercept))
    print("noncoding r_value: {}".format(nc_r_value))
    print("noncoding p_value: {}".format(nc_p_value))
    print("noncoding std_err: {}".format(nc_std_err))
    print()
    xvals = np.array(panel0.get_xlim())
    xvals[0] = 0
    yvals = nc_intercept + nc_slope * xvals
    panel0.plot(xvals, yvals, '--', color = "red", lw = 1)

    # LINEAR FIT FOR test points
    #if options.test_dir != None and len(tests_seqnames) > 1:
    #    xs = []
    #    ys = []
    #    for i in range(len(seqnames)):
    #        if seqnames[i] in tests_seqnames:
    #            xs = xs + [seqsize[seqnames[i]]]*len(data_lists[i])
    #            ys = ys + data_lists[i]
    #    slope, intercept, r_value, p_value, std_err = stats.linregress(xs,ys)
    #    print("slope: {}".format(slope))
    #    print("intercept: {}".format(intercept))
    #    print("r_value: {}".format(r_value))
    #    print("p_value: {}".format(p_value))
    #    print("std_err: {}".format(std_err))
    #    xvals = np.array(panel0.get_xlim())
    #    xvals[0] = 0
    #    yvals = intercept + slope * xvals
    #    panel0.plot(xvals, yvals, '--', color = "black", lw = 1)

    #panel0.set_ylim([len(data_lists), -1 ])
    #plt.gca().invert_yaxis()
    ymin  = results['ll_ratio'].min() * 1.1
    ymax = results['ll_ratio'].max() * 1.1
    #ax_width = abs(xmax - xmin)
    panel0.set_ylim(ymin, ymax)
    panel0.set_ylabel("Log-likelihood ratio")
    panel0.set_xlabel("Locus size")
    panel0.set_title("Codon Usage Log-likelihood Ratios")

    #plt.show()
    # Print image(s)
    print_images(
        base_output_name= kwargs["output_basename"] + "_violins_size",
        image_formats=kwargs["fileform"],
        dpi=kwargs["dpi"],
        no_timestamp = kwargs["no_timestamp"],
        transparent=kwargs["transparent"])

    # Now we calculate the percent of the LLs above the line
    print("seq\tper_above_coding_linear\tper_below_noncoding_linear")
    for i in range(len(seqnames)):
        thisname = seqnames[i]
        seqlen = seqsize[thisname]
        coding_val = c_intercept + c_slope * seqlen
        noncoding_val = nc_intercept + nc_slope * seqlen
        per_above_coding = len([val for val in data_lists[i] if val >= coding_val])/len(data_lists[i])*100
        per_below_ncoding = len([val for val in data_lists[i] if val <= noncoding_val])/len(data_lists[i])*100
        print("{}\t{}\t{}".format(thisname, per_above_coding, per_below_ncoding))
    print()

def plot_left_facing_bracket(panel, xmaxes, ymaxes, ymins, ax_width, text):
    """
    Author:
    - Darrin T Schultz (github@conchoecia)
    """

    nct_xmax = max(xmaxes)
    nct_ymax = max(ymaxes)
    nct_ymin = min(ymins)
    y_offset = 0.05

    topleftx = nct_xmax + (ax_width * 0.01)
    toplefty = nct_ymin - (y_offset * 2)
    toprightx = nct_xmax + (ax_width * 0.02)
    toprighty = nct_ymin - (y_offset * 2)
    bottomrightx = nct_xmax + (ax_width * 0.02)
    bottomrighty = nct_ymax + y_offset
    bottomleftx  = nct_xmax + (ax_width * 0.01)
    bottomlefty = nct_ymax  + y_offset
    #top part
    panel.plot([topleftx, toprightx],\
                [toplefty, toprighty],\
                color='k', linestyle='-',\
                linewidth = 0.5)
    #vertical part
    panel.plot([toprightx, bottomrightx],\
                [toprighty, bottomrighty],\
                color='k', linestyle='-',\
                linewidth = 0.5)
    #bottom part
    panel.plot([bottomleftx, bottomrightx],\
                [bottomlefty, bottomrighty],\
                color='k', linestyle='-',\
                linewidth = 0.5)
    tx = panel.text(bottomrightx + nct_xmax * 0.03,\
                     toprighty - ((toprighty - bottomrighty)/2),\
                    "noncoding",
                    verticalalignment='center',
                    horizontalalignment='left',
                    rotation = 90,
                    color='black', fontsize=10)

def plot_right_facing_bracket(panel, xmins, ymaxes, ymins, ax_width, text):
    """
    Author:
    - Darrin T Schultz (github@conchoecia)
    """

    nct_xmin = min(xmins)
    nct_ymax = max(ymaxes)
    nct_ymin = min(ymins)
    y_offset = 0.05

    topleftx = nct_xmin - (ax_width * 0.02)
    toplefty = nct_ymin - (y_offset * 2)
    toprightx = nct_xmin - (ax_width * 0.01)
    toprighty = nct_ymin - (y_offset * 2)
    bottomrightx = nct_xmin - (ax_width * 0.01)
    bottomrighty = nct_ymax + y_offset
    bottomleftx  = nct_xmin - (ax_width * 0.02)
    bottomlefty = nct_ymax + y_offset
    #top part
    panel.plot([topleftx, toprightx],\
                [toplefty, toprighty],\
                color='k', linestyle='-',\
                linewidth = 0.5)
    #vertical part
    panel.plot([bottomleftx, topleftx],\
                [toplefty, bottomlefty],\
                color='k', linestyle='-',\
                linewidth = 0.5)
    #bottom part
    panel.plot([bottomleftx, bottomrightx],\
                [bottomlefty, bottomrighty],\
                color='k', linestyle='-',\
                linewidth = 0.5)
    tx = panel.text(bottomleftx - (ax_width * 0.01),\
                     toprighty - ((toprighty - bottomrighty)/2),\
                    "coding",
                    verticalalignment='center',
                    horizontalalignment='right',
                    rotation = 90,
                    color='black', fontsize=10)

def run(args):
    dirichlet(args)
