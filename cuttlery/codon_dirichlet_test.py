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
conception: Jordan M Eizenga (and a tiny tiny bit DTS)

This program observes the nucleotide diversity of a protein-coding
locus in a population of sequences. Using this information it
generates other similar sequences with the same nucleotide diversity
and geneaology as the observed sequences. However, the generated
sequences have their mutation sites modeled from a randomly-chosen
observed sequence at randomly chosen sites. This simulation may
estimate the lower bounds of piN/piS for a neutrally evolving sequence
of similar base composition and geneology as the observed population
of protein-coding sequences.
"""

#
#@author: Jordan/DTS
#info about confusion matrix here: http://www.dataschool.io/simple-guide-to-confusion-matrix-terminology/


import argparse
import collections
import math
import itertools
import sys
from Bio import SeqIO
from builtins import range
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
from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
import time

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

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

#This class is used in argparse to expand the ~. This avoids errors caused on
# some systems.
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

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

def dirichlet_log_marginal_likelihood(counts_test, counts_background, pseudocounts):
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
    a string of length 3"""
    return [seq[i:i+3] for i in np.arange(0, len(seq), 3)]

# arguments are one test sequence and two lists of sequences representing the two hypothesis groups
# pseudocounts are added to each codon
# verbose flag toggles whether test will be interpreted for you on stdout
def codon_dirichlet_log_likelihood_ratio(test_seq, seq_list_1, seq_list_2,
                                         pseudocounts = 0.5, verbose = True):
    # check for data invariants
    assert len(test_seq) % 3 == 0
    for nt in test_seq:
        assert nt in nts
    for seq in seq_list_1:
        assert len(seq) % 3 == 0
        for nt in seq:
            assert nt in nts
    for seq in seq_list_2:
        assert len(seq) % 3 == 0
        for nt in seq:
            assert nt in nts

    # count up codons
    codon_counts_1 = collections.Counter()
    codon_counts_2 = collections.Counter()
    codon_counts_test = collections.Counter()

    for codon in gen_codons(test_seq):
        codon_counts_test[codon] += 1

    for seq in seq_list_1:
        for codon in gen_codons(seq):
            codon_counts_1[codon] += 1

    for seq in seq_list_2:
        for codon in gen_codons(seq):
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
    noncoding sequences to use in the analysis.  Instead of manually
    generating all three artificial reading frames, the function
    generates them automatically and returns in dict format.
    The keys are the sequence name and the values are a list of sequences.
    """
    seqdict = {}
    filelist = {os.path.splitext(os.path.basename(x))[0]:os.path.join(os.path.abspath(fasta_directory), x)
                for x in os.listdir(fasta_directory)
                if os.path.splitext(os.path.basename(x))[1] == ".fasta"}

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
    The keys are the sequence name and the values are a list of sequences"""
    #1.5 get a list of files from the directory we provided.
    # This is a dict object with key as
    seqdict = {}
    filelist = {os.path.splitext(os.path.basename(x))[0]:os.path.join(os.path.abspath(fasta_directory), x)
                for x in os.listdir(fasta_directory)
                if os.path.splitext(os.path.basename(x))[1] == ".fasta"}

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
    formatted as a python string. Just returns the original sequence otherwise
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
    """
    codons = [seq[i:i+3] for i in np.arange(0, len(seq), 3) if '-' not in seq[i:i+3]]
    return "".join([item for sublist in codons for item in sublist])

def remove_5p_starts(seq):
    """This function removes 5' start codons from sequences if they are present.
    Returns the original sequence otherwise.
    """
    codons = [seq[i:i+3] for i in np.arange(0, len(seq), 3)]
    if codons[0] in starts:
        new_codons = codons[1:]
    else:
        new_codons = codons
    return "".join([item for sublist in new_codons for item in sublist])

def remove_all_stops(seq):
    """This function removes all stop codons from a sequence.
    """
    codons = [seq[i:i+3] for i in np.arange(0, len(seq), 3) if seq[i:i+3] not in stops]
    return "".join([item for sublist in codons for item in sublist])

# estimates the power and false positive rate for the test based on a set of
# simulated positive and negative results
# takes the log likelihood of a test, a list of log likelihoods from leave-one-out
# positive simulations, and a list of log likelihoods leave-one-out negative
# simulations as inputs
def estimate_power_and_false_pos(log_likelihood_score, positive_simulations,
                                 negative_simulations):

    false_pos = np.mean(np.greater(negative_simulations, log_likelihood_score))
    power = np.mean(np.greater(positive_simulations, log_likelihood_score))
    return false_pos, power

def estimate_power_and_false_pos_from_distr(log_likelihood_scores,
                                            positive_simulations,
                                            negative_simulations):
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

def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")

def LOO_analysis_chunk(args):
    """randomly performs one leave-one-out analysis"""
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
            coding_seqlist = [np.random.choice(ncseqs_dict[genename])
                              for genename in other_gene_names]
            nc_seqlist = [np.random.choice(codingseqs_dict[genename])
                          for genename in codingseqs_dict]

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

def dirichlet(args):
    #First, read in the options
    global options
    options = args
    ## FOR INTERNAL CONSISTENCY, ANY SET_1 IS CODING and SET_2 is NONCODING!
    results_file = options.results_file
    if results_file == None:
        raise IOError("""You must specify a name for a results file. If you already
        have results and just want to plot your data, pass the existing file to --results_file""")
    chunksize = 60
    # If we've already done an analysis and the file is there already
    #  don't bother to do it again, but instead just plot
    if not os.path.exists(results_file):
        # check that the test_dir, coding dir, and noncoding dir exist and that
        # they have fasta files inside
        for dir_type, check_dir in (("noncoding_dir", options.noncoding_dir),
                          ("test_dir", options.test_dir),
                          ("coding_dir", options.coding_dir)):
            if not os.path.isdir(check_dir):
                raise IOError("""Your results file {} didn't exist so cuttlery dirichlet looked
                for the fasta files for analysis. We couldn't find the {}: {}.
                Please make sure that this directory exists.""".format(
                    results_file, dir_type, check_dir))
           ## In the future I should implement something to check that there are
           ##  some fasta files in this directory. Not sure if it is necessary right now though. 
        ##first generate a list of individuals of the sequence to test
        testseqs_dict = gen_codingseqs_dict(options.test_dir)
        #Then get a list of the known sequences
        codingseqs_dict = gen_codingseqs_dict(options.coding_dir)
        ncseqs_dict = gen_noncodingseq_dict(options.noncoding_dir)
        # then get the list of lists of noncoding sequences.
        # We first have to pass the argument of how many copies of each gene there
        #  are so that the program adds the sequences in groups, this determines how many
        #  sequences to add in each sublist.
        results_dict_list = []
        seqs_dicts_args = {'testseqs_dict': testseqs_dict,
                           'codingseqs_dict': codingseqs_dict,
                           'ncseqs_dict': ncseqs_dict,
                           'n_iterations': chunksize}
        # perform the test analyses of the unknown seqs
        num_simulations = int(options.numsims/chunksize) * len(testseqs_dict.keys())
        print("please wait, preparing test_seqs analyses")
        results = parallel_process([seqs_dicts_args for x in range(num_simulations)],
                                   unknown_seq_analysis_chunk, n_jobs = options.threads,
                                   use_kwargs = False, front_num=3)
        flat_results = [item for sublist in results for item in sublist]
        results_dict_list += flat_results

        num_gene_simulations = int(options.numsims/chunksize) * (len(codingseqs_dict.keys()) + len(ncseqs_dict.keys()))
        ## now perform the LOO analyses of all the known coding and nc seqs
        print("please wait, preparing loo analyses")
        results = parallel_process([seqs_dicts_args for x in range(num_gene_simulations)],
                                   LOO_analysis_chunk, n_jobs = options.threads,
                                   use_kwargs = False, front_num=3)
        flat_results = [item for sublist in results for item in sublist]

        results_dict_list += flat_results
        results_df = pd.DataFrame.from_dict(results_dict_list)
        print(results_df)
        results_df.to_csv(results_file, index=False)
    else:
        print("found {} so skipping analysis".format(results_file))
        results_df = pd.read_csv(results_file, index_col = False)

    if options.plot:
        #plot_results(results_df)
        plot_results_simple(results_df)


    # now output statistics about the data
    real_val = results_df.loc[results_df['analysis_type'] == 'LOO', 'real_val']
    observed_val = results_df.loc[results_df['analysis_type'] == 'LOO', 'observed']
    matrix = cmcalc(real_val, observed_val)
    print("confusion matrix")
    print (matrix)

    # now calculate the false positive rate and the power of the simulation for
    #  each test data type
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

def plot_results_simple(results):
    plt.style.use('BME163')
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

    noncodings = results[results['real_val'] == 'noncoding']
    codings = results[results['real_val'] == 'coding']
    tests = results[results['analysis_type'] == 'test']

    noncoding_seqnames = sorted(noncodings['seqname'].unique())
    coding_seqnames = sorted(codings['seqname'].unique())
    tests_seqnames = sorted(tests['seqname'].unique())
    all_seqnames = noncoding_seqnames + tests_seqnames + coding_seqnames

    data_lists = [list(results.loc[results['seqname'] == seqname, 'll_ratio']) for
                  seqname in all_seqnames]

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

    ## We need to plot each of the sections individually to accommodate the gaps
    ## Put all the colors into a list with the same indices as all_seqnames
    noncoding_colors = [autumn[noncoding_seqname] for noncoding_seqname in noncoding_seqnames]
    coding_colors = [winter[coding_seqname] for coding_seqname in coding_seqnames]
    tests_colors = [magma[tests_seqname] for tests_seqname in tests_seqnames]
    all_colors = noncoding_colors + tests_colors + coding_colors
    print(all_colors)

    plt.axvline(x=0, color=(0,0,0,0.5), linewidth=0.5)


    positions = []
    counter = -1
    for this_list in [noncoding_seqnames, tests_seqnames, coding_seqnames]:
        counter += 1
        positions = positions + list(np.arange(counter, counter+len(this_list), 1))
        counter += len(this_list)

    nc=panel0.violinplot(data_lists, \
                  positions=positions, \
                  widths=0.5, vert = False,
                  showmeans=False,showmedians=False,showextrema=False,
                  bw_method=0.05)
    for i in range(len(nc['bodies'])):
        nc['bodies'][i].set_facecolor(all_colors[i])
        nc['bodies'][i].set_alpha(0.75)
        nc['bodies'][i].set_edgecolor('black')

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
    plt.savefig("dirichlet_violins.png", dpi=600, transparent=True)

    #plt.savefig("dirichlet_histogram_{}.png".format(timestamp()), dpi=600, transparent=False)

def plot_left_facing_bracket(panel, xmaxes, ymaxes, ymins, ax_width, text):
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

def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")

def run(args):
    dirichlet(args)
