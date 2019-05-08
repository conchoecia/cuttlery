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
"""
title: cuttlery piNpiSsim
author: Darrin T Schultz - github@conchoecia
conception: Russell Corbett-Detig
thanks to: Jordan Eizenga (github@jeizenga)

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

#plotting stuff
# start with matplotlib since other libraries may load it
import matplotlib
#had to do this do get it to work on mac https://stackoverflow.com/questions/30138026
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib import cm
import matplotlib.patches as mplpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path


#math
import numpy as np
from numpy import cumsum
from numpy.random import rand
import scipy.stats as stats
import pandas as pd
import time

# functools because apply doesn't work with two arguments
from functools import partial
from itertools import permutations

# cuttlery stuff
from cuttlery.codonfunctions import fasta_dir_to_gene_filelist, fasta_path_to_codonseqs,\
    print_images, seqs_to_df, calculate_pi, seqfreqs, calculate_piN_piS

import urllib.request
import scipy.stats
import sys, os
import copy
from Bio import SeqIO
from Bio.codonalign import CodonSeq
import Bio.Data.CodonTable
from Bio.codonalign.codonalphabet import get_codon_alphabet
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import _get_pi

# Also import some stuff to parallelize the function
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

#set font to helvetica
global hfont
hfont = {'fontname':'Helvetica'}

from matplotlib import rc
rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}',    # helvetica font
        r'\usepackage{sansmath}',   # math-font matching  helvetica
        r'\sansmath'                # actually tell tex to use it!
        r'\usepackage{siunitx}',    # micro symbols
        r'\sisetup{detect-all}',    # force siunitx to use the fonts
        ]

# don't print out options
pd.options.mode.chained_assignment = None

# This program needs to do a few things
# 0. Select a codon usage table
# 1. Read in a fasta alignment
# 2. Determine all of the polymorphic sites in the alignment
# 3. Calculate pi, piN, and piS for the alignment


def random_sequence(seqs, mode = 'dominant'):
    """Generates a random consensus sequence given the probabilities of
       observing that base
    - The dominant mode will pick the most frequent allele and randomly
      select a minor allele if all minor alleles have the same
      probability (only for n<5)
    - The random mode selects a base given the probability of observing that
      base at that locus
    """
    df = seqs_to_df(seqs)
    consensus = df.iloc[:,np.random.choice(list(range(len(seqs))))]
    return consensus

# This isn't used in the current version of the script
def _num_mutations(row):
    unique = np.unique(row)
    num_seqs = len(row)
    num_mutations = len(unique) - 1
    if 'nummutations' in options.debug:
        if num_mutations > 0:
            print("mutations: {}".format(num_mutations))

def _edit_distance(seq, consensus):
    """Determines the edit distance between two sequences, gaps included"""
    edit_distance = np.sum(consensus != seq)
    if 'editdistance' in options.debug:
        print("consensus: ", "".join(consensus[1:10]))
        print("      seq: ", "".join(seq[1:10]))
        print("edit dist: ", edit_distance)
    return edit_distance

def _binary_comparison(seq, reference):
    """This outputs a binary column that, where True, indicates when a sequence
    and a reference do not share the same base"""
    return reference != seq

def _row_profile(row):
    """returns a list of row positions and the mutational profile.
    0 = is always used for the 0th position, then every other base
         in [G,A,T,C] is assigned an index
    """
    bases = ['G', 'A', 'T', 'C']
    zeroth = row[0]
    remaining = [x for x in bases if x != zeroth]
    replace_dict = {remaining[i]:i+1 for i in range(len(remaining))}
    replace_dict[zeroth] = 0
    #print(replace_dict)
    row.replace(replace_dict, inplace = True)
    #print(row)
    return row

def get_mutation_profile(consensus, seqs):
    """Input is a list of CodonSeq objects, output is a tuple.
    element 0: number of input seqs
    element 1: numpy array of num_seqs columns (N) by num_mutation sites (M).
      A 0 in a cell indicates that there shall be no mutation there, but a 1
      indicates that this position should be randomly mutated.

    For example, in the alignment:
    GATA
    GTTA
    GCCA
    0210

    ...the output of this function will be (3, [2, 1])

    To generate the mutation profile, we first compare the consensus to all of
    the seqs and determine which has the smallest edit distance. We will then
    use the lowest-edit distance sequence, not the consensus, to build the
    mutation profile. We do not use the consensus because it may contain alleles
    from multiple individuals and would obscure the phylogenetic signal.
    """
    df = seqs_to_df(seqs)
    cols = list(range(len(df.columns)))
    np.random.shuffle(cols)
    # now we shuffle the df to make the reference the 0th element
    newdf = df.iloc[:, cols]

    ####################
    # MUTATION PROFILE #
    ####################

    reference = newdf[newdf.columns[0]]
    # This returns a binary matrix showing which sites have mutations in
    #  which sequence relative to the reference
    mutation_matrix = newdf.apply(_row_profile, axis = 1)
    filtered_mut_matrix = mutation_matrix[mutation_matrix.apply(pd.Series.nunique, axis=1) > 1]
    # Now we trim down this matrix to only include rows where one of the elements
    # == True, aka where there is one mutation relative to the reference.
    filtered_mut_matrix.reset_index(inplace=True, drop = True)
    return filtered_mut_matrix

def mutate_consensus(consensus, mutation_profile, codon_table):
    """This mutates several sequences randomly according to the mutation_profile.
    This outputs a list of CodonSeq objects, the same format as the input for
     the calculate_pi() function, the seqfreqs() function, and in the
     piN_piS() function.
    """
    names= ["seq{}".format(i) for i in range(len(mutation_profile.columns))]
    names_orig = ["ori{}".format(i) for i in range(len(mutation_profile.columns))]
    num_seqs_to_generate = len(mutation_profile.columns)

    # to fix the problem of difficulty accessing codons, I will make a codon map
    # that maps the sequence index to a codon index, and then to a position
    # index within the codon. I can then use these indices to modify the codon

    # first, make sure that the consensus has no extra bases
    if len(consensus) % 3 != 0:
        raise Exception("The consensus isn't divisible by 3. Are there indels?")

    # the consensus sequence as codons
    #now, all of the mutable data will be in seqs_matrix
    seqs_matrix = [[list(consensus[x:x+3]) for x in range(0, len("".join(consensus)), 3)]
                   for n in range(num_seqs_to_generate)]
    seqs_matrix_orig = copy.deepcopy(seqs_matrix)
    ## the codon map turns the absolute index of a nucleotide position into a list of tuples.
    ##  where the 0th element of the tuple is ith codon index to get the codon
    ##   and the 1st element of the tuple is the position index in the codon
    codon_map = [( int(seq_index / 3 ), seq_index % 3 )
                 for seq_index in range(len(consensus))]
    # the mutation map is used to keep track of sites
    mutation_map = [[' '  for i in range(len("".join(consensus)))]
                    for n in range(num_seqs_to_generate)]
    # this randomly generates a list of sites in the sequence to mutate
    random_sites = list(range(len(consensus)))
    np.random.shuffle(random_sites)
    random_iterator = iter(random_sites)

    #this block makes a list of all possible mutations given the 0th element
    ref_base_to_possible_mutations = {}
    for this_base in ['G','A','T','C']:
        possible_bases = [x for x in ['G','A','T','C'] if x != this_base]
        possible_mutations = [dict(zip(list(range(1,4)), p)) for p in permutations(possible_bases)]
        for item in possible_mutations:
            item.update({0:this_base})
        ref_base_to_possible_mutations[this_base] = possible_mutations

    #print(codon_table.forward_table)
    num_sites_mutated = 0
    i = 0
    if 'mutate' in options.debug:
        print("num_mutations: {}".format(len(mutation_profile)))

    #only stop once we have mutated enough sites
    while num_sites_mutated < len(mutation_profile):
        #get a random site to try to mutate
        random_site_index = next(random_iterator)
        if 'mutate' in options.debug:
            print("mutation site: {}".format(random_site_index))
            print("num_sites_mitated: {}".format(num_sites_mutated))

        this_codon_ix = codon_map[random_site_index][0]
        this_codon_pos_ix = codon_map[random_site_index][1]
        this_mut_profile = list(mutation_profile.iloc[i,])
        #print(this_mut_profile)
        mutation_dict = {}
        successful_mutation = []
        ref_base = seqs_matrix[0][this_codon_ix][this_codon_pos_ix]
        original_codon = seqs_matrix[0][this_codon_ix]
        mutations_dict = np.random.choice(ref_base_to_possible_mutations[ref_base])
        all_new_codons = []
        for j in range(1, len(this_mut_profile)):
            # if our mutation profile shows that this should be mutated, change
            #if the mut profile is zero, then it matches the reference. Otherwise
            # it is a mutation relative to the reference.
            # first we need to determine which codon we're in
            # it's zero-based indexing so the rules for mod are 1 off
            done = False
            mutation_map[j][random_site_index] = '*'
            new_base = mutations_dict[this_mut_profile[j]]
            new_codon = copy.deepcopy(seqs_matrix[j][this_codon_ix])
            new_codon[this_codon_pos_ix] = new_base
            all_new_codons.append(new_codon)
            #seqs_matrix[j][this_codon_ix][this_codon_pos_ix] = new_base
            # if this new codon isn't in the stop_codons, proceed
            if "".join(new_codon) not in codon_table.stop_codons:
                successful_mutation.append(True)
                if 'mutate' in options.debug:
                       print("{} to {}".format(this_base, new_base))
            else:
                #we don't do anything since we've implemented a base-wise change
                successful_mutation.append(False)
                if 'mutate' in options.debug:
                       if "".join(seqs_matrix[j][this_codon_ix]) in codon_table.stop_codons:
                           print("accidentally a stop codon: {}".format(seqs_matrix[j][this_codon_ix]))

        # if all the elements in successful_mutation are true, we can safely add
        #  the new mutations to the original sequence and move to the next
        #  random site
        if 'mutate' in options.debug:
            print("original codon: ", end="")
            print(codon_map[random_site_index])
            print("new codons: ", end="")
            print(all_new_codons)
        if sum(successful_mutation) == len(successful_mutation):
            for j in range(1, len(this_mut_profile)):
                #print("changing this codon: {} to this {}".format(
                #    seqs_matrix[j][this_codon_ix], all_new_codons[j-1]))
                #if all_new_codons[j-1] in codon_table.stop_codons:
                #    print("There's a stop codon: {}".format(all_new_codons[j-1]))
                seqs_matrix[j][this_codon_ix] = all_new_codons[j-1]
            num_sites_mutated += 1
            i += 1

    #for codon in seqs_matrix:
    #    if codon not in codon_table.stop_codons:
    #        print("somehow stop in seq: {}".format(codon))
    flat_seqs = [[item for sublist in l for item in sublist] for l in seqs_matrix]
    # this is a nice diagnostic for a single mutation
    #if 'mutate' in options.debug:
    #    flat_ori = [[item for sublist in l for item in sublist] for l in seqs_matrix_orig]
    #    for i in range(len(flat_ori)):
    #        a = "".join(mutation_map[i])
    #        print("mut: {}".format(' '.join([a[i:i+3] for i in range(0, len(a), 3)])))
    #        a = "".join(flat_ori[i])
    #        print("ori: {}".format(' '.join([a[i:i+3] for i in range(0, len(a), 3)])))
    #        a = "".join(flat_seqs[i])
    #        print("new: {}\n".format(' '.join([a[i:i+3] for i in range(0, len(a), 3)])))

    codonseqs = [CodonSeq( "".join(row)) for row in flat_seqs]
    for i in range(len(codonseqs)):
        codonseqs[i].id = names[i]
    return codonseqs



def simulate_chunk(arg_dict):
    """This is the helper function for parallelism.
    This is the actual block of code that performs the mutations"""
    all_results = []
    for i in range(arg_dict['numSimulations']):
        # Fourth, generate a random (or somewhat random) consensus sequence
        consensus = random_sequence(arg_dict['codonseqs'], mode = arg_dict['mode']) 
        # Fifth, determine the mutation profile of the other sequences we will mutate
        mutation_profile = get_mutation_profile(consensus, arg_dict['codonseqs'])
        if 'mutate' in options.debug:
            print("num_mutations: {}".format(len(mutation_profile)))
            print("len mutation profile: {}".format(len(mutation_profile)))
            print(mutation_profile)
        # Sixth, Make a new list of sequences. The 0th will be the consensus, the remaining
        #  sequences will be mutated randomly according to the mutation profile
        mutated_seqs = mutate_consensus(consensus, mutation_profile,
                                        arg_dict['codon_table'])
        results = calculate_piN_piS(mutated_seqs, arg_dict['method'],
                                    arg_dict['codon_table'])
        results['pi'] = calculate_pi(mutated_seqs)
        results['seqname'] = arg_dict['genename']
        results['type'] = "simulation"
        all_results.append(results)
    return all_results

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

def pinpissim(args):
    #First, read in the options
    global options
    options = args
    print(options)
    if options.tt_options:
        print("Options for codon usage tables in --tt_code are:")
        for key in sorted(Bio.Data.CodonTable.generic_by_id):
            print("{}: {}".format(key, Bio.Data.CodonTable.generic_by_id[key].names))
        sys.exit()

    results_file = options.results_file
    # If we've already done an analysis and the file is there already
    #  don't bother to do it again, but instead just plot
    if not os.path.exists(results_file):
        #1.5 get a list of files from the directory we provided.
        # This is a dict object with filesnames as keys
        filelist = fasta_dir_to_gene_filelist(options.fasta_dir)

        #second, select the codon alphabet to use
        codon_alphabet = get_codon_alphabet(Bio.Data.CodonTable.generic_by_id[options.tt_code], gap_char="-")
        codon_table =    Bio.Data.CodonTable.generic_by_id[options.tt_code]
        all_results = []
        for genename in filelist:
            # Third, read in the sequences to mutate
            codonseqs = fasta_path_to_codonseqs(filelist[genename], codon_table, codon_alphabet)
            print(genename)

            method = options.method
            print("processing: {}".format(genename))
            #parallelized for loop goes here
            # first we calculate the real information about these data
            results = calculate_piN_piS(codonseqs, method, codon_table)
            results['pi'] = calculate_pi(codonseqs)
            results['seqname'] = genename
            results['type'] = "observed"
            all_results.append([results])

            #now we run a bunch of simulations and chop them up into blocks of 100
            numSimulations = options.numsims
            chunk_size = 1
            #if method == 'ML':
            #    numSimulations = int(numSimulations / 10)
            num_chunks = int(numSimulations/chunk_size)
            random_args = {'genename': genename, 'numSimulations': int(numSimulations/num_chunks),
                              'mode': 'dominant', 'codonseqs': codonseqs,
                              'codon_table': codon_table, 'method': method}

            results = parallel_process([random_args for x in range(num_chunks)],
                                       simulate_chunk, n_jobs = options.threads,
                                       use_kwargs = False, front_num=3)
            flat_results = [item for sublist in results for item in sublist]
            all_results.append(flat_results)

        all_flat = [item for sublist in all_results for item in sublist]
        results_df = pd.DataFrame.from_dict(all_flat)
        results_df.to_csv(results_file, index = False)
    else:
        print("\nFound {} so skipping analysis\n".format(results_file), file = sys.stderr)
        results_df = pd.read_csv(results_file)
    #print(" - Printing results", file = sys.stderr)
    #print_results(results_df)
    print(" - Plotting results", file = sys.stderr)
    plot_results(results_df, **vars(args))
    print(" - swarmplot of results", file = sys.stderr)
    piNpiS_swarmplot(results_df, **vars(args))
    print(" - Boxplot of results", file = sys.stderr)
    piNpiS_boxplot(results_df, **vars(args))


def print_results(results, **kwargs):
    """This dataframe prints out the results of the data depending on the type.
    Observed data:
      - pi, piN, piS, piN/piS, seqname
    Simulated data:
    """
    print("The observed data:")
    observed = results.query("type == 'observed'")
    print(observed)
    observed.to_csv("observed_sumstats_{}.csv".format(timestamp()), index = False)
    simulation = results.query("type == 'simulation'")
    seqnames = simulation['seqname'].unique()
    simulation_data = []
    for seqname in seqnames:
        these_results = simulation.query("seqname == '{}'".format(seqname))
        results_dict = {}
        results_dict['seqname']   = seqname

        results_dict["piN_mean"]   = np.mean(these_results['piN'])
        results_dict["piN_median"] = np.median(these_results['piN'])
        results_dict["piN_gmean"]  = scipy.stats.gmean(these_results['piN'])
        results_dict["piN_var"]    = np.var(these_results['piN'], ddof=1)
        results_dict["piN_stderr"] = scipy.stats.sem(these_results['piN'])

        results_dict["piS_mean"]   = np.mean(these_results['piS'])
        results_dict["piS_median"] = np.median(these_results['piS'])
        results_dict["piS_gmean"]  = scipy.stats.gmean(these_results['piS'])
        results_dict["piS_var"]    = np.var(these_results['piS'], ddof=1)
        results_dict["piS_stderr"] = scipy.stats.sem(these_results['piS'])

        results_dict["piNpiS_mean"]   = np.mean(these_results['piNpiS'])
        results_dict["piNpiS_median"] = np.median(these_results['piNpiS'])
        results_dict["piNpiS_gmean"]  = scipy.stats.gmean(these_results['piNpiS'])
        results_dict["piNpiS_var"]    = np.var(these_results['piNpiS'], ddof=1)
        results_dict["piNpiS_stderr"] = scipy.stats.sem(these_results['piNpiS'])
        simulation_data.append(results_dict)

    sim_results = pd.DataFrame.from_dict(simulation_data)
    sim_results.reindex_axis(sorted(sim_results.columns), axis=1)
    sim_results.sort_values(by='seqname', ascending = False, inplace=True)
    cols = list(sim_results)
    cols.insert(0, cols.pop(cols.index('seqname')))
    sim_results = sim_results.loc[:, cols]
    sim_results.reset_index(inplace=True, drop = True)

    if kwargs["output_basename"] is None:
        filename = "simulation_sumstats"
    else:
        filename = kwargs["output_basename"] + "_simulation_sumstats"
    if kwargs["no_timestamp"]:
        filename = "{}.csv".format(filename)
    else:
        filename = "{}_{}.csv".format(filename, timestamp())

    sim_results.to_csv(filename, index = False)
    print("Simulated results:")
    print(sim_results)
    print("")

def plot_results(results, **kwargs):

    sims = results[results['type'] == 'simulation']
    obs = results[results['type'] == 'observed']
    obs.sort_values(by="piNpiS", ascending = False, inplace=True)
    obs.reset_index(inplace=True, drop = True)

    seqnames = results['seqname'].unique()

    p_values = []
    for this_seqname in seqnames:
        observed_value = obs.loc[obs['seqname'] == this_seqname,'piNpiS']
        pi = obs.loc[obs['seqname'] == this_seqname,'pi']
        observations =  sims.loc[sims['seqname'] == this_seqname, ]
        num_observations = len(observations.query("piS > 0".format(float(observed_value))))
        min_piNpiS = np.min(observations['piNpiS'])
        num_ltet = len(observations.query("piNpiS <= {} and piS > 0".format(float(observed_value))))
        p_val = -1
        if num_observations > 0:
            p_val = num_ltet/num_observations
        p_values.append({'seqname': this_seqname, 'p_val': p_val, 'pi': float(pi), 'piNpiS': float(observed_value)})
    pval_df = pd.DataFrame.from_dict(p_values)

    if kwargs["output_basename"] is None:
        filename = "pval_sumstats"
    else:
        filename = kwargs["output_basename"] + "_sumstats"
    if kwargs["no_timestamp"]:
        filename = "{}.csv".format(filename)
    else:
        filename = "{}_{}.csv".format(filename, timestamp())

    pval_df.to_csv(filename, index = False)
    print("False positive rate of gene compared to simulated data:")
    print(pval_df)
    print("")


    #randomly sample because there are probably too many points to plot
    sims = sims.sample(frac=0.4, replace=False)

    cmap = cm.get_cmap('viridis')
    rgba = {seqnames[i]: cmap(i/len(seqnames))
            for i in range(len(seqnames)) }
    markers = {"observed": 'o', "simulation": "*"}
    sizes = {"observed": 3, "simulation": 1}

    #set the figure dimensions
    figWidth = 5
    figHeight = 5
    plt.figure(figsize=(figWidth,figHeight))

    #set the panel dimensions
    panelWidth = 4
    panelHeight = 2.5

    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2) + 0.25

    panel0=plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')

    #panel0.set_xlim([min(results['pi'])*0.9, max(results['pi'])*1.1])
    panel0.set_xlim([0, max(results['pi'])*1.1])
    #panel0.set_ylim([0, max(results['piNpiS']) * 0.5])
    panel0.set_ylim([0, 2])

    #panel0.set_yscale('log')
    panel0.scatter(sims['pi'], sims['piNpiS'],
                   marker = 'o',
                   alpha = 1,
                   #alpha = 0.01,
                   s = 3, #results['type'].apply(lambda x: sizes[x]),
                   c = sims['seqname'].apply(lambda x: rgba[x]))
    panel0.scatter(obs['pi'], obs['piNpiS'],
                   marker = 'X',
                   #alpha = 0.1,
                   s = 200, lw = 0,
                   c = obs['seqname'].apply(lambda x: rgba[x]))

    for i in range(0, 3):
        panel0.text(obs.loc[i,'pi'] - 0.00075, obs.loc[i, 'piNpiS'],
                    obs.loc[i, 'seqname'],
                    fontsize = 12,
                    ha='right', va='bottom',
                    color = 'black')

    panel0.set_ylabel("piN/piS")
    panel0.set_xlabel("pi")
    panel0.set_title("Simulated mutations pi and piN/piS")

    # Print image(s)
    if kwargs["output_basename"] is None:
        file_base = "piNpiS_picheck"
    else:
        file_base = kwargs["output_basename"] + "_picheck"
    print_images(
        base_output_name=file_base,
        image_formats=kwargs["fileform"],
        no_timestamp = kwargs["no_timestamp"],
        dpi=kwargs["dpi"],
        transparent=kwargs["transparent"])


def piNpiS_swarmplot(results, **kwargs):
    #print("kwargs is: \n", kwargs)
    #set the figure dimensions
    figWidth = 5
    figHeight = 6
    fig = plt.figure(figsize=(figWidth,figHeight))
    #set the panel dimensions
    panelWidth = 3
    panelHeight = 4
    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2) + 0.25
    panel0 =plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='off', labelleft='on', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panel0.spines['top'].set_visible(False)
    panel0.spines['right'].set_visible(False)
    panel0.spines['left'].set_visible(False)


    #panel0.set_ylim([0, max(results['piNpiS']) * 0.5])
    panel0.set_ylim([0, 2])
    sims = results[results['type'] == 'simulation']
    obs = results[results['type'] == 'observed']
    obs.sort_values(by="piNpiS", ascending = False, inplace=True)
    obs.reset_index(inplace=True, drop = True)
    seqnames = sorted(results['seqname'].unique(), reverse = True)

    cmap = cm.get_cmap('viridis')
    rgba = {seqnames[i]: cmap(i/len(seqnames))
            for i in range(len(seqnames)) }

    # This block mostly just prints the observed piNpiS values for the sequences
    for i in range(len(seqnames)):
        this_seqname = seqnames[i]
        observed_value = float(obs.loc[obs['seqname'] == this_seqname,'piNpiS'])
        top_y = i + 0.1 + 0.5
        bottom_y = i + 0.9 + 0.5
        panel0.plot([float(observed_value), float(observed_value)], [top_y, bottom_y],
                    c = 'black', lw = 1)
        observations =  sims.loc[sims['seqname'] == this_seqname, ]
        num_observations = len(observations)
        min_piNpiS = np.min(observations['piNpiS'])
        num_ltet = len(observations.query("piNpiS <= {} and piS > 0".format(float(observed_value))))
        p_val = num_ltet/num_observations
        #print("{} observed: {}, min: {}, #ltet: {}, #obs: {}".format(this_seqname,
        #       float(observed_value), min_piNpiS, num_ltet, num_observations))

    #randomly sample because there are probably too many points to plot
    sims = sims.query('piNpiS > 0')
    simulations = [sims.loc[sims['seqname'] == seqname, 'piNpiS'] for seqname in seqnames]

    bp=panel0.boxplot(simulations, \
                  positions=np.arange(1, len(seqnames) + 1, 1), \
                  patch_artist=True, widths=0.5, vert = False)
    for box in bp['boxes']:
        box.set(edgecolor=(0,0,0,0.5),facecolor=(0,0,0,0),linewidth=0)
    for whisker in bp['whiskers']:
        whisker.set(color=(0,0,0,0.25), linestyle='-',linewidth=1)
    for median in bp['medians']:
        median.set(color='red', linestyle='-',linewidth=2)
    for flier in bp['fliers']:
        flier.set(markersize=0)
    for cap in bp['caps']:
        cap.set(lw=0)

    #fig.set_size_inches(4, 6, forward=True)

    # this is here because I need to align the tick labels
    # https://stackoverflow.com/questions/15882249/
    plt.draw()
    panel0.set_yticklabels(seqnames)
    yax = panel0.get_yaxis()
    # find the maximum width of the label on the major ticks
    pad = max(T.label.get_window_extent().width for T in yax.majorTicks)
    yax.set_tick_params(pad=pad/3)


    sims = sims.sample(frac=0.05, replace=False)
    simulations = [sims.loc[sims['seqname'] == seqname, 'piNpiS'] for seqname in seqnames]
    masterPP = []
    for i in np.arange(0,len(seqnames),1):
        center = 1 + i
        left_bound = i
        right_bound = 2 + i
        step_size = 0.02
        cutoff = 0.01
        placed_points = []
        counter = 0
        seqname = seqnames[i]
        for y_value in simulations[i]:
            counter+=1
            if len(placed_points)==0:
                placed_points.append((center, y_value,rgba[seqname]))
            else:
                potential_x_position=[]
                for x_position in np.arange(left_bound,right_bound, step_size):
                    distances=[]
                    for placed_point in placed_points:
                        distance=((x_position-placed_point[0])**2+(y_value-placed_point[1])**2)**0.5
                        distances.append(distance)
                    if min(distances)>cutoff:
                        potential_x_position.append(x_position)
                if len(potential_x_position)>0:
                     best_x_position=sorted(potential_x_position,key=lambda x: np.absolute(x-center))[0]
                     placed_points.append((best_x_position,y_value, rgba[seqname]))
                else:
                     print('point not placed: ',y_value)
        masterPP += placed_points

    panel0.set_xlim([0, 3])
    panel0.set_xlabel(r"Efficiency of selection ($\pi N/\pi S$)")
    panel0.set_title("Observed and simulated $\pi N/\pi S$")


    for point in masterPP:
        panel0.plot(point[1],point[0],marker='o',ms=2,
                    mfc=point[2], mew=0,linewidth=0,
                    alpha=0.25)

    # Print image(s)
    if kwargs["output_basename"] is None:
        file_base = "piNpiS_swarmplot"
    else:
        file_base = kwargs["output_basename"] + "_swarmplot"
    print_images(
        base_output_name=file_base,
        image_formats=kwargs["fileform"],
        no_timestamp = kwargs["no_timestamp"],
        dpi=kwargs["dpi"],
        transparent=kwargs["transparent"])

def piNpiS_boxplot(results, **kwargs):
    print("kwargs is: \n", kwargs)
    #set the figure dimensions
    figWidth = 5
    figHeight = 6
    fig = plt.figure(figsize=(figWidth,figHeight))
    #set the panel dimensions
    panelWidth = 3
    panelHeight = 4
    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2) + 0.25
    panel0 =plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='off', labelleft='on', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panel0.spines['top'].set_visible(False)
    panel0.spines['right'].set_visible(False)
    panel0.spines['left'].set_visible(False)


    #panel0.set_ylim([0, max(results['piNpiS']) * 0.5])
    panel0.set_ylim([0, 2])
    sims = results[results['type'] == 'simulation']
    obs = results[results['type'] == 'observed']
    obs.sort_values(by="piNpiS", ascending = False, inplace=True)
    obs.reset_index(inplace=True, drop = True)
    seqnames = sorted(results['seqname'].unique(), reverse = True)

    # This block mostly just prints the observed piNpiS values for the sequences
    for i in range(len(seqnames)):
        this_seqname = seqnames[i]
        observed_value = float(obs.loc[obs['seqname'] == this_seqname,'piNpiS'])
        top_y = i + 0.1 + 0.5
        bottom_y = i + 0.9 + 0.5
        panel0.plot([float(observed_value), float(observed_value)], [top_y, bottom_y],
                    c =(0.88, 0.03, 0.005, 0.75), lw = 2)
        observations =  sims.loc[sims['seqname'] == this_seqname, ]
        num_observations = len(observations)
        min_piNpiS = np.min(observations['piNpiS'])
        num_ltet = len(observations.query("piNpiS <= {} and piS > 0".format(float(observed_value))))
        p_val = num_ltet/num_observations
        #print("{} observed: {}, min: {}, #ltet: {}, #obs: {}".format(this_seqname,
        #       float(observed_value), min_piNpiS, num_ltet, num_observations))

    #randomly sample because there are probably too many points to plot
    #sims = sims.query('piNpiS > 0')
    simulations = [sims.loc[sims['seqname'] == seqname, 'piNpiS'] for seqname in seqnames]

    bp=panel0.boxplot(simulations, \
                  positions=np.arange(1, len(seqnames) + 1, 1), \
                  patch_artist=True, widths=0.5, vert = False,
                  notch=True, whis=[5, 95])
    for box in bp['boxes']:
        box.set(#edgecolor=(0,0,0,0.5),
                facecolor="#d8d2d7",linewidth=1)
    #for whisker in bp['whiskers']:
    #    whisker.set(color=(0,0,0,0.25), linestyle='-',linewidth=1)
    for median in bp['medians']:
        median.set(color='#443968', linestyle='dotted',linewidth=1)
    #for flier in bp['fliers']:
    #    flier.set(markersize=0)
    #for cap in bp['caps']:
    #    cap.set(lw=0)

    #fig.set_size_inches(4, 6, forward=True)

    # this is here because I need to align the tick labels
    # https://stackoverflow.com/questions/15882249/
    plt.draw()
    panel0.set_yticklabels(seqnames)
    yax = panel0.get_yaxis()
    # find the maximum width of the label on the major ticks
    pad = max(T.label.get_window_extent().width for T in yax.majorTicks)
    yax.set_tick_params(pad=pad/3)

    panel0.set_xlabel(r"Efficiency of selection ($\pi N/\pi S$)")
    panel0.set_title("Observed and simulated $\pi N/\pi S$")

    # Print image(s)
    if kwargs["output_basename"] is None:
        file_base = "piNpiS_boxplot"
    else:
        file_base = kwargs["output_basename"] + "_boxplot_noxlim"
    print_images(
        base_output_name=file_base,
        image_formats=kwargs["fileform"],
        no_timestamp = kwargs["no_timestamp"],
        dpi=kwargs["dpi"],
        transparent=kwargs["transparent"])

    panel0.set_xlim([0, 4])
    # Print image(s)
    if kwargs["output_basename"] is None:
        file_base = "piNpiS_boxplot"
    else:
        file_base = kwargs["output_basename"] + "_boxplot"
    print_images(
        base_output_name=file_base,
        image_formats=kwargs["fileform"],
        no_timestamp = kwargs["no_timestamp"],
        dpi=kwargs["dpi"],
        transparent=kwargs["transparent"])

def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")

def run(args):
    pinpissim(args)
