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

# This program needs to do a few things
# 0. Select a codon usage table
# 1. Read in a fasta alignment
# 2. Determine all of the polymorphic sites in the alignment
# 3. Calculate pi, piN, and piS for the alignment

#This class is used in argparse to expand the ~. This avoids errors caused on
# some systems.
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))


def parse_arguments():
    """This method parses the arguments provided during command-line input
    to control how the kmers are counted. Please see the documentation for help.
    """

    #using Kevin's idea here to make te docstring the help file
    parser=argparse.ArgumentParser(description=__doc__,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outprefix",
                        type=str,
                        default=sys.stdout,
                        help="""the output will be saved here""")
    parser.add_argument("--tt_code",
                        type = int,
                        default = 4,
                        help="""Select which gene code to use. Default is
                        Standard""")
    parser.add_argument("--tt_options",
                        action = 'store_true',
                        default = False,
                        help="""Display the optional gene code names that you
                        may pass to the <--tt_code> argument in a subsequent
                        run""")
    parser.add_argument("--debug",
                        type = str,
                        nargs = '+',
                        default = ['NA'],
                        help = """Use this to print out special things while
                        debugging""")
    parser.add_argument("--fasta_dir",
                        action = FullPaths,
                        required = True,
                        help = """This is the directory where the fasta file
                        alignments are located. The filename will be used as the
                        figure label while plotting.""")
    parser.add_argument("--results_file",
                        required = True,
                        help = """The results will be written to this file as a
                        csv.""")
    parser.add_argument("--method",
                        default = 'NG86',
                        choices = ['NG86', 'LWL85', 'YN00', 'ML'],
                        help = """the method to use in the simulation""")
    parser.add_argument("--numsims",
                        required = True,
                        type = int,
                        help = """the results will be written to this file as a
                        csv.""")
    parser.add_argument("--threads",
                        required = False,
                        type = int,
                        help = """The num threads to attempt to use.""")

    args = parser.parse_args()
    return args

def seqfreqs(seqs):
    """Calculates the relative frequency of each sequence for calculating pi
    """
    if "seqfreqs" in options.debug:
        print("There are {} seqs".format(len(seqs)))
    x = []
    #this block calculates the frequencies of each sequence
    for i in range(len(seqs)):
        this_x = 0
        for j in range(len(seqs)):
            if str(seqs[i]) == str(seqs[j]):
                if "seqfreqs" in options.debug:
                    print("{} == {}".format(i, j))
                this_x += 1
        x.append(this_x/len(seqs))
    #print("done with these seqfreqs\n")
    if "seqfreqs" in options.debug:
        print("the frequencies are {}".format(x))
    return x

def seqs_to_df(seqs):
    """converts sequences to a pandas dataframe s.t. each column is one sequence
    and each row is one position in the sequence
    """
    df_list = []
    names_list = []
    for index in range(len(seqs)):
        names_list.append(seqs[index].id)
    newseqs = [[char for char in str(seq)] for seq in seqs]
    df = pd.DataFrame.from_items(zip(names_list, newseqs))
    return df

def calculate_pi(seqs):
    """
    there is some info here http://www.columbia.edu/cu/biology/courses/c3020/solutions-2.html
    calculates pi of a list of CodonSeq objects
    """
    df = seqs_to_df(seqs)
    #df.apply(rowwise_unique, axis=1)

    # this is called x due to tradition in the definition of pi
    x = seqfreqs(seqs)
    running_sum_pi = 0
    for i in range(len(seqs)):
        for j in range(i+1, len(seqs)):
            comp = df.copy().iloc[:,[i,j]]
            comp.replace('-', np.nan, inplace=True)
            comp.dropna(axis=0,how='any', inplace=True)
            #print(comp)
            num_difs = sum(comp.iloc[:,0] != comp.iloc[:,1])
            len_seqs = len(comp.iloc[:,0])
            this_sum = x[i] * x[j] * (num_difs/len_seqs)
            running_sum_pi += this_sum
            if 'pi' in options.debug:
                print("x[i] * x[j] * pi = {} * {} * {}/{}".format(x[i], x[j], num_difs, len_seqs))
    if 'pi' in options.debug:
        print("pi: {}".format(running_sum_pi))
    return running_sum_pi

def _weighted_base(row):
    """Returns a base using weighted probabilities of what is observed in this
    sequence column
    - input is a row of bases corresponding to one alignment column
    """
    unique, counts = np.unique(row, return_counts=True)
    weights = counts / len(row)
    maxweight = max(weights)
    choice = np.random.choice(unique, p=weights)
    return choice

def _max_or_rand_base(row):
    """Returns the base with maximum frequency, or if none exists,
    randomly selects one of the objects with maximum frequency

    - input is a row of bases corresponding to one alignment column
    """
    unique, counts = np.unique(row, return_counts=True)
    weights = counts / len(row)
    maxweight = max(weights)
    choice = np.random.choice(unique[np.where(weights == maxweight)])
    return choice

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
    # this is important and randomly chooses sites in the sequence to mutate randomly
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
            if 'mutate' in options.debug:
                   print("{} to {}".format(this_base, new_base))
            # if this new codon isn't in the stop_codons, proceed
            if 'mutate' in options.debug:
                   if "".join(seqs_matrix[j][this_codon_ix]) in codon_table.stop_codons:
                       print("accidentally a stop codon: {}".format(seqs_matrix[j][this_codon_ix]))
            if "".join(new_codon) not in codon_table.stop_codons:
                successful_mutation.append(True)
            else:
                #we don't do anything since we've implemented a base-wise change
                successful_mutation.append(False)
        # if all the elements in successful_mutation are true, we can safely add
        #  the new mutations to the original sequence and move to the next
        #  random site
        if 'mutate' in options.debug:
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

def calculate_piN_piS(codonseqs, method, codon_table):
    analysis = {"seqname": "", "piN": -1, "piS": -1, "piNpiS": -1, "pi": -1, "method":method}
    x = seqfreqs(codonseqs)
    if 'piNpiS' in options.debug:
        print("freqs are: {}".format(x))
        print("len codonseqs is: ", len(codonseqs))
    piN = 0
    piS = 0
    for i in range(len(codonseqs)):
        for j in range(i+1, len(codonseqs)):
            dN, dS = cal_dn_ds(codonseqs[i], codonseqs[j], codon_table=codon_table, method=method)
            piN = piN + (x[i] * x[j] * dN)
            piS = piS + (x[i] * x[j] * dS)
            if 'piNpiS' in options.debug:
                print("{0} dN{1}{2}={3} dS{1}{2}={4}".format(method, i, j, dN, dS))

    analysis['piN'] = piN
    analysis['piS'] = piS
    try:
        analysis['piNpiS'] = piN/piS
    except:
        analysis['piNpiS'] = 0
    if 'piNpiS' in options.debug:
        print ("{0} dN={1:.3f} dS={2:.3f} piN/piS = {3:.3f}".format(
            method, analysis['piN'], analysis['piS'], analysis['piNpiS']))


    return analysis

def remove_stop_from_end(codonseq, codon_table):
    """This checks to see if the last three bases of the sequence is a stop
       codon. If so, remove the last codon."""
    if str(codonseq[-3:]) in codon_table.stop_codons:
        #print("chopped", len(CodonSeq().from_seq(codonseq[:-3])))
        return CodonSeq().from_seq(codonseq[:-3])
    elif str(codonseq[-3:]) == '---':
        #print("chopped", len(CodonSeq().from_seq(codonseq[:-3])))
        return CodonSeq().from_seq(codonseq[:-3])
    else:
        #print("unchopped", len(codonseq))
        return codonseq

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

def main():
    #First, read in the options
    global options
    options = parse_arguments()
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
        # This is a dict object with key as
        filelist = {os.path.splitext(x)[0]:os.path.join(os.path.abspath(options.fasta_dir), x)
                    for x in os.listdir(options.fasta_dir)
                    if os.path.splitext(x)[1]}

        #second, select the codon alphabet to use
        codon_alphabet = get_codon_alphabet(Bio.Data.CodonTable.generic_by_id[options.tt_code], gap_char="-")
        codon_table =    Bio.Data.CodonTable.generic_by_id[options.tt_code]
        all_results = []
        for genename in filelist:
            # Third, read in the sequences to mutate
            codonseqs = []
            for record in SeqIO.parse(filelist[genename], "fasta"):
                this_CS = remove_stop_from_end(CodonSeq().from_seq(record.seq), codon_table)
                #this_CS = CodonSeq().from_seq(record.seq)
                this_CS.alphabet = codon_alphabet
                this_CS.id = record.id
                codonseqs.append(this_CS)
            print(genename)

            #methods = ['NG86', 'LWL85', 'YN00', 'ML']
            #methods = ['NG86', 'LWL85', 'ML']
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
        results_df.to_csv(results_file)
    else:
        print("found {} so skipping analysis".format(results_file))
        results_df = pd.read_csv(results_file)

    plot_results(results_df)

def plot_results(results):
    sims = results[results['type'] == 'simulation']
    obs = results[results['type'] == 'observed']
    obs.sort_values(by="piNpiS", ascending = False, inplace=True)
    obs.reset_index(inplace=True, drop = True)
    print(obs.head())

    seqnames = results['seqname'].unique()

    for this_seqname in seqnames:
        observed_value = obs.loc[obs['seqname'] == this_seqname,'piNpiS']
        observations =  sims.loc[sims['seqname'] == this_seqname, ]
        num_observations = len(observations)
        min_piNpiS = np.min(observations['piNpiS'])
        num_ltet = len(observations.query("piNpiS <= {} and piNpiS > 0".format(float(observed_value))))
        p_val = num_ltet/num_observations
        #print("{} observed: {}, min: {}, #ltet: {}, #obs: {}".format(this_seqname,
        #       float(observed_value), min_piNpiS, num_ltet, num_observations))
        print("pvalue of {}: {}".format(this_seqname, p_val))

    #randomly sample because there are probably too many points to plot
    sims = sims.sample(frac=0.4, replace=False)

    cmap = cm.get_cmap('viridis')
    rgba = {seqnames[i]: cmap(i/len(seqnames))
            for i in range(len(seqnames)) }
    markers = {"observed": 'o', "simulation": "*"}
    sizes = {"observed": 3, "simulation": 1}

    plt.style.use('BME163')

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
    #counts = generate_heat_map(panel0, results, purple1)

    plt.savefig("simulation_results_{}.png".format(timestamp()), dpi=600, transparent=False)

def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")

if __name__ == "__main__":
    sys.exit(main())
