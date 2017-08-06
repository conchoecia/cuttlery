#!/usr/bin/env python

#mathy stuff
import numpy as np
from numpy import cumsum
from numpy.random import rand
import scipy.stats as stats
import pandas as pd

# functools because apply doesn't work with two arguments
from functools import partial

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
from Bio.codonalign import CodonSeq
import Bio.Data.CodonTable
from Bio.codonalign.codonalphabet import get_codon_alphabet
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import _get_pi

# Also import some stufffff to parallelize the function
from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool



import progressbar
#pd.options.mode.chained_assignment = None

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
    parser.add_argument("--translation_table",
                        type=str,
                        default="ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt",
                        help="""The genetic code table to read from. The default
                        is from ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt""")
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


    args = parser.parse_args()
    return args

def seqfreqs(seqs):
    """Calculates the relative frequency of each sequence for calculating pi
    """
    x = []
    #this block calculates the frequencies of each sequence
    for i in range(len(seqs)):
        this_x = 0
        for j in range(len(seqs)):
            if np.array_equal(str(seqs[i]), str(seqs[j])):
                this_x += 1
        x.append(this_x/len(seqs))
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

def random_consensus(seqs, mode = 'dominant'):
    """Generates a random consensus sequence given the probabilities of
        observing that base
    - The dominant mode will pick the most frequent allele and randomly
      select a minor allele if all minor alleles have the same
      probability (only for n<5)
    - The random mode selects a base given the probability of observing that
      base at that locus
    """
    df = seqs_to_df(seqs)
    if mode == 'dominant':
        consensus = df.apply(_max_or_rand_base, axis = 1)
    else:
        consensus = df.apply(_weighted_base, axis = 1)
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
    edit_distance_with_consensus = df.apply(partial(_edit_distance, consensus), axis = 0)
    # sometimes there will be more than one seq with the minimum edit distance.
    #  For that reason we use a random function here to chose the final
    #   'reference' against which to determine the mutation profile
    min_edit_distance = min(edit_distance_with_consensus)
    # this is the somewhat-randomly selected reference sequence that will be
    #  used to build the mutation profile
    random_reference_index= np.random.choice(np.flatnonzero(edit_distance_with_consensus == min_edit_distance))
    if 'mutprof' in options.debug:
        print(random_reference_index)
    # Now we reshuffle the indices
    indices = list(range(0, len(df.columns.tolist())))
    if 'mutprof' in options.debug:
        print("unshuffl: ", indices)
    del indices[random_reference_index]
    indices = [random_reference_index] + indices
    cols = [df.columns.tolist()[i] for i in indices]
    if 'mutprof' in options.debug:
        print("shuffled: ", indices)
        print(df.columns.tolist())
    # now we shuffle the df to make the reference the 0th element
    newdf = df[cols]
    # print(newdf.columns.tolist()) #verified that shuffling works with print statement

    ####################
    # MUTATION PROFILE #
    ####################

    reference = newdf[newdf.columns[0]]
    # This returns a binary matrix showing which sites have mutations in
    #  which sequence relative to the reference
    mutation_matrix = newdf.apply(partial(_binary_comparison, reference))
    filtered_mut_matrix = mutation_matrix[mutation_matrix.sum(axis=1) > 0]
    # Now we trim down this matrix to only include rows where one of the elements
    # == True, aka where there is one mutation relative to the reference.
    filtered_mut_matrix.reset_index(inplace=True)
    filtered_mut_matrix.drop('index', 1, inplace=True)
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
    seqs_matrix_orig = [[list(consensus[x:x+3]) for x in range(0, len("".join(consensus)), 3)]
                for n in range(num_seqs_to_generate)]
    # this is important and randomly chooses sites in the sequence to mutate randomly
    random_sites = np.random.choice(list(range(len(consensus))), size = len(mutation_profile), replace = False)
    ## the codon map turns the absolute index of a nucleotide position into a list of tuples.
    ##  where the 0th element of the tuple is ith codon index to get the codon
    ##   and the 1st element of the tuple is the position index in the codon
    codon_map = [( int(seq_index / 3 ), seq_index % 3 )
                 for seq_index in range(len(consensus))]
    mutation_map = [[' '  for i in range(len("".join(consensus)))] for n in range(num_seqs_to_generate)]
    for i in range(len(random_sites)):
        random_site_index = random_sites[i]
        this_codon_ix = codon_map[random_site_index][0]
        this_codon_pos_ix = codon_map[random_site_index][1]
        this_mut_profile = list(mutation_profile.iloc[i,])
        new_row = []
        for j in range(len(this_mut_profile)):
            # if our mutation profile shows that this should be mutated, change
            this_base = seqs_matrix[j][this_codon_ix][this_codon_pos_ix]
            if this_mut_profile[j]:
                # first we need to determine which codon we're in
                # it's zero-based indexing so the rules for mod are 1 off
                done = False
                mutation_map[j][random_site_index] = '*'
                while not done:
                    new_base = np.random.choice([x for x in ['G','A','T','C'] if x != this_base])
                    seqs_matrix[j][this_codon_ix][this_codon_pos_ix] = new_base
                    if 'mutate' in options.debug:
                        print("{} to {}".format(this_base, new_base))
                    # if this new codon isn't in the stop_codons, proceed
                    if 'mutate' in options.debug:
                        if "".join(seqs_matrix[j][this_codon_ix]) in codon_table.stop_codons:
                            print("accidentally a stop codon: {}".format(seqs_matrix[j][this_codon_ix]))
                    if "".join(seqs_matrix[j][this_codon_ix]) not in codon_table.stop_codons:
                        done = True
            else:
                #we don't do anything since we've implemented a base-wise change
                new_base = this_base

    flat_seqs = [[item for sublist in l for item in sublist] for l in seqs_matrix]
    # this is a nice diagnostic for a single mutation
    if 'mutate' in options.debug:
        flat_ori = [[item for sublist in l for item in sublist] for l in seqs_matrix_orig]
        for i in range(len(flat_ori)):
            a = "".join(mutation_map[i])
            print("mut: {}".format(' '.join([a[i:i+3] for i in range(0, len(a), 3)])))
            a = "".join(flat_ori[i])
            print("ori: {}".format(' '.join([a[i:i+3] for i in range(0, len(a), 3)])))
            a = "".join(flat_seqs[i])
            print("new: {}\n".format(' '.join([a[i:i+3] for i in range(0, len(a), 3)])))

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
        print("chopped", len(CodonSeq().from_seq(codonseq[:-3])))
        return CodonSeq().from_seq(codonseq[:-3])
    elif str(codonseq[-3:]) == '---':
        print("chopped", len(CodonSeq().from_seq(codonseq[:-3])))
        return CodonSeq().from_seq(codonseq[:-3])
    else:
        print("unchopped", len(codonseq))
        return codonseq

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
            print(codonseqs)

            #methods = ['NG86', 'LWL85', 'YN00', 'ML']
            methods = ['NG86', 'LWL85', 'ML']
            for method in methods:
                # first we calculate the real information about these data
                results = calculate_piN_piS(codonseqs, method, codon_table)
                results['pi'] = calculate_pi(codonseqs)
                results['seqname'] = genename
                results['type'] = "observed"
                all_results.append(results)

                #now we run a bunch of simulations and chop them up into blocks of 100
                numSimulations = 500
                chunk_size = 100
                if method == 'ML':
                    numSimulations = int(numSimulations / 10)
                bar = progressbar.ProgressBar()
                random_args
                #Third, mutate this sequence 35 times
                print("processing: {}".format(genename))
                for i in bar(range(numSimulations)):
                    # Fourth, generate a random (or somewhat random) consensus sequence
                    consensus = random_consensus(codonseqs, mode = 'dominant')

                    # Fifth, determine the mutation profile of the other sequences we will mutate
                    mutation_profile = get_mutation_profile(consensus, codonseqs)

                    # Sixth, Make a new list of sequences. The 0th will be the consensus, the remaining
                    #  sequences will be mutated randomly according to the mutation profile
                    mutated_seqs = mutate_consensus(consensus, mutation_profile, codon_table)

                    results = calculate_piN_piS(mutated_seqs, method, codon_table)
                    results['pi'] = calculate_pi(mutated_seqs)
                    results['seqname'] = genename
                    results['type'] = "simulation"
                    all_results.append(results)
                results_df = pd.DataFrame.from_dict(all_results)
                results_df.to_csv(results_file)
    else:
        print("found {} so skipping analysis".format(results_file))
        results_df = pd.read_csv(results_file)



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
