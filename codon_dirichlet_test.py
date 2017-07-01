# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:54:48 2017

@author: Jordan
@contributor: DTS
"""

import collections
import math
import itertools
import sys
from Bio import SeqIO
from builtins import range
import progressbar
import numpy as np
from sklearn.metrics import confusion_matrix as cmcalc


# TODO
# 1. Make something that stores true positives vs false positives
#   - Use y_true = [2, 0, 2, 2, 0, 1]
#         y_pred = [0, 0, 2, 2, 0, 2] from scikit learn
#   - more info here: http://scikit-learn.org/stable/modules/generated/sklearn.metrics.confusion_matrix.html

global nts
global starts
global stops
nts = ['A', 'C', 'G', 'T']
starts = ['TTA', 'TTG','CTG',
               'ATT', 'ATC','ATA',
               'ATG','GTG']
stops = ['TAA','TAG']

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
    # it is necessary to remove the starts and stops to avoid the analysis just
    #  being based on their binary presence/absence
    #print(counts_test)
    #print(counts_background)
    precodons = ["".join(codon) for codon in itertools.product(nts, repeat = 3)]
    codons = [x for x in precodons if x not in startsstops]
    log_like = 0
    log_like += math.lgamma(sum(counts_test[codon] for codon in codons))
    log_like += math.lgamma(sum(pseudocounts + counts_background[codon] for codon in codons))
    log_like -= math.lgamma(sum(pseudocounts + counts_background[codon] + counts_test[codon] for codon in codons))
    for codon in codons:
        log_like -= math.lgamma(pseudocounts + counts_test[codon])
        log_like -= math.lgamma(pseudocounts + counts_background[codon])
        log_like += math.lgamma(pseudocounts + counts_background[codon] + counts_test[codon])
    return log_like

def remove_starts_stops(seq):
    seq1 = [test_seq[i:i+3] for i in
# arguments are one test sequence and two lists of sequences representing the two hypothesis groups
# pseudocounts are added to each codon
# verbose flag toggles whether test will be interpreted for you on stdout
def codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.5, verbose = True): 
    # check for data invariants
    assert len(test_seq) % 3 == 0
    for nt in test_seq:
        assert nt in nts
    for seq in seq_set_1:
        assert len(seq) % 3 == 0
        for nt in seq:
            assert nt in nts
    for seq in seq_set_2:
        assert len(seq) % 3 == 0
        for nt in seq:
            assert nt in nts

    # count up codons
    codon_counts_1 = collections.Counter()
    codon_counts_2 = collections.Counter()
    codon_counts_test = collections.Counter()

    for i in range(int(len(test_seq) / 3)):
        codon_counts_test[test_seq[i:i+3]] += 1

    for seq in seq_set_1:
        for i in range(int(len(seq) / 3)):
            codon_counts_1[seq[i:i+3]] += 1

    for seq in seq_set_2:
        for i in range(int(len(seq) / 3)):
            codon_counts_2[seq[i:i+3]] += 1

    # compute bayes factor
    log_likelihood_ratio = dirichlet_log_marginal_likelihood(codon_counts_test, codon_counts_1, pseudocounts) \
                           - dirichlet_log_marginal_likelihood(codon_counts_test, codon_counts_2, pseudocounts)

    # interpret as indicated
    if verbose:
        interpret_hypothesis_test(log_likelihood_ratio)

    return log_likelihood_ratio

def gen_noncodings(nc_fasta_filepath, seqs_per_set):
    """This function takes a filepath of the fasta file that contains the
    noncoding sequences to use in the analysis.  Instead of manually
    generating all three artificial reading frames, the function
    generates them automatically and returns in list format for
    codon_dirichlet_log_likelihood_ratio()
    """

    seqs = []
    for record in SeqIO.parse(nc_fasta_filepath, "fasta"):
        one_set = []
        # I'm not sure why I originally put this conditional here,
        #  but I will leave it
        if len(record.id) > 4:
            for i in range(3):
                seqs.append(record.seq[i:])
    for i in range(len(seqs)):
        mod = len(seqs[i]) % 3
        thisLen = len(seqs[i])
        seqs[i] = str(seqs[i][0:thisLen-mod])
    #This groups all the sequences together if the are in the same locus
    grouped_seqs = [seqs[i:i+(3*seqs_per_set)] for i in range(0,len(seqs),3*seqs_per_set)]
    return grouped_seqs

def gen_codings(c_fasta_filepath, seqs_per_set):
    seqs = []
    for record in SeqIO.parse(c_fasta_filepath, "fasta"):
        seqs.append(str(record.seq))
    #This groups all the sequences together if the are in the same locus
    grouped_seqs = [seqs[i:i+seqs_per_set] for i in range(0,len(seqs),seqs_per_set)]
    return grouped_seqs

def confusion_matrix(coding_sets, noncoding_sets, seqs_per_set):
    #first, test all the trues. Generate a list of all combinations of noncoding
    # and for each combination iterate through all of the groups of coding, and for each
    # set of coding remove one and test the other against everything else.
    # That is three for loops.

    real_val = []
    observed = []
    # This makes all the possible products of noncoding_sets
    possible_indices = [list(range(len(x))) for x in noncoding_sets]
    nc_comb_list = list(itertools.product(*possible_indices))
    bar = progressbar.ProgressBar()
    seq_sets = set()
    for bari in bar(range(len(nc_comb_list))):
    #for bari in bar(range(2)):
        each = nc_comb_list[bari]
        for i in range(len(coding_sets)):
            for j in range(seqs_per_set):
                # Each of these tests is a true positive. append 1 to real_val
                # using 1 for real seqs and 2 for false
                real_val.append(1)
                # the test seq is the ith gene and the jth individual
                test_seq = coding_sets[i][j]
                # the real set is all coding genes except for at the ith position
                seq_set_1 = [seq for x in coding_sets if x != coding_sets[i] for seq in x]
                # the false set is the noncoding genes in index each
                seq_set_2 = [noncoding_sets[k][each[k]]for k in range(len(each))]
                seq_sets.update(seq_set_1)
                #print (seq_set_2)
                log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.1, verbose = False)
                if log_likelihood_ratio < 0:
                    favored = 2
                else:
                    favored = 1
                observed.append(favored)

    # now test all the true negatives.
    # compile a list of tests
    # ([set index of test: subset index], {dict of nc_set index: index of read})
    test_list = []
    #print("NC combination")
    possible_indices = [list(range(len(noncoding_sets[0])))] * (len(noncoding_sets) - 1)
    nc_comb_list = list(itertools.product(*possible_indices))
    # set up the list of analyses in the tuple format
    # first we iterate through every set of sequences in the noncoding seq set
    for i in range(len(noncoding_sets)):
        # for every set of sequences, iterate through each seq
        for j in range(len(noncoding_sets[0])):
            # now make a dict of every other noncoding combo for this test
            for combo in nc_comb_list:
                nc_index_count = 0
                combo_index_count = 0
                test_dict = {}
                #print("len of noncoding sets: {}".format( len(noncoding_sets[0]) ))
                while nc_index_count < len(noncoding_sets):
                    if nc_index_count == i:
                        #if we're at the test position(i), skip over it in the dict
                        nc_index_count += 1
                    else:
                        #print("combo_index_count: {}".format(combo_index_count))
                        #print("nc_index_count: {}".format(nc_index_count))
                        #print("len combo: {}".format(len(combo)))
                        test_dict[nc_index_count] = combo[combo_index_count]
                        #increase both the index for the dict and the combo
                        nc_index_count += 1
                        combo_index_count += 1
                test_list.append(([i,j], test_dict))

    bar2 = progressbar.ProgressBar()
    # use all the coding sequences
    seq_set_1 = [seq for l in coding_sets for seq in l]
    for neg_index in bar2(range(len(test_list))):
        analysis = test_list[neg_index]
        i,j = analysis[0]
        test_dict = analysis[1]
        # These are all negatives, append 2
        real_val.append(2)
        # the test seq is the ith gene and the jth individual
        test_seq = noncoding_sets[i][j]
        # the false set is the noncoding genes in index each
        seq_set_2 = [noncoding_sets[key][test_dict[key]] for key in test_dict]
        log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.5, verbose = False)
        if log_likelihood_ratio < 0:
            favored = 2
        else:
            favored = 1
        observed.append(favored)

    real_val
    observed
    matrix = cmcalc(real_val, observed)
    print(matrix)

def main():
    #first generate a string of the sequence to test
    test_seq = str(next(SeqIO.parse("Bf201606_COIII.fasta", "fasta")).seq)
    #Then get a list of the known sequences
    seqs_per_group = 2
    seq_set_1 = gen_codings("c_seqs_bothind.fasta", seqs_per_group)
    # then get the list of lists of noncoding sequences.
    # We first have to pass the argument of how many copies of each gene there
    #  are so that the program adds the sequences in groups, this determines how many
    #  sequences to add in each sublist.
    seq_set_2 = gen_noncodings("nc_seqs_gt50_bothind.fasta", seqs_per_group)

    confusion_matrix(seq_set_1, seq_set_2, seqs_per_group)
    #codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.5, verbose = True)

if __name__ == "__main__":
    sys.exit(main())
