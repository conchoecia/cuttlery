# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:54:48 2017

@author: Jordan/DTS
#info about confusion matrix here: http://www.dataschool.io/simple-guide-to-confusion-matrix-terminology/
"""

import collections
import math
import itertools
import sys
from Bio import SeqIO
from builtins import range
import progressbar
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
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

#def remove_starts_stops(seq):
#    seq1 = [test_seq[i:i+3] for i in

def gen_codons(seq, codon_filtering):
    codons = [seq[i:i+3] for i in np.arange(0, len(seq), 3)]
    new_codons = codons
    if codon_filtering == 1:
        if codons[0] in starts:
            new_codons = new_codons[1:]
        if codons[-1] in stops:
            new_codons = new_codons[0:-1]
    elif codon_filtering == 2:
        if codons[0] in starts:
            new_codons = new_codons[1:]
        new_codons = [x for x in new_codons if x not in stops]
    elif codon_filtering == 3:
        new_codons = [x for x in codons if x not in starts + stops]
    return new_codons


# arguments are one test sequence and two lists of sequences representing the two hypothesis groups
# pseudocounts are added to each codon
# verbose flag toggles whether test will be interpreted for you on stdout
def codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.5, verbose = True, codon_filtering = True): 
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

    for codon in gen_codons(test_seq, codon_filtering):
        codon_counts_test[codon] += 1

    for seq in seq_set_1:
        for codon in gen_codons(seq, codon_filtering):
            codon_counts_1[codon] += 1

    for seq in seq_set_2:
        for codon in gen_codons(seq, codon_filtering):
            codon_counts_2[codon] += 1

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


def gen_test_seqlist(c_fasta_filepath):
    #this assumes that each test sequence has its own file
    seqs = []
    for record in SeqIO.parse(c_fasta_filepath, "fasta"):
        seqs.append(str(record.seq))
    return seqs

def gen_codings(c_fasta_filepath, seqs_per_set):
    seqs = []
    for record in SeqIO.parse(c_fasta_filepath, "fasta"):
        seqs.append(str(record.seq))
    #This groups all the sequences together if the are in the same locus
    grouped_seqs = [seqs[i:i+seqs_per_set] for i in range(0,len(seqs),seqs_per_set)]
    return grouped_seqs

def confusion_matrix(coding_sets, nc_seqsets, seqs_per_set, codon_filtering):
    #first, test all the trues. Generate a list of all combinations of noncoding
    # and for each combination iterate through all of the groups of coding, and for each
    # set of coding remove one and test the other against everything else.
    # That is three for loops.

    real_val = []
    observed = []
    llratio = []
    # This makes all the possible products of nc_seqsets
    possible_indices = [list(range(len(x))) for x in nc_seqsets]
    nc_comb_list = list(itertools.product(*possible_indices))
    bar = progressbar.ProgressBar()
    seq_sets = set()
    for bari in bar(range(len(nc_comb_list))):
    #for bari in bar(range(100)):
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
                seq_set_2 = [nc_seqsets[k][each[k]]for k in range(len(each))]
                seq_sets.update(seq_set_1)
                #print (seq_set_2)
                log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.1, verbose = False, codon_filtering = codon_filtering)
                llratio.append(log_likelihood_ratio)
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
    possible_indices = [list(range(len(nc_seqsets[0])))] * (len(nc_seqsets) - 1)
    nc_comb_list = list(itertools.product(*possible_indices))
    # set up the list of analyses in the tuple format
    # first we iterate through every set of sequences in the noncoding seq set
    for i in range(len(nc_seqsets)):
        # for every set of sequences, iterate through each seq
        for j in range(len(nc_seqsets[0])):
            # now make a dict of every other noncoding combo for this test
            for combo in nc_comb_list:
                nc_index_count = 0
                combo_index_count = 0
                test_dict = {}
                #print("len of noncoding sets: {}".format( len(nc_seqsets[0]) ))
                while nc_index_count < len(nc_seqsets):
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
    #for neg_index in bar2(range(2000)):
        analysis = test_list[neg_index]
        i,j = analysis[0]
        test_dict = analysis[1]
        # These are all negatives, append 2
        real_val.append(2)
        # the test seq is the ith gene and the jth individual
        test_seq = nc_seqsets[i][j]
        # the false set is the noncoding genes in index each
        seq_set_2 = [nc_seqsets[key][test_dict[key]] for key in test_dict]
        log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.5, verbose = False, codon_filtering = codon_filtering)
        llratio.append(log_likelihood_ratio)
        if log_likelihood_ratio < 0:
            favored = 2
        else:
            favored = 1
        observed.append(favored)

    matrix = cmcalc(real_val, observed)
    np.save("llratio", np.array(llratio))
    np.save("real_val", np.array(real_val))
    np.save("observed", np.array(observed))
    print(matrix)
    
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
    

def plot_ratios(test_seq_names):
    llratio_list = np.load("llratio.npy")
    true_value_list = np.load("real_val.npy")
    obs_value_list  = np.load("observed.npy")
    test_seq_sets = []
    for name in test_seq_names:
        test_seq_sets.append(np.load("{}_llratio.npy".format(name)))

    llratios_1 = np.array([llratio_list[i] for i in range(len(true_value_list)) if true_value_list[i] == 1])
    llratios_2 = np.array([llratio_list[i] for i in range(len(true_value_list)) if true_value_list[i] == 2])

    llratios_1weights = 100 * np.ones_like(llratios_1) / llratios_1.size
    llratios_2weights = 100 * np.ones_like(llratios_2) / llratios_2.size
    unk1_weights = 100 * np.ones_like(test_seq_sets[0]) / test_seq_sets[0].size
    nad2l_weights = 100 * np.ones_like(test_seq_sets[1]) / test_seq_sets[1].size

    fig, ax = plt.subplots()
    #ax.hist(llratios_1, bins='auto', weights=llratios_1weights, color='lightblue', alpha=0.5, normed=True)
    #ax.hist(llratios_2, bins='auto', weights=llratios_2weights, color='salmon', alpha=0.5, normed=True)
    ax.hist(llratios_1, bins='auto', label = 'coding', color='lightblue', alpha=0.5, normed=True)
    ax.hist(llratios_2, bins='auto', label = 'nc', color='salmon', alpha=0.5, normed=True)
    ax.hist(test_seq_sets[0], bins='auto', label = 'unk1', color='purple', alpha=0.5, normed=True)
    ax.hist(test_seq_sets[1], bins='auto', label = 'nd2l', color='gold', alpha=0.5, normed=True)

    #create legend
    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in ["lightblue", "salmon", "purple", "gold"]]
    labels= ["coding", "nc", "unk1", "nd2l"]
    ax.legend(handles, labels)
    ax.set(title='Histogram Comparison', ylabel='% of Dataset in Bin')
    ax.margins(0.05)
    ax.set_ylim(bottom=0)
    plt.show()

def test_unknown(test_seqlist, coding_seqsets, nc_seqsets, codon_filtering, test_gene_name):
    #for this test, we need to generate all possible combinations of nc set,
    # plus using the same coding set for each analysis. This is very similar to
    # testing the coding sequences above.
    llratio = []
    #use all the coding sequences
    seq_set_1 = [seq for l in coding_seqsets for seq in l]

    # This makes all the possible products of nc_seqsets
    possible_indices = [list(range(len(x))) for x in nc_seqsets]
    nc_comb_list = list(itertools.product(*possible_indices))
    bar = progressbar.ProgressBar()
    seq_sets = set()

    # like for the negative analysis above, make a list of analyses in a list of lists
    # each element of list:
    # [index of test seq in test_seqlist, [indices of nc_seqs]]

    analyses = []
    for i in range(len(test_seqlist)):
        for j in range(len(nc_comb_list)):
            analyses.append([i, nc_comb_list[j]])

    for an_index in bar(range(len(analyses))):
        # unpack the index of the test sequence and of the noncoding indices
        i, nc_indices = analyses[an_index]
        seq_set_2 = [nc_seqsets[k][nc_indices[k]]for k in range(len(nc_indices))]
        log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(test_seqlist[i], seq_set_1, seq_set_2, pseudocounts = 0.1, verbose = False, codon_filtering = codon_filtering)
        llratio.append(log_likelihood_ratio)
    np.save("{}_llratio".format(test_gene_name), np.array(llratio))


def main():
    seqs_per_group = 2
    #first generate a list of individuals of the sequence to test
    unk1 = gen_test_seqlist("unk1_bothind.fasta")
    nd2l = gen_test_seqlist("nd2l_bothind.fasta")
    test_seq_names = ["unk1", "nd2l"]
    #Then get a list of the known sequences
    seq_set_1 = gen_codings("c_seqs_bothind.fasta", seqs_per_group)
    # then get the list of lists of noncoding sequences.
    # We first have to pass the argument of how many copies of each gene there
    #  are so that the program adds the sequences in groups, this determines how many
    #  sequences to add in each sublist.
    seq_set_2 = gen_noncodings("nc_seqs_gt50_bothind.fasta", seqs_per_group)
    #test_unknown(unk1, seq_set_1, seq_set_2, 2, "unk1")
    #test_unknown(nd2l, seq_set_1, seq_set_2, 2, "nd2l")


    #print("confusion matrix without removing starts and stops (codon_filtering = 0)")
    #confusion_matrix(seq_set_1, seq_set_2, seqs_per_group, 0)
    #print("confusion matrix with removing starts and stops (codon_filtering = 1)")
    #confusion_matrix(seq_set_1, seq_set_2, seqs_per_group, 1)
    print("confusion matrix with removing starts from beginning and stops from the end and internally (codon_filtering = 2)")
    #confusion_matrix(seq_set_1, seq_set_2, seqs_per_group, 2)
    #print("confusion matrix with removing all starts and stops from all positions (codon_filtering = 3)")
    #confusion_matrix(seq_set_1, seq_set_2, seqs_per_group, 3)
    #codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.5, verbose = True)
    plot_ratios(test_seq_names)

if __name__ == "__main__":
    sys.exit(main())
