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
from multiprocessing import cpu_count
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool



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
global analyses
analyses = set()

def determine_pool_size(job_vector):
    """This function determines how large of a pool to make based on the
    system resources currently available and how many jobs there are to complete
    """
    available_threads = cpu_count()
    total_jobs = len(job_vector)
    threads_to_pass = total_jobs
    if total_jobs >= available_threads:
        threads_to_pass = available_threads
    if threads_to_pass > 90:
        threads_to_pass = 90
    print("There are {} threads available.\nThere are {} jobs:\nMaking a pool with {} threads".format(
        available_threads, total_jobs, threads_to_pass), file = sys.stderr)
    return threads_to_pass

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
    # the final list is organized such that each index of the grouped_seqs has a
    # single sequence, and each list at each index contains a sequence from one
    # individual
    grouped_seqs = [seqs[i:i+seqs_per_set] for i in range(0,len(seqs),seqs_per_set)]
    return grouped_seqs

def generate_remove_one_substitution(nontest_seqs_set, test_seqs_set,
                                     number_of_analyses):
    """This method returns a a list of tuples of combinations for
    remove-one tests.  The input format of a seqs set is a list of
    lists. Each sublist contains the same locus from multiple individuals.
    For example,
     seqs_set = [[gene1_ind1, gene1_ind2, gene1_ind3],
                 [gene2_ind1, gene2_ind2, gene2_ind3],
                 [gene3_ind1, gene3_ind3, gene3_ind3]]

    This analysis outputs a list of lists. Each list contains at each index:
      - index 0: The indices for the nontest sequences. a list of tuples.
                 each tuple's 0th index is the index of the the list in which the
                 sequence resides and the 1st index is the index of the gene in
                 the 0th list.
      - index 1: the (ith gene, jth individual) element used for the test sequence
      - index 2: a list of (i, j) elements that does not include the ith gene at
                  index 1 described above.
    """
    # all of the analyses will be stored in this list
    analysis_list = []
    num_nontest_seqs = len(nontest_seqs_set)
    num_nontest_individuals = len(nontest_seqs_set[0])
    num_test_seqs = len(test_seqs_set)
    num_individuals = len(test_seqs_set[0])
    # now we iterate through all of the combinations of nontest sequences
    for index in range(number_of_analyses):
        # randomly generate the test i,j
        random_i = np.random.randint(0,num_test_seqs)
        random_j = np.random.randint(0,num_individuals)
        nontest_indices = np.random.choice(num_nontest_individuals,
                                           num_nontest_seqs, replace = True)
        test_indices = np.random.choice(num_individuals,
                                        num_test_seqs - 1, replace = True)
        # In this next bit of code, I use final_index_count to keep
        #  track of how many genes we have iterated through. This is
        #  mainly used to skip over index i when we encounter it.
        final_index_count = 0
        test_indices_index_count = 0
        test_tuples = []
        while final_index_count < num_test_seqs:
            if final_index_count == random_i:
                #if we're at the test position(i), skip over it
                final_index_count += 1
            else:
                # the 0th index of final_index_count is the index of
                #  the gene in test_seqs_set, and
                #  test_indices[test_indices_index_count] is the individual's
                #  index in the gene seq set
                one_test_gene = (final_index_count, test_indices[test_indices_index_count])
                test_tuples.append(one_test_gene)
                #increase both the indices
                final_index_count += 1
                test_indices_index_count += 1
        this_analysis = (nontest_indices, (random_i, random_j), test_tuples)
        analysis_list.append(this_analysis)
    return analysis_list


def generate_seq_lists_from_analyses(analysis_list, nontest_seqs_set, test_seqs_set):
    nontest_seq_list = []
    test_seq = ""
    test_seq_list    = []
    # unpack everything
    nontest_indices, test_seq_indices, test_indices_tuples = analysis_list
    # first go through the nontest_indices
    for i in range(len(nontest_indices)):
        this_seq = nontest_seqs_set[i][nontest_indices[i]]
        nontest_seq_list.append(this_seq)
    # now get the test sequence
    i = test_seq_indices[0]
    j = test_seq_indices[1]
    test_seq = test_seqs_set[i][j]
    # now get all the test sequences
    for i,j in test_indices_tuples:
        this_seq = test_seqs_set[i][j]
        test_seq_list.append(this_seq)

    # now we return the nontest sequences, the test sequence, and the test sequences
    return (nontest_seq_list, test_seq, test_seq_list)


def confusion_matrix(coding_seqsets, nc_seqsets, seqs_per_set,
                     codon_filtering, number_of_analyses):
    #first, test all the trues. Generate a list of all combinations of noncoding
    # and for each combination iterate through all of the groups of coding, and for each
    # set of coding remove one and test the other against everything else.
    # That is three for loops.
    real_val = []
    observed = []
    llratio  = []
    test_dicts = [{"nontest_seqs_set": nc_seqsets,
                      "test_seqs_set": coding_seqsets,
                          "test_type": "coding"},
                  {"nontest_seqs_set": coding_seqsets,
                      "test_seqs_set": nc_seqsets,
                          "test_type": "noncoding"}]
    for this_dict in test_dicts:
        print("I'm doing some analyses")
        nontest_seqs_set = this_dict["nontest_seqs_set"]
        test_seqs_set =    this_dict["test_seqs_set"]
        bar = progressbar.ProgressBar()
        # Now we test all of the true sequences
        print("im coming up with the analyses")
        analyses = generate_remove_one_substitution(nontest_seqs_set, test_seqs_set,
                                                    number_of_analyses)

        for i in bar(range(len(analyses))):
            analysis_list = analyses[i]
            nontest_seqs, test_seq, test_seqs = generate_seq_lists_from_analyses(analysis_list, nontest_seqs_set, test_seqs_set)
            if this_dict["test_type"] == "coding":
                seq_set_1 = test_seqs
                seq_set_2 = nontest_seqs
                real_val.append(1)
            elif this_dict["test_type"] == "noncoding":
                seq_set_1 = nontest_seqs
                seq_set_2 = test_seqs
                real_val.append(2)
            log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(test_seq,
                                    seq_set_1, seq_set_2, pseudocounts = 0.1,
                                    verbose = False, codon_filtering = codon_filtering)
            if log_likelihood_ratio < 0:
                favored = 2
            else:
                favored = 1
            observed.append(favored)
            llratio.append(log_likelihood_ratio)

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
    plt.show()


def test_unknown(test_seqlist, coding_seqsets, nc_seqsets,
                 codon_filtering, test_gene_name, num_simulations):
    #for this test, we need to generate all possible combinations of nc set,
    # plus using the same coding set for each analysis. This is very similar to
    # testing the coding sequences above.
    llratio = []

    bar = progressbar.ProgressBar()
    for i in bar(list(range(num_simulations))):
        # pick a random test sequence
        seq_i = np.random.randint(0,len(test_seqlist))
        #use all the coding indices
        coding_seq_indices = np.random.choice(len(coding_seqsets[0]),
                                          len(coding_seqsets), replace = True)

        #These are all coding sequences
        seq_set_1 = [coding_seqsets[i][coding_seq_indices[i]]
                     for i in range(len(coding_seq_indices))]
        #use all the noncoding indices
        nc_seq_indices = np.random.choice(len(nc_seqsets[0]),
                                          len(nc_seqsets), replace = True)
        #These are all noncoding sequences
        seq_set_2 = [nc_seqsets[i][nc_seq_indices[i]] for i in range(len(nc_seq_indices))]
        log_likelihood_ratio = codon_dirichlet_log_likelihood_ratio(test_seqlist[seq_i],
                                    seq_set_1, seq_set_2, pseudocounts = 0.1,
                                    verbose = False, codon_filtering = codon_filtering)
        llratio.append(log_likelihood_ratio)
    np.save("{}_llratio".format(test_gene_name), np.array(llratio))


def main():
    # FOR INTERNAL CONSISTENCY, ANY SET_1 IS CODING and SET_2 is NONCODING!
    seqs_per_set = 3
    codon_filtering = 2
    num_analyses = 1111111
    #first generate a list of individuals of the sequence to test
    unk1 = gen_test_seqlist("unk1_threeind.fasta")
    nd2l = gen_test_seqlist("nd2l_threeind.fasta")
    test_seq_names = ["unk1", "nd2l"]
    #Then get a list of the known sequences
    coding_seqsets = gen_codings("bf_coding_threeindividuals.fasta", seqs_per_set)
    nc_seqsets = gen_noncodings("bf_noncoding_threeindividuals.fasta", seqs_per_set)
    # then get the list of lists of noncoding sequences.
    # We first have to pass the argument of how many copies of each gene there
    #  are so that the program adds the sequences in groups, this determines how many
    #  sequences to add in each sublist.
    test_unknown(unk1, coding_seqsets, nc_seqsets, 2, "unk1", num_analyses)
    test_unknown(nd2l, coding_seqsets, nc_seqsets, 2, "nd2l", num_analyses)


    #print("confusion matrix without removing starts and stops (codon_filtering = 0)")
    #confusion_matrix(seq_set_1, seq_set_2, seqs_per_set, 0)
    #print("confusion matrix with removing starts and stops (codon_filtering = 1)")
    #confusion_matrix(seq_set_1, seq_set_2, seqs_per_set, 1)
    print("confusion matrix with removing starts from beginning and stops from the end and internally (codon_filtering = 2)")
    confusion_matrix(coding_seqsets, nc_seqsets, seqs_per_set, codon_filtering, num_analyses)
    #print("confusion matrix with removing all starts and stops from all positions (codon_filtering = 3)")
    #confusion_matrix(seq_set_1, seq_set_2, seqs_per_set, 3)
    plot_ratios(test_seq_names)


if __name__ == "__main__":
    sys.exit(main())
