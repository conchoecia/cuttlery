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


nts = ['A', 'C', 'G', 'T']

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

    print("Bayes factor is " + str(math.exp(log_likelihood_ratio)))
    print("Hypothesis test favors hypothesis " + str(favored) + " and disfavors hypothesis " + str(unfavored))
    print("Level of support is '" + jeffreys_support + "' on Jeffreys' scale and '" + kass_support + "' on Kass and Raferty's scale")

def dirichlet_log_marginal_likelihood(counts_test, counts_background, pseudocounts):
    codons = ["".join(codon) for codon in itertools.product(nts, repeat = 3)]
    log_like = 0
    log_like += math.lgamma(sum(counts_test[codon] for codon in codons))
    log_like += math.lgamma(sum(pseudocounts + counts_background[codon] for codon in codons))
    log_like -= math.lgamma(sum(pseudocounts + counts_background[codon] + counts_test[codon] for codon in codons))
    for codon in codons:
        print(counts_test[codon])
        log_like -= math.lgamma(pseudocounts + counts_test[codon])
        log_like -= math.lgamma(pseudocounts + counts_background[codon])
        log_like += math.lgamma(pseudocounts + counts_background[codon] + counts_test[codon])
    return log_like


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

def gen_noncodings(nc_fasta_filepath):
    """This function takes a filepath of the fasta file that contains the
    noncoding sequences to use in the analysis.  Instead of manually
    generating all three artificial reading frames, the function
    generates them automatically and returns in list format for
    codon_dirichlet_log_likelihood_ratio()
    """
 
    seqs = []
    for record in SeqIO.parse(nc_fasta_filepath, "fasta"):
        if len(record.id) > 4:
            for i in range(3):
                seqs.append(record.seq[i:])
    for i in range(len(seqs)):
        mod = len(seqs[i]) % 3
        thisLen = len(seqs[i])
        seqs[i] = str(seqs[i][0:thisLen-mod])
    return seqs


seq_set_1 = []
null      = []

def main():
    #first generate a string of the sequence to test
    test_seq = str(next(SeqIO.parse("Bf201606_COIII.fasta", "fasta")).seq)
    #Then get a list of the known sequences
    seq_set_1 = []
    for record in SeqIO.parse("Bf201606_COIII.fasta", "fasta"):
        seq_set_1.append(str(record.seq))
    #then get the list of noncoding sequences
    seq_set_2 = gen_noncodings("Bf_nc_fwd.fasta")
    #print(seq_set_2)
    #for seq in seq_set_2:
    #    print(len(seq) % 3)
    codon_dirichlet_log_likelihood_ratio(test_seq, seq_set_1, seq_set_2, pseudocounts = 0.5, verbose = True)

if __name__ == "__main__":
    sys.exit(main())
