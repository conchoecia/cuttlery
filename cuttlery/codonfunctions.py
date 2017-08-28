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
filename: codon_functions.py
author: Darrin T. Schultz (github @conchoecia)
purpose: common functions shared by codon-related things in cuttlery
"""

# imports
import os
import pandas as pd
import numpy as np

# biopython imports
from Bio import SeqIO
from Bio.codonalign import CodonSeq
import Bio.Data.CodonTable
from Bio.codonalign.codonalphabet import get_codon_alphabet
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import _get_pi


def fasta_dir_to_gene_filelist(fasta_dir):
    abspath = os.path.abspath(fasta_dir)
    filelist = {os.path.splitext(x)[0]:os.path.join(abspath, x)
                for x in os.listdir(fasta_dir)
                if os.path.splitext(x)[1]}
    return filelist

def fasta_path_to_codonseqs(fasta_path, codon_table, codon_alphabet):
    codonseqs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        this_CS = remove_stop_from_end(CodonSeq().from_seq(record.seq), codon_table)
        this_CS.alphabet = codon_alphabet
        this_CS.id = record.id
        codonseqs.append(this_CS)
    return codonseqs

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

def codonseqs_sliced(codonseqs, windowsize):
    """
    Arguments:
      <codonseqs>  A list of Biopython CodonSeq() objects, each of the same length
      <windowsize> How many codons to put into each analysis window
    Returns:
      A list of lists of CodonSeq() objects, each CodonSeq() object being a
       codon slice of size <windowsize>

    First, make sure that all of the CodonSeq() objects have the same length,
      if not throw an exception. This will violate the assumptions of using the
      seqs to measure pi, piN, and piS.

    Downstream, the list of lists of CodonSeq() objects can be passed
     iteratively to Bio.codonalign.codonseq module's cal_dn_ds() and _get_pi()
     functions, or to the calculate_pi() and calculate_piN_piS() functions herein.
    """
    # TODO: CodonSeq object can't be initialized with an id class attribute
    # TODO: CodonSeq.get_codon can't be accessed with slices. How about turning it
    # Added try:except blocks to the _ng86 method
    #  into a slice-able attribute or a function that accepts a range?

    # First, verify that all of the CodonSeq() objects are the same length
    seqlens = [len(x) for x in codonseqs]
    if not len(set(seqlens)) == 1:
        raise IOError("""After removing the stop codons from the 3' ends of
        your sequences, not all of your alignments were the same length.
        Please make sure that your ORF alignments do not contain any
        5' or 3' overhangs and try again. Please submit an issue to
        https://github.com/conchoecia/cuttlery/issues if you would like the
        program to handle this problem differently.""")
    # now that we have verified that all of the sequences are the same length,
    #  determine how many codons there are in each sequence.
    num_codons = codonseqs[0].get_codon_num()
    # extrapolate the alphabet and id to use for all the sliced sequences
    codon_id = codonseqs[0].id
    codon_alphabet = codonseqs[0].alphabet
    print("num_codons: {}".format(num_codons))

    # Iteratively extract all of the windows of codons, and the result
    #  is that each element of slied_codonseqs is one of the codonseq slices
    sliced_codonseqs = []
    for i in range(0,num_codons, windowsize):
        this_codonseqs_list = []
        for j in range(len(codonseqs)):
            this_slice = ""
            # slicing doesn't work yet for CodonSeq.get_codon so we need to get
            #  codons iteratively
            for k in range(windowsize):
                this_slice += codonseqs[j].get_codon(i+k)
            #print("this_slice: {}".format(this_slice))
            this_CodonSeq = CodonSeq(this_slice, alphabet = codon_alphabet)
            this_CodonSeq.id = codon_id
            this_codonseqs_list.append(this_CodonSeq)
        sliced_codonseqs.append(this_codonseqs_list)

    #print("0th codon is: {}".format(codonseqs[0].get_codon(0)))
    return sliced_codonseqs, num_codons

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
            num_difs = sum(comp.iloc[:,0] != comp.iloc[:,1])
            len_seqs = len(comp.iloc[:,0])
            this_sum = x[i] * x[j] * (num_difs/len_seqs)
            running_sum_pi += this_sum
            #if 'pi' in options.debug:
            #    print("x[i] * x[j] * pi = {} * {} * {}/{}".format(x[i], x[j], num_difs, len_seqs))
    #if 'pi' in options.debug:
    #    print("pi: {}".format(running_sum_pi))
    return running_sum_pi

def seqfreqs(seqs):
    """Calculates the relative frequency of each sequence for calculating pi
    """
    #if "seqfreqs" in options.debug:
    #    print("There are {} seqs".format(len(seqs)))
    x = []
    #this block calculates the frequencies of each sequence
    for i in range(len(seqs)):
        this_x = 0
        for j in range(len(seqs)):
            if str(seqs[i]) == str(seqs[j]):
                #if "seqfreqs" in options.debug:
                #    print("{} == {}".format(i, j))
                this_x += 1
        x.append(this_x/len(seqs))
    #print("done with these seqfreqs\n")
    #if "seqfreqs" in options.debug:
    #    print("the frequencies are {}".format(x))
    return x

def calculate_piN_piS(codonseqs, method, codon_table):
    """
    takes a list of CodonSeq() objects and calculates piN, piS, pi, and piNpiS
    for them
    """
    analysis = {"seqname": "", "piN": -1, "piS": -1, "piNpiS": -1, "pi": -1, "method":method}
    x = seqfreqs(codonseqs)
    #if 'piNpiS' in options.debug:
    #    print("freqs are: {}".format(x))
    #    print("len codonseqs is: ", len(codonseqs))
    piN = 0
    piS = 0
    for i in range(len(codonseqs)):
        for j in range(i+1, len(codonseqs)):
            #print(codonseqs[i], codonseqs[j])
            dN, dS = cal_dn_ds(codonseqs[i], codonseqs[j], codon_table=codon_table, method=method)
            piN = piN + (x[i] * x[j] * dN)
            piS = piS + (x[i] * x[j] * dS)
            #if 'piNpiS' in options.debug:
            #    print("{0} dN{1}{2}={3} dS{1}{2}={4}".format(method, i, j, dN, dS))

    analysis['piN'] = piN
    analysis['piS'] = piS
    try:
        analysis['piNpiS'] = piN/piS
    except:
        analysis['piNpiS'] = 0
    #if 'piNpiS' in options.debug:
    #    print ("{0} dN={1:.3f} dS={2:.3f} piN/piS = {3:.3f}".format(
    #        method, analysis['piN'], analysis['piS'], analysis['piNpiS']))

    return analysis
