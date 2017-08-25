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
