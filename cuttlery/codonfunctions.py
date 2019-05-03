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
import os, time
import pandas as pd
import numpy as np

# biopython imports
from Bio import SeqIO
from Bio.codonalign import CodonSeq
import Bio.Data.CodonTable
from Bio.codonalign.codonalphabet import get_codon_alphabet
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import _get_pi
import matplotlib.pyplot as plt

def calculate_average_percent_difference(seqs):
    """
    there is some info here http://www.columbia.edu/cu/biology/courses/c3020/solutions-2.html
    calculates pi of a list of codonseq objects

    this method does not factor indels into the calculation. all sites
    with indels will be ignored.
    """
    df = seqs_to_df(seqs)
    #df.apply(rowwise_unique, axis=1)

    # this is called x due to tradition in the definition of pi
    x = seqfreqs(seqs)
    percent_difference = []
    for i in range(len(seqs)):
        for j in range(i+1, len(seqs)):
            comp = df.copy().iloc[:,[i,j]]
            comp.replace('-', np.nan, inplace=True)
            comp.dropna(axis=0,how='any', inplace=True)
            num_difs = sum(comp.iloc[:,0] != comp.iloc[:,1])
            len_seqs = len(comp.iloc[:,0])
            percent_difference.append(num_difs/len_seqs)
    return np.mean(percent_difference)

def calculate_pi(seqs):
    """
    there is some info here http://www.columbia.edu/cu/biology/courses/c3020/solutions-2.html
    calculates pi of a list of codonseq objects

    this method does not factor indels into the calculation. all sites
    with indels will be ignored.
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

def calculate_piN_piS(codonseqs, method, codon_table, het=False):
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
            if not het:
                dN, dS = cal_dn_ds(codonseqs[i], codonseqs[j], codon_table=codon_table, method=method)
                piN = piN + (x[i] * x[j] * dN)
                piS = piS + (x[i] * x[j] * dS)
                #if 'piNpiS' in options.debug:
                #    print("{0} dN{1}{2}={3} dS{1}{2}={4}".format(method, i, j, dN, dS))
            else:
                try:
                    dN, dS = cal_dn_ds(codonseqs[i], codonseqs[j], codon_table=codon_table, method=method)
                    piN = piN + (x[i] * x[j] * dN)
                    piS = piS + (x[i] * x[j] * dS)
                except:
                    pass

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
            #print("this_slice[i][j]:[{0}][{1}] {2}".format(i, j, this_slice))
            this_CodonSeq = CodonSeq(this_slice, alphabet = codon_alphabet)
            this_CodonSeq.id = codon_id
            this_codonseqs_list.append(this_CodonSeq)
            #print("here's a single codon seq: {}".format(this_CodonSeq))
        sliced_codonseqs.append(this_codonseqs_list)

    #print("0th codon is: {}".format(codonseqs[0].get_codon(0)))
    return sliced_codonseqs, num_codons

def fasta_dir_to_gene_filelist(fasta_dirs, debug = False):
    """
    This method looks through directories for alignment files and returns
     a dictionary of {<gene_name>: </path/to/algn/file/gene_name.fasta>}.

    Gene names that are returned in the dictionary object are taken from the
     filenames. To ensure that downstream plots are readable, use simple file
     names for the fasta alignment files.

    This accepts a single directory path as a string or an iterable
     of directory path strings."

    Parameters:
    <fasta_dirs> - a directory as a string or an iterable that contains
     directory paths as strings.

    Returns:
    - A dictionary of {<gene_name>: </path/to/algn/file/gene_name.fasta>}.
    """
    if debug:
        print("fasta_dirs", fasta_dirs)
    #First test if the input is a string, if so put it into a list and
    # process it in the for loop.
    if isinstance(fasta_dirs, str):
        process_list = [fasta_dirs]
        if debug:
            print("fasta_dirs is a string")
            print("process_list", process_list)
    else:
        try:
            process_list = iter(fasta_dirs)
        except:
            raise Exception("""The fasta directory that you provided, {},
            is neither a directory string nor an iterable that contains
            directories.""".format(fasta_dirs))
    # Make sure that the process list has stuff in it. It is possible
    #  that an empty list made it through.
    if len(fasta_dirs) == 0:
        raise IOError("""The list of fasta directories is length zero.
        For some reason the program received an empty list of directories
        of fasta files""")
    results_dict = {}
    for single_dir in process_list:
        if debug:
            print("single_dir is", single_dir)
        for root, dirs, files in os.walk(single_dir):
            for this_file in files:
                split_file = os.path.splitext(this_file)
                # should be a simple file with at least one period
                if len(split_file) > 1:
                    # only get the fasta files
                    filename  = split_file[0]
                    extension = split_file[-1]
                    if extension in [".fasta", ".fa"]:
                        full_path = os.path.abspath(os.path.join(root, this_file))
                        results_dict[filename] = full_path
    if len(results_dict) == 0:
        raise Exception("this search yielded no fasta files")
    return results_dict

def fasta_path_to_codonseqs(fasta_path, codon_table, codon_alphabet):
    codonseqs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        try:
            this_CS = remove_stop_from_end(CodonSeq().from_seq(record.seq), codon_table)
            this_CS.alphabet = codon_alphabet
            this_CS.id = record.id
            codonseqs.append(this_CS)
        except:
            pass
    return codonseqs

def fasta_path_to_seqs(fasta_path, codon_table=False, codon_alphabet=False):
    """This converts a fasta file into a list of Seq Objects
    """
    seqs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seqs.append(record)
    return seqs

def print_images(base_output_name, image_formats, dpi, transparent=False,
                 no_timestamp = False):
    file_base = base_output_name
    for fmt in image_formats:
        if no_timestamp:
            out_name = "{0}.{1}".format(file_base, fmt)
        else:
            out_name = "{0}_{1}.{2}".format(file_base, timestamp(), fmt)
        try:
            if fmt == 'png':
                plt.savefig(out_name, dpi=dpi, transparent=transparent)
            else:
                plt.savefig(out_name, format=fmt, transparent=transparent)
        except PermissionError:
            # thanks to https://github.com/wdecoster for the suggestion
            print("""You don't have permission to save cuttlery plots to this
            directory. Try changing the directory and running the script again!""")

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
    and each row is one position in the sequence.

    seqs can be a list of biopython seq objects, or codonseqs
    """
    df_list = []
    names_list = []
    for index in range(len(seqs)):
        names_list.append(seqs[index].id)
    if isinstance(seqs[0], CodonSeq):
        newseqs = [[char for char in str(thisseq)] for thisseq in seqs]
    else:
        newseqs = [[char for char in str(thisseq.seq)] for thisseq in seqs]
    df = pd.DataFrame.from_items(zip(names_list, newseqs))
    return df

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

def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")
