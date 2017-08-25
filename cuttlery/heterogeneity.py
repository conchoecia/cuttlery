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
#
#Created on Thu Jun  1 14:54:48 2017
#
#@author: Jordan/DTS
#info about confusion matrix here: http://www.dataschool.io/simple-guide-to-confusion-matrix-terminology/

#Python/System stuff
import sys

# Biopython stuff
from Bio import SeqIO
import Bio.Data.CodonTable
from Bio.codonalign.codonalphabet import get_codon_alphabet
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import _get_pi


# cuttlery stuff
from cuttlery.codonfunctions import fasta_dir_to_gene_filelist, fasta_path_to_codonseqs,\
    seqs_to_df, calculate_pi, seqfreqs, calculate_piN_piS


def timestamp():
    """
    Returns the current time in :samp:`YYYY-MM-DD HH:MM:SS` format.
    """
    return time.strftime("%Y%m%d_%H%M%S")

def run(args):
    heterogeneity(args)

def heterogeneity(options):
    print(options)
    # if the user wants to know which codon usage tables are available, print them
    if options.tt_options:
        print("Options for codon usage tables in --tt_code are:")
        for key in sorted(Bio.Data.CodonTable.generic_by_id):
            print("{}: {}".format(key, Bio.Data.CodonTable.generic_by_id[key].names))
        sys.exit()

    #1.5 get a list of files from the directory we provided.
    # This is a dict object with filesnames as keys
    filelist = fasta_dir_to_gene_filelist(options.fasta_dir)
    print(filelist)

    #second, select the codon alphabet to use
    codon_alphabet = get_codon_alphabet(Bio.Data.CodonTable.generic_by_id[options.tt_code], gap_char="-")
    codon_table =    Bio.Data.CodonTable.generic_by_id[options.tt_code]

    # now iterate through all genes in the file list and calculate piNpiS for
    #  each codon in the alignment.
    for genename in filelist:
        codonseqs = fasta_path_to_codonseqs(filelist[genename], codon_table, codon_alphabet)
        print(codonseqs)
