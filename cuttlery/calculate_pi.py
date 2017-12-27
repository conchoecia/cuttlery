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
title: cuttlery calculate-pi
author: Darrin T Schultz (github@conchoecia)

This program calculates the value of little pi, also known as
nucleotide diversity, of the sequences in a fasta alignment.

The most sensible way to use this program is to make one fasta
alignment with your favorite fasta alignment program and to save it as
a fasta alignment. Alternatively, you can make a fasta alignment and
save each sequence to a separate fasta file.

Then, calculate little pi using cuttlery with

    cuttlery calculate-pi --fasta_aln <fasta_aln1.fasta> <fasta_aln2.fasta> <et cetera>

The output will be printed to std out.
"""
from cuttlery.codonfunctions import fasta_path_to_seqs, calculate_pi,\
                                    calculate_average_percent_difference
import os

def pi(args):
    print(args)
    print("# cuttlery calculate-pi")
    for fastapath in args.fasta_aln:
        # Read in a fasta file and turn into codonseqs
        codonseqs = fasta_path_to_seqs(fastapath)
        # first we calculate the real information about these data
        pi = calculate_pi(codonseqs)
        print("   pi for {}: {:0.5f}".format(os.path.basename(fastapath), pi))
    for fastapath in args.fasta_aln:
        # Read in a fasta file and turn into codonseqs
        codonseqs = fasta_path_to_seqs(fastapath)
        # first we calculate the real information about these data
        avg_per_dif = calculate_average_percent_difference(codonseqs)
        print("   avg % difference for {}: {:0.5f}".format(os.path.basename(fastapath),
                                                           avg_per_dif))

def run(args):
    pi(args)
