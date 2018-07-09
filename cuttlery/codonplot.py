#!/usr/bin/env python3

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
title: cuttlery codonplot
author: Darrin T Schultz - github@conchoecia

This program plots the codon usage of a collection of genes as violin
plots.  This style of plot more accurately reflects the biological
variability in codon usage in different protein types. In addition,
one can plot single or multiple genes as dots to highlight their
position in the distribution of codon usage frequencies.
"""
from collections import Counter
import itertools
import operator
import os
import sys
import pandas as pd
pd.options.display.float_format = '{:,.3f}'.format
import urllib.request

# Biopython stuff
from Bio import SeqIO
import Bio.Data.CodonTable

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.colors as colors
from scipy.stats import mannwhitneyu as mu
from collections import Counter
from matplotlib import rcParams
import numpy as np
import sys
import os
import time

from matplotlib import rc
rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{tgheros}',    # helvetica font
        r'\usepackage{sansmath}',   # math-font matching  helvetica
        r'\sansmath'                # actually tell tex to use it!
        r'\usepackage{siunitx}',    # micro symbols
        r'\sisetup{detect-all}',    # force siunitx to use the fonts
        ]

# cuttlery stuff
from cuttlery.codonfunctions import fasta_dir_to_gene_filelist,\
                                    fasta_path_to_codonseqs, print_images
import cuttlery.rcparams as rc

def fasta_dir_to_codondf(fasta_dir, codon_table, coding=True):
    """This method takes a fasta directory as input and outputs a pandas
    dataframe that contains the codon counts, the codon frequencies per amino
    acid, the seqname"""
    # 1. Get the list of fasta files
    gene_filelist = fasta_dir_to_gene_filelist(fasta_dir)
    # 2. For each fasta record, remove the gaps, count codons, add gene name
    codon_data = []
    for seqname_key in gene_filelist:
        for record in SeqIO.parse(gene_filelist[seqname_key], "fasta"):
            # make all of the codons
            this_dict = {"".join(x): 0 for x in list(itertools.product('GATC', repeat = 3))}
            start_codon_dict = {"${}".format(x):0
                                for x in codon_table.start_codons if 'U' not in x}
            stop_codon_dict  = {"*{}".format(x):0
                                for x in codon_table.stop_codons if 'U' not in x}
            this_dict['seqname'] = seqname_key
            fasta_path = gene_filelist[seqname_key]
            this_seq = str(record.seq).upper().replace('-', '')
            # assert that the sequence is divisible by three if coding
            if coding:
                this_dict["seqtype"] = "coding"
                # assert that the sequence must be divisible by three
                if len(this_seq) % 3 != 0:
                    raise IOError("""The sequence {} in {} is a coding sequence, but
                    is not divisible by three after removing gaps. Please check the
                    sequences and try again""".format(seqname_key, fasta_path))
                # assert that the first codon is a start codon
                start_codon = "${}".format(this_seq[0:3])
                if start_codon not in start_codon_dict.keys():
                    raise IOError("The sequence {} in {} is coding, but does not begin with a start codon. Please fix the sequence and try again.""".format(seqname_key, fasta_path))
                stop_codon = "*{}".format(this_seq[-3::])
                if stop_codon not in stop_codon_dict.keys():
                    raise IOError("""The sequence {} in {} is coding, but does
                    not end with a stop codon. Please fix the sequence and
                    try again.""".format(seqname_key, fasta_path))
                start_codon_dict[start_codon] += 1
                stop_codon_dict[stop_codon] += 1
                codons = [this_seq[i:i+3] for i in range(3,len(this_seq) - 3,3)]
            else:
                # TODO: Implement codon count for noncoding sequences by generating
                this_dict["seqtype"] = "noncoding"
                codons = [this_seq[i:i+3] for i in range(0,len(this_seq),3)]
            #  all three frames and removing internal stop codons
            # add the codon counts to this_dict
            for codon in codons:
                this_dict[codon] += 1
            # now combine the codon counts, the start and stop counts
            super_dict = {}
            # https://stackoverflow.com/questions/3494906
            for d in [this_dict,start_codon_dict, stop_codon_dict]:
                super_dict.update(d)
            codon_data.append(super_dict)
    return pd.DataFrame.from_dict(codon_data)

def plot_codon_usage(freqs, codon_table, CTSAC, title, invert):

    # codon_coords has the codon as the key and the value is the x pos
    codon_coords = {}
    #set the figure dimensions
    figWidth = 7
    figHeight = 3
    plt.figure(figsize=(figWidth,figHeight))

    #set the panel dimensions
    panelWidth = 6
    panelHeight = 2

    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = ((figHeight - panelHeight)/2)

    panel0=plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='off', labelbottom='off',\
                       left='on', labelleft='on', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    # determine which label colors to use
    if invert:
        labelcolor = 'white'
        labelcolorlight = (1,1,1,0.75)
        panel0.spines['left'].set_color('white')
        panel0.yaxis.label.set_color('white')
        panel0.tick_params(axis='y', colors='white')
    else:
        labelcolor = 'black'
        labelcolorlight = (0,0,0,0.75)

    # This block sets the plot order for the codons based on alphabetic order
    #  of amino acids and of the codons
    aas = {aa:[] for aa in set(codon_table.forward_table[x]\
                            for x in codon_table.forward_table)}
    codons = [x for x in codon_table.forward_table if 'U' not in x]
    print("codons", codons)
    for codon in codons:
        aa = codon_table.forward_table[codon]
        aas[aa].append(codon)
    codon_plot_order = []
    for aa in sorted(aas):
        codons = sorted(aas[aa])
        for codon in codons:
            codon_plot_order.append(codon)
        print(aa, codons)

    # rectangle y starts and stops
    yStart = -0.25
    yStop = -0.05
    height = abs(yStop - yStart)
    horizontal_text_y_pos = -0.185


    # First, sort out the start codons and plot those as frequencies
    starts_df = freqs['starts'][list(sorted(freqs['starts'].index))]
    print(starts_df)
    starts_pos_0 = 0
    scaling_factor = 0.5
    starts_pos_end = scaling_factor * len(freqs['starts'].index)
    # plot the rectangles for the calculated frequencies
    width = starts_pos_end - starts_pos_0 - 0.5
    bar_thickness = 0.01
    for i in range(len(starts_df.index)):
        if starts_df[i] != 0:
            rect = mplpatches.Rectangle((starts_pos_0 + (i*scaling_factor) - 0.5,starts_df[i]-(bar_thickness/2)),
                                        scaling_factor,bar_thickness,\
                                    linewidth=0.00,\
                                    facecolor=labelcolor,\
                                    edgecolor=(0,0,0,0),\
                                    alpha=1)
            panel0.add_patch(rect)
    # plot a solid box around the starts
    rect = mplpatches.Rectangle((starts_pos_0 - 0.5, yStart+0.005),
                                len(starts_df.index)*scaling_factor,1 + abs(yStart)-0.01,\
                                linewidth=0.25,\
                                facecolor=(1,1,1,0),\
                                edgecolor=labelcolorlight)
    panel0.add_patch(rect)
    # now plot the labels for the start codons
    rotlabels = ["{0}\n{1}\n{2}".format(x[1],x[2],x[3]) for x in list(sorted(starts_df.index))]
    print("rotlabels", rotlabels)
    for i in range(len(rotlabels)):
        label = rotlabels[i]
        # Don't add one
        print("starts_pos", starts_pos_0)
        xpos = starts_pos_0 + (( i + scaling_factor) * scaling_factor ) - 0.5
        codon_coords[starts_df.index[i]] = xpos
        if starts_df[i] == 0:
            ypos = -0.06
        else:
            ypos = starts_df[i] - 0.02
        panel0.text(xpos, ypos, label,
                    horizontalalignment='center',
                    color = labelcolor,
                    verticalalignment='top',
                    size = 4)
    xpos = (starts_pos_end - starts_pos_0)/2 - 0.5
    panel0.text(xpos, horizontal_text_y_pos, "Starts",
                color=labelcolor,
                horizontalalignment='center',
                verticalalignment='top',
                size = 6)


    # Second sort the codons and plot those
    #print(freqs)
    codon_df = freqs['codons'][codon_plot_order]
    print(codon_df.columns.values)
    data = codon_df.transpose().values.tolist()
    #print(data)
    codon_pos_0 = starts_pos_end
    codon_pos_end = starts_pos_end + len(codon_plot_order)

    # colors mostly taken from https://biologywarakwarak.files.wordpress.com/2012/01/amino-acid-table-singlet-code.png
    colormap = {"F": '#c0d9ef',
                "L": '#e8f1df',
                "I": '#f7dedc',
                "M": '#ff2e00',
                "V": '#ffc239',
                "S": '#ffff54',
                "P": "#7fce66",
                "T": "#00ae60",
                "A": "#00aeec",
                "Y": "#006fbb",
                "*": "#ffffff",
                "H": "#ded9c5",
                "Q": "#ffc294",
                "N": "#b5a2c4",
                "K": "#968b5a",
                "D": "#00fc65",
                "E": "#00dcf0",
                "C": "#ff994e",
                "W": "#dc31e6",
                "R": "#d8d8d8",
                "G": "#abdce7"}
    AAcoords={}
    uniqueAAs = aas.keys()
    print("aas",sorted(uniqueAAs))

    AAcoords_start = codon_pos_0
    for aa in sorted(uniqueAAs):
        # add one to the end of these two because we start plotting at 1, not 0
        firstx = AAcoords_start
        lastx  = AAcoords_start + len(aas[aa])
        AAcoords_start = lastx
        AAcoords[aa] = (firstx, lastx)

    for key in sorted(AAcoords):
        xStart = AAcoords[key][0]-0.5
        xStop = AAcoords[key][1]-0.5
        print(key, xStart, xStop)
        width = xStop - xStart
        facec = colors.to_rgba(colormap[key])
        if key != '*':
            rectangle1=mplpatches.Rectangle((xStart,yStart),width,height,\
                                               linewidth=0.0,\
                                               facecolor=facec,\
                                               edgecolor=(0,0,0),\
                                               alpha=0.8)
            rectangle2=mplpatches.Rectangle((xStart,yStop),width,1+abs(yStop),\
                                               linewidth=0.00,\
                                               facecolor=facec,\
                                               edgecolor=(0,0,0),\
                                               alpha=0.2)
            panel0.add_patch(rectangle2)
        else:
            rectangle1=mplpatches.Rectangle((xStart,yStart),width,height,\
                                               linewidth=0.25,\
                                               facecolor=facec,\
                                               edgecolor=(0,0,0),\
                                               alpha=0.8)
        panel0.add_patch(rectangle1)
    #now add the text labels on top of the boxes
    labels = codon_plot_order
    rotlabels = ["{0}\n{1}\n{2}".format(x[0],x[1],x[2]) for x in labels]
    for i in range(len(rotlabels)):
        label = rotlabels[i]
        # Don't add one
        xpos = codon_pos_0 + i
        codon_coords[labels[i]] = xpos
        ypos = -0.06
        panel0.text(xpos, ypos, label,
                    color = labelcolor,
                    horizontalalignment='center',
                    verticalalignment='top',
                    size = 4)
    # Now plot the three-letter or one-letter code
    threeletcode = {'*': 'Sto',
                    'A': 'Ala',
                    'R': 'Arg',
                    'N': 'Asn',
                    'D': 'Asp',
                    'C': 'Cys',
                    'E': 'Glu',
                    'Q': 'Gln',
                    'G': 'Gly',
                    'H': 'His',
                    'I': 'Ile',
                    'L': 'Leu',
                    'K': 'Lys',
                    'M': 'Met',
                    'F': 'Phe',
                    'P': 'Pro',
                    'S': 'Ser',
                    'T': 'Thr',
                    'W': 'Trp',
                    'Y': 'Tyr',
                    'V': 'Val'}
    #now add the one-leter or three-letter codes
    for key in AAcoords:
        width = abs(AAcoords[key][1]-AAcoords[key][0])
        xpos = AAcoords[key][1] - (width/2) - 0.5
        ypos = -0.185
        if width > 2:
            label = threeletcode[key]
        else:
            label = key
        panel0.text(xpos, ypos, label,
                    color = labelcolor,
                    horizontalalignment='center',
                    verticalalignment='top',
                    size = 6)

    bp = panel0.violinplot(data, widths = 0.75,
                           positions = np.arange(codon_pos_0, codon_pos_end),
                           bw_method = 0.2, showextrema=False)
    for pc in bp['bodies']:
        if invert:
            pc.set_facecolor('white')
        else:
            pc.set_facecolor('black')
        pc.set_zorder(10)
        pc.set_alpha(0.9)
        pc.set_edgecolor(None)

    # Third, plot the stop codons as frequencies
    stops_df = freqs['stops'][list(sorted(freqs['stops'].index))]
    print(stops_df)
    stops_pos_0 = codon_pos_end
    stops_pos_end = stops_pos_0 + (len(freqs['stops'].index) * scaling_factor)
    # plot the rectangles for the calculated frequencies
    width = stops_pos_end - stops_pos_0 - 0.5
    bar_thickness = 0.01
    for i in range(len(stops_df.index)):
        if stops_df[i] != 0:
            rect = mplpatches.Rectangle((stops_pos_0 + (i*scaling_factor) - 0.5,
                                         stops_df[i]-(bar_thickness/2)),
                                    1 * scaling_factor,bar_thickness,\
                                    linewidth=0.00,\
                                    facecolor=labelcolor,\
                                    edgecolor=(0,0,0,0),\
                                    alpha=1)
            panel0.add_patch(rect)
    # plot a solid box around the starts
    rect = mplpatches.Rectangle((stops_pos_0 - 0.5, yStart+0.005),
                                len(stops_df.index) * scaling_factor,1 + abs(yStart)-0.01,\
                                linewidth=0.25,\
                                facecolor=(1,1,1,0),\
                                edgecolor=labelcolorlight)
    panel0.add_patch(rect)
    # now plot the labels for the start codons
    rotlabels = ["{0}\n{1}\n{2}".format(x[1],x[2],x[3]) for x in list(sorted(stops_df.index))]
    print("rotlabels", rotlabels)
    for i in range(len(rotlabels)):
        label = rotlabels[i]
        # Don't add one
        xpos = stops_pos_0 + (( i + scaling_factor) * scaling_factor ) -0.5
        codon_coords[stops_df.index[i]] = xpos
        if stops_df[i] == 0:
            ypos = -0.06
        else:
            ypos = stops_df[i] - 0.02
        panel0.text(xpos, ypos, label,
                    horizontalalignment='center',
                    color = labelcolor,
                    verticalalignment='top',
                    size = 4)
    xpos = stops_pos_0 + (stops_pos_end - stops_pos_0)/2 - 0.5
    panel0.text(xpos, horizontal_text_y_pos, "*",
                horizontalalignment='center',
                color=labelcolor,
                verticalalignment='top',
                size = 6)
    panel0.set_ylim([-0.25,1.1])
    panel0.set_yticks(np.arange(0, 1.25, 0.25))
    panel0.set_ylabel("Fraction", \
                      color = labelcolor,
                      fontsize = 8,
                      labelpad = 0.5)
    panel0.set_xlim([-0.5, stops_pos_end-0.45])
    for each in panel0.spines:
        panel0.spines[each].set_visible(False)

    panel0.set_title(title, color=labelcolor)

    return panel0, codon_coords

def codon_to_same_aa_codons(codon_table):
    all_codons = [x for x in codon_table.forward_table.keys() if 'U' not in x]
    codon_to_same_aa_dict = {codon: [x for x in all_codons
                                     if codon_table.forward_table[x] == \
                                     codon_table.forward_table[codon]]\
                             for codon in all_codons}
    return codon_to_same_aa_dict

def row_freqs(row, codon_table, CTSAC):
    """ got the idea from here: https://stackoverflow.com/questions/12182744
    """
    rowcopy = row.copy(deep=True)
    for i in range(len(row)):
        key = row.index[i]
        codons = CTSAC[key]
        aa_total = sum(rowcopy.filter(items = codons))
        if aa_total == 0:
            row[i] = 0
        else:
            floatval = row.iloc[i] / aa_total
            row[i] = row.iloc[i] / aa_total
    return row

def counts_to_freqs(counts, codon_table, CTSAC):
    freqs = {}
    # calculate the start and stop frequencies as actual frequencies
    #  since there is only one observation for each from each sequence
    # *** STARTS ***
    start_codons = ["${}".format(x) for x in codon_table.start_codons if 'U' not in x]
    start_codons_df = counts[start_codons]
    total_sum = int(start_codons_df.sum(axis=0).sum())
    freqs['starts'] = start_codons_df.sum(axis=0) / total_sum
    # *** STOPS ***
    stop_codons = ["*{}".format(x) for x in codon_table.stop_codons if 'U' not in x]
    stop_codons_df = counts[stop_codons]
    total_sum = int(stop_codons_df.sum(axis=0).sum())
    freqs['stops'] = stop_codons_df.sum(axis=0) / total_sum
    # *** CODONS ***
    cs = pd.DataFrame({x:codon_table.forward_table[x]\
                       for x in codon_table.forward_table.keys()\
                       if 'U' not in x},
                       index=[0])
    codons = cs.iloc[0]
    codons_df = counts[list(codons.index)].astype(float)
    # codon to aa values
    freqs['codons'] = codons_df.apply(row_freqs, args=(codon_table, CTSAC,), axis = 1)
    return freqs

def add_points(panel, codon_coords, freqs):
    print(freqs)
    observations = []
    x_values = []
    # first, get all the points for the start codons
    for i in range(len(freqs['starts'].index)):
        key = freqs['starts'].index[i]
        this_observation = [freqs['starts'][i]]
        this_x = codon_coords[key]
        observations.append(this_observation)
        x_values.append(this_x)
    # now, get all the points for the stop codons
    for i in range(len(freqs['stops'].index)):
        key = freqs['stops'].index[i]
        this_observation = [freqs['stops'][i]]
        this_x = codon_coords[key]
        observations.append(this_observation)
        x_values.append(this_x)
    # now, get the the points for the codons
    for i in range(len(freqs['codons'].columns.values)):
        key = freqs['codons'].columns.values[i]
        this_observation = list(freqs['codons'][key])
        this_x = codon_coords[key]
        observations.append(this_observation)
        x_values.append(this_x)

    #masterPP = _optimize_points(observations, x_values)
    # facec = colors.to_rgba("#00aeec", 0.6)
    facec = colors.to_rgba("#ff2e00", 0.5)
    for i in range(len(observations)):
        for y_value in observations[i]:
            p = panel.plot(x_values[i],y_value,marker='o',ms=3,mfc=facec, mew=0,
                           linewidth=0, zorder = 11)
    return panel

def _optimize_points(observations, centers):
    """This method returns a list of points to plot
    """
    masterPP = []
    for i in np.arange(0,len(observations),1):
        center = centers[i]
        left_bound= i
        right_bound=2 + i
        step_size=0.2
        cutoff=5
        placed_points=[]
        counter=0
        for y_value in observations[i]:
            counter+=1
            if len(placed_points)==0:
                placed_points.append((center, y_value))
            else:
                potential_x_position=[]
                left_search = np.arange(center, left_bound, -1 * step_size)
                right_search = np.arange(center,right_bound, step_size)
                xBinarySearch = [val for pair in zip(left_search, right_search) for val in pair][1:]
                for x_position in xBinarySearch:
                    distances=[]
                    binaryPoints = [tup for tup in placed_points if abs(y_value - tup[1]) < cutoff*2]
                    if len(binaryPoints) == 0:
                        potential_x_position.append(x_position)
                    else:
                        for placed_point in binaryPoints:
                            distance=((x_position-placed_point[0])**2+(y_value-placed_point[1])**2)**0.5
                            distances.append(distance)
                        if min(distances)>cutoff:
                            potential_x_position.append(x_position)
                if len(potential_x_position)>0:
                     best_x_position=sorted(potential_x_position,key=lambda x: np.absolute(x-center))[0]
                     placed_points.append((best_x_position,y_value))
                else:
                     print('point not placed: ',y_value)
        masterPP += placed_points
    return masterPP

def codonplot(args):
    rc.update_rcParams()
    options = args
    print(options)

    # if the user wants to know which codon usage tables are available, print them
    if options.tt_options:
        print("Options for codon usage tables in --tt_code are:")
        for key in sorted(Bio.Data.CodonTable.generic_by_id):
            print("{}: {}".format(key, Bio.Data.CodonTable.generic_by_id[key].names))
        sys.exit()

    if options.coding_fasta_dir == None:
        raise IOError("""You have not provided a directory of fasta alignment files
        or fasta sequence files to produce a codon plot. Please re-run the
        program and pass a directory to --coding_fasta_dir.""")

    codon_table = Bio.Data.CodonTable.generic_by_id[options.tt_code]

    # First get the codon usage of the real sequences
    results_file = "{}.csv".format(options.output_basename)
    if not os.path.exists(results_file):
        results_df = fasta_dir_to_codondf(options.coding_fasta_dir, codon_table,
                                          coding=True)
        #results_df.to_csv(results_file, index = False)
    else:
        print("\nFound {} so skipping analysis\n".format(results_file))
        results_df = pd.read_csv(results_file)

    # get all the frequency info

    CTSAC = codon_to_same_aa_codons(codon_table)
    freqs = counts_to_freqs(results_df, codon_table, CTSAC)
    panel, codon_coords = plot_codon_usage(freqs, codon_table, CTSAC, "Per-sequence codon usage distribution", options.invert)
    print(codon_coords)

    # now figure out the frequencies for the test sequences if they exist
    if options.test_fasta_dir:
        testseq_df = fasta_dir_to_codondf(options.test_fasta_dir, codon_table,
                                          coding=True)
        freqs = counts_to_freqs(testseq_df, codon_table, CTSAC)
        panel = add_points(panel, codon_coords, freqs)

    print_images(
        base_output_name=options.output_basename,
        image_formats=options.fileform,
        no_timestamp = options.no_timestamp,
        dpi=options.dpi,
        transparent=options.transparent)

def run(args):
    codonplot(args)
