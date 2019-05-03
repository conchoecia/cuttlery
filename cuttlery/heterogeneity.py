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
#
"""title: cuttlery heterogeneity
author: Darrin T Schultz - github@conchoecia
conception: Russell Corbett-Detig

This program plots synonymous and nonsynonymous mutations along the
length of a locus. The density of the mutations along the proteins'
sequences are represented by sticks and a density plot.
"""

#@author: DTS
#info about confusion matrix here: http://www.dataschool.io/simple-guide-to-confusion-matrix-terminology/

#Python/System stuff
import os, sys

# Biopython stuff
from Bio import SeqIO
import Bio.Data.CodonTable
from Bio.codonalign.codonalphabet import get_codon_alphabet
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import _get_pi

# pandas
import pandas as pd
import numpy as np

# cuttlery stuff
from cuttlery.codonfunctions import fasta_dir_to_gene_filelist,\
    fasta_path_to_codonseqs,\
    seqs_to_df, calculate_pi,\
    seqfreqs, calculate_piN_piS,\
    codonseqs_sliced, print_images

# plotting stuff
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib import rc
import seaborn as sns

#set font to helvetica
#global hfont
#hfont = {'fontname':'Helvetica'}

# set seaborn stuff
sns.set(rc={'text.usetex' : True})
sns.set_style("whitegrid", {'font.family': ['sans-serif'],
                            'font.sans-serif': ['Helvetica'],
                            'grid.color': '.95'})

rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',    # micro symbols
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{tgheros}',    # helvetica font
        r'\usepackage{sansmath}',   # math-font matching  helvetica
        r'\sansmath'                # actually tell tex to use it!
        r'\sisetup{detect-all}',    # force siunitx to use the fonts
        ]

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
    # print("This is the file list", filelist)

    #second, select the codon alphabet to use
    codon_alphabet = get_codon_alphabet(Bio.Data.CodonTable.generic_by_id[options.tt_code], gap_char="-")
    codon_table =    Bio.Data.CodonTable.generic_by_id[options.tt_code]

    # now iterate through all genes in the file list and calculate piNpiS for
    #  each codon in the alignment. We use windowsize=1 because we are determining
    #  if every site is nonsynonymous or synonymous mutation. See #notes for
    #  details on what excactly is going on since this bit is confusing.
    windowsize = 1
    # #notes all of the results of this program are saved to the results_file
    #  in tab-delimited format. The processing time is time-consuming, so this
    #  saves time when one just wants to plot.
    results_file = "{}.csv".format(options.output_basename)
    # If we've already done an analysis and the file is there already
    #  don't bother to do it again, but instead just plot
    if not os.path.exists(results_file):
        all_results = []
        for genename in filelist:
            print("looking at this gene", genename)
            codonseqs = fasta_path_to_codonseqs(filelist[genename], codon_table, codon_alphabet)
            #print("these are the codonseqs", codonseqs)
            sliced_codonseqs, num_codons = codonseqs_sliced(codonseqs, windowsize)
            #print("these are the sliced_codons. len={}", sliced_codonseqs, len(sliced_codonseqs))
            #print("there are {} codons".format(num_codons))
            # For each element in the list of CodonSeqs, calculate piN, piS, pi and
            #  save that info along with the codon index of the gene
            for i in range(len(sliced_codonseqs)):
                partial = sliced_codonseqs[i]
                results = calculate_piN_piS(partial, options.method,
                                            codon_table, het = True)
                #if results['piN'] > 0:
                #    print("piN > 0: {}".format(results['piN']))
                results['pi'] = calculate_pi(partial)
                results['seqname'] = genename
                results['type'] = "observed"
                results['slice_start'] = i * windowsize
                results['slice_stop'] = (i + 1) * windowsize
                all_results.append(results)

        results_df = pd.DataFrame.from_dict(all_results)
        # use boolean to mark where piN and piS do not == 0
        #results_df['piNsite'] = results_df.loc[results_df['piN'] != 0, 'slice_start']
        #results_df['piSsite'] = results_df.loc[results_df['piS'] != 0, 'slice_start']
        results_df["piNorS"] = "unassigned"
        for i,row in results_df.iterrows():
            if abs(row['piN']) > abs(row['piS']):
                results_df.loc[i, 'piNorS'] = 'piN'
            elif abs(row['piN']) < abs(row['piS']):
                results_df.loc[i, 'piNorS'] = 'piS'
            else:
                results_df.loc[i, 'piNorS'] = 'NA'
        #results_df.loc[results_df['piN'] != 0, 'piNorS'] = 'piN'
        #results_df.loc[results_df['piS'] != 0, 'piNorS'] = 'piS'
        results_df['piNSsite'] = results_df.query("piN != 0 or piS != 0")['slice_start']
        print(results_df)
        results_df.to_csv(results_file, index = False)
    else:
        print("\nFound {} so skipping analysis\n".format(results_file))
        results_df = pd.read_csv(results_file)

    genename_to_codonlen={seqname:max(list(results_df.loc[results_df["seqname"] == seqname, "slice_stop"])) \
                            for seqname in list(results_df["seqname"].unique())}
    #print("lens of seqs is:")
    #print(genename_to_codonlen)

    final_df = plot_results(results_df,
                 transparent=options.transparent,
                 dpi=options.dpi,
                 fileform = options.fileform,
                 no_timestamp = options.no_timestamp,
                 basename = options.output_basename,
                 TM_dir = options.TM_dir)

    final_df = plot_results_bars(results_df,
                 transparent=options.transparent,
                 dpi=options.dpi,
                 fileform = options.fileform,
                 no_timestamp = options.no_timestamp,
                 basename = "{}_bar".format(options.output_basename),
                 TM_dir = options.TM_dir)

    genename_to_site_type = { seqname:{thistype:len(final_df.loc[(final_df["TMtype"] == thistype) & (final_df["seqname"] == seqname),]) for thistype in ["inside", "TMhelix", "outside"]} for seqname in list(final_df["seqname"].unique())}
    #print("seqs sitecounts is:")
    #print(genename_to_site_type)


    #now we get some simple stats on how many sites we observed in each domain type
    columnnames = ["I_N_count", "I_S_count", "O_N_count", "O_S_count",
               "IO_N_count", "IO_S_count", "TM_N_count", "TM_S_count",
               "I_sites_count", "O_sites_count", "IO_sites_count", "TM_sites_count"]
    empty_array = np.array([np.repeat(-1, len(final_df["seqname"].unique()))]*len(columnnames)).T
    stats_df = pd.DataFrame(empty_array, index = list(final_df["seqname"].unique()), columns=columnnames)
    #print("here is the empty stats array")
    #print(stats_df)
    for genenamekey in list(final_df["seqname"].unique()):
        #inside#N
        size = len(final_df.loc[(final_df["seqname"] == genenamekey) & \
                          (final_df["piNorS"]  == "piN") & \
                          (final_df["TMtype"]  == "inside")])
        stats_df.loc[genenamekey, "I_N_count"] = size

        #inside#S
        size = len(final_df.loc[(final_df["seqname"] == genenamekey) & \
                          (final_df["piNorS"]  == "piS") & \
                          (final_df["TMtype"]  == "inside")])
        stats_df.loc[genenamekey, "I_S_count"] = size

        #outside#N
        size = len(final_df.loc[(final_df["seqname"] == genenamekey) & \
                          (final_df["piNorS"]  == "piN") & \
                          (final_df["TMtype"]  == "outside")])
        stats_df.loc[genenamekey, "O_N_count"] = size

        #outside#S
        size = len(final_df.loc[(final_df["seqname"] == genenamekey) & \
                          (final_df["piNorS"]  == "piS") & \
                          (final_df["TMtype"]  == "outside")])
        stats_df.loc[genenamekey, "O_S_count"] = size

        #TM#N
        size = len(final_df.loc[(final_df["seqname"] == genenamekey) & \
                          (final_df["piNorS"]  == "piN") & \
                          (final_df["TMtype"]  == "TMhelix")])
        stats_df.loc[genenamekey, "TM_N_count"] = size

        #TM#S
        size = len(final_df.loc[(final_df["seqname"] == genenamekey) & \
                          (final_df["piNorS"]  == "piS") & \
                          (final_df["TMtype"]  == "TMhelix")])
        stats_df.loc[genenamekey, "TM_S_count"] = size

        #IO#N
        size = int(stats_df.loc[genenamekey, "I_N_count"]) + int(stats_df.loc[genenamekey, "O_N_count"])
        stats_df.loc[genenamekey, "IO_N_count"] = size

        #IO#S
        size = int(stats_df.loc[genenamekey, "I_S_count"]) + int(stats_df.loc[genenamekey, "O_S_count"])
        stats_df.loc[genenamekey, "IO_S_count"] = size

        #I_sites_count
        size = genename_to_site_type[genenamekey]["inside"]
        stats_df.loc[genenamekey, "I_sites_count"] = size

        #O_sites_count
        size = genename_to_site_type[genenamekey]["outside"]
        stats_df.loc[genenamekey, "O_sites_count"] = size

        #IO_sites_count
        size = genename_to_site_type[genenamekey]["outside"] + genename_to_site_type[genenamekey]["inside"]
        stats_df.loc[genenamekey, "IO_sites_count"] = size

        #TM_sites_count
        size = genename_to_site_type[genenamekey]["TMhelix"]
        stats_df.loc[genenamekey, "TM_sites_count"] = size

    TM_NS_count = stats_df["TM_N_count"] + stats_df["TM_S_count"]
    IO_NS_count = stats_df["IO_N_count"] + stats_df["IO_S_count"]
    all_NS_count = TM_NS_count + IO_NS_count
    I_NS_count =  stats_df["I_N_count"]  + stats_df["I_S_count"]
    O_NS_count =  stats_df["O_N_count"]  + stats_df["O_S_count"]

    stats_df[ "per_of_all_O_sites_that_are_N"]  =  stats_df["O_N_count"]/stats_df["O_sites_count"]
    stats_df[  "per_of_all_O_sites_that_are_S"]  =  stats_df["O_S_count"]/stats_df["O_sites_count"]

    stats_df[ "per_of_all_I_sites_that_are_N"]  =  stats_df["I_N_count"]/stats_df["I_sites_count"]
    stats_df[  "per_of_all_I_sites_that_are_S"]  =  stats_df["I_S_count"]/stats_df["I_sites_count"]

    stats_df[ "per_of_all_IO_sites_that_are_N"]  =  stats_df["IO_N_count"]/stats_df["IO_sites_count"]
    stats_df[ "per_of_all_IO_sites_that_are_S"]  =  stats_df["IO_S_count"]/stats_df["IO_sites_count"]

    stats_df[ "per_of_all_TM_sites_that_are_N"]  =  stats_df["TM_N_count"]/stats_df["TM_sites_count"]
    stats_df[ "per_of_all_TM_sites_that_are_S"]  =  stats_df["TM_S_count"]/stats_df["TM_sites_count"]

    stats_df[ "total_sites"] = stats_df["IO_sites_count"] + stats_df["TM_sites_count"]

    stats_df[ "O_Ns_comprise_what_percent_of_locus"]  =  stats_df["O_N_count"]/stats_df["total_sites"]
    stats_df[ "O_Ss_comprise_what_percent_of_locus"]  =  stats_df["O_S_count"]/stats_df["total_sites"]
    stats_df[ "I_Ns_comprise_what_percent_of_locus"]  =  stats_df["I_N_count"]/stats_df["total_sites"]
    stats_df[ "I_Ss_comprise_what_percent_of_locus"]  =  stats_df["I_S_count"]/stats_df["total_sites"]
    stats_df[ "IO_Ns_comprise_what_percent_of_locus"]  =  stats_df["IO_N_count"]/stats_df["total_sites"]
    stats_df[ "IO_Ss_comprise_what_percent_of_locus"]  =  stats_df["IO_S_count"]/stats_df["total_sites"]
    stats_df[ "TM_Ns_comprise_what_percent_of_locus"]  =  stats_df["IO_N_count"]/stats_df["total_sites"]
    stats_df[ "TM_Ss_comprise_what_percent_of_locus"]  =  stats_df["IO_S_count"]/stats_df["total_sites"]

    stats_df.round(4)
    print("here is the stats array after filling")
    print(stats_df)
    stats_file = "{}.stats.csv".format(options.output_basename)
    stats_df.to_csv(stats_file, index = True)

def plot_results_bars(df, **kwargs):
    print("inside plot_results_bar")
    df["TMtype"] = "NA"
    print(df)
    width = 0.8
    delta = 0.05
    final_width = width - delta
    seqnames = sorted(df['seqname'].unique())
    # now that we have plotted the piN and piS distributions, plot a solid line
    #  indicating the whole length of the gene
    gene_lens = {seqname:max(df.loc[df['seqname'] == seqname, 'slice_stop'])
                 for seqname in seqnames}
    global_max_len = max(gene_lens.values())


    figWidth = 9
    figHeight = 3.5
    plt.figure(figsize=(figWidth,figHeight))
    #set the panel dimensions
    panelWidth = 8
    panelHeight = 2.5
    #find the margins to center the panel in figure
    leftMargin = (figWidth - panelWidth)/2
    bottomMargin = (figHeight - panelHeight)/2
    panel0 =plt.axes([leftMargin/figWidth, #left
                     bottomMargin/figHeight,    #bottom
                     panelWidth/figWidth,   #width
                     panelHeight/figHeight])     #height
    panel0.tick_params(axis='both',which='both',\
                       bottom='on', labelbottom='on',\
                       left='on', labelleft='on', \
                       right='off', labelright='off',\
                       top='off', labeltop='off')
    panel0.spines['top'].set_visible(False)
    panel0.spines['right'].set_visible(False)
    panel0.spines['left'].set_visible(False)


    piNcolor = "#74b6a1"
    piScolor = "#e79578"
    rectangle_patches = []
    add_first_patches = []
    bar_height = 0.4
    highlight_width = 6
    grad_num = 50
    grad_bar_width = highlight_width/grad_num
    grad_bar_pos = np.linspace(0, highlight_width, grad_num)
    alpha_index = np.linspace(0.75,0,grad_num)
    for i,row in df.iterrows():
        #plot the piN bars
        if row["piNorS"] != "NA":
            for moltype in ["piN", "piS"]:
                if row["piNorS"] == moltype:
                    thiscolor = "black"
                    bottom = -10
                    if moltype == "piN":
                        thiscolor= piNcolor
                        bottom = seqnames.index(row["seqname"])-bar_height
                    else:
                       thiscolor = piScolor
                       bottom = seqnames.index(row["seqname"])
                    left = row["slice_start"]
                    rectangle1=mplpatches.Rectangle((left,bottom),1,0.4,\
                                                    linewidth=0.0,\
                                                    facecolor='black',\
                                                    alpha=1,
                                                    edgecolor=(1,1,1))
                    rectangle_patches.append(rectangle1)
                    for j in range(0,grad_num):
                        thisalpha = alpha_index[j]
                        this_offset=grad_bar_pos[j]
                        thisleft = left - this_offset
                        if thisleft >= 0:
                            rectangle2=mplpatches.Rectangle((thisleft,bottom),grad_bar_width,0.4,\
                                                            linewidth=0.0,\
                                                            facecolor=thiscolor,\
                                                            alpha = thisalpha,
                                                            edgecolor=(1,1,1))
                            add_first_patches.append(rectangle2)
                        thisleft = left + 1 + this_offset
                        if thisleft < gene_lens[row["seqname"]]:
                            rectangle2=mplpatches.Rectangle((thisleft,bottom),grad_bar_width,0.4,\
                                                            linewidth=0.0,\
                                                            facecolor=thiscolor,\
                                                            alpha = thisalpha,
                                                            edgecolor=(1,1,1))
                            add_first_patches.append(rectangle2)

    #This just plots all the gene lengths as rectangles
    print("gene lens")
    print(gene_lens)
    left = 0
    bottom = 0
    height = 0
    width = 0
    for genenamekey in gene_lens:
        if kwargs["TM_dir"]:
            filepath = os.path.join(kwargs["TM_dir"],"{}_TM.txt".format(genenamekey))
            #This part plots the transmembrane domains
            if os.path.exists(filepath):
                left = 0
                height = 0.2
                bottom = seqnames.index(genenamekey)-(height/2)
                width = gene_lens[genenamekey]
                # white rectangle
                rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                            linewidth=0.0,\
                                            facecolor='white',\
                                            alpha=1,
                                            edgecolor=(1,1,1))
                rectangle_patches.append(rectangle1)

                TM_df = pd.read_csv(filepath, delimiter=',',
                                 comment='#', header=None,
                                 names = ["gene", "tool", "type",
                                          "start", "stop"])
                #print(df)
                for index, row in TM_df.iterrows():
                    thistype = row["type"]
                    height = 0.2
                    left = row["start"]
                    width = row["stop"] - row["start"] + 1
                    if thistype == "TMhelix":
                        bottom = seqnames.index(genenamekey)-(height/2)
                        color = "black"
                    if thistype == "outside":
                        height = height/4
                        bottom = seqnames.index(genenamekey)-(height*1.5)
                        color = "#ff00ff"
                    if thistype == "inside":
                        height = height/4
                        bottom = seqnames.index(genenamekey)+(height*0.6666)
                        color = "#3f3fff"
                    rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                                        linewidth=0.0,\
                                                        facecolor=color,\
                                                        alpha=1,
                                                        edgecolor=(1,1,1))
                    rectangle_patches.append(rectangle1)
                    #now add the info to the df to calculate stats
                    for i in range(row["start"], row["stop"]+1):
                        df.loc[(df["seqname"] == genenamekey) & (df["slice_start"] == i), "TMtype"] = row["type"] 
            else:
                message = "{} does not exist".format(filepath)
                raise(message)
        else:
            left = 0
            height = 0.1
            bottom = seqnames.index(genenamekey)-(height/2)
            width = gene_lens[genenamekey]
            rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                        linewidth=0.0,\
                                        facecolor='black',\
                                        alpha=1,
                                        edgecolor=(1,1,1))
            rectangle_patches.append(rectangle1)

    for patch in add_first_patches:
        patch.set_zorder(20)
        panel0.add_patch(patch)
    for patch in rectangle_patches:
        patch.set_zorder(20)
        panel0.add_patch(patch)

    panel0.set_xlim([0, global_max_len * 1.05])
    panel0.set_ylim([len(seqnames)-0.5, -0.5])
    panel0.set_yticks(list(range(len(seqnames))))
    panel0.set_yticklabels(seqnames)

    panel0.set_ylabel('')
    panel0.set_xlabel('Codon Number')
    panel0.set_title('Density of nucleotide diversity')


    # Print image(s)
    print_images(
        base_output_name=kwargs["basename"],
        image_formats=kwargs["fileform"],
        no_timestamp=kwargs["no_timestamp"],
        dpi=kwargs["dpi"],
        transparent=kwargs["transparent"])

    return df

def plot_results(df, **kwargs):
    print("inside plot_results")
    df["TMtype"] = "NA"
    print(df)
    inner = "sticks"
    width = 0.8
    delta = 0.05
    final_width = width - delta
    seqnames = sorted(df['seqname'].unique())
    fig = plt.figure(figsize=(8,4))
    ax = sns.violinplot(x="piNSsite", y="seqname",
                  hue="piNorS", hue_order = ["piN", "piS"],
                  data=df, palette="Set2", split=True,
                  scale="width", inner=inner,
                  scale_hue=False, bw=.1,
                  cut = 0, dodge = True,
                  width = final_width,
                  order = seqnames)

    ax.legend(loc='best')
    # took this bit of code to change the legend labels from here:
    #  https://stackoverflow.com/questions/45201514/
    new_labels = [r'$\mathrm{\pi}$ N', r'$\mathrm{\pi}$ S']
    for t, l in zip(ax.legend_.texts, new_labels): t.set_text(l)
    # this moves the violin plots away from each other a little bit.
    offset_violinplot_halves(ax, delta, final_width, inner, 'horizontal')

    # now that we have plotted the piN and piS distributions, plot a solid line
    #  indicating the whole length of the gene
    gene_lens = {seqname:max(df.loc[df['seqname'] == seqname, 'slice_stop'])
                 for seqname in seqnames}
    global_max_len = max(gene_lens.values())
    rectangle_patches = []
    #This just plots all the gene lengths as rectangles
    print("gene lens")
    print(gene_lens)
    left = 0
    bottom = 0
    height = 0
    width = 0
    for genenamekey in gene_lens:
        if kwargs["TM_dir"]:
            filepath = os.path.join(kwargs["TM_dir"],"{}_TM.txt".format(genenamekey))
            #This part plots the transmembrane domains
            if os.path.exists(filepath):
                left = 0
                height = 0.2
                bottom = seqnames.index(genenamekey)-(height/2)
                width = gene_lens[genenamekey]
                rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                            linewidth=0.0,\
                                            facecolor='white',\
                                            alpha=1,
                                            edgecolor=(1,1,1))
                rectangle_patches.append(rectangle1)

                TM_df = pd.read_csv(filepath, delimiter=',',
                                 comment='#', header=None,
                                 names = ["gene", "tool", "type",
                                          "start", "stop"])
                #print(df)
                for index, row in TM_df.iterrows():
                    thistype = row["type"]
                    height = 0.2
                    left = row["start"]
                    width = row["stop"] - row["start"] + 1
                    if thistype == "TMhelix":
                        bottom = seqnames.index(genenamekey)-(height/2)
                        color = "black"
                    if thistype == "outside":
                        height = height/4
                        bottom = seqnames.index(genenamekey)-(height*1.5)
                        color = "#ff00ff"
                    if thistype == "inside":
                        height = height/4
                        bottom = seqnames.index(genenamekey)+(height*0.6666)
                        color = "#3f3fff"
                    rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                                        linewidth=0.0,\
                                                        facecolor=color,\
                                                        alpha=1,
                                                        edgecolor=(1,1,1))
                    rectangle_patches.append(rectangle1)
                    #now add the info to the df to calculate stats
                    for i in range(row["start"], row["stop"]+1):
                        df.loc[(df["seqname"] == genenamekey) & (df["slice_start"] == i), "TMtype"] = row["type"] 
            else:
                message = "{} does not exist".format(filepath)
                raise(message)
        else:
            left = 0
            height = 0.1
            bottom = seqnames.index(genenamekey)-(height/2)
            width = gene_lens[genenamekey]
            rectangle1=mplpatches.Rectangle((left,bottom),width,height,\
                                        linewidth=0.0,\
                                        facecolor='black',\
                                        alpha=1,
                                        edgecolor=(1,1,1))
            rectangle_patches.append(rectangle1)

    for patch in rectangle_patches:
        patch.set_zorder(20)
        ax.add_patch(patch)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlim([0, global_max_len * 1.1])
    ax.set_ylabel('')
    ax.set_xlabel('Codon Number')
    ax.set_title('Density of nucleotide diversity')


    # Print image(s)
    print_images(
        base_output_name=kwargs["basename"],
        image_formats=kwargs["fileform"],
        no_timestamp=kwargs["no_timestamp"],
        dpi=kwargs["dpi"],
        transparent=kwargs["transparent"])

    return df

def offset_violinplot_halves(ax, delta, width, inner, direction):
    """
    This function offsets the halves of a violinplot to compare tails
    or to plot something else in between them. This is specifically designed
    for violinplots by Seaborn that use the option `split=True`.

    For lines, this works on the assumption that Seaborn plots everything with
     integers as the center.

    Args:
     <ax>    The axis that contains the violinplots.
     <delta> The amount of space to put between the two halves of the violinplot
     <width> The total width of the violinplot, as passed to sns.violinplot()
     <inner> The type of inner in the seaborn
     <direction> Orientation of violinplot. 'hotizontal' or 'vertical'.

    Returns:
     - NA, modifies the <ax> directly

    I submitted this code here: https://stackoverflow.com/questions/43357274/separate-halves-of-split-violinplot-to-compare-tail-data/45902085#45902085
    """
    # offset stuff
    if inner == 'sticks':
        lines = ax.get_lines()
        for line in lines:
            if direction == 'horizontal':
                data = line.get_ydata()
                if int(data[0] + 1)/int(data[1] + 1) < 1:
                    # type is top, move neg, direction backwards for horizontal
                    data -= delta
                else:
                    # type is bottom, move pos, direction backward for hori
                    data += delta
                line.set_ydata(data)
            elif direction == 'vertical':
                data = line.get_xdata()
                if int(data[0] + 1)/int(data[1] + 1) < 1:
                    # type is left, move neg
                    data -= delta
                else:
                    # type is left, move pos
                    data += delta
                line.set_xdata(data)

    for ii, item in enumerate(ax.collections):
        # axis contains PolyCollections and PathCollections
        if isinstance(item, matplotlib.collections.PolyCollection):
            # get path
            path, = item.get_paths()
            vertices = path.vertices
            half_type = _wedge_dir(vertices, direction)
            # shift x-coordinates of path
            if half_type in ['top','bottom']:
               if inner in ["sticks", None]:
                    if half_type == 'top': # -> up
                        vertices[:,1] -= delta
                    elif half_type == 'bottom': # -> down
                        vertices[:,1] += delta
            elif half_type in ['left', 'right']:
                if inner in ["sticks", None]:
                    if half_type == 'left': # -> left
                        vertices[:,0] -= delta
                    elif half_type == 'right': # -> down
                        vertices[:,0] += delta

def _wedge_dir(vertices, direction):
    """
    This has some function with turning the split violinplot horizontally or
    vertically.

    Args:
      <vertices>  The vertices from matplotlib.collections.PolyCollection
      <direction> Direction must be 'horizontal' or 'vertical' according to how
                   your plot is laid out.
    Returns:
      - a string in ['top', 'bottom', 'left', 'right'] that determines where the
         half of the violinplot is relative to the center.
    """
    if direction == 'horizontal':
        result = (direction, len(set(vertices[1:5,1])) == 1)
    elif direction == 'vertical':
        result = (direction, len(set(vertices[-3:-1,0])) == 1)
    outcome_key = {('horizontal', True): 'bottom',
                   ('horizontal', False): 'top',
                   ('vertical', True): 'left',
                   ('vertical', False): 'right'}
    # if the first couple x/y values after the start are the same, it
    #  is the input direction. If not, it is the opposite
    return outcome_key[result]
