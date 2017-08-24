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

# I took this code from https://github.com/arq5x/poretools/. Check it out. - DTS

import sys
import os.path
import argparse

#logger
import logging
logger = logging.getLogger('poretools')

# pauvre imports
import pauvre.version

#This class is used in argparse to expand the ~. This avoids errors caused on
# some systems.
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                os.path.abspath(os.path.expanduser(values)))

class FullPathsList(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
                [os.path.abspath(os.path.expanduser(value)) for value in values])

    commandDict={'dirichlet'    : parser_dirichlet.print_help,
                 'codonplot'    : parser_codonplot.print_help,
                 'piNpiSsim'    : parser_piNpiS.print_help,
                 'heterogeneity': parser_hetero.print_help}

def run_subtool(parser, args):
    if args.command == 'dirichlet':
        import cuttlery.codon_dirichlet_test as submodule
    elif args.command == 'codonplot':
        import cuttlery.codonplot as submodule
    elif args.command == 'piNpiSsim':
        import cuttlery.piNpiSsim as submodule
    elif args.command == 'heterogeneity':
        import cuttlery.heterogeneity as submodule
    # run the chosen submodule.
    submodule.run(args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")

def main():

    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='cuttlery', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed cuttlery version",
                        action="version",
                        version="%(prog)s " + str(cuttlery.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    #########################################
    # create the individual tool parsers
    #########################################

    #############
    # codonplot
    #############
    parser_codonplot = subparsers.add_parser('codonplot',
                                        help="""calculates codon usage from
                                        protein alignments, then prints and plots
                                        the results""")
    parser_codonplot.add_argument('-f', '--fastq',
                               metavar='FASTQ',
                               action=FullPaths,
                               help='The input FASTQ file.')
    parser_codonplot.add_argument('-n', '--no-transparent',
                               dest='TRANSPARENT',
                               action='store_false',
                               help="""Not the TV show. Specify this option if
                               you don't want a transparent background. Default
                               is on.""")
    parser_codonplot.add_argument('-t', '--title',
                               metavar='TITLE',
                               default='Read length vs mean quality',
                               help="""This sets the title for the whole plot.
                               Use --title "Crustacean's DNA read quality"
                               if you need single quote or apostrophe
                               inside title.""")
    parser_codonplot.add_argument('--maxlen',
                               metavar='MAXLEN',
                               type=int,
                               help="""This sets the max read length to plot.""")
    parser_codonplot.add_argument('-m', '--maxqual',
                               metavar='MAXQUAL',
                               type=int,
                               help="""This sets the max mean read quality
                               to plot.""")
    parser_codonplot.add_argument('--lengthbin',
                               metavar='LENGTHBIN',
                               type=int,
                               help="""This sets the bin size to use for length.""")
    parser_codonplot.add_argument('--qualbin',
                               metavar='QUALBIN',
                               type=float,
                               help="""This sets the bin size to use for quality""")
    parser_codonplot.add_argument('-y', '--add-yaxes',
                               dest='Y_AXES',
                               action='store_true',
                               help='Add Y-axes to both marginal histograms.')
    parser_codonplot.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
    parser_codonplot.add_argument('-o', '--output-base-name',
                               dest='BASENAME',
                               help='Specify a base name for the output file('
                                    's). The input file base name is the '
                                    'default.')
    parser_codonplot.add_argument('-d', '--dpi',
                               metavar='dpi',
                               default=600,
                               type=int,
                               help="""Change the dpi from the default 600
                               if you need it higher""")
    parser_codonplot.set_defaults(func=run_subtool)

    ############
    # dirichlet
    ############
    parser_dirichlet = subparsers.add_parser('dirichlet',
                                        help="""Determine if unknown sequences
                                        better fit coding or noncoding sequences""")
    parser_dirichlet.add_argument("--tt_code",
                        type = int,
                        default = 4,
                        help="""Select which gene code to use. Default is
                        Standard""")
    parser_dirichlet.add_argument("--tt_options",
                        action = 'store_true',
                        default = False,
                        help="""Display the optional gene code names that you
                        may pass to the <--tt_code> argument in a subsequent
                        run""")
    parser_dirichlet.add_argument("--plot",
                        action = 'store_true',
                        default = False,
                        help="""Plot the data file""")
    parser_dirichlet.add_argument("--coding_dir",
                        action = FullPaths,
                        required = True,
                        help = """The directory with known coding sequences.""")
    parser_dirichlet.add_argument("--test_dir",
                        action = FullPaths,
                        required = True,
                        help = """The directory with unknown sequences to
                        test.""")
    parser_dirichlet.add_argument("--noncoding_dir",
                        action = FullPaths,
                        required = True,
                        help = """The directory with noncoding sequences to
                        test.""")
    parser_dirichlet.add_argument("--numsims",
                        default = 10,
                        type = int,
                        help = """number of simulations per gene.""")
    parser_dirichlet.add_argument("--results_file",
                        action = FullPaths,
                        required = True,
                        help = """this saves the data to a csv file""")
    parser_dirichlet.add_argument("--threads",
                        type = int,
                        default = 2,
                        help="""The number of threads to use""")
    parser_dirichlet.set_defaults(func=run_subtool)

    #############
    # piNpiSsim
    #############
    parser_piNpiS = subparsers.add_parser('piNpiSsim',
                                        help="""Simulate neutral evolution of
                                        aligned sequences given observed
                                        alignments""")
    parser_piNpiS.add_argument("--outprefix",
                        type=str,
                        default=sys.stdout,
                        help="""the output will be saved here""")
    parser_piNpiS.add_argument("--tt_code",
                        type = int,
                        default = 4,
                        help="""Select which gene code to use. Default is
                        Standard""")
    parser_piNpiS.add_argument("--tt_options",
                        action = 'store_true',
                        default = False,
                        help="""Display the optional gene code names that you
                        may pass to the <--tt_code> argument in a subsequent
                        run""")
    parser_piNpiS.add_argument("--debug",
                        type = str,
                        nargs = '+',
                        default = ['NA'],
                        help = """Use this to print out special things while
                        debugging""")
    parser_piNpiS.add_argument("--fasta_dir",
                        action = FullPaths,
                        required = True,
                        help = """This is the directory where the fasta file
                        alignments are located. The filename will be used as the
                        figure label while plotting.""")
    parser_piNpiS.add_argument("--results_file",
                        required = True,
                        help = """The results will be written to this file as a
                        csv.""")
    parser_piNpiS.add_argument("--method",
                        default = 'NG86',
                        choices = ['NG86', 'LWL85', 'YN00', 'ML'],
                        help = """the method to use in the simulation""")
    parser_piNpiS.add_argument("--numsims",
                        required = True,
                        type = int,
                        help = """the results will be written to this file as a
                        csv.""")
    parser_piNpiS.add_argument("--threads",
                        required = False,
                        type = int,
                        help = """The num threads to attempt to use.""")
    parser_piNpiS.set_defaults(func=run_subtool)

    ###############
    # heterogenity
    ###############
    parser_hetero = subparsers.add_parser('heterogeneity',
                                        help="""Looks at site-wise heterogeneity
                                        in nonsynonymous and synonymous mutations
                                        to help infer which regions of a protein
                                        are under the heaviest selection""")
    parser_hetero.add_argument('--gff_paths',
                                metavar='gff_paths',
                                action=FullPathsList,
                                nargs = '+',
                                help="""The input filepath for the gff annotation
                                to plot""")
    parser_hetero.add_argument('--gff_labels',
                                metavar='gff_labels',
                                type = str,
                                nargs = '+',
                                help="""In case the gff names and sequence names
                                don't match, change the labels that will appear
                                over the text.""")
    parser_hetero.add_argument('--dpi',
                                metavar='dpi',
                                default=600,
                                type=int,
                                help="""Change the dpi from the default 600
                                if you need it higher""")
    parser_hetero.add_argument('--optimum_order',
                                action = 'store_true',
                                help="""If selected, this doesn't plot the
                                optimum arrangement of things as they are input
                                into gff_paths. Instead, it uses the first gff
                                file as the top-most sequence in the plot, and
                                reorganizes the remaining gff files to minimize
                                the number of intersections.""")
    parser_hetero.add_argument('--aln_dir',
                                metavar='aln_dir',
                                action=FullPaths,
                                help="""The directory where all the fasta
                                alignments are contained.""")
    parser_hetero.add_argument('--stop_codons',
                                action='store_true',
                                default = True,
                                help="""Performs some internal corrections if
                                the gff annotation includes the stop
                                codons in the coding sequences.""")
    parser_hetero.add_argument('--center_on',
                                type = str,
                                default = None,
                                help="""centers the plot around the gene that
                                you pass as an argument""")
    parser_hetero.add_argument('--start_with_aligned_genes',
                                action='store_true',
                                default = False,
                                help="""Minimizes the number of intersections
                                but only selects combos where the first gene in
                                each sequence is aligned.""")
    parser_hetero.add_argument('--fileform',
                               dest='fileform',
                               metavar='STRING',
                               choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help='Which output format would you like? Def.=png')
    parser_hetero.add_argument('-o', '--output-base-name',
                               dest='BASENAME',
                               help='Specify a base name for the output file('
                                's). The input file base name is the '
                                'default.')
    parser_hetero.add_argument('-n', '--no-transparent',
                               dest='TRANSPARENT',
                               action='store_false',
                               help="""Not the TV show. Specify this option if
                               you don't want a transparent background. Default
                               is on.""")
    parser_hetero.add_argument('--sandwich',
                                action='store_true',
                                default = False,
                                help="""Put an additional copy of the first gff
                                file on the bottom of the plot for comparison.""")

    parser_hetero.set_defaults(func=run_subtool)


    #######################################################
    # parse the args and call the selected function
    #######################################################

    args = parser.parse_args()

    # If there were no args, print the help function
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # If there were no args, but someone selected a program,
    #  print the program's help.
    commandDict={'dirichlet'    : parser_dirichlet.print_help,
                 'codonplot'    : parser_codonplot.print_help,
                 'piNpiSsim'    : parser_piNpiS.print_help,
                 'heterogeneity': parser_hetero.print_help}

    if len(sys.argv)==2:
        commandDict[args.command]()

        sys.exit(1)

    if args.quiet:
        logger.setLevel(logging.ERROR)

    try:
        args.func(parser, args)
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

if __name__ == "__main__":
    main()
