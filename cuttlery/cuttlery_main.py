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

"""
CUTTlery: Codon Usage Table Tools and Bayesian determination of ORF protein coding-ness. 
Darrin Schultz (github@conchoecia) and Jordan Eizenga (github@jeizenga)
"""

import sys
import os.path
import argparse

# cuttlery imports
import cuttlery.version

class InterpretColors(argparse.Action):
    """Interpret what type of color argument the user passed and modify to be
    usable by matplotlib

    Author:
    - Darrin T Schultz (github@conchoecia)
    """
    def __call__(self, parser, namespace, values, option_string=None):
        # if it is a string that starts with #, it is hex, pass as string
        if values[0] == "#":
            setattr(namespace, self.dest, values)
        # if it is a tuple or list, convert to a tuple
        elif values[0] in ["(", "["]:
            from ast import literal_eval as make_tuple
            color_tuple = make_tuple(values)
            setattr(namespace, self.dest, color_tuple)
        elif values.isalpha():
            setattr(namespace, self.dest, values)
        else:
            raise Exception("""Color methods only accept color names as hex strings
            '#f4dc41', float RGBA tuple-types represented as strings
            '(0.957, 0.863, 0.016, 1.0)', or as color names as strings that are
            interpretable by matplotlib, ie 'red'.""")


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

def run_subtool(parser, args):
    if args.command == 'calculate-pi':
        import cuttlery.calculate_pi as submodule
    elif args.command == 'codonplot':
        import cuttlery.codonplot as submodule
    elif args.command == 'dirichlet':
        import cuttlery.codon_dirichlet_test as submodule
    elif args.command == 'heterogeneity':
        import cuttlery.heterogeneity as submodule
    elif args.command == 'piNpiSsim':
        import cuttlery.piNpiSsim as submodule
    # run the chosen submodule.
    submodule.run(args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def main():
    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='cuttlery',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=__doc__)
    parser.add_argument("-v", "--version", help="Installed cuttlery version",
                        action="version",
                        version="%(prog)s " + str(cuttlery.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    #########################################
    # create the individual tool parsers
    #########################################


    ###############
    # calculate-pi
    ###############
    parser_calcpi = subparsers.add_parser('calculate-pi',
                                        help="""Calculates the pi value of a
                                        fasta alignment.""")
    parser_calcpi.add_argument("--fasta_aln",
                        action = FullPathsList,
                        nargs = '+',
                        help = """Calculates the pi value of a fasta alignment.
                        Can take an alignment of multiple sequences in one fasta
                        file, or split over many fasta files. Please note that
                        this program does not produce fasta alignments, just
                        analyzes them!""")
    parser_calcpi.set_defaults(func=run_subtool)

    #############
    # codonplot
    #############
    parser_codonplot = subparsers.add_parser('codonplot',
                                        help="""Plots the codon usage of a
                                        collection of genes as a violinplot.""")
    parser_codonplot.add_argument("--coding_fasta_dir",
                        action = FullPaths,
                        help = """The directory that contains the fasta files
                        (alignments or normal sequences) of known coding sequences
                        sequences""")
    parser_codonplot.add_argument('--dpi',
                                metavar='dpi',
                                default=600,
                                type=int,
                                help="""Change the dpi from the default 600
                                if you need it higher""")
    parser_codonplot.add_argument('--fileform',
                               dest='fileform',
                               choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help="""Which output format would you like?
                               Default is png. Select multiple options by putting
                               a space between them: --fileform png pdf jpg""")
    parser_codonplot.add_argument("--invert",
                        default = False,
                        action = 'store_true',
                        help="""Inverts some of the colors to make the plot look
                        better on dark backgrounds while preserving alpha.""")
    parser_codonplot.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_codonplot.add_argument("-o", "--output_basename",
                        default = "codonplot",
                        help = """The results will be written to this file as a
                        csv.""")
    parser_codonplot.add_argument('-T', '--transparent',
                               action='store_false',
                               default = True,
                               help="""Specifies if
                               you want a transparent background. Default
                               is on.""")
    parser_codonplot.add_argument("--test_fasta_dir",
                        action = FullPaths,
                        help = """If there are ORFs but their status as
                        protein-coding is dubious, use this to plot frequencies
                        of individuals reads over the histogram as dots.""")
    parser_codonplot.add_argument("--tt_code",
                        type = int,
                        default = 1,
                        help="""Select which gene code to use. Default is
                        Standard""")
    parser_codonplot.add_argument("--tt_options",
                        action = 'store_true',
                        default = False,
                        help="""Display the optional gene code names that you
                        may pass to the <--tt_code> argument in a subsequent
                        run""")
    parser_codonplot.set_defaults(func=run_subtool)

    ############
    # dirichlet
    ############
    parser_dirichlet = subparsers.add_parser('dirichlet',
                                        help="""Determine if unknown sequences
                                        better fit coding or noncoding sequences.""")
    parser_dirichlet.add_argument("--coding_dir",
                        action = FullPaths,
                        help = """The directory with known coding sequences.""")
    parser_dirichlet.add_argument("--color_nc",
                        default = None,
                        action = InterpretColors,
                        type=str,
                        help = """Select a color for all of the noncoding seqs
                        in the violinplot. Can accept a hex string like:
                        --color_nc '#f4dc41'
                        or a tuple/list type surrounded with quotation marks
                        '(R, G, B, A)' values scaled
                        between 0.0-1.0 like:
                        --color_nc '(0.957, 0.863, 0.016, 1.0)', or a color
                        name that is interpretable by matplotlib, like 'red'. """)
    parser_dirichlet.add_argument("--color_coding",
                        default = None,
                        action = InterpretColors,
                        type=str,
                        help = """Select a color for all of the coding seqs
                        in the violinplot. Can accept a hex string like:
                        --color_nc '#f4dc41'
                        or a tuple/list type surrounded with quotation marks
                        '(R, G, B, A)' values scaled
                        between 0.0-1.0 like:
                        --color_nc '(0.957, 0.863, 0.016, 1.0)', or a color
                        name that is interpretable by matplotlib, like 'red'. """)
    parser_dirichlet.add_argument("--color_test",
                        default = None,
                        action = InterpretColors,
                        type=str,
                        help = """Select a color for all of the test seqs
                        in the violinplot. Can accept a hex string like:
                        --color_nc '#f4dc41'
                        or a tuple/list type surrounded with quotation marks
                        '(R, G, B, A)' values scaled
                        between 0.0-1.0 like:
                        --color_nc '(0.957, 0.863, 0.016, 1.0)', or a color
                        name that is interpretable by matplotlib, like 'red'.""")
    parser_dirichlet.add_argument('--dpi',
                        metavar='dpi',
                        default=600,
                        type=int,
                        help="""Change the dpi from the default 600
                        if you need it higher.""")
    parser_dirichlet.add_argument('--fileform',
                        dest='fileform',
                        choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                 'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                 'svg', 'svgz', 'tif', 'tiff'],
                        default=['png'],
                        nargs='+',
                        help="""Which output format would you like?
                        Default is png. Select multiple options by putting
                        a space between them: --fileform png pdf jpg""")
    parser_dirichlet.add_argument("--noncoding_dir",
                        action = FullPaths,
                        help = """The directory with noncoding sequences to
                        test.""")
    parser_dirichlet.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_dirichlet.add_argument("--numsims",
                        default = 10,
                        type = int,
                        help = """number of simulations per gene.""")
    parser_dirichlet.add_argument('-o', '--output-basename',
                        default = "dirichlet",
                        help="""Specify a base name for the output file('
                        's). The input file base name is the '
                        'dirichlet.' Do not include any trailing periods,
                        characters, or file extensions in the basename.""")
    parser_dirichlet.add_argument("--results_file",
                        action = FullPaths,
                        help = """this saves the data to a csv file""")
    parser_dirichlet.add_argument('-s', '--sort_type',
                        choices=['alpha','meaninc','meandec',
                                 'medinc','meddec'],
                        default='alpha',
                        help="""How to sort the genes in the violinplot.
                        Select one option.
                        - 'alpha': alphabetical by gene name
                        - 'meaninc': increasing mean of log likelihood
                        - 'meandec': decreasing mean of log likelihood
                        - 'medinc' : increasing median of log likelihood
                        - 'meddec' : decreasing median of log likelihood""")
    parser_dirichlet.add_argument("--test_dir",
                        action = FullPaths,
                        help = """The directory with unknown sequences to
                        test.""")
    parser_dirichlet.add_argument("-@", "--threads",
                        type = int,
                        default = 2,
                        help="""The number of threads to use""")
    parser_dirichlet.add_argument('-T', '--transparent',
                        action='store_false',
                        default = True,
                        help="""Specifies if you want a transparent background.
                        Default is on.""")
    parser_dirichlet.set_defaults(func=run_subtool)

    ###############
    # heterogenity
    ###############
    parser_hetero = subparsers.add_parser('heterogeneity',
                                        help="""Looks at site-wise heterogeneity
                                        in nonsynonymous and synonymous mutations
                                        to help infer which regions of a protein
                                        are under the heaviest selection.""")
    parser_hetero.add_argument('--dpi',
                        metavar='dpi',
                        default=600,
                        type=int,
                        help="""Change the dpi from the default 600
                        if you need it higher.""")
    parser_hetero.add_argument("--fasta_dir",
                        action = FullPathsList,
                        nargs = '+',
                        help = """This is the directory where the fasta file
                        alignments are located. Looks recursively so files can
                        be in any subdirectory. The filenames of the fasta files
                        will be used as the figure label while plotting.""")
    parser_hetero.add_argument("--TM_dir",
                        action = FullPaths,
                        help = """This is the directory where the TMHMM
                        files result files are located.""")
    parser_hetero.add_argument('--fileform',
                        dest='fileform',
                        choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                 'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                 'svg', 'svgz', 'tif', 'tiff'],
                        default=['png'],
                        nargs='+',
                        help="""Which output format would you like?
                        Default is png. Select multiple options by putting
                        a space between them: --fileform png pdf jpg""")
    parser_hetero.add_argument("--method",
                        default = 'NG86',
                        choices = ['NG86', 'LWL85', 'YN00', 'ML'],
                        help = """This is the method to use in the simulation to
                        calculate the nonsynonymous and synonymous mutation
                        sites.""")
    parser_hetero.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_hetero.add_argument('-o', '--output_basename',
                        default = "heterogeneity_output",
                        help='Specify a base name for the output file('
                        's). The input file base name is the '
                        'default.')
    parser_hetero.add_argument("--tt_code",
                        type = int,
                        default = 1,
                        help="""Select which gene code to use. Default is
                        Standard.""")
    parser_hetero.add_argument("--tt_options",
                        action = 'store_true',
                        default = False,
                        help="""Display the optional gene code names that you
                        may pass to the <--tt_code> argument in a subsequent
                        run""")
    parser_hetero.add_argument('-T', '--transparent',
                        action='store_false',
                        default = True,
                        help="""Specifies if you want a transparent background.
                        Default is on.""")
    parser_hetero.set_defaults(func=run_subtool)

    #############
    # piNpiSsim
    #############
    parser_piNpiS = subparsers.add_parser('piNpiSsim',
                                        help="""Simulate neutral evolution of
                                        aligned sequences given observed
                                        nucleotide diversity.""")
    parser_piNpiS.add_argument("--debug",
                        type = str,
                        nargs = '+',
                        default = ['NA'],
                        help = """Use this to print out special things while
                        debugging""")
    parser_piNpiS.add_argument('--dpi',
                                metavar='dpi',
                                default=600,
                                type=int,
                                help="""Change the dpi from the default 600
                                if you need it higher""")
    parser_piNpiS.add_argument("--fasta_dir",
                        action = FullPathsList,
                        nargs = '+',
                        help = """This is the directory where the fasta file
                        alignments are located. The filename will be used as the
                        figure label while plotting.""")
    parser_piNpiS.add_argument('--fileform',
                               dest='fileform',
                               choices=['png','pdf', 'eps', 'jpeg', 'jpg',
                                        'pdf', 'pgf', 'ps', 'raw', 'rgba',
                                        'svg', 'svgz', 'tif', 'tiff'],
                               default=['png'],
                               nargs='+',
                               help="""Which output format would you like?
                               Default is png. Select multiple options by putting
                               a space between them: --fileform png pdf jpg""")
    parser_piNpiS.add_argument("--results_file",
                        help = """The results will be written to this file as a
                        csv.""")
    parser_piNpiS.add_argument("--method",
                        default = 'NG86',
                        choices = ['NG86', 'LWL85', 'YN00', 'ML'],
                        help = """the method to use in the simulation""")
    parser_piNpiS.add_argument("--numsims",
                        type = int,
                        help = """the results will be written to this file as a
                        csv.""")
    parser_piNpiS.add_argument("--no_timestamp",
                                  action = 'store_true',
                                  help="""Turn off time stamps in the filename
                                  output.""")
    parser_piNpiS.add_argument('-o', '--output-basename',
                        help="""Specify a base name for the output file('
                        's). The input file base name is the '
                        'dirichlet.' Do not include any trailing periods,
                        characters, or file extensions in the basename.""")
    parser_piNpiS.add_argument('-T', '--transparent',
                               action='store_false',
                               default = True,
                               help="""Specifies if
                               you want a transparent background. Default
                               is on.""")
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
    parser_piNpiS.add_argument("-@","--threads",
                        type = int,
                        help = """The num threads to attempt to use.""")
    parser_piNpiS.set_defaults(func=run_subtool)


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
                 'heterogeneity': parser_hetero.print_help,
                 'calculate-pi' : parser_calcpi.print_help}

    if len(sys.argv)==2:
        #if sys.argv[1] == "dirichlet":
        #    import cuttlery.codon_dirichlet_test as temp_name
        #else:
        #    exec("import cuttlery.{} as temp_name".format(
        #        args.command.replace("-", "_")), globals())
        #print(temp_name.__doc__)
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
