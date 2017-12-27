# `cuttlery`
Codon Usage Table Tools for python

- `cuttlery` has several tools

## Modules

### `cuttlery calculate-pi`

This program calculates the value of little pi, also known as
nucleotide diversity, of the sequences in a fasta alignment.

The most sensible way to use this program is to make one fasta
alignment with your favorite fasta alignment program and to save it as
a fasta alignment. Alternatively, you can make a fasta alignment and
save each sequence to a separate fasta file.

Then, calculate little pi using cuttlery with

    cuttlery calculate-pi --fasta_aln <fasta_aln1.fasta> <fasta_aln2.fasta> <et cetera>

The output will be printed to std out.

### `cuttlery codonplot`

This program plots the codon usage of a collection of genes as violin
plots.  This style of plot more accurately reflects the biological
variability in codon usage in different protein types. In addition,
one can plot single or multiple genes as dots to highlight their
position in the distribution of codon usage frequencies.


### `cuttlery dirichlet`

In short, this test asks if the trinucleotide frequency of a test ORF
more closely matches the trinucleotides frequencies of coding or
noncoding DNA from the same species.

### `cuttlery heterogeneity`

This program plots synonymous and nonsynonymous mutations along the
length of a locus. The density of the mutations along the proteins'
sequences are represented by sticks and a density plot.

### `cuttlery piNpiSsim`

This program observes the nucleotide diversity of a protein-coding
locus in a population of sequences. Using this information it
generates other similar sequences with the same nucleotide diversity
and geneaology as the observed sequences. However, the generated
sequences have their mutation sites modeled from a randomly-chosen
observed sequence at randomly chosen sites. This simulation may
estimate the lower bounds of piN/piS for a neutrally evolving sequence
of similar base composition and geneology as the observed population
of protein-coding sequences.

## Authors

Darrin T Schultz ([github@conchoecia](https://github.com/conchoecia))
Jordan M Eizenga ([github@jeizenga](https://github.com/jeizenga))

Special thanks to
[Russell Corbett-Detig](https://corbett-detig-lab.soe.ucsc.edu/) for
conceiving of the tests programmed in `cuttlery piNpiSsim` and
`cuttlery heterogeneity`. ([github@russcd](https://github.com/russcd))
