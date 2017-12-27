# cuttlery
Codon usage table tool for python

- `cuttlery` has several tools
  - `cuttlery calculate-pi`
  - `cuttlery codonplot`
    - This program plots the codon usage of a collection of genes as
      violin plots.  This style of plot more accurately reflects the
      biological variability in codon usage in different protein
      types. In addition, one can plot single or multiple genes as
      dots to highlight their position in the distribution of codon
      usage frequencies.
  - `cuttlery dirichlet`
  - `cuttlery heterogeneity`
    - This program plots synonymous and nonsynonymous mutations along
      the length of a locus. The density of the mutations along the
      proteins' sequences are represented by sticks and a density
      plot.
  - `cuttlery piNpiSsim`
    - This program observes the nucleotide diversity of a
      protein-coding locus in a population of sequences. Using this
      information it generates other similar sequences with the same
      nucleotide diversity and geneaology as the observed
      sequences. However, the generated sequences have their mutation
      sites modeled from a randomly-chosen observed sequence at
      randomly chosen sites. This simulation may estimate the lower
      bounds of piN/piS for a neutrally evolving sequence of similar
      base composition and geneology as the observed population of
      protein-coding sequences.

