# cutlery
Codon usage table tool for python

- Currently, this repo only contains the `codon_dirichlet_test.py` file and some test files
- There are two unknown fasta files currently
 - `NAD2L.fasta`
 - `unk1_2013Bf.fasta`
 - I also added one known gene as a test (COIII)
   - `Bf201606_COIII.fasta`
   - The fasta file of known genes without COIII
     for the test is: `Bf_all_noCOIII.fasta`
- There is one set of known genes
 - `knownCDS.fasta`
- There are two sets of null hypothesis fasta files
 - `Bf_nc_rev.fasta` is the noncoding regions on the antisense strand (probably not transcribed)
 - `Bf_nc_fwd.fasta` is the noncoding regions on the sense strand (transcribed)
