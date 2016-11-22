This is a perl script for detecting circular plasmids contigs/scaffolds given paired-end reads.

It requires two additional software in order to work properly:

ncbi-blast+
bowtie2


The usage of this script is:
Usage: perl check_circular.pl (input fasta file) (reads pair1) (reads pair2) (output fasta with circular seqs) (output fasta without circular seqs)

Briefly, the script takes in five parameters, the first three indicates the contig/scaffold and two paired-end reads files, and the latter two indicates the output files.

Please contact me Yu-Wei Wu at ywwei@lbl.gov if you encounter problems running this software. Thank you.

