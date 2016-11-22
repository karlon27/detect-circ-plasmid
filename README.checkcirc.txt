This is a perl script for detecting circular plasmids contigs/scaffolds given paired-end reads.

It requires two additional software in order to work properly:

ncbi-blast+
bowtie2


The usage of this script is:
Usage: perl check_circular.pl (input fasta file) (reads pair1) (reads pair2) (output fasta with circular seqs) (output fasta without circular seqs)

Briefly, the script takes in five parameters, the first three indicates the contig/scaffold and two paired-end reads files, and the latter two indicates the output files.

Please contact me Yu-Wei Wu at ywwei@lbl.gov if you encounter problems running this software. Thank you.


License
============
Copyright (C) 2016 Yu-Wei Wu
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
