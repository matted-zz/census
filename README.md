Census 0.8
==============

Census is a tool to estimate the complexity of sequencing libraries
from read count samples.  See the wiki at
https://github.com/matted/census/wiki for more details, including
guides to interpret the results of Census.

Installation:
==

Several Python packages are required.  On a Ubuntu-like system, these
commands will get the appropriate dependencies:

    sudo apt-get install python python-setuptools python-numpy python-scipy python-pylab
    sudo easy_install pysam

Quick usage:
==

Census operates in two phases, a read duplicate count generation step
and an estimation step.

    ./bam_to_histo.py dummy.bed input.bam | ./calculate_libsize.py -

The default is to use paired-end information to improve the accuracy
of duplicate detection.  Since this won't work for single-end reads,
those experiments must be analyzed with the "-s" option passed to
bam_to_histo.py.

The reads in the input bam must be coordinate-sorted.  The input bed
serves a dual purpose: it gives regions that should be filtered out in
duplicate detection, and only the chromosomes appearing in the bed
file will be used to create duplicates.  This allows for quick
filtering of mitochondrial reads and other sources that do not carry
the same assumptions as the rest of the genome.

Filtering regions for hg19 described by the Pritchard lab (Pickrell et
al., Bioinformatics 2011) are included in the repository (downloaded
from http://eqtl.uchicago.edu/Masking/).

Extended usage and options:
==

Histogram generator usage: 

    usage: bam_to_histo.py [-h] [-v] [-s] [-q MAPQ] [-d MINDIST] [-r REGEXP]
                       excluded_regions.bed sorted_reads.bam

    Histogram generator for Census library complexity package.

    positional arguments:
      excluded_regions.bed
      sorted_reads.bam

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      -s, --single_ended    Include only single-ended reads, instead of only
                            paired-end reads where both ends map.
      -q MAPQ, --mapq MAPQ  Minimum read mapping quality for a read or read pair
                            to be included. Default is 1.
      -d MINDIST, --mindist MINDIST
                            Maximum distance in reported flowcell coordinates for
                            reads to be considered optical duplicates. Default is
                            100.
      -r REGEXP, --regexp REGEXP
                            Regular expression for finding flowcell coordinates
                            from read names. Default is [\w\.
                            ]+:([\d]):([\d]+):([\d]+):([\d]+).*

Library complexity estimation usage:

    usage: calculate_libsize.py [-h] [-v] [-l MINCOUNT] [-r MAXCOUNT]
                            [-s SUBSAMPLE]
                            count_histogram.txt

    Census, library complexity estimation.

    positional arguments:
      count_histogram.txt   File for duplicate count histogram, or - for stdin.

    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
      -l MINCOUNT, --mincount MINCOUNT
                            Minimum duplicate count to use in estimation. Default
                            is 1.
      -r MAXCOUNT, --maxcount MAXCOUNT
                            Maximum duplicate count to use in estimation. Default
                            is 10.
      -s SUBSAMPLE, --subsample SUBSAMPLE
                            Fraction of counts to use (float), useful for testing.
                            Default is 1 (no downsampling).
