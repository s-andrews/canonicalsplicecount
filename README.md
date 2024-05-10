Splice Site Quantitation
========================

This program implements a method for quantitating RNA-Seq data to produce a raw count matrix suitable for feeding into downstream analyses with programs such as DESeq2 or EdgeR.

What does the program do?
-------------------------

Rather than count all of the reads in your data this quantitation is based off the observed splice sites seen when matching the data to the reference genome.  When you quantitate your BAM files you provide a GTF file with the details of your gene models and the program will count
only reads which cross a known and annotated splice junction in the gene.


Why do it this way?
-------------------

In some circumstances quantitating your data this way can give cleaner results than quantitation methods which use all of your reads.  Cases where this could be useful would be:

1. If you have DNA contamination in your samples
2. If you have variable amounts of unprocessed RNA mixed in with mature RNA
3. If you have chimeric sequences which are producing artefactual trans-splice constructs


Why not always do it this way?
------------------------------

This type of quantitation is not always going to be appropriate.  You will get lower overall counts than using all of your reads, since only reads which cross at least one splice boundary will be counted.  You will also not get coverage of any unspliced genes in your genome.  By default we only
quantitate protein coding genes, since spliced non-coding genes will generally be removed by nonsense-mediated-decay, so you will only quantitate this subset of the data.


Installation
============

You can download the latest version either from the Releases section of this site and un-tarring it, or by cloning the whole repository.

You will need a python3 (>=3.7) installation, with the pysam package installed.  If you want to create a virtual environment to run this program then you can do:

```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

Running the program
===================

To run the program you will need one or more BAM files from a spliced alignment program (Hisat2 or STAR for example) and a GTF file of gene annotations.  Since the program requires an exact match between your observed splice and the position of the intron we recommend that you provide the same GTF file to your aligner when running the alignment so it knows exactly where to position the intron in cases where the position is ambiguous.


This would be a typical run of the program (split over multiple lines here for reasability but you can put it all on the same line):

```
canonical_splice_count.py \
--gtf Mus_musculus.GRCm38.97.gtf \
--outfile count_matrix.txt \
--statsfile count_stats.txt \
*.bam
```

This will quantitate all BAM files in the current folder using the introns from the gtf file to filter the reads.  The BAM file doesn't need to be sorted, and it doesn't matter if the library is directional or not.

Full options
------------

You can see the complete options for the program by running ```canonical_splice_count.py --help```


```
usage: canonical_splice_count.py [-h] [--quiet] [--version] [--mapq MAPQ] [--biotype BIOTYPE] --outfile OUTFILE --statsfile STATSFILE --gtf GTF bamfiles [bamfiles ...]

Quanitate RNA-Seq BAM files based on canonical splice sites

positional arguments:
  bamfiles               The BAM files to analyse

options:
  -h, --help             show this help message and exit
  --quiet                Supress all but essential messages
  --version              show program's version number and exit
  --mapq MAPQ            Minimum MAPQ value to use (default 20)
  --biotype BIOTYPE      Gene biotype to use (default protein_coding)
  --outfile OUTFILE      The output file name for the assembled counts
  --statsfile STATSFILE  The output file for quantitation statistics
  --gtf GTF              The GTF file from which to read canonical splice sites
```

Output
======

The program generates two output files, both consisting of tab delimited text.

Stats file
----------

The stats file contains the following fields:

1. The total number of reads
2. The number of unmapped reads
3. The number of reads with MAPQ values below the threshold set (default 20)
4. The number of reads contaning no splice junctions
5. The number of reads containing unrecognised splice junctions (novel, or in genes not being quantitated)
6. The number of reads containing only valid splice junctions


Count matrix file
-----------------

The count matrix file contains for all of the data quantitated in a single run.  It has one line per gene and the colums represent the BAM files in the order which they were supplied to the program.  Values are raw counts.

```
Gene                    file1_GRCm38_hisat2.bam    file2_GRCm38_hisat2.bam
ENSMUSG00000000001      1842                       1803
ENSMUSG00000000003      0                          0
ENSMUSG00000000028      4034                       3204
ENSMUSG00000000037      4236                       6355
ENSMUSG00000000049      8                          0
ENSMUSG00000000056      650                        437
ENSMUSG00000000058      0                          0
ENSMUSG00000000078      62                         12
ENSMUSG00000000085      824                        1582
ENSMUSG00000000088      23277                      15153
ENSMUSG00000000093      0                          0
ENSMUSG00000000094      1420                       1135
[etc]
```














