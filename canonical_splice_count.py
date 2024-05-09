#!/usr/bin/env python3

import pysam
import argparse
import gzip

VERSION = "1.0"

def main():
    global options
    options = read_options()

    genes,valid_introns = read_canonical_splices(options.gtf)

    breakpoint()

    counts = []

    for bamfile in options.bamfiles:
        counts.append(quantitate_bam_file(bamfile, genes, valid_introns))

    write_output(counts)

def quantitate_bam_file(file,genes, introns):
    bam = pysam.AlignmentFile(file, "rb")

    for read in bam.fetch(until_eof=True):
        if read.mapping_quality < options.mapq:
            continue

def read_canonical_splices(file):

    genes = set()

    exons = {}

    if file.lower().endswith(".gz"):
        fh = gzip.open(file, "rt", encoding="utf8")

    else:
        fh = open(file,"rt",encoding="utf8")


    for line in fh:
        if line.startswith("#"):
            continue

        sections = line.strip().split("\t")

        if sections[2] == "gene":

            biotype = None
            gene_id = None
            # We need to extract the biotype
            tags = sections[8].strip().split(";")
            for tag in tags:
                if not tag:
                    continue
                t,value = tag.split(maxsplit=1)
                if t=="gene_biotype":
                    biotype = value.replace("\"","")
                    
                if t=="gene_id":
                    gene_id = value.replace("\"","")
            
            if biotype is not None and biotype == options.biotype:
                genes.add(gene_id)

            continue


        if not sections[2] == "exon":
            continue

        chrom = sections[0]
        start = int(sections[3])
        end = int(sections[4])

        if start > end:
            temp = start
            start = end
            end = temp

        gene_id = None
        transcript_id = None

        tags = sections[8].split(";")
        for tag in tags:
            if not tag.strip():
                continue

            t,value = tag.strip().split(maxsplit=1)
            if t=="transcript_id":
                transcript_id = value.replace("\"","")
                
            if t=="gene_id":
                gene_id = value.replace("\"","")

        if not gene_id in genes:
            # We don't want genes we haven't validated as having the correct biotype
            continue

        if not gene_id in exons:
            exons[gene_id] = {}

        if not transcript_id in exons[gene_id]:
            exons[gene_id][transcript_id] = []

        exons[gene_id][transcript_id].append((chrom,start,end))


    fh.close()
    
    # We now need to get the full list of introns for each gene and a full list of genes
    final_gene_list = []

    for gene in genes:
        if gene in exons:
            final_gene_list.append(gene)

    final_gene_list.sort()

    final_introns = {}


    for gene in final_gene_list:
        for transcript in exons[gene]:
            these_exons = exons[gene][transcript]
            these_exons.sort(key=lambda x: x[1]) # Sort by start position

            for i in range(len(these_exons)-1):
                intron = these_exons[i][0]+":"+str(these_exons[i][2]+1)+"-"+str(these_exons[i+1][1]-1)

                final_introns[intron] = gene

    return final_gene_list, final_introns


def read_options():
    parser = argparse.ArgumentParser(description="Quanitate RNA-Seq BAM files based on canonical splice sites")

    parser.add_argument('--quiet', dest="quiet", action='store_true', default=False, help="Supress all but essential messages")
    parser.add_argument('--version', action='version', version=f"Canonical Splice Count v{VERSION}")
    parser.add_argument('--mapq', type=int, default=20, help="Minimum MAPQ value to use (default 20)")
    parser.add_argument('--biotype', type=str, default="protein_coding", help="Gene biotype to use (default protein_coding)")
    parser.add_argument('--outfile', type=str, help="The output file name for the assembled counts", required=True)    
    parser.add_argument('--gtf', type=str, help="The GTF file from which to read canonical splice sites", required=True)
    parser.add_argument("bamfiles", nargs="+", type=str, help="The BAM files to analyse")

    options = parser.parse_args()

    return options

if __name__ == "__main__":
    main()
