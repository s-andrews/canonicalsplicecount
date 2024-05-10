#!/usr/bin/env python3

import pysam
import argparse
import gzip
import sys

VERSION = "1.0"

def main():
    global options
    options = read_options()

    genes,valid_introns = read_canonical_splices(options.gtf)

    counts = []

    statsfile = open(options.statsfile,"wt", encoding="utf8")
    stats_header = ["File","Total_Seqs","Unmapped","Low_MAPQ","Unspliced","Non_Canonical","Valid"]
    print("\t".join(stats_header),file=statsfile)

    for bamfile in options.bamfiles:
        counts.append(quantitate_bam_file(bamfile, genes, valid_introns, statsfile))

    statsfile.close()

    write_output(options.bamfiles, genes, counts, options.outfile)

def write_output(bamfiles, genes, counts,outfile):
    if not options.quiet:
        print("Writing results to",outfile, file=sys.stderr)

    with open(outfile,"wt",encoding="utf8") as out:
        # Header
        headerline = ["Gene"]
        headerline.extend(bamfiles)
        print("\t".join(headerline), file=out)

        for gene in genes:
            dataline = [gene]
            for x in counts:
                dataline.append(str(x[gene]))

            print("\t".join(dataline), file=out)



def quantitate_bam_file(file,genes, introns, statsfile):

    if not options.quiet:
        print("Quantitating",file, file=sys.stderr)

    bam = pysam.AlignmentFile(file, "rb")

    total_reads = 0
    unmapped = 0
    low_mapq = 0
    unspliced = 0
    non_canonical = 0
    canonical = 0

    counts = {}

    for gene in genes:
        counts[gene] = 0


    for read in bam.fetch(until_eof=True):
        total_reads += 1

        if total_reads % 1000000 == 0:
            if not options.quiet:
                print("Read",int(total_reads/1000000),"million reads", file=sys.stderr)

        if not read.is_mapped:
            unmapped += 1
            continue

        if read.mapping_quality < options.mapq:
            low_mapq += 1
            continue

        cigar_tuples = read.cigartuples
        if len(cigar_tuples) == 1:
            unspliced += 1
            continue

        # Now we work our way through the tuples, to find
        # splice sites.

        chromosome = read.reference_name
        current_position = read.reference_start

        found_gene = None

        for opcode,oplength in cigar_tuples:
            # We go through the different codes
            if opcode == 0:
                # Match to the reference
                current_position += oplength

            elif opcode == 1:
                # Insertion
                pass

            elif opcode == 2:
                # Deletion
                current_position += oplength
            
            elif opcode == 3:
                # A splice site.
                splice_start = current_position+1
                splice_end = current_position+oplength
                current_position += oplength

                # Let's see if this is a valid splice
                splice_string = f"{chromosome}:{splice_start}-{splice_end}"
                if splice_string in introns:
                    if found_gene is None:
                        found_gene = introns[splice_string]
                    else:
                        if introns[splice_string] != found_gene:
                            # We've got a mix of genes, so we'll bin this
                            found_gene = "NON_CANONICAL"
                            break
                else:
                    # It's a splice but not as we know it
                    found_gene = "NON_CANONICAL"
                    break

        if found_gene is None:
            # Turns out there wasn't a splice in here after all
            unspliced += 1
        elif found_gene == "NON_CANONICAL":
            non_canonical +=1
        else:
            # We can add to the count for this gene
            canonical += 1
            counts[found_gene] += 1

    # Print out the stats
    statsline = [file, total_reads, unmapped, low_mapq, unspliced, non_canonical, canonical]

    print("\t".join([str(x) for x in statsline]), file=statsfile)

    return counts

def read_canonical_splices(file):

    if not options.quiet:
        print("Reading introns from",file, file=sys.stderr)

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


    if not options.quiet:
        print("Found",len(final_gene_list),"genes to quantitate", file=sys.stderr)

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
    parser.add_argument('--statsfile', type=str, help="The output file for quantitation statistics", required=True)    
    parser.add_argument('--gtf', type=str, help="The GTF file from which to read canonical splice sites", required=True)
    parser.add_argument("bamfiles", nargs="+", type=str, help="The BAM files to analyse")

    options = parser.parse_args()

    return options

if __name__ == "__main__":
    main()
