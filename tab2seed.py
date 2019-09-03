"""
Convert a tab file from PATRIC to a the SEED files that we need for PhiSpy

We need the following files:

1. assigned_functions - a tab separated list of FIG ID and function
2. contigs - the fasta DNA sequence. Note we may download this separately
3. genome - the name of the genome -- may not be required
4. taxonomy - the taxonomy of the genome -- may not be required
5. taxonomy_id - the tax id. -- also may not be required
6. Features/peg/tbl - the tbl that has id,contig_start_stop, [alt ids]
7. Features/rna/tbl - the RNA genes


The files that PhiSpy opens are:

a. dir/contigs
b. dir/Features/peg/tbl
c. dir/assigned_functions
d. dir/Features/rna/tbl


"""

import os
import sys
import argparse

def parse_tab(filename, outputdir):
    """
    Parse a patric tab file
    :param filename: the file to parse
    :return: ummm
    """

    if not (os.path.exists(os.path.join(outputdir, "Features"))):
        os.mkdir(os.path.join(outputdir, "Features"))
    if not (os.path.exists(os.path.join(outputdir, "Features/peg"))):
        os.mkdir(os.path.join(outputdir, "Features/peg"))
    if not (os.path.exists(os.path.join(outputdir, "Features/rna"))):
        os.mkdir(os.path.join(outputdir, "Features/rna"))

    peg = open(os.path.join(outputdir, "Features/peg/tbl"), 'w')
    rna = open(os.path.join(outputdir, "Features/rna/tbl"), 'w')
    asf = open(os.path.join(outputdir, "assigned_functions"), 'w')

    wrote_genome = False

    with open(filename, 'r') as fin:
        for l in fin:
            if l.startswith('genome_id'):
                continue

            # genome_id	genome_name	accession	annotation	feature_type	patric_id	refseq_locus_tag	alt_locus_tag
            # uniprotkb_accession	start	end	strand	na_length	gene	product	figfam_id	plfam_id	pgfam_id
            # go	ec	pathway
            l = l.replace("\n", "") # this is a hack because I can't figure out how to do chomp
            p = l.split("\t")

            if not wrote_genome:
                with open(os.path.join(outputdir, "GENOME"), 'w') as gout:
                    gout.write("{}\n".format(p[1]))
                wrote_genome = True

            gid, name, acc, who, ftype, fid, refseq_locus, alt, uni, start, stop, strand, length, gene, prod, ffid, plid, pgid, go, ec, pw = p

            if start > stop:
                (start, stop) = (stop, start)

            if "CDS" in p[4]:
                peg.write("{}\t{}_{}_{}\n".format(fid, acc, start, stop))
                asf.write("{}\t{}\n".format(fid, prod))
            elif "rna" in p[4].lower():
                rna.write("{}\t{}_{}_{}\n".format(fid, acc, start, stop))
    peg.close()
    rna.close()
    asf.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert a patric tab file to a minimal seed directory")
    parser.add_argument('-f', help='The patric tab file', required=True)
    parser.add_argument('-o', help='output directory', required=True)
    parser.add_argument('-v', help='verbose output', action="store_true")
    args = parser.parse_args()

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    parse_tab(args.f, args.o)