
"""
A module to write the output in different formats
"""

import os
import sys
__author__ = 'Rob Edwards'


def write_gff3(output_dir, pp, verbose=False):
    """
    Write GFF3 code. This was adapted from code contribued by [Jose Francisco Sanchez-Herrero]
    (https://github.com/JFsanchezherrero/)

    :param output_dir: the location to write to
    :param pp: the list of prophage objects
    :param verbose: more output
    """

    ## GFF output
    out_GFF = open(os.path.join(output_dir, "prophage.gff3"), "w")
    out_GFF.write("##gff-version 3")
    out_GFF.write("\n")
    
    ## loop through pp results
    for i in pp:
        ## list
        list_write = ["pp" + str(i), str(pp[i]['contig']), str(pp[i]['start']), str(pp[i]['stop'])]

        ## GFF  
        ## strand is not known...
        out_GFF.write(str(pp[i]['contig']) + 
                      '\tPhiSpy\tprophage_region\t' +
                      str(pp[i]['start']) + '\t' + 
                      str(pp[i]['stop']) + '\t.\t.' + 
                      ' \t.\tID=pp' + str(i))
        out_GFF.write('\n')

        ## is there att sequence?
        if 'att' not in pp[i]:
            pp[i]['att']=""
        else:
            ## attL 
            out_GFF.write(str(pp[i]['contig']) +
                          '\tPhiSpy\tattL\t' +
                          str(pp[i]['att'][0]) + '\t' + 
                          str(pp[i]['att'][1]) + '\t.\t.' + 
                          ' \t.\tID=pp' + str(i))
            out_GFF.write('\n')

            ## attR 
            out_GFF.write(str(pp[i]['contig']) +
                          '\tPhiSpy\tattR\t' +
                          str(pp[i]['att'][2]) + '\t' +
                          str(pp[i]['att'][3]) + '\t.\t.' +
                          ' \t.\tID=pp' + str(i))
            out_GFF.write('\n')
        
    out_GFF.close()

def write_phage_and_bact(output_dir, pp, dna):
    # print(pp)
    print('writing bacterial and phage DNA')
    phage_out = open(os.path.join(output_dir, "phage.fasta"), "w")
    bacteria_out = open(os.path.join(output_dir, "bacteria.fasta"), "w")
    
    contig_to_phage = {}
    for i in pp:
        contig = pp[i]['contig']  # contig name
        # print(contig)
        if contig not in contig_to_phage:
            contig_to_phage[contig] = set()
        contig_to_phage[contig].add(i)
    # print('dictionary of contigs with phages')
    # print(contig_to_phage)

    for contig in dna:
        if contig not in contig_to_phage:
            bacteria_out.write(f">{contig}\n{dna[contig]}\n")
            continue
        pps1 = list(contig_to_phage[contig])
        pps = sorted(pps1, key=lambda k: pp[k]['start'])
        bactstart = 0
        for ppnum in pps:
            pphagestart = pp[ppnum]['start']
            pphagestop = pp[ppnum]['stop']

            phageseq = dna[contig][pphagestart:pphagestop]
            phage_out.write(f">{contig}_{pphagestart}_{pphagestop} [pp {ppnum}]\n{phageseq}\n")

            dnaseq = dna[contig][bactstart:pphagestart - 1]
            nphageseq = str.replace(phageseq, 'A', 'N').replace('T', 'N').replace('C', 'N').replace('G', 'N')
            bacteria_out.write(f">{contig}_{i}_{pphagestart}\n{dnaseq}{nphageseq}\n")

            bactstart = pphagestop + 1
        dnaseq = dna[contig][bactstart:]
        bacteria_out.write(f">{contig}_{i}_{len(dna[contig])}\n{dnaseq}\n")

    phage_out.close()
    bacteria_out.close()
