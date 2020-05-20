
"""
A module to write the output in different formats
"""

import os
import sys
import PhiSpyModules.version as version
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

def write_genbank(infile, record, output_directory, pp):
    """
    Write prophages and their potential attachment sites in updated input GenBank file.
    :param infile: path to input file
    :param record: SeqRecord generator of input file
    :param output_dir: the location to write to
    :param pp: the list of prophage objects
    """

    prophage_feature_type = 'misc_feature' # / prophage_region
    outfile = os.path.join(output_directory, os.path.basename(infile))
    for i in pp:
        record.get_entry(pp[i]['contig']).add_feature(SeqFeature(
                    location = FeatureLocation(pp[i]['start'], pp[i]['stop']),
                    type = prophage_feature_type,
                    strand = 1,
                    qualifiers = OrderedDict({'note': f'prophage region pp{i} identified with PhiSpy v{version.__version__}'})))
        record.get_entry(pp[i]['contig']).add_feature(SeqFeature(
                    location = FeatureLocation(int(pp[i]['att'][0]), int(pp[i]['att'][1])) + FeatureLocation(int(pp[i]['att'][2]), int(pp[i]['att'][3])),
                    type = 'repeat_region',
                    strand = 1,
                    qualifiers = OrderedDict({'note': f'prophage region pp{i} potential attachment sites'})))

    SeqIO.write(record, outfile, 'genbank')


def write_phage_and_bact(output_dir, pp, dna):
    # sys.stderr.write(pp)
    sys.stderr.write('writing bacterial and phage DNA')
    phage_out = open(os.path.join(output_dir, "phage.fasta"), "w")
    bacteria_out = open(os.path.join(output_dir, "bacteria.fasta"), "w")
    
    contig_to_phage = {}
    for i in pp:
        contig = pp[i]['contig']  # contig name
        # sys.stderr.write(contig)
        if contig not in contig_to_phage:
            contig_to_phage[contig] = set()
        contig_to_phage[contig].add(i)
    # sys.stderr.write('dictionary of contigs with phages')
    # sys.stderr.write(contig_to_phage)

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

def write_prophage_tbl(outputdir, pp):
    """
    Create a prophage_tbl file from our pp dictionary
    :param outputdir: the directory to write to
    :param pp: array of pp dictionaries
    :return: None
    """

    with open(os.path.join(outputdir, "prophage.tbl"), 'w') as out:
        for i in pp:
            out.write("pp_" + str(i) + "\t" + str(pp[i]['contig']) + "_" + str(pp[i]['start']) + "_" + str(pp[i]['stop']) + "\n")

def write_prophage_tsv(outputdir, pp):
    """
    Create a tsv with headers for this data. Issue #28 item 2
    :param outputdir: the directory to write to
    :param pp: array of pp dictionaries
    :return: None
    """
    with open(os.path.join(outputdir, "prophage.tsv"), 'w') as out:
        out.write("Prophage number\tContig\tStart\tStop\n")
        for i in pp:
            out.write("pp_" + str(i) + "\t" + str(pp[i]['contig']) + "\t" + str(pp[i]['start']) + "\t" + str(pp[i]['stop']) + "\n")




def prophage_measurements_to_tbl(inputf, outputf):
    try:
        f = open(inputf, 'r')
        fw = open(outputf, 'w')
    except:
        print('Cant open', inputf, ' or ', outputf)
        return
    pp = {}
    ppindx = 0
    prev_contig = None
    inphage = False
    header = f.readline()
    for line in f:
        temp = line.strip().split("\t")
        if int(temp[9]) > 0:
            newphage = False
            if temp[2] != prev_contig:
                newphage = True
            if not inphage:
                newphage = True
            if newphage:
                ppindx += 1
                pp[ppindx] = {}
                pp[ppindx]['contig'] = temp[2]
                pp[ppindx]['start'] = min(int(temp[3]), int(temp[4]))
                pp[ppindx]['stop'] = max(int(temp[3]), int(temp[4]))
            else:
                pp[ppindx]['stop'] = max(int(temp[3]), int(temp[4]))
            inphage = True
            prev_contig = temp[2]
        else:
            inphage = False
    for i in pp:
        fw.write("pp_" + str(i) + "\t" + pp[i]['contig'] + "_" + str(pp[i]['start']) + "_" + str(pp[i]['stop']) + "\n")
    f.close()
    fw.close()
