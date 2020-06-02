
"""
A module to write the output in different formats
"""

import os
import gzip


from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import OrderedDict

import PhiSpyModules.version as version
from .formatting import message
from .helper_functions import is_gzip_file
import logging

__author__ = 'Rob Edwards'

def log_and_message(msg, c="WHITE", stderr=False, stdout=False, quiet=False, loglevel="INFO"):
    """
    Write the message to both the log and an output stream. By default we will just log the message in the
    log.

    Note that we also adhere to the quiet option of self, and will only write to the log
    Set either stderr or stdout to true to write to those streams too (but you will need to reset quiet if appropriate)

    :param self: The arg parser object
    :param msg: the message to write
    :param c: the color to write the message
    :param stderr: write the message to stderr
    :param stdout: write the message to stdout
    :param loglevel: the logging level. See https://docs.python.org/3/library/logging.html#levels for a list
    :return:
    """

    if loglevel == 'CRITICAL':
        logging.getLogger('PhiSpy').critical(msg.strip())
    elif loglevel == 'ERROR':
        logging.getLogger('PhiSpy').error(msg.strip())
    elif loglevel == 'WARNING':
        logging.getLogger('PhiSpy').warning(msg.strip())
    elif loglevel == 'DEBUG':
        logging.getLogger('PhiSpy').debug(msg.strip())
    else:
        logging.getLogger('PhiSpy').info(msg.strip())

    if not quiet:
        if stderr:
            message(msg,c, "stderr")
        if stdout:
            message(msg, c, "stdout")


def write_gff3(self, pp):
    """
    Write GFF3 code. This was adapted from code contribued by [Jose Francisco Sanchez-Herrero]
    (https://github.com/JFsanchezherrero/)

    :param self: the data object
    :param pp: array of pp dictionaries
    :return: None
    """

    log_and_message("Writing GFF3 output file", c="GREEN", stderr=True, quiet=self.quiet)
    # GFF output
    out_gff = open(os.path.join(self.output_dir, self.file_prefix + "prophage.gff3"), "w")
    out_gff.write("##gff-version 3")
    out_gff.write("\n")
    
    # loop through pp results
    for i in pp:
        # GFF
        # strand is not known...
        out_gff.write(str(pp[i]['contig']) +
                      '\tPhiSpy\tprophage_region\t' +
                      str(pp[i]['start']) + '\t' + 
                      str(pp[i]['stop']) + '\t.\t.' + 
                      ' \t.\tID=pp' + str(i))
        out_gff.write('\n')

        if 'att' not in pp[i]:
            pp[i]['att'] = ""
        else:
            # attL
            out_gff.write(str(pp[i]['contig']) +
                          '\tPhiSpy\tattL\t' +
                          str(pp[i]['att'][0]) + '\t' + 
                          str(pp[i]['att'][1]) + '\t.\t.' + 
                          ' \t.\tID=pp' + str(i))
            out_gff.write('\n')

            # attR
            out_gff.write(str(pp[i]['contig']) +
                          '\tPhiSpy\tattR\t' +
                          str(pp[i]['att'][2]) + '\t' +
                          str(pp[i]['att'][3]) + '\t.\t.' +
                          ' \t.\tID=pp' + str(i))
            out_gff.write('\n')
        
    out_gff.close()


def write_genbank(self, pp):
    """
    Write prophages and their potential attachment sites in updated input GenBank file.
    :param self: the data object
    :param pp: array of pp dictionaries
    :return: None
    """

    log_and_message("Writing GenBank output file", c="GREEN", stderr=True, quiet=self.quiet)
    prophage_feature_type = 'misc_feature'  # / prophage_region
    outfile = os.path.join(self.output_dir, self.file_prefix + os.path.basename(self.infile))
    for i in pp:
        self.record.get_entry(pp[i]['contig']).append_feature(SeqFeature(
                    location=FeatureLocation(pp[i]['start'], pp[i]['stop']),
                    type=prophage_feature_type,
                    strand=1,
                    qualifiers=OrderedDict(
                        {'note': f'prophage region pp{i} identified with PhiSpy v{version.__version__}'}
                    )))
        if 'atts' in pp[i]:
            self.record.get_entry(pp[i]['contig']).append_feature(SeqFeature(
                        location=FeatureLocation(int(pp[i]['att'][0]), int(pp[i]['att'][1])) +
                                 FeatureLocation(int(pp[i]['att'][2]), int(pp[i]['att'][3])),
                        type='repeat_region',
                        strand=1,
                        qualifiers=OrderedDict({'note': f'prophage region pp{i} potential attachment sites'})))

    # are we writing a gzip file
    if is_gzip_file(self.infile):
        handle = gzip.open(outfile, 'wt')
    else:
        handle = open(outfile, 'w')

    SeqIO.write(self.record, handle, 'genbank')


def write_phage_and_bact(self, pp, dna):
    """
    Separate out the phage and bacterial fractions into fasta files
    :param self: the data object
    :param pp: array of pp dictionaries
    :param dna: the DNA sequence object
    :return: None
    """
    log_and_message('Writing bacterial and phage DNA as fasta', c="GREEN", stderr=True, quiet=self.quiet)
    phage_out = open(os.path.join(self.output_dir, self.file_prefix + "phage.fasta"), "w")
    phage_genbank = open(os.path.join(self.output_dir, self.file_prefix + "phage.gbk"), "w")
    bacteria_out = open(os.path.join(self.output_dir, self.file_prefix + "bacteria.fasta"), "w")
    bacteria_genbank = open(os.path.join(self.output_dir, self.file_prefix + "bacteria.gbk"), "w")
    
    contig_to_phage = {}
    for i in pp:
        contig = pp[i]['contig']  # contig name
        if contig not in contig_to_phage:
            contig_to_phage[contig] = set()
        contig_to_phage[contig].add(i)

    for contig in dna:
        if contig not in contig_to_phage:
            bacteria_out.write(f">{contig}\n{dna[contig]}\n")
            SeqIO.write(self.record.get_entry(contig), bacteria_genbank, "genbank")
            continue
        pps1 = list(contig_to_phage[contig])
        pps = sorted(pps1, key=lambda k: pp[k]['start'])
        bactstart = 0
        dnaseq = ""
        hostcounter = 0
        for ppnum in pps:
            pphagestart = pp[ppnum]['start']
            pphagestop = pp[ppnum]['stop']

            phageseq = dna[contig][pphagestart:pphagestop]
            phage_out.write(f">{contig}_{pphagestart}_{pphagestop} [pp {ppnum}]\n{phageseq}\n")
            pp_gbk = self.record.get_entry(contig)[pphagestart:pphagestop]
            # set the locus tag
            pp_gbk.name += f"_PP{ppnum}"
            # set the Accession
            pp_gbk.id += f"_PP{ppnum}"
            #set the definition line
            pp_gbk.description += f" prophage PP{ppnum} on {contig} from {pphagestart} to {pphagestop}"
            SeqIO.write(pp_gbk, phage_genbank, "genbank")

            dnaseq += dna[contig][bactstart:pphagestart - 1]
            dnaseq += "N" * len(phageseq)
            host_gbk = self.record.get_entry(contig)[bactstart:pphagestart - 1]
            hostcounter += 1
            host_gbk.name += f"_region_{hostcounter}"
            host_gbk.id += f"_region_{hostcounter}"
            host_gbk.description += f" region {hostcounter} on {contig} from {bactstart} to {pphagestart - 1}"
            SeqIO.write(host_gbk, bacteria_genbank, "genbank")
            bactstart = pphagestop + 1

        dnaseq += dna[contig][bactstart:]
        host_gbk = self.record.get_entry(contig)[bactstart:]
        hostcounter += 1
        host_gbk.name += f"_region_{hostcounter}"
        host_gbk.id += f"_region_{hostcounter}"
        host_gbk.description += f" region {hostcounter} on {contig} from {bactstart} onwards"
        bacteria_out.write(f">{contig} [phage regions replaced with N]\n{dnaseq}\n")
        SeqIO.write(host_gbk, bacteria_genbank, "genbank")

    phage_out.close()
    phage_genbank.close()
    bacteria_out.close()
    bacteria_genbank.close()


def write_prophage_coordinates(self, pp):
    """
    Write the coordinates and other details about the prophages
    :param self: the data object
    :param pp: array of pp dictionaries
    :return: None
    """

    log_and_message("Writing prophage_coordinates output file", c="GREEN", stderr=True, quiet=self.quiet)
    with open(os.path.join(self.output_dir, self.file_prefix + 'prophage_coordinates.tsv'), 'w') as out:
        for i in pp:
            if 'atts' not in pp[i]:
                pp[i]['atts'] = ""
            locs = [
                "pp" + str(i),
                pp[i]['contig'],
                pp[i]['start'],
                pp[i]['stop'],
                pp[i]['atts']
            ]
            out.write("\t".join(map(str, locs)) + "\n")


def write_prophage_tbl(self, pp):
    """
    Create a prophage_tbl file from our pp dictionary
    :param self: the data object
    :param pp: array of pp dictionaries
    :return: None
    """

    log_and_message("Writing prophage.tbl output file", c="GREEN", stderr=True, quiet=self.quiet)
    with open(os.path.join(self.output_dir, self.file_prefix + "prophage.tbl"), 'w') as out:
        for i in pp:
            locs = [
                "pp_" + str(i),
                pp[i]['contig'],
                pp[i]['start'],
                pp[i]['stop']
            ]
            out.write(locs[0] + "\t" + "_".join(map(str, locs[1:])) + "\n")


def write_prophage_tsv(self, pp):
    """
    Create a tsv with headers for this data. Issue #28 item 2
    :param self: the data object
    :param pp: array of pp dictionaries
    :return: None
    """

    log_and_message("Writing prophage.tsv output file", c="GREEN", stderr=True, quiet=self.quiet)
    with open(os.path.join(self.output_dir, self.file_prefix + "prophage.tsv"), 'w') as out:
        out.write("Prophage number\tContig\tStart\tStop\n")
        for i in pp:
            locs = [
                "pp_" + str(i),
                pp[i]['contig'],
                pp[i]['start'],
                pp[i]['stop']
            ]
            out.write("\t".join(map(str, locs)) + "\n")


def write_test_data(self, pp):
    """
    Write the testing measurements
    :param self: the data object
    :param measurements: the measurements
    :return: None
    """
    log_and_message("Writing test_data output file", c="GREEN", stderr=True, quiet=self.quiet)
    with open(os.path.join(self.output_dir, self.file_prefix + "test_data.tsv"), 'w') as out:
        if self.phmms:
            out.write('identifier\torf_length_med\tshannon_slope\tat_skew\tgc_skew\tmax_direction\tphmms\tstatus\n')
        else:
            out.write('identifier\torf_length_med\tshannon_slope\tat_skew\tgc_skew\tmax_direction\tstatus\n')

        for i, d in enumerate(self.test_data):
            out.write("\t".join(map(str, [self.initial_tbl[i][0]] + d)) + "\n")


def prophage_measurements_to_tbl(inputf, outputf):
    try:
        f = open(inputf, 'r')
        fw = open(outputf, 'w')
    except IOError as e:
        log_and_message(f"There was an error trying to convert the measurements to a tbl: {e}", c="RED", stderr=True)
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


