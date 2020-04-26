"""
A module to write the output in different formats
"""

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import OrderedDict
import os
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