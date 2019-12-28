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

