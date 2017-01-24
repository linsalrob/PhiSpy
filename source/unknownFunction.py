import re
import string
import os

def unknown_func(x):

    x_lower = x.lower()

    if (
       (len(x) == 0) or
       #('hypoth' in x_lower) or
       ('conserved protein' in x_lower) or
       ('gene product' in x_lower) or
       ('interpro' in x_lower) or
       ('uncharacterized' in x_lower) or
       ('pseudogene' in x_lower) or
       ('similar to' in x_lower) or
       ('similarity' in x_lower) or
       ('glimmer' in x_lower) or
       ('unknown' in x_lower) or
       ('complete' in x_lower) or
       ('ensang' in x_lower) or
       ('unnamed' in x_lower) or
       ('Expressed' in x_lower) or
       ('similar to' in x_lower) or
       (' identi' in x_lower) or
       ('ortholog of' in x_lower) or
       ('structural feature' in x_lower) or
       ('cds_' in x_lower) or
       ('predicted by Psort' in x) or
       ('AGR_' in x) or
       ('EG:' in x) or
       ('RIKEN' in x) or
       re.search('lmo\d+ protein', x_lower) or
       re.search('lmo\d+protein', x_lower) or
       re.search('B[sl][lr]\d', x_lower) or
       re.search('^U\d', x) or
       re.search('[a-zA-Z]{2,3}\|', x) or
       re.search('orf\d+', x_lower) or
       re.match('orf[^_]', x_lower) or
       re.match('predicted', x_lower) or
       re.match('bh\d+', x_lower) or
       re.match('y[a-z]{2,4}\\b', x) or
       re.match('[a-z]{2,3}\d+[^:\+\-0-9]', x_lower) ):
        
        return 1
    else:
        return 0

def add_unknown_function_initial_tbl(infile,outfile):
    try:
        f = open(infile,'r')
        fw = open(outfile,'w')
    except:
        return 0
    
    flag = 0
    for line in f:
        if flag == 0:
            fw.write(line)
            flag = 1
            continue

        line = line.strip()
        temp = re.split('\t',line)
        i = 0
        while i<8:
            fw.write(temp[i]+'\t')
            i = i + 1

        x = unknown_func(temp[1])
        if x == 0:
            fw.write(temp[8]+'\n')
        else:
            fw.write('0.5\n')
 
    f.close()
    fw.close()
    return 1

def consider_unknown(output_dir):
    x = add_unknown_function_initial_tbl(output_dir+'initial_tbl.txt',output_dir+'initial_tbl_2.txt')
    if (x == 1):
        cmd2 = "rm " + output_dir+'initial_tbl.txt'
        os.system(cmd2)

        cmd2 = "mv "+output_dir+'initial_tbl_2.txt '+output_dir+'initial_tbl.txt'
        os.system(cmd2)
