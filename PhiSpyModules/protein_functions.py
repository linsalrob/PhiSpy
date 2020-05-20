import re
import string
import os

def is_phage_func(func):
    func = func.lower()
    func = func.replace('-',' ')
    func = func.replace(',',' ')
    a = re.split(' ',func)
    if (
            'phage'         in a or
            'lysin'         in a or
            'endolysin'     in a or
            'holin'         in a or
            'capsid'        in a or
            'tail'          in a or
            'bacteriophage' in a or
            'prophage'      in a or
            'portal'        in a or
            'terminase'     in a or
            'tapemeasure'   in a or
            'baseplate'     in a or
            'virion'        in a or
            'antirepressor' in a or
            'excisionase'   in a or
            'mobile element protein' == func or
            re.search(r"\b%s\b" % "tape measure", func) or
            re.search(r"\b%s\b" % "Cro-like repressor", func) or
            re.search(r"\b%s\b" % "CI-like repressor", func) or
            re.search(r"\b%s\b" % "rIIA lysis", func) or
            re.search(r"\b%s\b" % "rI lysis", func) or
            re.search(r"\b%s\b" % "rIIB lysis", func) or
            re.search(r"\b%s\b" % "base plate", func)
    ):
        return True
    return False



def is_unknown_func(x):
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
       ('expressed' in x_lower) or
       ('similar to' in x_lower) or
       (' identi' in x_lower) or
       ('ortholog of' in x_lower) or
       ('structural feature' in x_lower) or
       ('cds_' in x_lower) or
       ('predicted by psort' in x_lower) or
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
        return True
    return False

def add_unknown_function_initial_tbl(infile,outfile):
    try:
        f = open(infile,'r')
        fw = open(outfile,'w')
    except:
        print('ERROR: Cannot open initial_tbl.tsv add_unknown_function_initial_tbl.')
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
        if is_unknown_func(temp[1]):
            fw.write('0.5\n')
        else:
            fw.write(temp[8]+'\n')

    f.close()
    fw.close()
    return 1

def consider_unknown(output_dir):
    it = os.path.join(output_dir, 'initial_tbl.tsv')
    it2 = os.path.join(output_dir, 'initial_tbl_2.tsv')
    x = add_unknown_function_initial_tbl(it, it2)
    if (x == 1):
        cmd2 = "rm " + it
        os.system(cmd2)
        cmd2 = "mv " + it2 + ' ' + it
        os.system(cmd2)
