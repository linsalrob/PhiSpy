import re
import sys
import string
import os
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
#from numpy import array
import numpy as np
from argparse import Namespace

def find_training_genome(trainingFlag, INSTALLATION_DIR):
    try:
        f = open(os.path.join(INSTALLATION_DIR, 'data/trainingGenome_list.txt'), 'r')
    except:
        print('cannot open ' + os.path.join(INSTALLATION_DIR, 'data/trainingGenome_list.txt'))
        return ''

    for line in f:
        temp = re.split('\t',line.strip())
        if int(temp[0]) == trainingFlag:
            f.close()
            return temp[1].strip()
    return ''

def call_randomforest(**kwargs):
    output_dir = kwargs['output_dir']
    trainingFile = kwargs['training_set']
    bin_path = os.path.join(os.path.dirname(os.path.dirname(os.path.relpath(__file__))),'bin')
    infile = os.path.join(output_dir, "testSet.txt")
    outfile = os.path.join(output_dir, "classify.tsv")
    train_data = np.genfromtxt(fname=trainingFile, delimiter="\t", skip_header=1, filling_values=1) # why not fill missing values with 0?
    test_data = np.genfromtxt(fname=infile, delimiter="\t", skip_header=1, filling_values=1)
    # Przemek's comment
    # by default 10 until version 0.22 where default is 100
    # number of estimators also implies the precision of probabilities, generally 1/n_estimators
    # in R's randomForest it's 500 and the usage note regarding number of trees to grow says:
    # "This should not be set to too small a number, to ensure that every input row gets predicted at least a few times."
    clf = RandomForestClassifier(n_estimators = kwargs['randomforest_trees'])
    clf.fit(train_data[:, :-1], train_data[:, -1].astype('int'))
    np.savetxt(outfile, clf.predict_proba(test_data)[:,1])

    # cmd = "Rscript " + bin_path + "/randomForest.r " + trainingFile + " " + infile + " " + outfile
    # os.system(cmd)

def my_sort(orf_list):
     n = len(orf_list)
     i = 1
     while( i <= n ):
          j = i + 1
          while( j < n ):
               flag = 0
               #direction for both
               if( orf_list[i]['start'] < orf_list[i]['stop'] ):
                    dir_i = 1
               else:
                    dir_i = -1
               if( orf_list[j]['start'] < orf_list[j]['stop'] ):
                    dir_j = 1
               else:
                    dir_j = -1

               #check whether swap need or not
               if dir_i == dir_j:
                    if orf_list[i]['start']>orf_list[j]['start']:
                         flag = 1
               else:
                    if dir_i == 1:
                         if orf_list[i]['start']>orf_list[j]['stop']:
                              flag = 1
                    else:
                         if orf_list[i]['stop']>orf_list[j]['start']:
                              flag = 1
               #swap
               if flag == 1:
                    temp = orf_list[i]
                    orf_list[i] = orf_list[j]
                    orf_list[j] = temp
               j = j+1
          i = i+1
     return orf_list

def find_mean(all_len):
     sum = 0.0
     for i in all_len:
          sum = sum + i
     return float(sum)/len(all_len)

def calc_pp(func):
    x = 0
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
           x = 1
    elif ('unknown' in func) or ('hypothetical' in func):
           x = 0.5
    else:
           x = 0
    if 'recombinase' in func or 'integrase' in func:
           x = 1.5
    if ('phage' in func) and ('shock' in func):
           x = 0
    if "dna binding domain" in func:
           x = 0
    return x

def calc_function_3files(organism):
    my_func = {}
    x = 0 #no need it for computation.. just a flag in "except"
    try:
        f_fun = open(organism+'/proposed_non_ff_functions','r')
        for line in f_fun:
            temp = re.split('\t',line.strip())
            if len(temp)>=2:
                my_func[temp[0]] = temp[1]
        f_fun.close()
    except:
        x = x + 1
    try:
        f_fun = open(organism+'/proposed_functions','r')
        for line in f_fun:
            temp = re.split('\t',line.strip())
            if len(temp)>=2:
                my_func[temp[0]] = temp[1]
        f_fun.close()
    except:
        x = x + 1
    try:
        f_fun = open(organism+'/assigned_functions','r')
        for line in f_fun:
            temp = re.split('\t',line.strip())
            if len(temp)>=2:
                my_func[temp[0]] = temp[1]
        f_fun.close()
    except:
        x = x + 1

    return my_func

def input_bactpp(**kwargs):
    print(kwargs)
    exit()
   
    #bact_file = organism+'/Features/peg/tbl'
    #try:
    #    fh = open(bact_file,'r')
    #except:
    #    print('cant open file- assigned functions/tbl file:',organism)
    #    return {}
    #my_func = calc_function_3files(organism)
    all_orf_list = {}
    for i in fh:
        temp = re.split('\t',i.strip())
        temp1 = re.split('_',temp[1])
        if ',' in temp[1]:
            ttemp = re.split(',',temp[1])
            temp[1] = ttemp[len(ttemp)-1]
        temp1 = re.split('_',temp[1])
        contig = temp[1][:temp[1][:temp[1].rfind('_')].rfind('_')]
        start = int(temp1[len(temp1)-2])
        stop = int(temp1[len(temp1)-1])
        #save info for sorting orf
        if contig in all_orf_list:
            x = len(all_orf_list[contig]) + 1
        else:
            x = 1
            all_orf_list[contig]={}
        all_orf_list[contig][x]={}
        all_orf_list[contig][x]['fig'] = temp[0]
        all_orf_list[contig][x]['contig'] = str(contig)
        all_orf_list[contig][x]['start'] = start
        all_orf_list[contig][x]['stop'] = stop
        if temp[0] in my_func:
            all_orf_list[contig][x]['function'] = my_func[temp[0]]
            all_orf_list[contig][x]['pp'] = calc_pp(my_func[temp[0]].lower(),INSTALLATION_DIR)
        else:
            all_orf_list[contig][x]['function'] = "-"
            all_orf_list[contig][x]['pp'] = 0.5
    fh.close()
    all = {}
    index = 1
    for mycontig in all_orf_list:
        orf_list = my_sort(all_orf_list[mycontig])
        i = 1
        while i <= len(orf_list):
            all[index] = {}
            all[index]['fig'] = orf_list[i]['fig']
            all[index]['function'] = orf_list[i]['function']
            all[index]['contig'] = orf_list[i]['contig']
            all[index]['start'] = orf_list[i]['start']
            all[index]['stop'] = orf_list[i]['stop']
            all[index]['rank'] = 0.0
            all[index]['status'] = 0
            all[index]['pp'] = orf_list[i]['pp']
            i = i+1
            index = index+1
    return all

def make_initial_tbl(**kwargs): #organismPath, output_dir, window, INSTALLATION_DIR):
    self = Namespace(**kwargs)
    x = []
    for entry in self.record:
        for feature in entry.get_features('CDS'):
            all = {}
            all['fig'] = feature.id
            all['function'] = feature.function
            all['contig'] = entry.id
            all['start'] = feature.start
            all['stop'] = feature.stop
            all['rank'] = 0.0
            all['status'] = 0
            all['pp'] = calc_pp(feature.function)
            x.append(all)
    try:
        infile = open(os.path.join(self.output_dir, 'classify.tsv'), 'r')
        outfile = open(os.path.join(self.output_dir, 'initial_tbl.tsv'), 'w')
    except:
        sys.exit('ERROR: Cannot open classify.tsv in make_initial_tbl')
    #x = input_bactpp(**kwargs)
    j = 0
    ranks = [[] for n in range(len(x))]
    for line in infile:
        val = float(line.strip())
        for k in range(j-int(self.window_size/2), j+int(self.window_size/2)):
            if k < 0 or k >= len(x) or j >= len(x) or x[k]['contig'] != x[j]['contig']:
                continue
            ranks[k].append(val)
        j += 1
    infile.close()
    #calculate threshold
    y = []
    j = 0
    while j < len(x):
        x[j]['rank'] = sum(ranks[j]) / len(ranks[j]) 
        x[j]['extra'] =  ranks[j] 
        y.append(x[j]['rank'])
        #y.append([x[j]['rank']])
        j = j+1
    #threshold = max(y)/2
    y2 = np.array(y).reshape(-1, 1)
    km = KMeans(n_clusters = 2)
    km.fit(y2)
    centers = km.cluster_centers_
    threshold = max(centers[0][0], centers[1][0])
    """
    Note added by Rob:
    At this point we have the classifications for each ORF and we want to take a sliding window and decide where the phage should
    start. We have two calculations for a threshold for the rank: either the kmeans centers and finding things above the larger center
    or just a plain threshold.
    
    """
    j = 0
    outfile.write('fig_no\tfunction\tcontig\tstart\tstop\tposition\trank\tmy_status\tpp\tFinal_status\tstart of attL\tend of attL\tstart of attR\tend of attR\tsequence of attL\tsequence of attR\tReason for att site\n')
    while j < len(x):
        if x[j]['rank'] > threshold:
            x[j]['status'] = 1
        outfile.write(str(x[j]['fig']))
        outfile.write('\t')
        outfile.write(str(x[j]['function']))
        outfile.write('\t')
        outfile.write(str(x[j]['contig']))
        outfile.write('\t')
        outfile.write(str(x[j]['start']))
        outfile.write('\t')
        outfile.write(str(x[j]['stop']))
        outfile.write('\t')
        outfile.write(str(j))
        outfile.write('\t')
        outfile.write(str(x[j]['rank']))
        outfile.write('\t')
        outfile.write(str(x[j]['status']))
        outfile.write('\t')
        outfile.write(str(x[j]['pp']))
        outfile.write('\n')
        j = j+1
    outfile.close()
