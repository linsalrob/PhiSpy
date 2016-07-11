import re
import string
import os

def find_training_genome(trainingFlag,INSTALLATION_DIR):
    try:
        f = open(INSTALLATION_DIR+"data/trainingGenome_list.txt","r")
    except:
        print'cannot open '+INSTALLATION_DIR+'data/trainingGenome_list.txt'
        return ''

    for line in f:
        line = line.strip()
        temp = re.split('\t',line)
        if int(temp[0]) == trainingFlag:
            f.close()
            return temp[1].strip()
    return ''

def call_randomForest_generic(a,trainingFlag,wrtfile,INSTALLATION_DIR):
    print 'Using training flag: ', trainingFlag
    x = find_training_genome(trainingFlag,INSTALLATION_DIR)
    if len(x)<2:
        return
    cmd = "Rscript "+INSTALLATION_DIR+"source/randomForest.r "+INSTALLATION_DIR+"data/trainingSet/"+x+" "+a+" "+wrtfile+" "+str(trainingFlag)  
    os.system(cmd)
        
def my_sort(orf_list):
     n = len(orf_list)
     i = 1
     while(i<=n):
          j = i+1
          while( j<=n):
               flag = 0
               #direction for both
               if orf_list[i]['start']<orf_list[i]['stop']:
                    dir_i = 1
               else:
                    dir_i = -1
               if orf_list[j]['start']<orf_list[j]['stop']:
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

def phage_function(c,INSTALLATION_DIR):
    m = open(INSTALLATION_DIR+"data/phage_functions.txt",'r')
    for line in m:
        t = re.split('\t',line.strip())
        if len(t) == 4:
            a = t[3].upper().strip()
            a = a.replace('-',' ')
            if a in c or c in a:
                if 'hypothetical' in c:
                    continue
                m.close()
                return 1
    m.close()
    return 0
            

def calc_pp(func,INSTALLATION_DIR):
    x = 0
    func = func.replace('-',' ')
    func = func.replace(',',' ')
    a = re.split(' ',func)
    if ('phage' in a) or ('lysin' in a) or ('holin' in a) or ('capsid' in a) or ('tail' in a) or ('bacteriophage' in a) or ('prophage' in a) or ('portal' in a) or phage_function(func.lower(),INSTALLATION_DIR) == 1:
        x = 1
    else:
        if 'unknown' in func:
            x = 0.5
        else:
            x = 0
    if 'recombinase' in a or 'integrase' in a:
        x = 1.5
    if ('phage' in a) and ('shock' in a):
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

def input_bactpp(organism,INSTALLATION_DIR):
    bact_file = organism+'/Features/peg/tbl'     
    try:
        fh = open(bact_file,'r')
    except:
        print 'cant open file- assigned functions/tbl file:',organism
        return {}
 
    my_func = calc_function_3files(organism)
    
    all_orf_list = {}
    for i in fh:
        temp = re.split('\t',i.strip())
        temp1 = re.split('_',temp[1])

        if ',' in temp[1]:
            ttemp = re.split(',',temp[1])
            temp[1] = ttemp[len(ttemp)-1]
        temp1 = re.split('_',temp[1])

        #contig = temp1[len(temp1)-3]
        contig = temp[1][:temp[1][:temp[1].rfind('_')].rfind('_')]
            
        start = int(temp1[len(temp1)-2])
        stop = int(temp1[len(temp1)-1])

              
        #save info for sorting orf
        if contig in all_orf_list:
            x = len(all_orf_list[contig])+1
        else:
            x = 1
            all_orf_list[contig]={}

        all_orf_list[contig][x]={}
        all_orf_list[contig][x]['start'] = start
        all_orf_list[contig][x]['stop'] = stop
        if temp[0] in my_func:
            all_orf_list[contig][x]['fig'] = temp[0]+'\t'+my_func[temp[0]]+'\t'+contig+'\t'+str(start)+'\t'+str(stop)
            all_orf_list[contig][x]['pp'] = calc_pp(my_func[temp[0]].lower(),INSTALLATION_DIR)
        else:
            all_orf_list[contig][x]['fig'] = temp[0]+'\t-'+'\t'+contig+'\t'+str(start)+'\t'+str(stop)
            all_orf_list[contig][x]['pp'] = 0.5
    fh.close()
    
    all = {}
    index = 1
    for mycontig in all_orf_list:
        orf_list = my_sort(all_orf_list[mycontig])
        i = 1
        while i < len(orf_list):
            all[index] = {}
            all[index]['fig'] = orf_list[i]['fig']
            all[index]['status'] = 0
            all[index]['rank'] = 0.0
            all[index]['pp'] = orf_list[i]['pp']
            i = i+1
            index = index+1
    return all

def make_tbl(alist,input_classify,output_tblfile,window,INSTALLATION_DIR):
    for l in alist:
        orgf = l+'/Features/peg/tbl'
        try:
        #fr = open(orgf,'r')
            frm = open(input_classify,'r')
            fw = open(output_tblfile,'w')
        except:
            print 'cant open',input_classify
            continue
        x = input_bactpp(l,INSTALLATION_DIR)
        j = 1
        for line in frm:
            temp1 = re.split('\t',line.strip())
            val = float(temp1[0])
            k= j
            while  k < (j + window) and k < len(x):
                x[k]['rank'] = x[k]['rank'] + val
                k = k+1
            j = j+1
        frm.close()
            
        #claculate threshold
        y = []
        j = 1
        while j < len(x):
            y.append(x[j]['rank']/window)
            j = j+1
        y.sort()
        #print len(y)
        threshold = (y[len(y)-1])/2
        #print threshold,y[len(y)-1],len(y)
            
        j = 1
        fw.write('fig_no\tfunction\tcontig\tstart\tstop\tposition\trank\tmy_status\tpp\tFinal_status\tstart of attL\tend of attL\tstart of attR\tend of attR\tsequence of attL\tsequence of attR\n')
        while j < len(x):
            x[j]['rank'] = x[j]['rank']/window
            if x[j]['rank'] > threshold:
                x[j]['status'] = 1
            fw.write(str(x[j]['fig'])+'\t'+str(j)+'\t'+str(x[j]['rank'])+'\t'+str(x[j]['status'])+'\t'+str(x[j]['pp'])+'\n')
            j = j+1
        fw.close()
        

##########################################################################

def call_classificaton(organismPath,output_dir,trainingFlag,INSTALLATION_DIR):
    call_randomForest_generic(output_dir+'testSet.txt',trainingFlag,output_dir+'classify.txt',INSTALLATION_DIR)        
    make_tbl([organismPath],output_dir+'classify.txt',output_dir+'initial_tbl.txt',40,INSTALLATION_DIR)            

    cmd2 = "rm "+output_dir+'testSet.txt'
    os.system(cmd2)
    cmd2 = "rm "+output_dir+'classify.txt'
    os.system(cmd2)
