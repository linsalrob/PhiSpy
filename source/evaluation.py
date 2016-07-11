import os
import re
import math
import string

def read_contig(organism): 
     try:
          f_dna = open(organism+'/contigs','r')
     except:
          print 'cant open contig file ',organism
          return ''   

     dna = {}
     seq = ''
     name = ''
     for i in f_dna:
          if i[0]=='>':
               if len(seq)>10:
                    dna[name]=seq
               name = i.strip()
               if ' ' in name: #12149.1
                    temp = re.split(' ',name)
                    name = temp[0]
               '''
               if '_' in name:
                    temp = re.split('_',name)
                    name = temp[len(temp)-1]
               else:
               '''
               name = name[1:len(name)]

               seq = ''
          else:
               seq = seq+i.strip()


     dna[name]=seq
     f_dna.close()
     return dna

def find_repeat(fn,st,INSTALLATION_DIR):
    if len(fn) == 0:
        return {}
    #write dna in file
    try:
        fw = open("tempRepeatDNA.txt","w")
    except:
        print 'cant open make_dna_file'
        return {}
    
    fw.write('>hello\n'+fn)
    fw.close()

    # call repeat finder
    try:
         cmd1 =INSTALLATION_DIR+"source/repeatFinder -f tempRepeatDNA.txt" 
         #print cmd1
         os.system(cmd1)
    except:
         print 'repeat finder did not work for ',len(fn)
         return {}

    #read repeats
    try:
        f = open("repeat.txt","r")
    except:
        print 'cant open repeat.txt'
        return {}

    rep = {}
    index = 1
    for line in f:
        '''
        if 'note' in line:
            continue
        temp = re.split('join',line.strip())
        temp1 = temp[1][1:len(temp[1])-1]
        temp = re.split(',',temp1)
        t1 = re.split('\.\.',temp[0])
        t2 = re.split('\.\.',temp[1])
        if math.fabs(int(t1[0]) - int(t2[1])) > 10000:
            rep[index] = {}
            rep[index]['s1'] = int(t1[0])+st
            rep[index]['s2'] = int(t2[0])+st
            rep[index]['e1'] = int(t1[1])+st
            rep[index]['e2'] = int(t2[1])+st
            index = index + 1
        '''
        temp = re.split('\t',line.strip())
        if math.fabs(int(temp[0]) - int(temp[3])) > 10000:
            rep[index] = {}
            rep[index]['s1'] = int(temp[0])+st
            rep[index]['s2'] = int(temp[2])+st
            rep[index]['e1'] = int(temp[1])+st
            rep[index]['e2'] = int(temp[3])+st
            index = index + 1

    f.close()
    cmd1 ="rm repeat.txt" 
    os.system(cmd1)
    cmd1 ="rm tempRepeatDNA.txt"
    os.system(cmd1)
    #print cmd1
    return rep

def check_intg(prophage_sta,prophage_sto,rep,integ):
    for m in integ:
        if integ[m]['start'] < (prophage_sta+(prophage_sto-prophage_sta)/2):
            l = integ[m]['start']
        else:
            l = integ[m]['stop']
        if (l - rep['s1'] <= 500 and l - rep['s1'] > 0) or (l - rep['e1'] <= 500 and l - rep['e1'] > 0) :
            return 1
        if (l - rep['s2'] <= 500 and l - rep['s2'] > 0) or (l - rep['e2'] <= 500 and l - rep['e2'] > 0) :
            return 1
    return 0 

def find_smallest(a,b):
    mm = 1000000
    for i in a:
        for j in b:
            if math.fabs(i-j)<mm:
                mm = math.fabs(i-j)
    return mm

def find_rna(prophage_start,prophage_stop,repeat_list,org,cont,integrs):
    try:
        f = open(org+'/Features/rna/tbl','r')
    except:
        print 'cant open',org
        return '0_0'
    my_start = 1000000
    start_end = 0
    end_start = 0
    my_end = 1000000
    mydiff = 1000000
    for line in f:
        temp = re.split('\t',line.strip())
        if len(temp)<3:
            continue
        if 'trna' in temp[2].lower() or 'tmrna' in temp[2].lower():
            if ',' in temp[1]:
                 ttemp = re.split(',',temp[1])
                 temp[1] = ttemp[len(ttemp)-1]
            temp1 = re.split('_',temp[1])

            #contig = temp1[len(temp1)-3]
            contig = temp[1][:temp[1][:temp[1].rfind('_')].rfind('_')]
            
            start = int(temp1[len(temp1)-2])
            stop = int(temp1[len(temp1)-1])

            
            if cont == contig:
                i = 1
                while i <= len( repeat_list):
                    a = find_smallest([start,stop],[repeat_list[i]['s1'],repeat_list[i]['s2'],repeat_list[i]['e1'],repeat_list[i]['e2']])
                    
                    if check_intg(prophage_start,prophage_stop,repeat_list[i],integrs) == 1:
                        if (math.fabs(repeat_list[i]['s1'] - repeat_list[i]['e2'])>math.fabs(my_start-my_end)) or mydiff == 1000000: 
                            my_start = repeat_list[i]['s1']
                            my_end = repeat_list[i]['e2']
                            start_end = repeat_list[i]['e1']
                            end_start = repeat_list[i]['s2']
                            mydiff = a
                    
                    if (a <= 500 and a < mydiff):
                        my_start = repeat_list[i]['s1']
                        my_end = repeat_list[i]['e2']
                        start_end = repeat_list[i]['e1']
                        end_start = repeat_list[i]['s2'] 
                        mydiff = a
                        
                    i = i+1
    f.close()
    if mydiff == 1000000:
        return '0_0' #'null'
    return str(my_start)+'_'+ str(my_end)+'_'+ str(start_end)+'_'+ str(end_start) 

def check_pp(contig,start,stop,pp):
    if start>stop:
        t = start
        start =stop
        stop = t
    
    j= 1
    while j <= len(pp):
        if contig == pp[j]['contig']:
          if pp[j]['start']<=start and pp[j]['stop']>=stop:
              return j
        j = j+1
    return 0


def check_phage_word_Start(sjcontig,a,b,c):
     j = 0
     tot = 0
     for i in c:
          start = c[i]['start']
          stop = c[i]['stop']
          if start>stop:
               t = start
               start = stop
               stop = t
          if a <= start and stop <= b  and c[i]['contig'] == sjcontig:
               if c[i]['pp']>0.5:
                    j = j +1
               tot = tot + 1
     #print j,tot,a,b
     if tot < 4 * j:
          return a
     else:
          return b


def check_phage_word_End(sjcontig,a,b,c):
     j = 0
     tot = 0
     for i in c:
          start = c[i]['start']
          stop = c[i]['stop']
          if start>stop:
               t = start
               start = stop
               stop = t
          if a <= start and stop <= b  and c[i]['contig'] == sjcontig:
               if c[i]['pp']>0.5:
                    j = j +1
               tot = tot + 1
     #print j,tot,a,b
     if tot < 4 * j:
          return b
     else:
          return a
               
def final_check_phage_word(sjcontig,a,b,c):
     j = 0
     tot = 0
     for i in c:
          start = c[i]['start']
          stop = c[i]['stop']
          if start>stop:
               t = start
               start = stop
               stop = t
          if a <= start and stop <= b and c[i]['contig'] == sjcontig:
               if c[i]['pp']>0:
                    j = j +1
               tot = tot + 1
     #print j, tot,a,b
     if j>5 and tot< 2*j:
          return str(a)+'_'+str(b)
     else:
          return '0_0'

def clarification_byPhageWord(sjcontig,bef_start,bef_stop,aft_start,aft_stop,myfunction):
     
     if aft_start == 0 and aft_stop == 0:
          return '0_0'
     if bef_start <= aft_start:
          s = check_phage_word_Start(sjcontig,bef_start,aft_start,myfunction)
     else:
          s = check_phage_word_Start(sjcontig,aft_start,bef_start,myfunction)
     
     if bef_stop <= aft_stop:
          e = check_phage_word_End(sjcontig,bef_stop,aft_stop,myfunction)
     else:
          e = check_phage_word_End(sjcontig,aft_stop,bef_stop,myfunction)
          
     se = final_check_phage_word(sjcontig,s,e,myfunction)
     
     #print bef_start,bef_stop,aft_start,aft_stop,se,'!!!!!!!!!!!!!!!!!!!!!!!!!!'
     return se


def fixing_start_end(file,of,org,INSTALLATION_DIR):
    try:
        f = open(file,'r')
    except:
        print 'cant open',file
        return

    #make all predicted pp list
    pp = {}
    i = 1
    flag = 0
    intg = {}
    intg_index = 1
    function = {}
    index_function = 1
    for line in f:
        temp = re.split('\t',line.strip())
        if temp[1]=='function':
            continue
        me = int(temp[7])
        #for interagse
        if float(temp[8]) == 1.5: 
            if int(temp[3]) < int(temp[4]):
                intg[intg_index] = {} 
                intg[intg_index]['start']=int(temp[3])
                intg[intg_index]['stop']=int(temp[4])
            else:
                intg[intg_index] = {}
                intg[intg_index]['start']=int(temp[4])
                intg[intg_index]['stop']=int(temp[3])
            intg_index = intg_index + 1
        ########
        if me == 1 and flag == 0:
            flag = 1
            pp[i]={}
            if int(temp[3]) < int(temp[4]):
                pp[i]['start'] = int(temp[3])
            else:
                pp[i]['start'] = int(temp[4])
            pp[i]['contig'] = temp[2]
        
        if me == 0 and flag == 1: 
            flag = 0
            if int(temp[3])<int(temp[4]):
                pp[i]['stop'] = int(temp[3])
            else:
                pp[i]['stop'] = int(temp[4]) 
            i = i+1

        function[index_function]={}
        function[index_function]['start'] = int(temp[3])
        function[index_function]['stop'] = int(temp[4])
        function[index_function]['pp'] = float(temp[8])
        function[index_function]['contig'] = temp[2]
        index_function = index_function + 1
    f.close()
    #print 'total integrase =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=',len(intg),pp
        
    # find start end for all pp using repeat finder
    dna = read_contig(org)
    for i in pp:
        start = pp[i]['start'] -2000
        if 'stop' in pp[i]:
             stop = pp[i]['stop'] + 2000
        else:
             stop = function[len(function)-1]['stop']

        if stop - start <10000:
            pp[i]['start'] = 0
            pp[i]['stop'] = 0
            continue
        if stop - start >200000:
            pp[i]['start'] = 0
            pp[i]['stop'] = 0
            continue

        repeat_list = find_repeat(dna[pp[i]['contig']][start:stop],start,INSTALLATION_DIR)
        s_e = find_rna(start,stop,repeat_list,org,pp[i]['contig'],intg)
        if s_e != 'null':
            t = re.split('_',s_e)
            
            if s_e == '0_0':
                 s_e1 = clarification_byPhageWord(pp[i]['contig'],pp[i]['start'],pp[i]['stop'],pp[i]['start'],pp[i]['stop'],function)
            else:
                 s_e1 = clarification_byPhageWord(pp[i]['contig'],pp[i]['start'],pp[i]['stop'],float(t[0]),float(t[1]),function)

            t1 = re.split('_',s_e1)
            pp[i]['start'] = float(t1[0])
            pp[i]['stop'] = float(t1[1])
            
            if float(t[0])!= 0 and pp[i]['start'] == float(t[0]) and pp[i]['stop'] == float(t[1]):
                 if int(t[0])<int(t[2]):
                      temps1 = t[0]
                      tempe1 = t[2]
                 else:
                      temps1 = t[2]
                      tempe1 = t[0]
                 if int(t[1])>int(t[3]):
                      temps2 = t[3]
                      tempe2 = t[1]
                 else:
                      temps2 = t[1]
                      tempe2 = t[3]

                 pp[i]['att'] = t[0] +'\t'+t[2]+'\t'+t[3]+'\t'+t[1]+'\t'+dna[pp[i]['contig']][int(temps1)-1:int(tempe1)]+'\t'+dna[pp[i]['contig']][int(temps2)-1:int(tempe2)]
                 
    # fix start end for all pp
    try:
        f = open(file,'r')
        fw = open(of,'w')
    except:
        print 'cant open',file
        return

    for line in f:
        temp = re.split('\t',line.strip())
        if temp[1]=='function':
            fw.write(line)
            continue

        me = check_pp(temp[2],int(temp[3]),int(temp[4]),pp)
        if me == 0:
             fw.write(line.strip()+'\t0'+'\n')
        else:
             fw.write(line.strip()+'\t1')
             if 'att' in pp[me]:
                  fw.write('\t'+pp[me]['att'])
             fw.write('\n')
    f.close()
    fw.close()

############################## added in new version  ##################################################
def fixing_false_negetive(temp_outf, outf, threshold_for_FN, phageWindowSize):
     try:
          f = open(temp_outf,'r')
     except:
          print 'cant open',file
          return

     pp_change = []
 
     fn_end = 0
     fn_start = 0
     count_fn = 0
     print "Threshold for fn is ", threshold_for_FN
     for line in f:
        temp = re.split('\t',line.strip())
        if temp[1]=='function':
            continue

        me = int(temp[9])
        pp = float(temp[8])

        if me > 1:
            print "Checking me: ", me , " and pp: ", pp

        if count_fn == 0:
            if me == 0 and pp >= 1 and count_fn == 0:
                count_fn =  1
                fn_start = int(temp[5])
        else:
            if me >= 1 or pp >= 1:
                count_fn = count_fn + 1 # we are in a run of prophage genes
                fn_end = int(temp[5])
            else:
                # we are not in a run of prophage genes, or in the midst of one, but we'll go to the last
                # gene we've seen
                if (int(temp[5])-fn_start) > 2 * phageWindowSize:
                     if count_fn > threshold_for_FN:
                          while fn_start <= fn_end:
                               pp_change.append(fn_start)
                               fn_start = fn_start + 1
                     count_fn = 0
               
     
     
     f.close()

     try:
          ft = open(temp_outf,'r')
          fw = open(outf,'w')
     except:
          print 'cant open',file
          return

     for line in ft:
          temp = re.split('\t',line.strip())
          if temp[1]=='function':
               fw.write(line)
               continue

          position = int(temp[5])
          if position in pp_change:
               line = line.strip()
               line = line[0:len(line)-1] + '1' + '\n'
          fw.write(line)
     ft.close()
     fw.close()
          
def make_prophage_tbl(input,output):
     try:
          f = open(input,'r')
          fw = open(output,'w')
     except:
          print 'Cant open',input

     pp = {}
     ppindx = 1
     flag = 0
     total_phage_gene = 0
     prev_contig = ''
     
     for line in f:
          if flag == 0:
               flag = 1
               continue
          temp = re.split('\t',line.strip())
          if int(temp[9]) == 1:
               if total_phage_gene == 0 or prev_contig != temp[2]:
                    id_temp = temp[0][:len(temp[0])-temp[0][::-1].find('.')-4]
                    id_temp = id_temp + 'pp.' + str(ppindx)
                    loc_temp = temp[2]+'_'
                    if int(temp[3])>int(temp[4]):
                         loc_temp = loc_temp + temp[4] + '_' 
                    else:
                         loc_temp = loc_temp + temp[3] + '_' 

                    pp[ppindx] = id_temp + '\t' + loc_temp
                    ppindx = ppindx + 1
                    prev_contig = temp[2]
                    total_phage_gene = 0

               total_phage_gene = total_phage_gene+1
               
               if int(temp[3])>int(temp[4]):
                    pp[ppindx-1] = pp[ppindx-1][:len(pp[ppindx-1])-pp[ppindx-1][::-1].find('_')] +temp[3] 
               else:
                    pp[ppindx-1] = pp[ppindx-1][:len(pp[ppindx-1])-pp[ppindx-1][::-1].find('_')] +temp[4] 
          else:
               total_phage_gene = 0

     for i in pp:
          fw.write(pp[i]+'\n')

     f.close()
     fw.close()

                    

################################################################################

def call_start_end_fix(output_dir, organismPath, INSTALLATION_DIR, threshold_for_FN, phageWindowSize):
     fixing_start_end(output_dir+'initial_tbl.txt',output_dir+'prophage_tbl_temp.txt',organismPath,INSTALLATION_DIR)
     cmd2 = "rm "+output_dir+'initial_tbl.txt'
     #os.system(cmd2) # if we comment this line out we keep the initial table file, and then we can just reevaluate
     fixing_false_negetive(output_dir+'prophage_tbl_temp.txt', output_dir+'prophage_tbl.txt', threshold_for_FN, phageWindowSize)
     cmd2 = "rm "+output_dir+'prophage_tbl_temp.txt'
     os.system(cmd2)
     #another outputfile for prophage list
     make_prophage_tbl(output_dir+'prophage_tbl.txt',output_dir+'prophage.tbl')
     
