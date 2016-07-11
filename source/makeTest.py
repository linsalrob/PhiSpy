import re
import math
import string
def read_contig(organism): 
     try:
          #f_dna = open('/home/sajia/Organisms/'+organism+'/contigs','r')
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

def complement(gene):
    x = ''
    gene = gene.upper()
    for i in gene:
        if i == 'A':
            x = x+'T'
        else:
            if i == 'T':
                x = x+'A'
            else:
                if i == 'C':
                    x = x+'G'
                else:
                    if i == 'G':
                        x = x+'C'
                    else: 
                         x = x + i
    return x

def t_skewf(region):
     region = region.upper()
     tot = 0.0
     x = 0
     for i in region:
          if i == 'A':
               tot = tot +1
          if i == 'T':
               x = x + 1
               tot = tot + 1

     if tot == 0:
          return 0
     return float(x)/tot

def a_skewf(region):
     region = region.upper()
     tot = 0.0
     x = 0
     for i in region:
          if i == 'A':
               x = x + 1
               tot = tot +1
          if i == 'T':
               tot = tot + 1
   
     if tot == 0:
          return 0     
     return float(x)/tot

def g_skewf(region):
     region = region.upper()
     tot = 0.0
     x = 0
     for i in region:
          if i == 'G':
               x = x + 1
               tot = tot +1
          if i == 'C':
               tot = tot + 1
     
     if tot == 0:
          return 0
     return float(x)/tot

def c_skewf(region):
     region = region.upper()
     tot = 0.0
     x = 0
     for i in region:
          if i == 'G':
               tot = tot +1
          if i == 'C':
               x = x + 1
               tot = tot + 1
     
     if tot == 0:
          return 0
     return float(x)/tot


def find_mean(all_len):
     
     sum = 0.0
     for i in all_len:
          sum = sum + i
     if len(all_len) == 0:
          return 0
     return float(sum)/len(all_len)

def find_all_median(x):
     all_len = []
     for i in x:
          all_len.append(abs(x[i]['start']-x[i]['stop']))
     return find_median(all_len)

def find_median(all_len):
     n = len(all_len)/2
     all_len.sort()
     if len(all_len) == n*2:
          return (all_len[n]+all_len[n-1])/float(2)
     else:
          return all_len[n]

############################################# shannon ################################

kmers={}

def input_shannon_mer(input):
    global kmers
    try:
        f = open(input,'r')
    except:
        print 'cant open file 1'
        return 0
    
    for line in f:
        line = line.strip()
        kmers[line] = 0
    return 1

def initialize_shannon_mer():
    global kmers
    for i in kmers:
        kmers[i] = 0

def shannon_restricted(seq,mer,flag,total): #if flag == 1 then all ORF done so calculate shannon stuff; 
    	
    global kmers
    #initialize_shannon_mer()

    seq = seq.strip().upper()
    #total = 0;

    pos = 0
    while pos <=len(seq)-mer:
       	substr = seq[pos:pos+mer]
       	pos = pos+mer

       	if substr in kmers:
            kmers[substr] = kmers[substr] + 1
        total = total +1

    if flag == 1:
         if total == 0:
              return 0
         H = 0.0
         found_total = 0.0
         for i in kmers:
              p = float(kmers[i])/total
              if p >0:
                   H = H + p * (math.log(p)/math.log(2))
                   found_total = found_total+ kmers[i]
         
         H = -H
         if H <= 0:
              return 0
         freq_found = found_total/float(total)
         myslope = freq_found/H
         return myslope
    else:
         return total

########################################

def find_avg_length(orf_list):
     x = []
     for i in orf_list:
          x.append(abs(orf_list[i]['start']-orf_list[i]['stop']))
     return find_mean(x)


def find_avg_atgc_skew(orf_list,mycontig,dna):

     a_skew = []
     c_skew = []
     g_skew = []
     t_skew = []

     for i in orf_list:
          start = orf_list[i]['start']
          stop = orf_list[i]['stop']

          if start<stop:
               bact = dna[mycontig][start-1:stop]
          else:
               bact = dna[mycontig][stop-1:start]
               bact = bact[::-1]
               bact = complement(bact)
                             
          if len(bact)<3:
               continue

          a_skew.append(a_skewf(bact))
          t_skew.append(t_skewf(bact))
          g_skew.append(g_skewf(bact))
          c_skew.append(c_skewf(bact))
     a = find_mean(a_skew)
     t = find_mean(t_skew)
     g = find_mean(g_skew)
     c = find_mean(c_skew)
     at = math.fabs(a-t)
     gc = math.fabs(g-c)
     
     return str(at)+'_'+str(gc)

######################################################################################

def make_set_test(organism_list,outfile,window,mer):
  
    for organism in organism_list:

     #for each genome which has prophages
      
         bact_file = organism+'/Features/peg/tbl'
         
         #for host
         try:
              fh = open(bact_file,'r')
         except:
              print 'cant open file',organism
              continue
 
         dna = read_contig(organism)
         try:
              fw = open(outfile,'w')
         except:
              print 'cant open file for write:', outfile
              continue

         fw.write('orf_length_med\tshannon_slope\tat_skew\tgc_skew\tmax_direction\n')
    

         #open host/bact dna file which has a contig
         all_orf_list = {}
         for i in fh:
              temp = re.split('\t',i.strip())
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
         fh.close()
    
         
         for mycontig in all_orf_list:
               
              orf_list = my_sort(all_orf_list[mycontig])
              ######################
              #avg_length = find_avg_length(orf_list)
              all_median = find_all_median(orf_list)
              avg_at_skew =find_avg_atgc_skew(orf_list,mycontig,dna)
              temp = re.split('_',avg_at_skew)
              avg_at_skew = float(temp[0])
              avg_gc_skew = float(temp[1])
              
              #####################
              
              i = 1
              while i<len(orf_list)-window +1:
                   
                   #initialize
                   initialize_shannon_mer()
                   total = 0
                   length = []
                   direction = []
                   a_skew = []
                   t_skew = []
                   g_skew = []
                   c_skew = []
                   pp = 0
                   
                   j = i
                   while j < (i+window):
                        start = orf_list[j]['start']
                        stop = orf_list[j]['stop']

                        if start<stop:
                             bact = dna[mycontig][start-1:stop]
                             direction.append(1) # direction
                        else:
                             bact = dna[mycontig][stop-1:start]
                             bact = bact[::-1]
                             bact = complement(bact)
                             direction.append(-1) # direction
                             
                        if len(bact)<3:
                             print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                             j = j+1
                             continue

                        
                        #at skew
                        a_skew.append(a_skewf(bact))
                        t_skew.append(t_skewf(bact))
                        g_skew.append(g_skewf(bact))
                        c_skew.append(c_skewf(bact))
                   
                        #length
                        length.append(len(bact))
                        
                        #shannon
                        total = shannon_restricted(bact,mer,0,total)
                        
                        j = j+1
                   
                   # write in file for one window
              
                   mylength = find_median(length)-all_median #find_mean(length)
                   fileWriteStr =  str(mylength)+'\t'

                   myshannon = shannon_restricted('',mer,1,total)
                   fileWriteStr = fileWriteStr  + str(myshannon)+'\t'
                   
                   a = find_mean(a_skew)
                   t = find_mean(t_skew)
                   g = find_mean(g_skew)
                   c = find_mean(c_skew)
                   at = math.fabs(a-t)/avg_at_skew
                   gc = math.fabs(g-c)/avg_gc_skew
                   fileWriteStr = fileWriteStr +str(at)+'\t'
                   fileWriteStr = fileWriteStr + str(gc)+'\t'

                   #orf direction
                   orf = []
                   x = 0
                   flag = 0
                   for ii in direction:
                        if ii == 1:
                             if flag == 0:
                                  x = x+1
                             else:
                                  orf.append(x)
                                  x = 1
                                  flag = 0
                        else:
                             if flag == 1:
                                  x = x+1
                             else:
                                  if flag < 1 and x>0:
                                       orf.append(x)
                                  x = 1
                                  flag = 1
                   orf.append(x)
                   orf.sort()
                   if len(orf) == 1:
                        fileWriteStr = fileWriteStr +str(orf[len(orf)-1])+'\n'
                   else:
                        fileWriteStr = fileWriteStr +str(orf[len(orf)-1]+orf[len(orf)-2])+'\n'         
                   fw.write(fileWriteStr)
                   i = i + 1                             
         fw.close()

##################### function call #################################


def call_make_test_set(organismPath,output_dir,INSTALLATION_DIR):
    
     mer = 12
     window = 40
     a = input_shannon_mer(INSTALLATION_DIR+'data/mer_ORF_list.txt')
     

     if a>0:
          make_set_test([organismPath],output_dir+'testSet.txt',window,mer) 
     else:
          print 'error: Shannon does not work'
          return 0
    
     # to check whether the output file has any data for further analysis. For shorter genome (less that 40 genes phiSpy will not work)
     try:
          f = open(output_dir+'testSet.txt','r')
          i = 0
          for line in f:
               i = i + 1
               if (i>1):
                    break
          f.close()
          if i <= 1:
               return 0
     except:
          return 0
     return 1

     
     
