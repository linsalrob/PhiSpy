import re
import math
import string
import pprint
import sys
import array


class ShannonScore:

    def __init__(self, INSTALLATION_DIR):
        # Create a hash of the kmers that points to the index of an array that holds the value
        self._key_to_index = {}
        self._values = array.array('i')
        self.total = 0
        try:
            infile = open(INSTALLATION_DIR + 'data/mer_ORF_list.txt', 'r')
        except:
            sys.exit('ERROR: Cannot open data/mer_ORF_list.txt')
        for line in infile:
            line = line.strip()
            self._values.append(0)
            self._key_to_index[line] = len(self._values) - 1

    def reset(self):
        self.total = 0
        self._values = array.array('i', '\x00' * self._values.itemsize * len(self._values))

    def addValue(self, seq):
        mer = 12
        seq = seq.strip().upper()
        pos = 0
        while (pos <= (len(seq) - mer)):
            substr = seq[pos:pos + mer]
            pos = pos + mer
            if substr in self._key_to_index:
                self._values[self._key_to_index[substr]] += 1
            self.total += 1

    def getSlope(self):
        if self.total == 0:
            return 0
        H = 0.0
        found_total = 0.0
        for i in self._key_to_index:
            p = float(self._values[self._key_to_index[i]]) / self.total
            if (p > 0):
                H = H + p * (math.log(p) / math.log(2))
                found_total = found_total + self._values[self._key_to_index[i]]

        H = -H
        if H <= 0:
            return 0
        freq_found = found_total / float(self.total)
        myslope = freq_found / H
        return myslope


def read_contig(organismPath):
    try:
        f_dna = open(organismPath + '/contigs', 'r')
    except:
        print('cant open contig file ', organismPath)
        return ''

    dna = {}
    seq = ''
    name = ''
    for i in f_dna:
        if i[0] == '>':
            if len(seq) > 10:
                dna[name] = seq
            name = i.strip()
            if ' ' in name:
                temp = re.split(' ', name)
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
            seq = seq + i.strip()

    dna[name] = seq
    f_dna.close()
    return dna


def my_sort(orf_list):
    n = len(orf_list)
    i = 0
    while (i < n):
        j = i + 1
        while (j < n):
            flag = 0
            # direction for both
            if orf_list[i]['start'] < orf_list[i]['stop']:
                dir_i = 1
            else:
                dir_i = -1
            if orf_list[j]['start'] < orf_list[j]['stop']:
                dir_j = 1
            else:
                dir_j = -1

            # check whether swap need or not
            if dir_i == dir_j:
                if orf_list[i]['start'] > orf_list[j]['start']:
                    flag = 1
            else:
                if dir_i == 1:
                    if orf_list[i]['start'] > orf_list[j]['stop']:
                        flag = 1
                else:
                    if orf_list[i]['stop'] > orf_list[j]['start']:
                        flag = 1

            # swap
            if flag == 1:
                temp = orf_list[i]
                orf_list[i] = orf_list[j]
                orf_list[j] = temp
            j += 1
        i += 1
    return orf_list


def complement(gene):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = gene.translate(complements)[::-1]
    return rcseq

def find_all_median(x):
    all_len = []
    for i in x:
        all_len.append((abs(x[i]['start'] - x[i]['stop'])) + 1)
    return find_median(all_len)


def find_median(all_len):
    n = len(all_len) / 2
    all_len.sort()
    if len(all_len) == n * 2:
        return (all_len[n] + all_len[n - 1]) / float(2)
    else:
        return all_len[n]


def find_avg_length(orf_list):
    x = []
    for i in orf_list:
        x.append(abs(orf_list[i]['start'] - orf_list[i]['stop']))
    return sum(x) / len(x)


def find_atgc_skew(seq):
    seq = seq.upper()
    total_at = 0.0
    total_gc = 0.0

    #          A    G     C   T
    scores = {
        'A': [1.0, 0.0, 0.0, 0.0],
        'T': [0.0, 0.0, 0.0, 1.0],
        'G': [0.0, 1.0, 0.0, 0.0],
        'C': [0.0, 0.0, 1.0, 0.0],
        'R': [0.5, 0.5, 0.0, 0.0], #ag
        'Y': [0.0, 0.0, 0.5, 0.5], #ct
        'S': [0.0, 0.5, 0.5, 0.0], #gc
        'W': [0.5, 0.0, 0.0, 0.5], #at
        'K': [0.0, 0.5, 0.0, 0.5], #gt
        'M': [0.5, 0.0, 0.5, 0.0], #ac
        'B': [0.0, 0.3, 0.3, 0.3], #cgt
        'D': [0.3, 0.3, 0.0, 0.3], #agt
        'H': [0.3, 0.0, 0.3, 0.3], #act
        'V': [0.3, 0.3, 0.3, 0.0], #acg
        'N': [0.25, 0.25, 0.25, 0.25], #acgt
    }

    counts = [0, 0, 0, 0]

    for base in seq:
        if base not in scores:
            sys.stderr.write("ERROR: found base " + base + " that is not in the iupac code. Skipped\n")
            continue
        for i,j in enumerate(scores[base]):
            counts[i] += j

    total_at = counts[0] + counts[3]
    total_gc = counts[1] + counts[2]
    if (total_at * total_gc) == 0:
        sys.exit("a total of zero")
    return float(counts[0]) / total_at, float(counts[3]) / total_at, float(counts[1]) / total_gc, float(counts[2]) / total_gc


def find_avg_atgc_skew(orf_list, mycontig, dna):
    a_skew = []
    t_skew = []
    g_skew = []
    c_skew = []
    for i in orf_list:
        start = orf_list[i]['start']
        stop = orf_list[i]['stop']

        if start < stop:
            bact = dna[mycontig][start - 1:stop]
        else:
            bact = dna[mycontig][stop - 1:start]
            bact = bact[::-1]
            bact = complement(bact)

        if len(bact) < 3:
            continue

        xa, xt, xg, xc = find_atgc_skew(bact)
        a_skew.append(xa)
        t_skew.append(xt)
        g_skew.append(xg)
        c_skew.append(xc)

    a = sum(a_skew) / len(a_skew)
    t = sum(t_skew) / len(t_skew)
    g = sum(g_skew) / len(g_skew)
    c = sum(c_skew) / len(c_skew)
    at = math.fabs(a - t)
    gc = math.fabs(g - c)
    return at, gc


######################################################################################

def make_set_test(organismPath, output_dir, window, INSTALLATION_DIR):
    my_shannon_scores = ShannonScore(INSTALLATION_DIR)
    all_orf_list = {}
    try:
        infile = open(organismPath + '/Features/peg/tbl', 'r')
    except:
        sys.exit('ERROR: Cannot open file ' + organismPath + '/Features/peg/tbl')

    dna = read_contig(organismPath)

    # open host/bact dna file which has a contig
    for line in infile:
        temp = re.split('\t', line.strip())
        if ',' in temp[1]:
            ttemp = re.split(',', temp[1])
            temp[1] = ttemp[len(ttemp) - 1]
        temp1 = re.split('_', temp[1])

        contig = temp[1][:temp[1][:temp[1].rfind('_')].rfind('_')]

        start = int(temp1[len(temp1) - 2])
        stop = int(temp1[len(temp1) - 1])

        # save info for sorting orf
        if contig in all_orf_list:
            x = len(all_orf_list[contig])
        else:
            x = 0
            all_orf_list[contig] = {}

        all_orf_list[contig][x] = {}
        all_orf_list[contig][x]['start'] = start
        all_orf_list[contig][x]['stop'] = stop
        all_orf_list[contig][x]['peg'] = temp[0]
    infile.close()

    try:
        outfile = open(output_dir + 'testSet.txt', 'w')
    except:
        sys.exit('ERROR: Cannot open file for writing: testSet.txt')

    outfile.write('orf_length_med\tshannon_slope\tat_skew\tgc_skew\tmax_direction\n')
    for mycontig in all_orf_list:
        orf_list = my_sort(all_orf_list[mycontig])
        ######################
        # avg_length = find_avg_length(orf_list)
        all_median = find_all_median(orf_list)
        avg_at_skew, avg_gc_skew = find_avg_atgc_skew(orf_list, mycontig, dna)
        #####################
        i = 0
        # while i<len(orf_list)-window +1:
        while (i < len(orf_list)):
            # initialize
            my_shannon_scores.reset()
            length = []
            direction = []
            a_skew = []
            t_skew = []
            g_skew = []
            c_skew = []

            for j in range(i - int(window / 2), i + int(window / 2)):
                if ((j < 0) or (j >= len(orf_list))):
                    continue
                start = orf_list[j]['start']
                stop = orf_list[j]['stop']

                if start < stop:
                    bact = dna[mycontig][start - 1:stop]
                    direction.append(1)  # direction
                else:
                    bact = dna[mycontig][stop - 1:start]
                    bact = bact[::-1]
                    bact = complement(bact)
                    direction.append(-1)  # direction

                if len(bact) < 3:
                    print('Short Protein Found!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    j += 1
                    continue

                    # at skew
                xa, xt, xg, xc = find_atgc_skew(bact)
                a_skew.append(xa)
                t_skew.append(xt)
                g_skew.append(xg)
                c_skew.append(xc)

                # length
                length.append(len(bact))

                # shannon
                my_shannon_scores.addValue(bact)
                j += 1

                # write in file for one window
                mylength = find_median(length) - all_median  # find_mean(length)
                fileWriteStr = ''
                fileWriteStr += str(mylength) + '\t'
                fileWriteStr += str(my_shannon_scores.getSlope()) + '\t'
                a = sum(a_skew) / len(a_skew)
                t = sum(t_skew) / len(t_skew)
                c = sum(c_skew) / len(c_skew)
                g = sum(g_skew) / len(g_skew)
                at = math.fabs(a - t) / avg_at_skew if avg_at_skew else 0
                gc = math.fabs(g - c) / avg_gc_skew if avg_gc_skew else 0
                fileWriteStr += str(at) + '\t'
                fileWriteStr += str(gc) + '\t'

                # orf direction
                orf = []
                x = 0
                flag = 0
                for ii in direction:
                    if (ii == 1):
                        if (flag == 0):
                            x += 1
                        else:
                            orf.append(x)
                            x = 1
                            flag = 0
                    else:
                        if (flag == 1):
                            x += 1
                        else:
                            if (flag < 1 and x > 0):
                                orf.append(x)
                            x = 1
                            flag = 1
                orf.append(x)
                orf.sort()
                if len(orf) == 1:
                    fileWriteStr += str(orf[len(orf) - 1]) + '\n'
                else:
                    fileWriteStr += str(orf[len(orf) - 1] + orf[len(orf) - 2]) + '\n'
                outfile.write(fileWriteStr)
                i += 1
    outfile.close()


##################### function call #################################


def call_make_test_set(organismPath, output_dir, INSTALLATION_DIR):
    window = 40

    make_set_test(organismPath, output_dir, window, INSTALLATION_DIR)

    # Check whether the output file has data. For shorter genomes (less that 40 genes) phiSpy will not work)
    num_lines = sum(1 for line in open(output_dir + 'testSet.txt', 'r'))
    if (num_lines > 0):
        return 1
    else:
        return 0
