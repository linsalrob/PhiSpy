import os
import sys
import math
import sys
import pkg_resources
from argparse import Namespace

from .log_and_message import log_and_message
from .errors import NoBasesCounted

class ShannonScore:
    def __init__(self, kmers_type):
        # Create a hash of the kmers that points to the index of an array that holds the value
        self._kmers = {}
        self._kmers_phage = []
        self._kmers_all = []
        self._kmers_type = kmers_type
        kmers_file = 'data/phage_kmers_' + self._kmers_type + '_wohost.txt'
        if not pkg_resources.resource_exists:
            log_and_message(f"ERROR: Kmers file {kmers_file} not found", "RED", stderr=True, quiet=self.quiet)
            sys.exit(13)

        for line in pkg_resources.resource_stream('PhiSpyModules', kmers_file):
            line = line.decode().strip()
            self._kmers[line] = ''

    def reset(self):
        self._kmers_phage = []
        self._kmers_all = []

    def addValue(self, seq):
        mer = 12
        pos = 0
        self._kmers_phage.append([])
        self._kmers_all.append(0)
        kmers = self.kmerize_orf(seq, mer, self._kmers_type)
        for kmer in kmers:
            self._kmers_all[-1] += 1
            try:
                self._kmers[kmer]
                self._kmers_phage[-1].append(kmer)
            except KeyError:
                continue

    def getSlope(self, start, stop):
        total = sum(self._kmers_all[start : stop])
        found_total = sum([len(x) for x in self._kmers_phage[start : stop]])
        if total == 0:
            return 0
        H = 0.0
        window_kmers = {}
        for kl in self._kmers_phage[start : stop]:
            for k in kl:
                try:
                    window_kmers[k] += 1.0
                except KeyError:
                    window_kmers[k] = 1.0
        for i in window_kmers.values():
            p = i/total
            H = H + p * (math.log(p)/math.log(2))
        if H >= 0:
            return 0
        freq_found = found_total / float(total)
        myslope = - freq_found / H
        return myslope

    def kmerize_orf(self, orf, k, t):

        kmers = []
        if t == 'simple':
            stop = len(orf) - (len(orf) % k)
            for i in range(0, stop, k):
                kmers.append(orf[i : i + k])
        elif t == 'all':
            for j in range(0, k):
                stop = len(orf) - ((len(orf) - j) % k)
                for i in range(j, stop, k):
                    kmers.append(orf[i : i + k])
        elif t == 'codon':
            for j in range(0, k, 3):
                stop = len(orf) - ((len(orf) - j) % k)
                for i in range(j, stop, k):
                    kmers.append(orf[i : i + k])

        return kmers

def my_sort(orf_list):
    n = len(orf_list)
    i = 0
    while i < n:
        j = i + 1
        while j < n:
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

def find_all_median(orf_list):
    all_len = []
    for i in orf_list:
        all_len.append((abs(i['start'] - i['stop'])) + 1)
    return find_median(all_len)

def find_median(all_len):
    n = int(round(len(all_len) / 2))
    all_len.sort()
    if len(all_len) == n * 2:
        return (all_len[n] + all_len[n - 1]) / float(2)
    else:
        return all_len[n]

def find_atgc_skew(seq):
    total_at = 0.0
    total_gc = 0.0
    a = 0
    t = 0
    c = 0
    g = 0
    for base in seq:
        if base == 'A':
            a += 1
            total_at += 1
        elif base == 'T':
            t += 1
            total_at += 1
        elif base == 'G':
            g += 1
            total_gc += 1
        elif base == 'C':
            c += 1
            total_gc += 1
        elif base == 'R':
            a += 0.5
            total_at += 0.5
            g += 0.5
            total_gc += 0.5
        elif base == 'Y':
            c += 0.5
            total_gc += 0.5
            t += 0.5
            total_at += 0.5
        elif base == 'S':
            g += 0.5
            total_gc += 0.5
            c += 0.5
            total_gc += 0.5
        elif base == 'W':
            a += 0.5
            total_at += 0.5
            t += 0.5
            total_at += 0.5
        elif base == 'K':
            g += 0.5
            total_gc += 0.5
            t += 0.5
            total_at += 0.5
        elif base == 'M':
            c += 0.5
            total_gc += 0.5
            a += 0.5
            total_at += 0.5
        elif base == 'B':
            c += 0.3
            total_gc += 0.3
            g += 0.3
            total_gc += 0.3
            t += 0.3
            total_at += 0.3
        elif base == 'D':
            a += 0.3
            total_at += 0.3
            g += 0.3
            total_gc += 0.3
            t += 0.3
            total_at += 0.3
        elif base == 'H':
            a += 0.3
            total_at += 0.3
            c += 0.3
            total_gc += 0.3
            t += 0.3
            total_at += 0.3
        elif base == 'V':
            a += 0.3
            total_at += 0.3
            c += 0.3
            total_gc += 0.3
            g += 0.3
            total_gc += 0.3
        elif base == 'N':
            a += 0.25
            total_at += 0.25
            t += 0.25
            total_at += 0.25
            g += 0.25
            total_gc += 0.25
            c += 0.25
            total_gc += 0.25
        else:
            log_and_message(f"A non nucleotide base ({base}) was found in {seq}", "RED", stderr=True)
            sys.exit(25)
    if(total_at * total_gc) == 0:
        raise NoBasesCounted("a total of zero total_at*total_gc")
    return float(a)/total_at, float(t)/total_at, float(g)/total_gc, float(c)/total_gc

def find_avg_atgc_skew(orf_list, mycontig, dna):
    a_skew = []
    t_skew = []
    g_skew = []
    c_skew = []
    for i in orf_list:
        start = i['start']
        stop = i['stop']
        try:
            if start < stop:
                bact = dna[mycontig][start - 1:stop]
                xa, xt, xg, xc = find_atgc_skew(bact)
            else:
                bact = dna[mycontig][stop - 1:start]
                xt, xa, xc, xg = find_atgc_skew(bact)
        except NoBasesCounted as e:
            msg = e.message + "\n"
            msg += f"No bases were counted for orf {i} from {start} to {stop}\n"
            msg += "This error is usually thrown with an exceptionally short ORF that is only a "
            msg += " few bases. You should check this ORF and confirm it is real!\n"
            log_and_message(msg, "RED", stderr=True)
            sys.exit(26)
        if len(bact) < 3:
            continue
        a_skew.append(xa)
        t_skew.append(xt)
        g_skew.append(xg)
        c_skew.append(xc)
    return a_skew, t_skew, g_skew, c_skew

def reverse_complement(seq):

    rcd = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
           'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M',
           'S': 'S', 'W': 'W',
           'H': 'D', 'D': 'H',
           'B': 'V', 'V': 'B',
           'N': 'N'}
    rcseq = ''
    for base in seq[::-1]:
        try:
            rcseq += rcd[base]
        except KeyError:
            rcseq += base

    return rcseq

def measure_features(**kwargs):
    """
    This is a combination of the old make_test_set // make_train_set
    but now returns an array of the data and so reduces redundancy
    :param kwargs: The keyworded arguments object
    :return: an array of arrays. Each element contains median orf length, shannon slope,
            at_skew, gc_skew, max_direction, phmms
    """

    self = Namespace(**kwargs)
    data = []  # the array of arrays we will return
    my_shannon_scores = ShannonScore(self.kmers_type)
    all_orf_list = {}
    dna = {}
    window = self.window_size
    for entry in self.record:
        dna[entry.id] = str(entry.seq)
        for feature in entry.get_features('CDS'):
            orf_list = all_orf_list.get(entry.id, [])
            is_phage = int(feature.qualifiers['is_phage'][0]) if 'is_phage' in feature.qualifiers else 0
            orf_list.append(
                   {'start' : feature.start,
                    'stop'  : feature.stop,
                    'phmm'  : sum([-math.log10(x) if x != 0 else 500 for x in feature.phmm])/100,
                    'peg'   : 'peg',
                    'is_phage': is_phage  # note that we include this here, even if we are generating test data
                   }
            )
            all_orf_list[entry.id] = orf_list

    if len(all_orf_list) == 0:
        log_and_message(f"There were no ORFs predicted in {self.infile}. Please annotate the genome and try again",
                        "RED", stderr=True, quiet=self.quiet)
        sys.exit(40)

    for mycontig in all_orf_list:
        orf_list = all_orf_list[mycontig]
        if not orf_list:
            # an empty list of orfs
            log_and_message("No ORFs were found in {mycontig}", "YELLOW", stderr=True, quiet=self.quiet)
            continue
        all_median = find_all_median(orf_list)
        lengths = []
        directions = []
        phmms = []
        for i in orf_list:
            lengths.append(abs(i['start'] - i['stop']) + 1) # find_all_median can be deleted now
            directions.append(1 if i['start'] < i['stop'] else -1)
            phmms.append(i['phmm'])
            if i['start'] < i['stop']:
                seq = dna[mycontig][i['start'] - 1 : i['stop']].upper()
            else:
                seq = reverse_complement(dna[mycontig][i['stop'] - 1 : i['start']].upper())
            my_shannon_scores.addValue(seq)

        ga_skew, gt_skew, gg_skew, gc_skew = find_avg_atgc_skew(orf_list, mycontig, dna)
        a = sum(ga_skew) / len(ga_skew)
        t = sum(gt_skew) / len(gt_skew)
        g = sum(gg_skew) / len(gg_skew)
        c = sum(gc_skew) / len(gc_skew)
        avg_at_skew, avg_gc_skew = math.fabs(a - t), math.fabs(g - c)

        for i, orf_data in enumerate(orf_list):
            # orf_length_med\tshannon_slope\tat_skew\tgc_skew\tmax_direction\tphmms
            this_orf = []  # the data for this particular orf that gets added to test_data
            # initialize
            j_start = i - int(window / 2)
            j_stop = i + int(window / 2)
            if j_start < 0:
                j_start = 0
            elif j_stop >= len(orf_list):
                j_stop = len(orf_list)
            # at and gc skews
            ja_skew = ga_skew[j_start:j_stop]
            jt_skew = gt_skew[j_start:j_stop]
            jc_skew = gc_skew[j_start:j_stop]
            jg_skew = gg_skew[j_start:j_stop]
            ja = sum(ja_skew) / len(ja_skew)
            jt = sum(jt_skew) / len(jt_skew)
            jc = sum(jc_skew) / len(jc_skew)
            jg = sum(jg_skew) / len(jg_skew)
            jat = math.fabs(ja - jt) / avg_at_skew if avg_at_skew else 0
            jgc = math.fabs(jg - jc) / avg_gc_skew if avg_gc_skew else 0
            my_length = find_median(lengths[j_start:j_stop]) - all_median
            # orf direction
            orf = []
            x = 0
            flag = 0
            for ii in directions[j_start:j_stop]:
                if ii == 1:
                    if flag == 0:
                        x += 1
                    else:
                        orf.append(x)
                        x = 1
                        flag = 0
                else:
                    if flag == 1:
                        x += 1
                    else:
                        if flag < 1 and x > 0:
                            orf.append(x)
                        x = 1
                        flag = 1
            orf.append(x)
            orf.sort()
            this_orf.append(my_length)
            if self.expand_slope:
                s = my_shannon_scores.getSlope(j_start, j_stop)
                this_orf.append(s * s)
            else:
                this_orf.append(my_shannon_scores.getSlope(j_start, j_stop))
            this_orf.append(jat)
            this_orf.append(jgc)
            this_orf.append(
                orf[len(orf) - 1] if len(orf) == 1 else orf[len(orf) - 1] + orf[len(orf) - 2]
            )
            this_orf.append(sum(phmms[j_start:j_stop]))
            assert(len(this_orf) == 6)  # confirm I added everything!
            if 'making_training_set' in self:
                this_orf.append('1' if orf_data['is_phage'] else '0')
            data.append(this_orf)
        my_shannon_scores.reset()
    return data

def make_set_train(**kwargs):
    kwargs['making_training_set'] = True
    self = Namespace(**kwargs)

    if self.output_dir:
        of = os.path.join(self.output_dir, self.make_training_data)
    else:
        of = self.make_training_data

    with open(of,'w') as outfile:
        outfile.write('orf_length_med\tshannon_slope\tat_skew\tgc_skew\tmax_direction\tphmms\tstatus\n')
        for d in measure_features(**kwargs):
            outfile.write("\t".join(map(str, d)) + "\n")
