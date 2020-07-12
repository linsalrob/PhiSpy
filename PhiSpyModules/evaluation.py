import os
import re
import math
import sys
from argparse import Namespace

from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature

from .log_and_message import log_and_message
import PhiSpyModules.version as version

import PhiSpyRepeatFinder


def find_repeat(self, contig, fn, st, ppno, extra_dna):
    """
    Find repeats in the DNA sequence
    :param self: The data object
    :param contig: the name of the contig we are searching on
    :param fn: the nuclotide sequence to search
    :param st: the start to find repeats at
    :param ppno: the prophage number
    :param extra_dna: the extra dna that flanks the sequence
    :return: a list of repeat regions
    """

    if len(fn) == 0:
        log_and_message("Len sequence is 0 so ignoring\n", c="RED", stderr=True, loglevel="WARNING")
        return {}

    rep = {}
    index = 0

    # with open(os.path.join(output_dir, "repeat_finding"), 'a') as rptout:
    #     rptout.write(f">pp{ppno} {st}\n{fn}\n")

    try:
        # set the False parameter to True to enable debugging of repeat finder
        repeats = PhiSpyRepeatFinder.repeatFinder(fn, 3, self.min_repeat_len, ppno, False)
    except Exception as e:
        log_and_message(f"There was an error running repeatfinder for {fn}:{e}\n", c="RED", stderr=True,
                        loglevel="WARNING")
        return {}

    for r in repeats:
        if (r['first_start'] < (3 * extra_dna)) and (r['second_start'] > (len(fn) - (3 * extra_dna))):
            # check that start is always less than end
            # This always causes an off by one error, so we have to increment our ends
            if r['first_end'] < r['first_start']:
                [r['first_start'], r['first_end']] = [r['first_end'] + 1, r['first_start'] + 1]
            if r['second_end'] < r['second_start']:
                [r['second_start'], r['second_end']] = [r['second_end'] + 1, r['second_start'] + 1]

            rep[index] = {}
            rep[index]['s1'] = r['first_start'] + st
            rep[index]['e1'] = r['first_end'] + st
            rep[index]['s2'] = r['second_start'] + st
            rep[index]['e2'] = r['second_end'] + st
            if self.include_all_repeats:
                replen = max(rep[index]['e1'] - rep[index]['s1'], rep[index]['e2'] - rep[index]['s2'])
                r1loc = FeatureLocation(rep[index]['s1'], rep[index]['e1'], strand=+1)
                r2loc = FeatureLocation(rep[index]['s2'], rep[index]['e2'], strand=+1)
                rptloc = CompoundLocation([r1loc, r2loc])
                rptsf = SeqFeature(rptloc,type="repeat_region",
                                   qualifiers={'note':f"{replen}bp repeat identified by PhiSpy v{version.__version__}"})
                self.record.get_entry(contig).features.append(rptsf)
            index += 1

    return rep


def check_intg(prophage_sta, prophage_sto, rep, integ, con):
    for m in integ:
        if integ[m]['contig'] != con:
            continue
        if integ[m]['start'] < (prophage_sta + (prophage_sto - prophage_sta) / 2):
            l = integ[m]['start']
        else:
            l = integ[m]['stop']
        if (l - rep['s1'] <= 500 and l - rep['s1'] > 0) or (l - rep['e1'] <= 500 and l - rep['e1'] > 0):
            return 1
        if (l - rep['s2'] <= 500 and l - rep['s2'] > 0) or (l - rep['e2'] <= 500 and l - rep['e2'] > 0):
            return 1
    return 0


def find_smallest(a, b):
    mm = 1000000
    for i in a:
        for j in b:
            if math.fabs(i - j) < mm:
                mm = math.fabs(i - j)
    return mm


def find_rna(prophage_start, prophage_stop, repeat_list, record, cont, integers):
    my_start = 1000000
    start_end = 0
    end_start = 0
    my_end = 1000000
    mydiff = 1000000
    for feature in record.get_entry(cont).get_features('tRNA'):
        i = 0
        while i < len(repeat_list):
            a = find_smallest([feature.start, feature.stop], 
                              [repeat_list[i]['s1'], repeat_list[i]['s2'],
                               repeat_list[i]['e1'], repeat_list[i]['e2']])
            if check_intg(prophage_start, prophage_stop, repeat_list[i], integers, cont) == 1:
                if (math.fabs(repeat_list[i]['s1'] - repeat_list[i]['e2']) > math.fabs(
                        my_start - my_end)) or mydiff == 1000000:
                    my_start = repeat_list[i]['s1']
                    my_end = repeat_list[i]['e2']
                    start_end = repeat_list[i]['e1']
                    end_start = repeat_list[i]['s2']
                    mydiff = a
            if a <= 500 and a < mydiff:
                my_start = repeat_list[i]['s1']
                my_end = repeat_list[i]['e2']
                start_end = repeat_list[i]['e1']
                end_start = repeat_list[i]['s2']
                mydiff = a
            i += 1
    if mydiff == 1000000:
        return '0_0'  # 'null'
    return str(my_start) + '_' + str(my_end)+'_' + str(start_end) + '_' + str(end_start)


def check_pp(contig, start, stop, pp):
    if start > stop:
        (start, stop) = (stop, start)

    for j in pp:
        if contig == pp[j]['contig']:
            if pp[j]['start'] <= start and pp[j]['stop'] >= stop:
                return j
    return 0


def check_phage_word_start(sjcontig, a, b, c):
    j = 0
    tot = 0
    for i in c:
        start = c[i]['start']
        stop = c[i]['stop']
        if start > stop:
            t = start
            start = stop
            stop = t
        if a <= start and stop <= b and c[i]['contig'] == sjcontig:
            if c[i]['pp'] > 0.5:
                j += 1
            tot += 1
    if tot < 4 * j:
        return a
    else:
        return b


def count_phage_word_ends(sjcontig, a, b, c, maxpp):
    j = 0
    tot = 0
    for i in c:
        start = c[i]['start']
        stop = c[i]['stop']
        if start > stop:
            t = start
            start = stop
            stop = t
        if a <= start and stop <= b and c[i]['contig'] == sjcontig:
            if c[i]['pp'] > maxpp:
                j = j + 1
            tot = tot + 1
    return j, tot


def check_phage_word_end(sjcontig, a, b, c):
    j, tot = count_phage_word_ends(sjcontig, a, b, c, 0.5)
    if tot < 4 * j:
        return b
    else:
        return a


def final_check_phage_word(sjcontig, a, b, c):
    j, tot = count_phage_word_ends(sjcontig, a, b, c, 0)
    if j > 5 and tot < 2 * j:
        return str(a) + '_' + str(b)
    else:
        return '0_0'


def clarification_by_phage_word(sjcontig, bef_start, bef_stop, aft_start, aft_stop, genome):
    if aft_start == 0 and aft_stop == 0:
        return '0_0'
    if bef_start <= aft_start:
        s = check_phage_word_start(sjcontig, bef_start, aft_start, genome)
    else:
        s = check_phage_word_start(sjcontig, aft_start, bef_start, genome)

    if bef_stop <= aft_stop:
        e = check_phage_word_end(sjcontig, bef_stop, aft_stop, genome)
    else:
        e = check_phage_word_end(sjcontig, aft_stop, bef_stop, genome)
    se = final_check_phage_word(sjcontig, s, e, genome)
    return se


def fixing_start_end(**kwargs):
    self = Namespace(**kwargs)
    # make all predicted pp list
    log_and_message("Checking prophages we might have found", c="GREEN", stderr=True, quiet=self.quiet)
    pp = {}
    i = 0
    flag = 0
    intg = {}
    intg_index = 1
    genome = {}
    index = 1
    distance_from_last_prophage = 1000
    last_pp = ["" for x in range(9)]

    for this_pp in self.initial_tbl:
        distance_from_last_prophage += 1

        # Find location of all prophage regions
        # columns 7 and 8 are my_status and pp
        if float(this_pp[7]) >= 1 or float(this_pp[8]) >= 1:
            # This is a prophage region. Is it a new prophage?
            new_prophage = False
            # check the sequences are on the same contig. If not, definitely a new prophage
            if this_pp[2] != last_pp[2]:
                new_prophage = True
                distance_from_last_prophage = 1000
            if flag == 0 and distance_from_last_prophage > self.nonprophage_genegaps:
                # we need at least 10 non phage genes betweeen prophages.
                new_prophage = True
            if new_prophage:
                i += 1
                pp[i] = {}
                pp[i]['contig'] = this_pp[2]
                pp[i]['start'] = min(int(this_pp[3]), int(this_pp[4]))
                pp[i]['stop'] = max(int(this_pp[3]), int(this_pp[4]))
                pp[i]['num genes'] = 1
                pp[i]['annotated_as_pp'] = False
                pp[i]['phage_genes'] = 0
                flag = 1
                last_pp = this_pp
            else:
                pp[i]['stop'] = max(pp[i]['stop'], int(this_pp[3]), int(this_pp[4]))
                pp[i]['num genes'] += 1
            distance_from_last_prophage = 0
        else:
            flag = 0
            # Find location of integrases
        if float(this_pp[8]) == 1.5:
            intg[intg_index] = {}
            intg[intg_index]['start'] = min(int(this_pp[3]), int(this_pp[4]))
            intg[intg_index]['stop'] = max(int(this_pp[3]), int(this_pp[4]))
            intg[intg_index]['contig'] = str(this_pp[2])
            intg_index += 1

        if this_pp[8] >=1:
            pp[i]['phage_genes'] += 1 # number of genes annotated as phage genes

        genome[index] = {}
        genome[index]['start'] = int(this_pp[3])
        genome[index]['stop'] = int(this_pp[4])
        genome[index]['pp'] = float(this_pp[8])
        genome[index]['contig'] = this_pp[2]
        genome[index]['function'] = this_pp[1]
        genome[index]['rank'] = float(this_pp[6])
        index += 1

    #######################################################################################
    #                                                                                     #
    # Filter the potential prophages based on how many potential prophage genes are       #
    # in the window, and tell us which ones were dropped.                                 #
    #                                                                                     #
    #######################################################################################
    temppp = {}
    j = 1
    prophagesummary = []
    for i in pp:
        if pp[i]['num genes'] >= self.number and pp[i]['phage_genes'] >= self.phage_genes:
            temppp[j] = pp[i]
            j += 1
            prophagesummary.append([pp[i]['contig'], pp[i]['start'], pp[i]['stop'], pp[i]['num genes'], "Kept"])
        elif pp[i]['num genes'] >= self.number and pp[i]['phage_genes'] > 0:
            prophagesummary.append([pp[i]['contig'], pp[i]['start'], pp[i]['stop'], pp[i]['num genes'],
                                    f"Dropped. Only {pp[i]['phage_genes']} gene(s) were identified as phage genes"])
        elif pp[i]['num genes'] >= self.number:
            prophagesummary.append([pp[i]['contig'], pp[i]['start'], pp[i]['stop'], pp[i]['num genes'],
                                    "Dropped. No genes were identified as phage genes"])
        else:
            prophagesummary.append([pp[i]['contig'], pp[i]['start'], pp[i]['stop'], pp[i]['num genes'],
                                    "Dropped. Not enough genes"])
    # print a list of all prophages and the number of genes
    if prophagesummary:
        log_and_message('Potential prophages (sorted highest to lowest)', stderr=True, quiet=self.quiet)
        log_and_message('Contig\tStart\tStop\tNumber of potential genes\tStatus', stderr=True, quiet=self.quiet)
        for p in sorted(prophagesummary, key=lambda x: (x[3], x[0], x[1]), reverse=True):
            log_and_message("\t".join(map(str, p)), stderr=True, quiet=self.quiet)
    pp = temppp
    # End filtering
    # find start end for all pp using repeat finder
    dna = {entry.id: str(entry.seq) for entry in self.record}
    for i in pp:
        log_and_message(f"PROPHAGE: {i} Contig: {pp[i]['contig']} Start: {pp[i]['start']} Stop: {pp[i]['stop']}",
                        c="PINK", stderr=True, quiet=self.quiet)
        start = pp[i]['start'] - self.extra_dna
        if start < 1:
            start = 1
        if 'stop' in pp[i]:
            stop = pp[i]['stop'] + self.extra_dna
        else:
            stop = genome[len(genome) - 1]['stop']
        if stop > len(dna[pp[i]['contig']]):
            stop = len(dna[pp[i]['contig']])
        if stop - start > 200000:
            log_and_message(f"Not checking repeats for pp {i} because it is too big: {stop - start} bp",
                            c="PINK", stderr=True, quiet=self.quiet)
            continue
        repeat_list = find_repeat(self, pp[i]['contig'], dna[pp[i]['contig']][start:stop], start, i, self.extra_dna)
        s_e = find_rna(start, stop, repeat_list, self.record, pp[i]['contig'], intg)
        if s_e != 'null':
            t = re.split('_', s_e)
            if s_e == '0_0':
                s_e1 = clarification_by_phage_word(pp[i]['contig'], pp[i]['start'], pp[i]['stop'], pp[i]['start'],
                                                   pp[i]['stop'], genome)
            else:
                s_e1 = clarification_by_phage_word(pp[i]['contig'], pp[i]['start'], pp[i]['stop'], float(t[0]),
                                                   float(t[1]), genome)
            t1 = re.split('_', s_e1)
            # we only accept this if it increases the phage length
            if float(t1[0]) > 0 and float(t1[0]) < pp[i]['start'] and float(t1[1]) > 0 and float(t1[1]) > pp[i]['stop']:
                pp[i]['start'] = float(t1[0])
                pp[i]['stop'] = float(t1[1])

            if (float(t[0]) != 0) and (pp[i]['start'] == float(t[0])) and (pp[i]['stop'] == float(t[1])):
                # 1 if (float(t[0])!= 0) and (float(t1[0]) == float(t[0])) and (float(t1[1]) == float(t[1])):
                temps1 = min(t[0], t[2])
                tempe1 = max(t[0], t[2])
                temps2 = min(t[1], t[3])
                tempe2 = max(t[1], t[3])
                pp[i]['att'] = [
                    t[0],
                    t[2], 
                    t[3],
                    t[1], 
                    dna[pp[i]['contig']][int(temps1) - 1:int(tempe1)],
                    dna[pp[i]['contig']][int(temps2) - 1:int(tempe2)],
                    'Repeat exactly at the end'
                ]
                # string representation of the above
                pp[i]['atts'] = "\t".join(map(str, pp[i]['att']))
            else:
                # this approach will just append the longest repeat
                longestrep = 0
                bestrep = None
                samelenrep = 0
                for idx in repeat_list:
                    lengthrep = math.fabs(repeat_list[idx]['e1'] - repeat_list[idx]['s1'])
                    if lengthrep > longestrep and lengthrep < 150:
                        longestrep = lengthrep
                        bestrep = repeat_list[idx]
                        samelenrep = 1
                    elif lengthrep == longestrep:
                        samelenrep += 1
                if bestrep:
                    attLseq = dna[pp[i]['contig']][int(bestrep['s1']) - 1:int(bestrep['e1']) - 1]
                    attRseq = dna[pp[i]['contig']][int(bestrep['s2']) - 1:int(bestrep['e2']) - 1]
                    if len(attLseq) == 0:
                        msg = f"The attL sequence had no length from {int(bestrep['s1']) - 1} to "
                        msg += f"{int(bestrep['e1']) - 1} on contig {pp[i]['contig']} "
                        msg += f"(length: {len(dna[pp[i]['contig']])} )"
                        log_and_message(msg, c="YELLOW", stderr=True, quiet=self.quiet)
                    if len(attRseq) == 0:
                        msg = f"The attL sequence had no length from {int(bestrep['s2']) - 1} to "
                        msg += f"{int(bestrep['e2']) - 1} on contig {pp[i]['contig']} "
                        msg += f"(length: {len(dna[pp[i]['contig']])} )"
                        log_and_message(msg, c="YELLOW", stderr=True, quiet=self.quiet)
                    pp[i]['att'] = [
                        bestrep['s1'],
                        bestrep['e1'],
                        bestrep['s2'],
                        bestrep['e2'],
                        attLseq,
                        attRseq,
                        f"Longest Repeat flanking phage and within {self.extra_dna} bp"
                    ]
                    pp[i]['atts'] = "\t".join(map(str, pp[i]['att']))
                    if samelenrep > 1:
                        msg=f"There were {samelenrep} repeats with the same length as the best. One chosen somewhat randomly!"
                        log_and_message(msg, c="YELLOW", stderr=True, quiet=self.quiet)
        if 'atts' not in pp[i]:
            pp[i]['atts'] = "No potential att site found"
    return pp





