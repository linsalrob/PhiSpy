import os
import re
import math
import sys
from argparse import Namespace

from .writers import write_phage_and_bact, write_gff3, write_prophage_tbl
from .writers import write_prophage_tsv, write_genbank, write_prophage_coordinates
from .formatting import message

import PhiSpyRepeatFinder


def find_repeat(fn, st, ppno, extra_dna, output_dir):
    if len(fn) == 0:
        message("Len sequence is 0 so ignoring\n", "RED", 'stderr')
        return {}

    rep = {}
    index = 0

    # with open(os.path.join(output_dir, "repeat_finding"), 'a') as rptout:
    #     rptout.write(f">pp{ppno} {st}\n{fn}\n")

    try:
        # set the False parameter to True to enable debugging of repeat finder
        repeats = PhiSpyRepeatFinder.repeatFinder(fn, 3, ppno, False)
    except Exception as e:
        message(f"There was an error running repeatfinder for {fn}:{e}\n", "RED", 'stderr')
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
    message("Checking prophages we might have found\n", "GREEN", 'stderr')
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
        if pp[i]['num genes'] >= self.number:
            temppp[j] = pp[i]
            j += 1
            prophagesummary.append([pp[i]['contig'], pp[i]['start'], pp[i]['stop'], pp[i]['num genes'], "Kept"])
        else:
            prophagesummary.append([pp[i]['contig'], pp[i]['start'], pp[i]['stop'], pp[i]['num genes'],
                                    "Dropped. Not enough genes"])
    # print a list of all prophages and the number of genes
    if prophagesummary:
        sys.stderr.write('Potential prophages (sorted highest to lowest)\n')
        sys.stderr.write('Contig\tStart\tStop\tNumber of potential genes\tStatus\n')
        for p in sorted(prophagesummary, key=lambda x: x[3], reverse=True):
            sys.stderr.write("\t".join(map(str, p)) + "\n")
    pp = temppp
    sys.stderr.write("\n")
    # End filtering
    # find start end for all pp using repeat finder
    dna = {entry.id: str(entry.seq) for entry in self.record}
    extra_dna = 2000
    for i in pp:
        message(f"PROPHAGE: {i} Contig: {pp[i]['contig']} Start: {pp[i]['start']} Stop: {pp[i]['stop']}", "PINK", 'stderr')
        start = pp[i]['start'] - extra_dna
        if start < 1:
            start = 1
        if 'stop' in pp[i]:
            stop = pp[i]['stop'] + extra_dna
        else:
            stop = genome[len(genome) - 1]['stop']
        if stop > len(dna[pp[i]['contig']]):
            stop = len(dna[pp[i]['contig']])
        if stop - start > 200000:
            message(f"Not checking repeats for pp {i} because it is too big: {stop - start} bp", "PINK", 'stderr')
            continue
        repeat_list = find_repeat(dna[pp[i]['contig']][start:stop], start, i, extra_dna, self.output_dir)
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
                        message(msg, "YELLOW", 'stderr')
                    if len(attRseq) == 0:
                        msg = f"The attL sequence had no length from {int(bestrep['s2']) - 1} to "
                        msg += f"{int(bestrep['e2']) - 1} on contig {pp[i]['contig']} "
                        msg += f"(length: {len(dna[pp[i]['contig']])} )"
                        message(msg, "YELLOW", 'stderr')
                    pp[i]['att'] = [
                        bestrep['s1'],
                        bestrep['e1'],
                        bestrep['s2'],
                        bestrep['e2'],
                        attLseq,
                        attRseq,
                        "Longest Repeat flanking phage and within " + str(extra_dna) + " bp"
                    ]
                    pp[i]['atts'] = "\t".join(map(str, pp[i]['att']))
                    if samelenrep > 1:
                        message(f"There were {samelenrep} repeats with the same length as the best. " +
                                         "One chosen somewhat randomly!\n", "YELLOW", 'stderr')
    # fix start end for all pp
    try:
        outfile = open(os.path.join(self.output_dir, self.file_prefix + 'prophage_information.tsv'), 'w')
    except IOError as e:
        message(f"There was an error opening the files: {e}\n", "RED", 'stderr')
        sys.exit(-1)

    outfile.write("identifier\tfunction\tcontig\tstart\tstop\tposition\trank\tmy_status\tpp\tFinal_status\t")
    outfile.write("start of attL\tend of attL\tstart of attR\tend of attR\tsequence of attL\tsequence of attR\t")
    outfile.write("Reason for att site\n")

    for this_pp in self.initial_tbl:
        me = check_pp(this_pp[2], int(this_pp[3]), int(this_pp[4]), pp)
        if me == 0:
            outfile.write("\t".join(map(str, this_pp + [0])) + '\n')
        else:
            outfile.write("\t".join(map(str, this_pp)) + '\t' + str(me) + '\n')

    """
    now we need to decide which files to keep
    It is based on this code:
        Code | File
        --- | ---
        1 | prophage_coordinates.tsv 
        2 | GenBank format output 
        4 | prophage and bacterial sequences  
        8 | prophage_information.tsv  
        16 | prophage.tsv  
        32 | GFF3 format  
        64 | prophage.tbl
    As explained in the README. 
    """

    oc = self.output_choice
    if oc >= 64:
        # write the prophage location table
        write_prophage_tbl(self.output_dir, pp, self.file_prefix)
        oc -= 64
    if oc >= 32:
        # write the prophage in GFF3 format
        write_gff3(self.output_dir, pp, self.file_prefix)
        oc -= 32
    if oc >= 16:
        # write a tsv file of this data
        write_prophage_tsv(self.output_dir, pp, self.file_prefix)
        oc -= 16
    if oc < 8:
        # we want to KEEP the prophage_information.tsv that we created earlier
        # if oc >= 8, or if we get to this point and oc is 7 or less we want to
        # delete it.
        os.remove(os.path.join(self.output_dir, self.file_prefix + 'prophage_information.tsv'))
    if oc == 8:
        oc = 0
    if oc >= 4:
        # separate out the bacteria and phage as fasta files
        write_phage_and_bact(self.output_dir, pp, dna, self.file_prefix)
        oc -= 4
    if oc >= 2:
        # update input GenBank file and incorporate prophage regions
        write_genbank(self.infile, self.record, self.output_dir, pp, self.file_prefix)
    if oc >= 1:
        # print the prophage coordinates:
        write_prophage_coordinates(self.output_dir, pp, self.file_prefix)





