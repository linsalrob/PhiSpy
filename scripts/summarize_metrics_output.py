"""
Summarize the output from the snakefile that compares all metrics to each other.

This will make a table with "codes" for each combination of metrics, as well
as output the results table.

"""

import os
import sys
import re
from roblib import bcolors
from itertools import combinations
__author__ = 'Rob Edwards'


def all_metrics():
    metrics = ['none', 'noannotation', 'pg0', 'orf_length_med', 'shannon_slope', 'at_skew', 'gc_skew', 'max_direction', 'phmms']
    all_metrics = []
    for i in range(1,len(metrics)+1):
        for p in combinations(metrics, i):
            all_metrics.append("-".join(p))
    return all_metrics


def all_genomes():
    phispydir = "/home3/redwards/GitHubs/PhiSpy/test_genbank_files"
    genomes = []
    for f in os.listdir(phispydir):
        if f.endswith('.gb.gz'):
            genomes.append(f.replace('.gb.gz', ''))
    return genomes


def get_data(tptnf):

    if not os.path.exists(tptnf):
        sys.stderr.write(f"FATAL: {tptnf} does not exist\n")
        sys.exit(-1)

    sys.stderr.write(f"Parsing {tptnf}\n")
    
    are = re.compile('Accuracy:\s+([\d\.]+)\s')
    pre = re.compile('Precision:\s+([\d\.]+)\s')
    rre = re.compile('Recall:\s+([\d\.]+)\s')
    sre = re.compile('Specificity:\s+([\d\.]+)\s')
    fre = re.compile('f1 score:\s+([\d\.]+)\s')

    res = ["", "", "", "", ""]
    with open(tptnf, 'r') as f:
        for l in f:
            if m:= are.search(l):
                res[0] = float(m.group(1))
            elif m:= pre.search(l):
                res[1] = float(m.group(1))
            elif m:= rre.search(l):
                res[2] = float(m.group(1))
            elif m:= sre.search(l):
                res[3] = float(m.group(1))
            elif m:= fre.search(l):
                res[4] = float(m.group(1))

    return res




if __name__ == "__main__":
    #ms = ['xxxx']
    #ms += all_metrics()
    ms = all_metrics()

    metric_printed = set()

    with open("metric_codes.tsv", 'w') as mc, open("phispy_metrics_tptn.tsv", "w") as out:
            for i, j in enumerate(ms):
                if  j == 'xxxx':
                    msn = "No metrics"
                    pmd = "phispy_no_metrics"
                else:
                    msn = " ".join(j.split("-"))
                    pmd = "phispy_metrics"

                if i not in metric_printed:
                    mc.write(f"{i}\t{msn}\n")
                    metric_printed.add(i)
                for g in all_genomes():
                    res = [g, i]
                    if not os.path.exists(os.path.join(pmd, f"{g}.phispy.{j}", "tptn.txt")):
                        continue
                    res += get_data(os.path.join(pmd, f"{g}.phispy.{j}", "tptn.txt"))

                    out.write("\t".join(map(str, res)))
                    out.write("\n")


