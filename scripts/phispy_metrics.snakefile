########################################################################
#                                                                      #
#                   PhiSpy Metrics Snakefile                           #
#                                                                      #
# Snakemake code to compare all the possible metrics against each      #
# other, and to run all combinations of those metrics.                 #
#                                                                      #
# You should be able to run this with:                                 #
#   snakemake -s phispy_metrics.snakefile -j 10                        #
#                                                                      #
# But I really recommend running this across a cluster                 #
# as it will perform a lot of pairwise combinations!                   #
#                                                                      #
# Rob Edwards, 2020                                                    #
#                                                                      #
########################################################################

# These are the only two things you need to change:

# Where is PhiSpy installed?
phispydir = "/home3/redwards/GitHubs/PhiSpy/"

# What is the output directory
outdir = 'phispy_metrics'


import os
import sys
from itertools import combinations 

metrics = ['none', 'pg0', 'orf_length_med', 'shannon_slope', 'at_skew', 'gc_skew', 'max_direction', 'phmms']

def all_metrics():
    all_metrics = []
    for i in range(1,len(metrics)+1):
        for p in combinations(metrics, i):
            all_metrics.append("-".join(p))
    return all_metrics


def get_params(wildcards):
    r = []
    for m in wildcards.metric.split("-"):
        if m == 'phmms':
            r.append(' --phmms /home3/redwards/VOGs/VOGs.hmm ')
        elif m == "pg0":
            r.append(' --phage_genes 0 ')
        else:
            r.append(f" --metrics {m} ")
    return " ".join(r)


def base_filename(wildcards):
    m = "-".join(wildcards.metric)
    return f"{wildcards.sample}.{m}"


GENOMES, = glob_wildcards(os.path.join(phispydir, "test_genbank_files", '{genome}.gb.gz'))
METRICS = all_metrics()

rule all:
    input:
        expand(os.path.join(outdir, "{genome}.phispy.{metric}", "tptn.txt"),
               genome=GENOMES, metric=METRICS)

rule run_phispy:
    input:
        g = os.path.join(phispydir, "test_genbank_files", "{genome}.gb.gz")
    params:
        o = os.path.join(outdir, "{genome}.phispy.{metric}"),
        m = get_params
    benchmark:
        os.path.join(outdir, "benchmarks", "{genome}.phispy.{metric}")
    output:
        temporary(os.path.join(outdir, "{genome}.phispy.{metric}", "bacteria.fasta")),
        temporary(os.path.join(outdir, "{genome}.phispy.{metric}", "bacteria.gbk")),
        temporary(os.path.join(outdir, "{genome}.phispy.{metric}", "phage.fasta")),
        os.path.join(outdir, "{genome}.phispy.{metric}", "phage.gbk"),
        os.path.join(outdir, "{genome}.phispy.{metric}", "phispy.log"),
    shell:
        """
        PhiSpy.py {params.m} --color -o {params.o} --output_choice 256 {input.g}
        """

rule count_tp_tn:
    input:
        gen = os.path.join(phispydir, "test_genbank_files", "{genome}.gb.gz"),
        phg = os.path.join(outdir, "{genome}.phispy.{metric}", "phage.gbk")
    output:
        tp = os.path.join(outdir, "{genome}.phispy.{metric}", "tptn.txt")
    params:
        comp = os.path.join(phispydir, "scripts", "compare_predictions_to_phages.py")
    shell:
        """
        python3 {params.comp} -t {input.gen} -p {input.phg} > {output.tp}
        """
