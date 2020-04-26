#!/usr/bin/python3
__author__ = 'Przemek Decewicz'

from argparse import Namespace
from Bio import SeqIO
from glob import glob
from os import path, remove
from subprocess import call

def read_genbank(infile, infile_type, target_type):
    """
    This function reads the GenBank file and extracts protein and/or nucleotide sequences in a manner that will allow easier mapping of results against it's genome.
    """

    if infile_type == 'genbank':
        genome = { 'contigs': list(SeqIO.parse(infile, 'genbank')),
                   'proteins': [],
                   'info': {}}
    elif infile_type == 'contigs':
        genome = { 'contigs': list(SeqIO.parse(infile, 'fasta')),
                   'proteins': [],
                   'info': {}}

    if target_type == 'proteins':
        gpid_cnt = 1

        for r in genome['contigs']:
            rid = r.id
            genome['info'][rid] = {'length': len(r.seq),
                                   'description': r.description,
                                   'features': {}}

            for f in r.features:
                try:
                    genome['info'][rid]['features'][f.type] += 1
                except KeyError:
                    genome['info'][rid]['features'][f.type] = 1

                if f.type == 'CDS':
                    gpid = 'pid%07d' % gpid_cnt
                    gpid_cnt += 1

                    try:
                        pid = f.qualifiers['protein_id'][0]
                    except KeyError:
                        pid = gpid
                        f.qualifiers['protein_id'] = [gpid]

                    try:
                        product = f.qualifiers['product'][0]
                    except KeyError:
                        product = 'unknown product of protein %s' % pid

                    try:
                        aa = f.qualifiers['translation'][0]
                    except KeyError:
                        if f.location.strand == 1:
                            aa = r.seq[f.location.start:f.location.end].translate(table = 11, to_stop = True)
                        else:
                            aa = r.seq[f.location.start:f.location.end].reverse_complement().translate(table = 11, to_stop = True)

                    header = [rid, gpid, pid, str(f.location)] # product was ommited
                    genome['proteins'].append('>%s\n%s' % ('|'.join(header), aa))

        print('  Read %i contigs.' % len(genome['contigs']))
        print('  Read %i proteins.' % len(genome['proteins']))

    elif target_type == 'contigs':
        print('  Read %i contigs.' % len(genome['contigs']))

    return genome


def write_search_input(genome, target_type, outdir):
    """
    Writes protein FASTA file.
    """

    if target_type == 'proteins':
        outfile = path.join(outdir, 'prots.fasta')
        with open(outfile, 'w') as outf:
            outf.write('\n'.join(genome['proteins']))
    elif target_type == 'contigs':
        outfile = path.join(outdir, 'cont.fasta')
        with open(outfile, 'w') as outf:
            outf.write('\n'.join(genome['contigs']))

        target = path.join(outdir, 'cont_prot.fasta')
        cmd = ['transeq', '-sequence', outfile, '-outseq', target, '-clean', '-table', '11', '-frame', '6']
        print(' '.join(cmd))
        call(cmd)
        outfile = target

    return outfile


def search_target(searchfile, target_type, phmms, outdir, threads):
    """
    Performs the search of reference databases against target sequences (searchfile) using hmmsearch.
    """

    results = {'hsout': {}, # unused
               'hstbl': {}}
    print('  ----> %s' % path.basename(phmms))

    hstbl = path.join(outdir, path.basename(phmms) + '.tbl')
    hsout = path.join(outdir, path.basename(phmms) + '.out')

    cmd = ['hmmsearch', '--cpu', str(threads), '-E', '1e-10', '--domE', '1e-5', '--noali', '--tblout', hstbl, '-o', hsout, phmms, searchfile]
    # print(' '.join(cmd))
    call(cmd)

    # Read hmmsearch results
    # print('  Reading output file.')
    # results['hsout'].update(parse_hmmsearch_output(hsout, target_type))

    print('  Reading hmmsearch tbl output files.')
    results['hstbl'].update(parse_hmmersearch_tbl_output(hstbl))

    return results


def update_genbank(records, results, outfile, color):
    """
    Goes through each record and its proteins for which HMMs were assigned and adds them into annotation.
    """
    # results {sid: {pid: [hmms]}}
    for r in records:
        for f in r.features:
            if f.type != 'CDS': continue
            pid = f.qualifiers['protein_id'][0]
            try:
                phmms = results['hstbl'][r.id][pid]
            except KeyError:
                phmms = []
            for phmm in phmms:
                try:
                    f.qualifiers['phmm'].append(phmm)
                except KeyError:
                    f.qualifiers['phmm'] = [phmm]
                if color:
                    f.qualifiers['color'] = ['6']


    SeqIO.write(records, outfile, 'genbank')

    return records


def parse_hmmersearch_tbl_output(infile):
    """
    This function parses hmmsearch tbl output to further map HMM profiles identifiers into GenBank file.
    """

    results = {}    # {sid: {pid: [hmms]}}
    with open(infile) as inf:
        for line in inf:
            if line.startswith('#'): continue
            line = line.split()
            sid, gpid, pid = line[0].split('|')[:3]
            hit = f'{line[2]}:{line[4]}'
            try:
                results[sid][pid].append(hit)
            except KeyError:
                try:
                    results[sid][pid] = [hit]
                except KeyError:
                    results[sid] = {pid: [hit]}
    return results


def search_phmms(**kwargs):

    """
    Wraps phmms search.
    """

    self = Namespace(**kwargs)
    input_type = 'genbank' # alternatively contigs
    target = 'proteins' # alternatively contigs
    outfile = path.join(self.output_dir, path.basename(self.infile))


    if not self.skip_search:
        # 1 - Read GenBank file and extract required information
        print('  Reading input file: %s' % self.infile)
        genome = read_genbank(self.infile, input_type, target)


        # 2 - Prepare DNA sequence for profile searching
        print('  Writing target data.')
        searchfile = write_search_input(genome, target, self.output_dir)


        # 3 - Perform searches against prepared reference datafiles
        print('  Performing search against reference datasets.')
        results = search_target(searchfile, target, self.phmms, self.output_dir, self.threads)


        # 4 - Print the summary
        print('  Search summary.')
        for sid, pids in results['hstbl'].items():
            print(f'   - {sid}: found HMMs for {len(pids)}/{len(genome["proteins"])} proteins.')


        # 5 - Update GenBank file
        print('  Updating GenBank file.')
        genome['contigs'] = update_genbank(genome['contigs'], results, outfile, self.color)

    if not self.keep:
        if path.isfile(path.join(self.output_dir, 'prots.fasta')): remove(path.join(self.output_dir, 'prots.fasta'))
        for hmm in glob(path.join(self.output_dir, '*.hmm.*')):
            remove(hmm)

    print('  Done!')

    return outfile
