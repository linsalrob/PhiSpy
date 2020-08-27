#!/usr/bin/python3
__author__ = 'Przemek Decewicz'

import gzip
import sys
import pkg_resources
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from Bio import SeqIO
from compare_predictions_to_phages import genbank_seqio
from glob import glob
from os import makedirs, path
from PhiSpyModules import log_and_message
from re import sub
from numpy import arange
from subprocess import call

DATA_DIR = pkg_resources.resource_filename('PhiSpyModules', 'data')
TEST_DIR = pkg_resources.resource_filename('PhiSpyModules', 'data/testSets')
if not path.isdir(TEST_DIR): makedirs(TEST_DIR)

def read_genbank(gbkfile, full_analysis=False):
    """
    Parses GenBank file's CDSs and groups them into host or phage groups based on the '/is_phage' qualifier.
    :param gbkfile: path to GenBank file
    :param full_analysis: wether to read just the taxonomy (False) or all CDSs
    :return infile_data: dictionary with taxonomy, bact/phage CDSs lists
    """

    log_and_message(f"Reading {gbkfile}.", stderr=True)

    infile_data = {
        'taxonomy': [],
        'bact_cds': [],
        'phage_cds': []
    }

    tax_present = False
    for record in genbank_seqio(gbkfile):
        # check records until taxonomy information is provided
        if not tax_present:
            if len(record.annotations['taxonomy']) == 0:
                infile_data['taxonomy'] = ['Bacteria']
            else:
                infile_data['taxonomy'] = check_taxa_names(record.annotations['taxonomy'])
                tax_present = True

        # get bacteria and phage CDSs nucleotide sequences to make kmers
        if full_analysis:
            for f in record.features:
                if f.type == 'CDS':
                    dna = str(f.extract(record).seq)
                    try:
                        status = f.qualifiers['is_phage'][0]
                        if status == '1':
                            infile_data['phage_cds'].append(dna)
                        else:
                            infile_data['bact_cds'].append(dna)
                    except KeyError:
                        infile_data['bact_cds'].append(dna)
    if not tax_present:
        log_and_message(f"- WARNING! Taxonomy was missing!!! Assigning to Bacteria.", c="RED", stderr=True)
    if full_analysis:
        log_and_message(f"- Bact CDSs: {len(infile_data['bact_cds'])}", stderr=True)
        log_and_message(f"- Phage CDSs: {len(infile_data['phage_cds'])}", stderr=True)

    return infile_data


def check_taxa_names(taxonomy):
    """
    Checks wether taxa names contain any illegal characters.
    :param taxonomy: list of taxonomic levels names
    :return taxonomy: corrected Taxonomy
    """

    ILLEGAL_CHARACTERS = {
        " ": "-", # e.g. Mycobacteium tuberculosis comples
        "/": "-and-", # e.g. Rhizobium/Agrobacterium group
    }

    for i in range(len(taxonomy)):
        for ic, lc in ILLEGAL_CHARACTERS.items():
            taxonomy[i] = sub(ic, lc, taxonomy[i])

    return taxonomy


def read_groups(groups_file, training_data):
    """
    Reads tab-delimited input file with the path to input file to use for training in the first column
    or just the file name (if --indir provided, it will be later added to this file path)
    and the name of the group to put it afterwards while creating trainingGenome_list.txt.
    :param groups_file: path to the groups file
    :param training_data: dictionary with groups and infiles
    :return training_data: updated dictionary with genomes and infiles
    """

    trained_groups = len(training_data['groups'])

    with open(groups_file) as inf:
        for line in inf:
            line = line.strip().split('\t')
            training_data['genomes'].add(line[0])
            try:
                training_data['groups'][line[1]].add(line[0])
            except KeyError:
                training_data['groups'][line[1]] = set([line[0]])

    log_and_message(f"Read {len(training_data['groups']) - trained_groups} new groups.", stderr=True)

    return training_data


def read_kmers(kmerfile):
    """
    Simply reads input file with kmers, ony by line.
    :param kmerfile: path to file with kmers to read
    :return kmers: set of kmers
    """

    with gzip.open(kmerfile, 'rt') as inf:
        kmers = {x.strip() for x in inf}

    return kmers


def read_training_genomes_list(training_data):
    """
    Read trainingGenome_list.txt in PhiSpy's data directory and check currently set
    groups and already trained genomes.
    :param training_data: dictionary storing groups, genomes and taxonomy information of processed genomes
    :return training_data: updated dictionary
    """

    with open(path.join(DATA_DIR,'trainingGenome_list.txt')) as infile:
        infile.readline() # skip the first line with Generic test set
        for line in infile: #pkg_resources.resource_stream('PhiSpyModules', 'data/trainingGenome_list.txt'):
            n, group, genomes, genomes_number = line.strip().split('\t') #decode().strip().split('\t')
            group = group.rsplit('.txt', 1)[0][9:] #remove .txt and trainSet_
            genomes = set(genomes.split(';')) if ';' in genomes else set([genomes])
            training_data['groups'][group] = genomes
            training_data['genomes'].update(genomes)

    log_and_message(f"Read {len(training_data['genomes'])} genomes assigned to {len(training_data['groups'])} groups.", stderr=True)

    return training_data


def kmerize_orf(orf, k, t):
    """
    Creates kmers of certain size and type from provided sequence.
    :param orf: nucleotide sequence
    :param k: size of a kmer
    :param t: type of a kmer, i.e. 'simple', 'all' or 'codon'
    :return kmers: list of kmers
    """

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


def prepare_taxa_groups(training_data):
    """
    Reads all input files / genomes and creates a groups file based on all genomes taxonomy.
    :param training_data: currently trained genomes and set groups
    :return training_data
    """

    current_groups = len(training_data['groups'])

    for file_name in training_data['genomes']:
        for tax in training_data['taxonomy'][file_name]:
            try:
                training_data['groups'][tax].add(file_name)
            except KeyError:
                training_data['groups'][tax] = set([file_name])

    log_and_message(f"Created {len(training_data['groups']) - current_groups} new groups based on taxonomy.",stderr=True)

    return training_data


def write_kmers_file(file_name, bact_orfs_list, phage_orfs_list, kmer_size, kmers_type):
    """
    Calculates host/phage kmers from input file and writes phage-specific kmers to file.
    """

    MIN_RATIO = 1.0
    # host_kmers = set()
    # phage_kmers = set()
    kmers_dict = {} # {kmer: [# in host, # in phage]}
    kmers_total_count = 0
    kmers_host_unique_count = 0
    kmers_phage_unique_count = 0
    kmers_ratios = {}
    kmers_ratios_stats = {round(x, 2): 0 for x in arange(0, 1.01, 0.01)}
    kmers_ratios_file = path.join(TEST_DIR, f'{file_name}.kmers_ratios.txt')
    kmers_phage_file = path.join(TEST_DIR, f'{file_name}.kmers_phage.gz')
    kmers_host_file = path.join(TEST_DIR, f'{file_name}.kmers_host.gz')

    # kmrize CDSs
    for orf in bact_orfs_list:
        for kmer in kmerize_orf(orf, kmer_size, kmers_type):
            try:
                kmers_dict[kmer][0] += 1
            except KeyError:
                kmers_dict[kmer] = [1, 0]
                kmers_host_unique_count += 1
            kmers_total_count += 1
        # host_kmers.update(kmerize_orf(i, kmer_size, kmers_type))
    for orf in phage_orfs_list:
        for kmer in kmerize_orf(orf, kmer_size, kmers_type):
            try:
                kmers_dict[kmer][1] += 1
            except KeyError:
                kmers_dict[kmer] = [0, 1]
                kmers_phage_unique_count += 1
            kmers_total_count += 1
        # phage_kmers.update(kmerize_orf(i, kmer_size, kmers_type))

    log_and_message(f"Analyzed {kmers_total_count} kmers.", stderr=True)
    log_and_message(f"Identified {len(kmers_dict)} unique kmers.", stderr=True)
    log_and_message(f"- Bact unique {kmers_host_unique_count} ({kmers_host_unique_count / len(kmers_dict) * 100:.2f}%).", stderr=True)
    log_and_message(f"- Phage unique {kmers_phage_unique_count} ({kmers_phage_unique_count / len(kmers_dict) * 100:.2f}%).", stderr=True)

    ##################
    # the below part could be simplified if ratios will not be considered
    # it just slightly extends the calculations
    ##################

    # calculate ratios
    for kmer, freqs in kmers_dict.items():
        ratio = round(freqs[1] / sum(freqs), 2)
        kmers_ratios[(ratio, kmer)] = freqs
        kmers_ratios_stats[ratio] += 1
    del kmers_dict

    # write ratios_stats
    log_and_message(f"Writing kmers ratios stats.", stderr=True)
    with open(kmers_ratios_file, 'w') as outf:
        outf.write("Ratio\tNumber of kmers\tPerc of such kmers\tCumulative perc\n")
        tot = 0
        tot_perc = 0
        for x in reversed(arange(0, 1.01, 0.01)):
            x = round(x, 2)
            tot += kmers_ratios_stats[x]
            perc = kmers_ratios_stats[x]/len(kmers_ratios) * 100
            tot_perc = tot / len(kmers_ratios) * 100
            outf.write(f"{x}\t{kmers_ratios_stats[x]}\t{perc:.3f}%\t{tot_perc:.3f}%\n")

    # write unique phage kmers
    log_and_message(f"Writing kmers into phage and host kmers files.", stderr=True)
    cnt = 0
    with gzip.open(kmers_phage_file, 'wt') as out_phage:
        with gzip.open(kmers_host_file, 'wt') as out_host:
            for ratio, kmer in kmers_ratios.keys():
                if ratio >= MIN_RATIO:
                    cnt += 1
                    out_phage.write(f"{kmer}\n")
                else:
                    out_host.write(f"{kmer}\n")
    log_and_message(f"- Wrote {cnt} kmers with ratios >= {MIN_RATIO}.", stderr=True)

    return


def write_training_sets(training_data):
    """
    Writes trainSets based on provided groups.
    :param training_data: dictionary with genomes and groups used for making trainging sets
    """

    training_data['groups']['genericAll'] = training_data['genomes']
    for i, group in enumerate(sorted(training_data['groups']), 1):
        header = False
        log_and_message(f"[{i}/{len(training_data['groups'])}] Writing trainSet for {group} ({len(training_data['groups'][group])} genome(s)).", c="PINK", stderr=True)
        with open(path.join(DATA_DIR, f"trainSet_{group}.txt"), 'w') as outf:
            for genome in training_data['groups'][group]:
                with open(path.join(TEST_DIR, f"{genome}.testSet")) as inf:
                    if header: inf.readline()
                    outf.write(inf.read())
                    header = True
    training_data['groups'].pop('genericAll')


def write_training_genome_list(training_data):
    """
    Writes trainingGenome_list.txt file with currently available training sets.
    :param training_data: dictionary with genomes and groups used for making trainging sets
    """

    with open(path.join(DATA_DIR, 'trainingGenome_list.txt'), 'w') as outf:
        outf.write(f"0\ttestSet_genericAll.txt\tGeneric Test Set\t{len(training_data['genomes'])}\n")
        for i, group in enumerate(sorted(training_data['groups']), 1):
            outf.write(f"{i}\t")
            outf.write(f"trainSet_{group}.txt\t")
            outf.write(f"{';'.join(training_data['groups'][group])}\t")
            outf.write(f"{len(training_data['groups'][group])}\n")


def print_groups(groups):
    """
    Prints groups with their genomes.
    :param groups: a dictionary with group name as key and genomes list as value
    """

    for group, genomes in sorted(groups.items()):
        gg = '\n- '.join(genomes)
        log_and_message(f"{group}", c="PINK", stderr=True)
        log_and_message(f"- {gg}", stderr=True)


def get_file_path(file_name, infiles, indir):
    """
    Return path to input file.
    :param file_name: GenBank file names
    :param infiles: list of input files from user's input directory
    :param indir: user's input directory
    :return infile: path to file of interest
    """

    if file_name in infiles:
        # if indicated within directory provided by user
        log_and_message(f"File in user's input directory.", stderr=True)
        infile = path.join(indir, file_name)
    else:
        log_and_message(f"File {file_name} not present in user's input directory. If retraing with PhiSpy's default training genomes consider using its test_genbank_files directory.\nIf you want to use just your own data, run the script with --absolute_retrain flag.\nQuiting.", c="RED", stderr=True)
        exit(2)

    return infile


def main():
    args = ArgumentParser(prog = 'make_training_sets.py',
                          description = 'Automates making new or extending current PhiSpy\'s training sets. By default these will be created in PhiSpyModules/data directory so keep that in mind preparing groups file. ',
                          epilog = 'Example usage:\npython3 scripts/make_training_sets.py -d test_genbank_files -g test_genbank_files/groups.txt --retrain --use_taxonomy --phmms pVOGs.hmm --threads 4',
                          formatter_class = RawDescriptionHelpFormatter)

    args.add_argument('-d', '--indir',
                      type = str,
                      help = 'Path to input directory with GenBank file(s) for training. This path will be added to file names in groups file.')

    args.add_argument('-g', '--groups',
                      type = str,
                      help = 'Path to file two tab-delimited columns: file name and group name. If not provided each file will have its own training set.')

    args.add_argument('--use_taxonomy',
                      action = 'store_true',
                      help = 'If set, taxonomy from input files will be used to update or create new groups. This is performed after reading groups file.')

    args.add_argument('-k', '--kmer_size',
                      type = int,
                      help = 'The size of required kmers. For codon approach use multiplicity of 3. [Default: 12]',
                      default = 12)

    args.add_argument('-t', '--kmers_type',
                      type = str,
                      help = 'Approach for creating kmers. Options are: simple (just slicing the sequence from the first position), all (all possible kmers), codon (all possible kmers made with step of 3 nts to get kmers corresponding translated aas). [Default: all]',
                      default = 'all')

    args.add_argument('--phmms',
                      type = str,
                      help = 'Phage HMM profile database (like pVOGs) will be mapped against the genome of interest and used as additional feature to identify prophages.')

    args.add_argument('--threads',
                      type = str,
                      help = 'Number of threads to use while searching with phmms.',
                      default = '4')

    args.add_argument('--retrain',
                      action = 'store_true',
                      help = 'Set if any of reference files previously used for training has changed, e.g. prophage protein features indication was modified.')

    args.add_argument('--absolute_retrain',
                      action = 'store_true',
                      help = 'If set, ignores trainingGenome_list file and PhiSpy\'s default reference genomes. This option allows to train PhiSpy with files provided just by the user.')

    if len(sys.argv[1:]) == 0:
        args.print_usage()
        args.exit()

    try:
        args = args.parse_args()
    except:
        args.exit()

    if not args.indir:
        log_and_message(f"You have to provide input directory --indir.", c="RED",
                        stderr=True, stdout=False)
        sys.exit(2)


    log_and_message("Reading input directory", c="GREEN", stderr=True)
    infiles = glob(path.join(args.indir, r'*.gb'))
    infiles += glob(path.join(args.indir, r'*.gb[kf]'))
    infiles += glob(path.join(args.indir, r'*.gbff'))
    infiles += glob(path.join(args.indir, r'*.gb.gz'))
    infiles += glob(path.join(args.indir, r'*.gb[kf].gz'))
    infiles += glob(path.join(args.indir, r'*.gbff.gz'))
    infiles = {path.basename(infile) for infile in infiles}
    log_and_message(f"Read {len(infiles)} GenBank files from input directory.", stderr=True)


    # read currently available genomes - either by reading trainingGenome_list.txt
    log_and_message("Checking currently available training sets.", c="GREEN", stderr=True, stdout=False)
    training_data = {
        'groups': {},
        'genomes': set(),
        'taxonomy': {}
    }
    not_trained = set()
    if args.absolute_retrain:
        log_and_message(f"Ignoring PhiSpy's trainingGenome_list.txt file and default test GenBank files.Files provided with --indir and/or --groups will overwrite current reference sets.", stderr=True)
        args.retrain = True
    elif pkg_resources.resource_exists('PhiSpyModules', 'data/trainingGenome_list.txt'):
        training_data = read_training_genomes_list(training_data)
    else:
        log_and_message(f"trainingGenome_list.txt is missing.", c="RED", stderr=True)
    log_and_message(f"{len(training_data['groups'])} groups based on trainingGenome_list file:", c="GREEN", stderr=True)
    old_training_groups = training_data['groups'].copy()
    print_groups(training_data['groups'])


    log_and_message(f"Checking which genomes are considered as new.", c="GREEN", stderr=True)
    for infile in infiles:
        file_name = path.basename(infile)
        if file_name not in training_data['genomes']:
            not_trained.add(infile)
            training_data['genomes'].add(infile)
            log_and_message(f"- {file_name}", c="YELLOW", stderr=True)
    log_and_message(f"In total there are {len(not_trained)} new genomes.", stderr=True)


    # check what new groups were requested
    if args.groups:
        log_and_message(f"Reading provided groups file.", c="GREEN", stderr=True)
        training_data = read_groups(args.groups, training_data)
        log_and_message(f"{len(training_data['groups'])} currently considered groups:", c="GREEN", stderr=True)


    # make kmers files or read taxonomy information from input file if needed
    # Comment:
    # In general, all files need to be retrained if:
    # (i) there's at least one new reference file or
    # (ii) annotation of CDSs with is_phage qualifier has changed in any file
    # The change in CDSs marked as phage CDSs triggers the change of phage-specific
    # kmers set and therefore changed the Shannon Score statistics.
    log_and_message(f"Checking which genomes need to be read/retrained.", c="GREEN", stderr=True)
    for i, file_name in enumerate(sorted(training_data['genomes']), 1):
        full_analysis = False
        if file_name in not_trained:
            log_and_message(f"This file has not been used for training yet: {file_name}", c="RED", stderr=True)
            full_analysis = True
        elif args.retrain:
            log_and_message(f"Retraining file upon user's request.", c="PINK", stderr=True)
            full_analysis = True
        elif not pkg_resources.resource_exists('PhiSpyModules', f"data/testSets/{file_name}.kmers_phage.gz") or \
             not pkg_resources.resource_exists('PhiSpyModules', f"data/testSets/{file_name}.kmers_host.gz") or \
             not pkg_resources.resource_exists('PhiSpyModules', f"data/testSets/{file_name}.testSet"):
            log_and_message(f"Training files missing for: {file_name}", c="RED", stderr=True)
            full_analysis = True
        # full_analysis = file_name in not_trained or args.retrain
        if full_analysis or args.use_taxonomy:
            log_and_message(f"[{i}/{len(training_data['genomes'])}] Reading {file_name}.", c="YELLOW", stderr=True)
            infile = get_file_path(file_name, infiles, args.indir)
            infile_data = read_genbank(infile, full_analysis)
            training_data['taxonomy'][file_name] = infile_data['taxonomy']
            if full_analysis:
                write_kmers_file(file_name, infile_data['bact_cds'], infile_data['phage_cds'], args.kmer_size, args.kmers_type)
            else:
                log_and_message("No need to write kmers files.", stderr=True)
        else:
            log_and_message(f"[{i}/{len(training_data['genomes'])}] Skipping. {file_name} already analyzed.", c="YELLOW", stderr=True)


    # use taxonomy information to create/update groups
    if args.use_taxonomy:
        log_and_message(f"Using taxonomy from input files to create new or update current groups.", c="GREEN", stderr=True)
        training_data = prepare_taxa_groups(training_data)
        log_and_message(f"{len(training_data['groups'])} considered groups, including taxonomy-based ones:", c="GREEN", stderr=True)
        print_groups(training_data['groups'])


    # Check wether there's a point to go further
    if len(not_trained) == 0 and training_data['groups'] == old_training_groups:
        log_and_message(f"There are 0 new genomes and groups. Quiting.", c="GREEN", stderr=True)
        exit(1)


    log_and_message("Making phage unique kmers file from all considered genomes: phage_kmers_all_wohost.txt", c="GREEN", stderr=True)
    phage_kmers_all_wohost_file = path.join(DATA_DIR, 'phage_kmers_all_wohost.txt')
    phage_kmers = set()
    for i, file_name in enumerate(sorted(training_data['genomes']), 1):
        log_and_message(f"[{i}/{len(training_data['genomes'])}] Reading kmers from {file_name}.", c="YELLOW", stderr=True)
        kmers_file = path.join(TEST_DIR, file_name)
        phage_kmers.update(read_kmers(kmers_file + '.kmers_phage.gz'))
        phage_kmers.difference_update(read_kmers(kmers_file + '.kmers_host.gz'))
    log_and_message(f"Writing {len(phage_kmers)} into {phage_kmers_all_wohost_file}.", stderr=True)
    with open(phage_kmers_all_wohost_file, 'w') as outf:
        outf.write("\n".join(phage_kmers))


    retrain_all = False
    if len(not_trained) > 0:
        log_and_message(f"{len(not_trained)} file{' has' if len(not_trained) == 1 else 's have'} not been analyzed.\nMaking testSets for all genomes using new kmers file.", c="GREEN", stderr=True)
        retrain_all = True
    elif args.retrain:
        log_and_message(f"Making testSets for all genomes upon user's request.", c="GREEN", stderr=True)
        retrain_all = True
    else:
        log_and_message(f"There's no need to make testSets for any genome.", c="GREEN", stderr=True)

    if retrain_all:
        for i, file_name in enumerate(sorted(training_data['genomes']), 1):
            log_and_message(f"[{i}/{len(training_data['genomes'])}] Making testSet for {file_name}.", c="YELLOW", stderr=True)
            infile = get_file_path(file_name, infiles, args.indir)

            cmd = ['PhiSpy.py', infile,
                '-o', TEST_DIR,
                '-m', file_name + '.testSet',
                # '--kmer_size', str(args.kmer_size), # TODO not supported by PhiSpy yet
                '--kmers_type', args.kmers_type]
            if args.phmms: cmd.extend(['--phmms', args.phmms, '--threads', args.threads])
            log_and_message(f"PhiSpy command: {' '.join(cmd)}", stderr=True)
            log_and_message(f"{'PhiSpy start':=^30s}", c="PINK", stderr=True)
            call(cmd)
            log_and_message(f"{'PhiSpy stop':=^30s}", c="PINK", stderr=True)


    log_and_message("Updating training sets based on new groups.", c="GREEN", stderr=True)
    write_training_sets(training_data)


    # update trainingGenome_list file - this will act as a new groups file
    # for genomes available in PhiSpy's data directory
    log_and_message("Writing updated trainingGenome_list.", c="GREEN", stderr=True)
    write_training_genome_list(training_data)

    log_and_message("Done!", c="GREEN", stderr=True)

if __name__ == '__main__':
    main()
