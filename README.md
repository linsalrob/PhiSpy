[![DOI](https://www.zenodo.org/badge/60999054.svg)](https://www.zenodo.org/badge/latestdoi/60999054)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/PhiSpy)
[![Build Status](https://travis-ci.org/linsalrob/PhiSpy.svg?branch=master&label=Travis%20Build)](https://travis-ci.org/linsalrob/PhiSpy)
[![PyPi](https://img.shields.io/pypi/pyversions/phispy.svg?style=flat-square&label=PyPi%20Versions)](https://pypi.org/project/PhiSpy/)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/phispy.svg?style=flat-square&label=BioConda%20install)](https://anaconda.org/bioconda/phispy)
[![Downloads](https://img.shields.io/github/downloads/linsalrob/PhiSpy/total?style=flat-square)](https://github.com/linsalrob/PhiSpy/releases)


# What is PhiSpy?

PhiSpy identifies prophages in Bacterial (and probably Archaeal) genomes. Given an annotated genome it will use several approaches to identify the most likely prophage regions.

Initial versions of PhiSpy were written by 

Sajia Akhter (sajia@stanford.edu)
[Edwards Bioinformatics Lab](http://edwards.sdsu.edu/research/)

Improvements, bug fixes, and other changes were made by

Katelyn McNair
[Edwards Bioinformatics Lab](http://edwards.sdsu.edu/research/)
and Przemyslaw Decewicz 
[DEMB at the University of Warsaw](http://ddlemb.com/)



# Installation

## Conda

The easiest way to install for all users is to use `bioconda`.

```bash
conda install -c bioconda phispy
```


## PIP

`python-pip` requires a C++ compiler and the Python header files. You should be able to install it like this:


```bash
sudo apt install -y build-essential python3-dev python3-pip
python3 -m pip install --user PhiSpy
```
This will install `PhiSpy.py` in `~/.local/bin` which should be in your `$PATH` but might not be (see [this](https://bugs.launchpad.net/ubuntu/+source/bash/+bug/1588562) detailed discussion). See the tips and tricks below for a solution to this.


## Advanced Users

For advanced users, you can clone the git repository and use that (though `pip` is the recommended install method). 
```bash
git clone https://github.com/linsalrob/PhiSpy.git
cd PhiSpy`
python3 setup.py install --user --record installed_files.txt
```

Note that we recommend using --record to save a list of all the files that were installed by `PhiSpy`. If you ever want to uninstall it, or to remove everything to reinstall e.g. from `pip`, you can simply use the contents of that file:

```
cat installed_files.txt | xargs rm -f
```



If you have root and you want to install globally, you can change the setup command.

```bash
git clone https://github.com/linsalrob/PhiSpy.git
cd PhiSpy`
python3 setup.py install
```

For ease of use, you may wish to add the location of PhiSpy.py to your $PATH.

## Software Requirements

PhiSpy requires following programs to be installed in the system. Most of these are likely already on your system or will be installed using the mechanisms above.

1. `Python` - version 3.4 or later
2. `Biopython` - version 1.58 or later 
3. `gcc` - GNU project C and C++ compiler - version 4.4.1 or later
4. The `Python.h` header file. This is included in `python3-dev` that is available on most systems.

# Testing PhiSpy.py

Download the [Streptococcus pyogenes M1 genome](https://raw.githubusercontent.com/linsalrob/PhiSpy/master/tests/Streptococcus_pyogenes_M1_GAS.gb) 

```bash
curl -Lo Streptococcus_pyogenes_M1_GAS.gb https://bit.ly/37qFArb
PhiSpy.py -o Streptococcus.phages Streptococcus_pyogenes_M1_GAS.gb
```

or to run it with the `Streptococcus` training set:

```bash
PhiSpy.py -o Streptococcus.phages -t data/trainSet_160490.61.txt Streptococcus_pyogenes_M1_GAS.gb
```

This uses the `GenBank` format file for *Streptococcus pyogenes* M1 GAS that we provide in the [tests/](tests/) directory, and we use the training set for *S. pyogenes* M1 GAS that we have pre-calculated. This quickly identifies the four prophages in this genome, runs the repeat finder on all of them, and outputs the answers.

You will find the output files from this query in `output_directory`.


# Running PhiSpy.py

The simplest command is:

```bash
PhiSpy.py genbank_file -o output_directory
```

where:
- `genbank file`: The input DNA sequence file in GenBank format.
- `output directory`: The output directory is the directory where the final output file will be created.

If you have new genome, we recommend annotating it using the [RAST server](http://rast.nmpdr.org/rast.cgi) or [PROKKA](https://github.com/tseemann/prokka).
 
After annotation, you can download the genome directory from the server.

When also considering the signal from HMM profile search:
```bash
PhiSpy.py genbank_file -o output_directory --phmms hmm_db --threads 4 --color
```
where:
- `hmm_db`: reference HMM profiles database to search with genome-encoded proteins (at the moment)

Training sets were searched with [pVOG database](http://dmk-brain.ecn.uiowa.edu/pVOGs) HMM profiles: [AllvogHMMprofiles.tar.gz](http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz).
To use it:
```bash
wget http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz
tar -zxvf AllvogHMMprofiles.tar.gz
cat AllvogHMMprofiles/* > pVOGs.hmm
```
Then use `pVOGs.hmm` as `hmm_db`.

Since extra step before the regular processing of PhiSpy is performed, input `genbank file` is updated and saved in `output_directory/phmms_search`.
When `--color` flag is used, additional qualifier `/color` will be added in the updated GenBank file so that the user could easily distinguished proteins with hits to `hmm_db` while viewing the file in [Artemis](https://www.sanger.ac.uk/science/tools/artemis)

When running PhiSpy again on the same input data and with `--phmms` option you can skip the search step by `--skip_search` flag.

# Help

For the help menu use the `-h` option:
```bash
python PhiSpy.py -h
```

# Output Files

There are 3 output files, located in output directory.

1. prophage.tbl: This file has two columns separated by tabs [id, location]. 
The id is in the format: pp_number, where number is a sequential number of the prophage (starting at 1). 
Location is be in the format: contig_start_stop that encompasses the prophage.
 
2. prophage_tbl.tsv: This is a tab seperated file. The file contains all the genes of the genome. The tenth colum represents the status of a gene. If this column is 1 then the gene is a phage like gene; otherwise it is a bacterial gene. 

This file has 16 columns:(i) fig_no: the id of each gene; (ii) function: function of the gene;	(iii) contig; (iv) start: start location of the gene; (v) stop: end location of the gene; (vi) position: a sequential number of the gene (starting at 1); (vii)	rank: rank of each gene provided by random forest; (viii) my_status: status of each gene based on random forest; (ix) pp: classification of each gene based on their function; (x) Final_status: the status of each gene. For prophages, this column has the number of the prophage as listed in prophage.tbl above; If the column contains a 0 we believe that it is a bacterial gene. If we can detect the _att_ sites, the additional columns will be: (xi) start of _attL_; (xii) end of _attL_; (xiii) start of _attR_; (xiv) end of _attR_; (xv) sequence of _attL_; (xvi) sequence of _attR_.

3. prophage_coordinates.tsv: This file has the prophage ID, contig, start, stop, and potential _att_ sites identified for the phages.

# Example Data

* _Streptococcus pyogenes_ M1 GAS which has a single genome contig. The genome contains four prophages.

To analyze this data, you can use:

```
PhiSpy.py -o output_directory -t data/trainSet_160490.61.txt tests/Streptococcus_pyogenes_M1_GAS.gb
```

And you should get a prophage table that has this information (for example, take a look at `output_directory/prophage.tbl`).

| Prophage number | Contig | Start | Stop |
| --- | --- | --- | --- | 
pp_1 | NC_002737 | 529631 | 569288
pp_2 | NC_002737 | 778642 | 820599
pp_3 | NC_002737 | 1192630 | 1222549
pp_4 | NC_002737 | 1775862 | 1782822


# Tips, Tricks, and Errors

If you are feeling lazy, you actually only need to use `sudo apt install -y python3-pip; python3 -m pip install phispy` since python3-pip requires `build-essential` and `python3-dev`! 

If you try `PhiSpy.py -v` and get an error like this:

```bash
$ PhiSpy.py -v
-bash: PhiSpy.py: command not found
```

Then you can either use the full path:

```bash
~/.local/bin/PhiSpy.py -v
```

or add that location to your `$PATH`:

```bash
echo "export PATH=\$HOME/.local/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc
PhiSpy.py -v
```

# Making your own training sets

If within reference datasets, close relatives to bacteria of your interest are missing, you can make your own training sets by providing at least a single genome in which you indicate prophage proteins. This is done by adding a new qualifier GenBank annotation for each CDS feature within a prophage region: `/is_phage="1"`. This allows PhiSpy to distinguish the signal from bacterial/phage regions and make a training set to use afterwards during classification with random forest algorithm. 

To make a training set out of your files use `make_training_sets.py` script:

```bash
python3 scripts/make_training_sets.py -d input_directory -o output_directory -k kmer_size -t kmers_type -g groups_file --retrain --phmms hmm_db --color --threads 4
```
where:
- `input_directory`: a directory where GenBank files to use for training are stored.
- `output_directory`: a directory where are temporary and final training sets will be written.
- `kmer_size`: is the size of kmers that will be produces. By default it's 12. If changed, remember to also change that parameter while running PhiSpy with produced training sets.
- `kmers_type`: type of generated kmers. By default 'all' means generating kmers by 1 nt. If changed, remember to also change that parameter while running PhiSpy with produced training sets.
- `groups_file`: a file mapping GenBank file names with extension and the name of group they will make; each file can be assigned to more than one group - take a look at how the reference data grouping file was constructed: `tests/groups.txt`.
Beside the flags that allow training with phmm signal, there's also a `--retrain` flag. When used, it overwrites all the training sets in the `output_directory` that will be produced while training. That includes: `phage_kmers_all_wohost.txt`, `trainSets_genericAll.txt` and `trainingGenome_list.txt`. The same will happen when `trainingGenome_list.txt` is missing in `output_directory`.
If `--retrain` is not set, the script extends the `trainingGenome_list.txt`, adds new files to `output_directory` (overwrites the ones with the same group name) and updates `phage_kmers_all_wohost.txt`. 

The default reference data is trained by the following command:
```bash
python3 scripts/make_training_sets.py -d tests -o PhiSpyModules/data -g tests/groups.txt --retrain --phmms pVOGs.hmm --color --threads 4
```
You can modify/update the `test` directory and `groups.txt` file for your needs.

Within the `output_directory` you will find a `trainSets` directory with a single trainSet file for each genome and (if requested) `phmms_search` directory with updated input GenBank files and temp files from the last hmmsearch.

## Preparing GenBank files
- it is recommended to mark prophage proteins even from prophage remnants/disrupted regions composed of a few proteins with `is_phage="1"` to minimize the loss of good signal, kmers in particular,
- don't use too many genomes (e.g. a 100) as you may end up with a small set of phage-specific kmers,
- try to pich several genomes with different prophages too increase the diversity.