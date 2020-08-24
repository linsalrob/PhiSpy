[![Edwards Lab](https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4)](https://edwards.sdsu.edu/research)
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
cd PhiSpy
python3 setup.py install --user --record installed_files.txt
```

Note that we recommend using --record to save a list of all the files that were installed by `PhiSpy`. If you ever want to uninstall it, or to remove everything to reinstall e.g. from `pip`, you can simply use the contents of that file:

```
cat installed_files.txt | xargs rm -f
```



If you have root and you want to install globally, you can change the setup command.

```bash
git clone https://github.com/linsalrob/PhiSpy.git
cd PhiSpy
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

## Download more testing data

You can also download all the genomes in [tests/](tests). These are not installed with PhiSpy if you use pip/conda, but will be if you clone the repository.
Please note that these are stored on [git lfs](https://git-lfs.github.com/), and so if you notice an error that the files are small and don't ungzip, you may need to (i) install
`git lfs` and (ii) use `git lfs fetch` to update this data.

# Running PhiSpy.py

The simplest command is:

```bash
PhiSpy.py genbank_file -o output_directory
```

where:
- `genbank file`: The input DNA sequence file in GenBank format.
- `output directory`: The output directory is the directory where the final output file will be created.

If you have new genome, we recommend annotating it using the [RAST server](http://rast.nmpdr.org/rast.cgi) or [PROKKA](https://github.com/tseemann/prokka). RAST has a server that allows you to upload and download the genome (and can handle lots of genomes), while PROKKA is stand-alone software.

### phage_genes

By default, `PhiSpy.py` uses *strict* mode, where we look for two or more genes that are likely to be a phage in each prophage region. If you increase the value of `--phage_genes` that will reduce the number of prophages that are predicted. Conversely, if you reduce this, or set it to `0` we will overcall mobile elements. 

When `--phage_genes` is set to `0`, `PhiSpy.py` will identify other mobile elements like plasmids, integrons, and pathogenicity islands. Somewhat unexpectedly, it will also identify the ribosomal RNA operons as likely being mobile since they are unlike the host's backbone!

### color

If you add the `--color` flag, we will color the CDS based on their function. The colors are primarily used in [artemis](https://sanger-pathogens.github.io/Artemis/) for visualizing phage regions.

### file name prefixes

By default the outputs from `PhiSpy.py` have standard names. If you supply a file name prefix it will be prepended to all the file so that you can run `PhiSpy.py` on multiple genomes and have the outputs in the same directory without overwriting each other.

### gzip support

`PhiSpy.py` natively supports both reading and writing files in `gzip` format. If you provide a `gzipped` input file, we will write a `gzipped` output file.

### HMM Searches

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

Since extra step before the regular processing of PhiSpy is performed, input `genbank file` is updated and saved in `output_directory`.
When `--color` flag is used, additional qualifier `/color` will be added in the updated GenBank file so that the user could easily distinguished proteins with hits to `hmm_db` while viewing the file in [Artemis](https://www.sanger.ac.uk/science/tools/artemis)

When running PhiSpy again on the same input data and with `--phmms` option you can skip the search step by `--skip_search` flag.

Another database that maybe of interest is the [VOGdb](http://vogdb.org/) database. You can download all their VOGs, and the press them into a compiled format for `hmmer`:

```bash
curl -LO http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz
mkdir vog
tar -C vog -xf vog.hmm.tar.gz
cat vog/* > VOGs.hmms
hmmpress VOGs.hmms
```

### Metrics

We use several different metrics to predict regions that are prophages, and there are some optional metrics you can add. The default set of metrics are:

 - `orf_length_med`: median ORF length
 - `shannon_slope`: the slope of Shannon's diversity of _k_-mers across the window under consideration. You can also expand this with the `--expand_slope` option.
 - `at_skew`: the normalized AT skew across the window under consideration
 - `gc_skew`: the normalized GC skew across the window under consideration
 - `max_direction`: The maximum number of genes in the same direction

You can also add

 - `phmms`: The [phmm](#HMM-Searches) search results
 - `phage_genes`: The number of genes that must be annotated as phage in the region
 - `nonprophage_genegaps` : The maximum number of non-phage genes between two phage-like regions that will enable them to be merged

# Help

For the help menu use the `-h` option:
```bash
python PhiSpy.py -h
```

# Output Files

`PhiSpy` has the option of creating multiple output files with the prophage data:

1. **prophage_coordinates.tsv** (code: 1)

This is the coordinates of each prophage identified in the genome, and their _att_ sites (if found) in tab
separated text format. 

The columns of the file are:
  - 1. Prophage number
  - 2. The contig upon which the prophage resides
  - 3. The start location of the prophage
  - 4. The stop location of the prophage
If we can detect the _att_ sites, the additional columns are:
  - 11. start of _attL_;
  - 12. end of _attL_;
  - 13. start of _attR_;
  - 14. end of _attR_;
  - 15. sequence of _attL_;
  - 16. sequence of _attR_;
  - 17. The explanation of why this _att_ site was chosen for this prophage.

2. **GenBank format output** (code: 2)

We provide a duplicate GenBank record that is the same as the input record, but we have inserted the
prophage information, including _att_ sites into the record.

If the original GenBank file was provided in `gzip` format this file will also be created in gzip format.

3. **prophage and bacterial sequences**  (code: 4)

`PhiSpy` can automatically separate the DNA sequences into prophage and bacterial components. If this output is chosen, we generate both fasta and
GenBank format outputs: 
 - _GenBank files_: Two files are made, one for the bacteria and one for the phages. Each contains the appropriate fragments of the genome annotated as in the original.
 - _fasta files_: Two files are made, the first contains the entire genome, but the prophage regions have been masked with `N`s. We explicitly chose this format for a
few reasons: (i) it is trivial to convert this format into separate contigs without the Ns but it is more complex to go from separate contigs
back to a single joined contig; (ii) when read mapping against the genome, understanding that reads map either side of a prophage maybe
important; (iii) when looking at insertion points this allows you to visualize the where the prophage was lying.
   
4. **prophage_information.tsv**  (code: 8)
 
 This is a tab separated file, and is the key file to assess prophages in genomes (see [assessing predictions](#assessing-predictions), below). The file contains all the genes of the genome, one per line.
 The tenth colum represents the status of a gene. If this column is 0 then we consider this a bacterial gene. 
 If it is non-zero it is probably a phage gene, and the higher the score the more likely we believe it is a phage
  gene. This is the raw data that we use to identify the prophages in your genome.

This file has 16 columns:
  - 1. The id of each gene; 
  - 2. function: function of the gene (or `product` from a GenBank file);
  - 3. contig;
  - 4. start: start location of the gene;
  - 5. stop: end location of the gene; 
  - 6. position: a sequential number of the gene (starting at 1);
  - 7. rank: rank of each gene provided by random forest; 
  - 8. my_status: status of each gene based on random forest;
  - 9. pp: classification of each gene based on their function;
  - 10. Final_status: the status of each gene. For prophages, this column has the number of the prophage as listed in prophage.tbl above; 
   If the column contains a 0 we believe that it is a bacterial gene. Otherwise we believe that it is possibly a phage gene.
   
If we can detect the _att_ sites, the additional columns are:
  - 11. start of _attL_;
  - 12. end of _attL_;
  - 13. start of _attR_;
  - 14. end of _attR_;
  - 15. sequence of _attL_;
  - 16. sequence of _attR_;


5. **prophage.tsv**  (code: 16)

This is a simpler version of the _prophage_coordinates.tsv_ file that only has prophage number, contig, start, and stop.

6. **GFF3 format**  (code: 32)

This is the prophage information suitable for insertion into a [GFF3](https://m.ensembl.org/info/website/upload/gff3.html). 
This is a legacy file format, however, since GFF3 is no longer widely supported, this only has the prophage coordinates. 
Please post an issue on GitHub if more complete GFF3 files are required.

7. **prophage.tbl**  (code: 64)

This file has two columns separated by tabs [prophage_number, location]. This is a also a 
legacy file that is not generated by default. The prophage number is a sequential number of the 
prophage (starting at 1), and the location is in the format: contig_start_stop that encompasses the prophage.

8. **test data** (code: 128) 

This file has the data used in the random forest. The columns are:
 - Identifier
 - Median ORF length
 - Shannon slope
 - Adjusted AT skew
 - Adjusted GC skew
 - The maxiumum number of ORFs in the same direction
 - PHMM matches
 - Status

The numbers are averaged across a window of size specified by `--window_size`

## Choosing which output files are created.

We have provided the option (`--output_choice`) to choose which output files are created. Each file above has a code associated with it, 
and to include that file add up the codes:

Code | File
--- | ---
1 | prophage_coordinates.tsv 
2 | GenBank format output 
4 | prophage and bacterial sequences
8 | prophage_information.tsv
16 | prophage.tsv
32 | GFF3 format
64 | prophage.tbl 
128 | test data used in the random forest


So for example, if you want to get `GenBank format output` (2) and `prophage_information.tsv` (8), then 
enter an `--output_choice` of 10.

The default is 3: you will get both the `prophage_coordinates.tsv` and `GenBank format output` files.

If you want _all_ files output, use `--output_choice 255`.

# Example Data

* _Streptococcus pyogenes_ M1 GAS which has a single genome contig. The genome contains four prophages.

To analyze this data, you can use:

```
PhiSpy.py -o output_directory -t data/trainSet_160490.61.txt tests/Streptococcus_pyogenes_M1_GAS.gb.gz 
```

And you should get a prophage table that has this information (for example, take a look at `output_directory/prophage.tbl`).

| Prophage number | Contig | Start | Stop |
| --- | --- | --- | --- | 
pp_1 | NC_002737 | 529631 | 569288
pp_2 | NC_002737 | 778642 | 820599
pp_3 | NC_002737 | 1192630 | 1222549
pp_4 | NC_002737 | 1775862 | 1782822

# Assessing predictions

As with any software, it is critical that you assess the output from `phispy` to see if it actually makes sense! We start be ensuring we have the `prophage_information.tsv` file output (this is not output by default, and requires adding 8 to the `--output-choice` flag).

That is a tab-separated text file that you can import into Microsoft Excel, LibreOffice Calc, Google Sheets, or your favorite spreadsheet viewing program. 

There are a few columns that you should pay attention to:
- _position_ (the 6<sup>th</sup> column) is the position of the gene in the genome. If you sort by this column you will always return the genome to the original order.
- _Final status_ (the 10<sup>th</sup> column) is whether this region is predicted to be a prophage or not. The number is the prophage number. If the entry is 0 it is not a prophage.
- _pp_ and _my status_ (the 8<sup>th</sup> and 9<sup>th</sup> columns) are interim indicators about whether this gene is potentially part of a phage.

We recommend:
1. Freeze the first row of the spreadsheet so you can see the column headers
2. Sort the spreadsheet by the _my status_ column and color any row red where the value in this column is greater than 0
3. Sort the spreadsheet by the _final status_ column and color those rows identified as a prophage green.
4. Sort the spreadsheet by the _position_ column.

Now all the prophages are colored green, while all the potential prophage genes that are not included as part of a prophage are colored red. You can easily review those non-prophage regions and determine whether _you_ think they should be included in prophages. Note that in most cases you can adjust the `phispy` parameters to include regions you think are prophages.

**Note:** Ensure that while you are reviewing the results, you pay particular attention to the _contig_ column. In partial genomes, contig breaks are very often located in prophages. This is usual because prophages often contain sequences that are repeated around the genome. We have an [open issue](https://github.com/linsalrob/PhiSpy/issues/33) open issue to try and resolve this in a meaningful way.

# Interactive PhiSpy

We have created a [jupyter notebook](https://github.com/linsalrob/PhiSpy/blob/master/jupyter_notebooks/PhiSpy.ipynb) 
example where you can run `PhiSpy` to test the effect of the different parameters
on your prophage predictions. Change the name of the genbank file to point to your genome, and 
change the values in `parameters` and see how the prophage predictions vary!

# Tips, Tricks, and Errors

If you are feeling lazy, you actually only need to use `sudo apt install -y python3-pip; python3 -m pip install phispy` 
since python3-pip requires `build-essential` and `python3-dev`! 

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

# Exit (error) codes

We use a few different error codes to signify things that we could not compute. So far, we have:

Exit Code | Meaning | Suggested solution
--- | --- | ---
2 | No input file provided | We need a file to work with!
3 | No output directory provided | We need somewhere to write the results to!
10 | No training sets available | This should be in the default install. Please check your installation
11 | The specific training set is not available | Check the argument passed to the `--training_set` parameter
13 | No kmers file found | This should be in the default install. Please check your installation
20 | IO Error | There was an error reading your input file.
25 | Non nucleotide base found | Check for a non-standard base in your sequence
26 | An ORF with no bases | This is probably a really short ORF and should be deleted.
30 | No contigs | We filter contigs by length, and so try adjusting the `--min_contig_size` parameter, though the default is 5,000 bp and you will need some adjacent genes!
40 | No ORFs in your genbank file | Please annotate your genome, e.g. using [RAST](http://rast.nmpdr.org/) or [PROKKA](https://github.com/tseemann/prokka)
41 | Less than 100 ORFs are in your annotated genome. This is not enough to find a prophage | Please annotate your genome, e.g. using [RAST](http://rast.nmpdr.org/) or [PROKKA](https://github.com/tseemann/prokka)

# Making your own training sets

If within reference datasets, close relatives to bacteria of your interest are missing, you can make your own training sets by providing at least a single genome in which you indicate prophage proteins. This is done by adding a new qualifier GenBank annotation for each CDS feature within a prophage region: `/is_phage="1"`. This allows PhiSpy to distinguish the signal from bacterial/phage regions and make a training set to use afterwards during classification with random forest algorithm. 

To make a training set out of your files use the `PhiSpy.py` option `-m`:

```bash
PhiSpy.py -o output_directory -k kmer_size -t kmers_type -g groups_file --retrain --phmms hmm_db --color --threads 4 genome.gb.gz
```
where:
- `output_directory`: a directory where are temporary and final training sets will be written.
- `kmer_size`: is the size of kmers that will be produces. By default it's 12. If changed, remember to also change that parameter while running PhiSpy with produced training sets.
- `kmers_type`: type of generated kmers. By default 'all' means generating kmers by 1 nt. If changed, remember to also change that parameter while running PhiSpy with produced training sets.
- `groups_file`: a file mapping GenBank file names with extension and the name of group they will make; each file can be assigned to more than one group - take a look at how the reference data grouping file was constructed: `tests/groups.txt`.
Beside the flags that allow training with phmm signal, there's also a `--retrain` flag. When used, it overwrites all the training sets in the `output_directory` that will be produced while training. That includes: `phage_kmers_all_wohost.txt`, `trainSets_genericAll.txt` and `trainingGenome_list.txt`. The same will happen when `trainingGenome_list.txt` is missing in `output_directory`.
- `genome.gb.gz` is a gzip compressed GenBank file that has the `/is_phage="1"` flags set.
If `--retrain` is not set, the script extends the `trainingGenome_list.txt`, adds new files to `output_directory` (overwrites the ones with the same group name) and updates `phage_kmers_all_wohost.txt`. 


## Preparing GenBank files
- it is recommended to mark prophage proteins even from prophage remnants/disrupted regions composed of a few proteins with `/is_phage="1"` to minimize the loss of good signal, kmers in particular,
- don't use too many genomes (e.g. a 100) as you may end up with a small set of phage-specific kmers,
- try to pick several genomes with different prophages to increase the diversity.
