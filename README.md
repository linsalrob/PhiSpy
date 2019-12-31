[![DOI](https://www.zenodo.org/badge/60999054.svg)](https://www.zenodo.org/badge/latestdoi/60999054)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub language count](https://img.shields.io/github/languages/count/linsalrob/PhiSpy)
[![Build Status](https://travis-ci.org/linsalrob/PhiSpy.svg?branch=master)](https://travis-ci.org/linsalrob/PhiSpy)
[![PyPi](https://img.shields.io/pypi/pyversions/phispy)](https://pypi.org/project/PhiSpy/)
![Conda](https://img.shields.io/conda/pn/bioconda/phispy)


# What is PhiSpy?

PhiSpy identifies prophages in Bacterial (and probably Archaeal) genomes. Given an annotated genome it will use several approaches to identify the most likely prophage regions.

Initial versions of PhiSpy were written by 

Sajia Akhter (sajia@stanford.edu)
[Edwards Bioinformatics Lab](http://edwards.sdsu.edu/research/)

Improvements, bug fixes, and other changes were made by

Katelyn McNair
[Edwards Bioinformatics Lab](http://edwards.sdsu.edu/research/)
and Przemyslaw Decewicz 
[University of Warsaw](http://en.uw.edu.pl/)



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
python3 setup.py install --user
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
PhiSpy.py -f genbank_file -o output_directory
```

where:
- `genbank file`: The input DNA sequence file in GenBank format.
- `output directory`: The output directory is the directory where the final output file will be created.

If you have new genome, we recommend annotating it using the [RAST server](http://rast.nmpdr.org/rast.cgi) or [PROKKA](https://github.com/tseemann/prokka).
 
After annotation, you can download the genome directory from the server.


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

