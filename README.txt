[DOI] [License: MIT] [GitHub language count] [Build Status]



INTRODUCTION


PhiSpy is a computer program written in C++ and Python to identify
prophages in a complete bacterial genome sequences.

Initial versions of PhiSpy were written by

Sajia Akhter (sajia@stanford.edu) PhD Student Edwards Bioinformatics Lab
(http://edwards.sdsu.edu/labsite/) Computational Science Research Center
(http://www.csrc.sdsu.edu/csrc/) San Diego State University
(http://www.sdsu.edu/)

Improvements, bug fixes, and other changes were made by

Katelyn McNair Edwards Bioinformatics Lab San Diego State University

Przemyslaw Decewicz University of Warsaw



SYSTEM REQUIREMENTS


The program should run on Mac and Unix platforms, although it was not
tested in all platforms.



SOFTWARE REQUIREMENTS


PhiSpy requires following programs to be installed in the system. NOTE:
You can ignore this if you’re using the singularity container method of
installation.

1.  Python - version 3.4 or later
2.  Biopython - version 1.58 or later
3.  gcc - GNU project C and C++ compiler - version 4.4.1 or later



INSTALLATION


For most users, this will create a local installation for you

    git clone https://github.com/linsalrob/PhiSpy.git
    cd PhiSpy`
    python3 setup.py install --user

If you have root and you want to install globally, you can change the
setup command.

    git clone https://github.com/linsalrob/PhiSpy.git
    cd PhiSpy`
    python3 setup.py install

For ease of use, you may wish to add the location of PhiSpy.py to your
$PATH.

TO TEST THE PROGRAM

Change to the install location and run PhiSpy with an example genome:

    % cd PhiSpy`
    % python3 PhiSpy.py   -o output_directory -t data/trainSet_160490.61.txt tests/Streptococcus_pyogenes_M1_GAS.gb

We use the GenBank format file for _Streptococcus pyogenes_ M1 GAS that
we hvae provided in the tests/ directory, and we use the training set
for _S. pyogenes_ M1 GAS that we have pre-calculated. This quickly
identifies the four prophages in this genome, runs the repeat finder on
all of them, and outputs the answers.

You will find the output files from this query in output_directory.



TO RUN PHISPY


The simplest command is:

    % python3 PhiSpy.py -f genbank_file -o output_directory

where: - genbank file: The input DNA sequence file in GenBank format. -
output directory: The output directory is the directory where the final
output file will be created.

If you have new genome, we recommend annotating it using the RAST server
or PROKKA.

After annotation, you can download the genome directory from the server.



HELP


For the help menu use the -h option:

    % python PhiSpy.py -h



OUTPUT FILES


There are 3 output files, located in output directory.

1.  prophage.tbl: This file has two columns separated by tabs [id,
    location]. The id is in the format: pp_number, where number is a
    sequential number of the prophage (starting at 1). Location is be in
    the format: contig_start_stop that encompasses the prophage.

2.  prophage_tbl.tsv: This is a tab seperated file. The file contains
    all the genes of the genome. The tenth colum represents the status
    of a gene. If this column is 1 then the gene is a phage like gene;
    otherwise it is a bacterial gene.

This file has 16 columns:(i) fig_no: the id of each gene; (ii) function:
function of the gene; (iii) contig; (iv) start: start location of the
gene; (v) stop: end location of the gene; (vi) position: a sequential
number of the gene (starting at 1); (vii) rank: rank of each gene
provided by random forest; (viii) my_status: status of each gene based
on random forest; (ix) pp: classification of each gene based on their
function; (x) Final_status: the status of each gene. For prophages, this
column has the number of the prophage as listed in prophage.tbl above;
If the column contains a 0 we believe that it is a bacterial gene. If we
can detect the _att_ sites, the additional columns will be: (xi) start
of _attL_; (xii) end of _attL_; (xiii) start of _attR_; (xiv) end of
_attR_; (xv) sequence of _attL_; (xvi) sequence of _attR_.

3.  prophage_coordinates.tsv: This file has the prophage ID, contig,
    start, stop, and potential _att_ sites identified for the phages.



EXAMPLE DATA


We have provided two different example data sets.

-   _Streptococcus pyogenes_ M1 GAS which has a single genome contig.
    The genome contains four prophages.

To analyze this data, you can use:

    python3 PhiSpy.py -o output_directory -t data/trainSet_160490.61.txt tests/Streptococcus_pyogenes_M1_GAS.gb

And you should get a prophage table that has this information (for
example, take a look at output_directory/prophage.tbl).

  Prophage number   Contig      Start     Stop
  ----------------- ----------- --------- ---------
  pp_1              NC_002737   529631    569288
  pp_2              NC_002737   778642    820599
  pp_3              NC_002737   1192630   1222549
  pp_4              NC_002737   1775862   1782822
