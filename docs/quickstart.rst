Quick Start
===========

This guide will help you get started with PhiSpy quickly.

Testing PhiSpy
--------------

To test that PhiSpy is working correctly, download the *Streptococcus pyogenes* M1 genome:

.. code-block:: bash

   curl -Lo Streptococcus_pyogenes_M1_GAS.gb https://bit.ly/37qFArb
   PhiSpy.py -o Streptococcus.phages Streptococcus_pyogenes_M1_GAS.gb

This should identify four prophages in the genome.

To use a specific training set:

.. code-block:: bash

   PhiSpy.py -o Streptococcus.phages -t data/trainSet_160490.61.txt Streptococcus_pyogenes_M1_GAS.gb

Basic Usage
-----------

The simplest command to run PhiSpy is:

.. code-block:: bash

   PhiSpy.py genbank_file -o output_directory

Where:

- ``genbank_file``: The input DNA sequence file in GenBank format
- ``output_directory``: The directory where output files will be created

Preparing Your Genome
---------------------

PhiSpy requires an annotated genome in GenBank format. If you have a new genome that needs annotation, we recommend:

1. **RAST** (http://rast.nmpdr.org/rast.cgi) - A web server that allows you to upload and download annotated genomes
2. **PROKKA** (https://github.com/tseemann/prokka) - Stand-alone annotation software

Example Workflow
----------------

1. **Annotate your genome** using RAST or PROKKA
2. **Run PhiSpy** with default parameters:

   .. code-block:: bash

      PhiSpy.py my_genome.gb -o phispy_results

3. **Review the output** files in the ``phispy_results`` directory
4. **Adjust parameters** if needed and re-run

Common Use Cases
----------------

Find prophages with default settings:

.. code-block:: bash

   PhiSpy.py genome.gb -o results

Find prophages and color CDSs for Artemis visualization:

.. code-block:: bash

   PhiSpy.py genome.gb -o results --color

Find prophages using HMM profiles:

.. code-block:: bash

   PhiSpy.py genome.gb -o results --phmms pVOGs.hmm --threads 4

Get all output files:

.. code-block:: bash

   PhiSpy.py genome.gb -o results --output_choice 512

Next Steps
----------

- Read the `Usage Guide <usage.html>`_ for detailed information about parameters
- Learn about `Output Files <output.html>`_ to understand the results
- See `Training Sets <training.html>`_ to create custom training data
