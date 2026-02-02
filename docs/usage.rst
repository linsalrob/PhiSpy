Usage Guide
===========

This guide provides detailed information about using PhiSpy effectively.

Basic Usage
-----------

The simplest command is:

.. code-block:: bash

   PhiSpy.py genbank_file -o output_directory

Where:

- ``genbank_file``: Input DNA sequence file in GenBank format
- ``output_directory``: Directory where output files will be created

Important Parameters
--------------------

phage_genes
^^^^^^^^^^^

The ``--phage_genes`` parameter controls how strict PhiSpy is in calling prophages.

By default, PhiSpy uses *strict* mode, looking for 1 or more genes that are likely to be phage genes in each prophage region.

**Increasing the value** (e.g., ``--phage_genes 5``) will:

- Reduce the number of prophages predicted
- Make predictions more specific but less sensitive
- Only report regions with strong phage signals

**Setting to 0** (``--phage_genes 0``) will:

- Identify other mobile elements (plasmids, integrons, pathogenicity islands)
- Also identify ribosomal RNA operons (they're unlike the host backbone)
- Generate many false positives

Example::

   PhiSpy.py genome.gb -o results --phage_genes 5

Training Sets
^^^^^^^^^^^^^

Training sets improve prediction accuracy by providing information about prophages in related organisms.

List available training sets::

   PhiSpy.py --list short

Use a specific training set::

   PhiSpy.py genome.gb -o results -t data/trainSet_Streptococcus.txt

The default training set (``trainSet_genericAll.txt``) works well for most genomes.

File Name Prefixes
^^^^^^^^^^^^^^^^^^

When analyzing multiple genomes, use prefixes to avoid overwriting outputs::

   PhiSpy.py genome1.gb -o results -p genome1_
   PhiSpy.py genome2.gb -o results -p genome2_

All output files for genome1 will have the prefix ``genome1_``.

Gzip Support
^^^^^^^^^^^^

PhiSpy natively supports gzip format for both input and output:

- If you provide a gzipped input file, PhiSpy will write gzipped output files
- No need to manually decompress/compress files

Example::

   PhiSpy.py genome.gb.gz -o results

HMM Searches
------------

PhiSpy can use HMM profile searches to improve prophage detection.

Using pVOGs Database
^^^^^^^^^^^^^^^^^^^^

Download and prepare the pVOGs database:

.. code-block:: bash

   wget http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz
   tar -zxvf AllvogHMMprofiles.tar.gz
   cat AllvogHMMprofiles/* > pVOGs.hmm

Run PhiSpy with pVOGs::

   PhiSpy.py genome.gb -o results --phmms pVOGs.hmm --threads 4 --color

Using VOGdb Database
^^^^^^^^^^^^^^^^^^^^

Download and prepare VOGdb:

.. code-block:: bash

   curl -LO http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz
   mkdir vog
   tar -C vog -xf vog.hmm.tar.gz
   cat vog/* > VOGs.hmms
   hmmpress VOGs.hmms

Run PhiSpy with VOGdb::

   PhiSpy.py genome.gb -o results --phmms VOGs.hmms --threads 4

HMM Search Features
^^^^^^^^^^^^^^^^^^^

When using ``--phmms``:

- The input GenBank file is updated and saved in the output directory
- With ``--color`` flag, proteins with HMM hits are colored for Artemis visualization
- Use ``--skip_search`` to skip the search step when re-running on the same data

Metrics
-------

PhiSpy uses several metrics to identify prophages:

Default Metrics
^^^^^^^^^^^^^^^

When no ``--metrics`` flag is provided, all metrics are used:

- **orf_length_med**: Median ORF length
- **shannon_slope**: Slope of Shannon's diversity of k-mers
- **at_skew**: Normalized AT skew
- **gc_skew**: Normalized GC skew
- **max_direction**: Maximum number of genes in the same direction

Specifying Metrics
^^^^^^^^^^^^^^^^^^

You can choose specific metrics:

Single metric::

   PhiSpy.py genome.gb -o results --metrics shannon_slope

Multiple metrics (method 1)::

   PhiSpy.py genome.gb -o results --metrics shannon_slope gc_skew

Multiple metrics (method 2)::

   PhiSpy.py genome.gb -o results --metrics shannon_slope --metrics gc_skew

The metrics used are recorded in the log file.

Additional Options
^^^^^^^^^^^^^^^^^^

**expand_slope**: Improves Shannon score calculations::

   PhiSpy.py genome.gb -o results --expand_slope

**kmers_type**: Controls k-mer generation method::

   PhiSpy.py genome.gb -o results --kmers_type codon

Advanced Usage
--------------

Color Output for Artemis
^^^^^^^^^^^^^^^^^^^^^^^^^

The ``--color`` flag adds color qualifiers to the GenBank output::

   PhiSpy.py genome.gb -o results --color --phmms pVOGs.hmm

Open the resulting GenBank file in `Artemis <https://sanger-pathogens.github.io/Artemis/>`_ to see colored CDS features.

Choosing Output Files
^^^^^^^^^^^^^^^^^^^^^^

Control which files are generated using ``--output_choice``:

Minimal output (coordinates only)::

   PhiSpy.py genome.gb -o results --output_choice 1

Default output (coordinates + GenBank)::

   PhiSpy.py genome.gb -o results --output_choice 3

Include prophage information file::

   PhiSpy.py genome.gb -o results --output_choice 11

All outputs::

   PhiSpy.py genome.gb -o results --output_choice 512

See the `Output Files <output.html>`_ page for details.

Interactive PhiSpy
------------------

A `Jupyter notebook <https://github.com/linsalrob/PhiSpy/blob/master/jupyter_notebooks/PhiSpy.ipynb>`_ is available for interactive exploration.

Use it to:

- Test different parameters interactively
- Visualize how parameters affect predictions
- Explore your genome's prophage content

Change the genome file path and parameter values to see how predictions vary.

Performance Optimization
------------------------

Multi-threading
^^^^^^^^^^^^^^^

Use the ``--threads`` parameter to speed up HMM searches and random forest::

   PhiSpy.py genome.gb -o results --phmms pVOGs.hmm --threads 8

Memory Usage
^^^^^^^^^^^^

For large genomes or many genomes:

- Process one genome at a time
- Use appropriate ``--min_contig_size`` to filter small contigs
- Consider using a high-performance computing cluster

Best Practices
--------------

1. **Always annotate your genome** first using RAST or PROKKA
2. **Start with default parameters** to get a baseline
3. **Examine the output** carefully (see `Assessing Predictions <output.html#assessing-predictions>`_)
4. **Adjust parameters** based on results:
   
   - Too many prophages? Increase ``--phage_genes``
   - Too few prophages? Decrease ``--phage_genes`` or use ``--metrics``
   - Related organism available? Use appropriate ``--training_set``

5. **Use HMM searches** when available for better accuracy
6. **Review predictions** in context of genome biology

Common Workflows
----------------

Standard Workflow
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # 1. Annotate genome (using PROKKA as example)
   prokka --outdir prokka_results genome.fasta
   
   # 2. Run PhiSpy with default settings
   PhiSpy.py prokka_results/genome.gbk -o phispy_results
   
   # 3. Review output files
   less phispy_results/prophage_coordinates.tsv
   
   # 4. If needed, adjust and re-run
   PhiSpy.py prokka_results/genome.gbk -o phispy_results_v2 --phage_genes 3

High-Quality Workflow
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Get HMM database
   wget http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz
   tar -zxvf AllvogHMMprofiles.tar.gz
   cat AllvogHMMprofiles/* > pVOGs.hmm
   
   # Run PhiSpy with HMM and appropriate training set
   PhiSpy.py genome.gb -o results \
       --phmms pVOGs.hmm \
       --threads 8 \
       --color \
       --output_choice 11 \
       -t data/trainSet_YourOrganism.txt

Batch Processing
^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Process multiple genomes
   for genome in genomes/*.gb; do
       name=$(basename "$genome" .gb)
       PhiSpy.py "$genome" -o results -p "${name}_" --threads 4
   done
