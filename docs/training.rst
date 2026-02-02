Custom Training Sets
====================

PhiSpy can be trained on your own data to improve predictions for specific organisms or clades.

Why Create Training Sets?
--------------------------

Create custom training sets when:

- Your organism is distantly related to reference datasets
- You have manually curated prophage data
- You want organism-specific predictions
- Default training sets don't work well

Requirements
------------

To create training sets, you need:

1. **Annotated GenBank files** with prophage regions marked
2. **Prophage markers**: CDS features in prophage regions must have ``/is_phage="1"`` qualifier
3. **At least one genome** (more is better for diversity)

Marking Prophage Features
--------------------------

Use the ``mark_prophage_features.py`` script to add ``/is_phage="1"`` qualifiers to prophage proteins.

Input Format
^^^^^^^^^^^^

The script accepts a tab-delimited file with:

1. Path to GenBank file
2. Replicon ID
3. Prophage start coordinate
4. Prophage end coordinate

**Example input file (prophages.txt):**

.. code-block:: text

   genome1.gb    NC_002737    529631    569288
   genome1.gb    NC_002737    778642    820599
   genome2.gb    NC_003028    123456    145678

Usage
^^^^^

.. code-block:: bash

   mark_prophage_features.py -i prophages.txt -o marked_genomes/

This updates the GenBank files with ``/is_phage="1"`` qualifiers for all CDS features within the specified prophage regions.

Creating Training Sets
-----------------------

Use ``make_training_sets.py`` to generate training sets from marked GenBank files.

Basic Usage
^^^^^^^^^^^

.. code-block:: bash

   make_training_sets.py -d marked_genomes/ -g groups.txt

The script:

1. Reads marked GenBank files
2. Generates phage-specific and bacteria-specific k-mer sets
3. Calculates training features
4. Creates training set files

Parameters
^^^^^^^^^^

.. option:: -d INDIR, --indir INDIR

   Path to directory containing marked GenBank files for training.

.. option:: -g GROUPS, --groups GROUPS

   Path to groups file mapping GenBank files to training set names.
   
   If not provided, each file gets its own training set.

.. option:: --use_taxonomy

   Use taxonomy information from GenBank files to create groups.
   
   Files without taxonomy are assigned to "Bacteria" group.

.. option:: -k KMER_SIZE, --kmer_size KMER_SIZE

   Size of k-mers to generate. Default: 12
   
   For codon approach, use multiples of 3.

.. option:: -t KMERS_TYPE, --kmers_type KMERS_TYPE

   K-mer generation method:
   
   - ``simple``: Slice sequence from first position
   - ``all``: All possible k-mers (default)
   - ``codon``: K-mers with step of 3 nucleotides

.. option:: --phmms PHMMS

   HMM profile database for additional feature calculation.

.. option:: --threads THREADS

   Number of threads for HMM searches.

.. option:: --retrain

   Trigger complete reanalysis when reference files have changed.
   
   Use when you've modified ``/is_phage="1"`` qualifiers.

.. option:: --absolute_retrain

   Ignore PhiSpy's default reference genomes.
   
   Train only on your provided data.

Groups File Format
^^^^^^^^^^^^^^^^^^

The groups file has two tab-separated columns:

1. GenBank filename (with extension)
2. Group name

**Example groups.txt:**

.. code-block:: text

   genome1.gb    Streptococcus
   genome2.gb    Streptococcus
   genome3.gb    Lactobacillus
   genome4.gb    Streptococcus
   genome5.gb    Lactobacillus

- Multiple files can belong to the same group
- Files can be assigned to multiple groups (list on separate lines)

Complete Example
----------------

Step 1: Prepare Prophage Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a file listing prophages:

.. code-block:: bash

   cat > prophages.txt << EOF
   genomes/strep1.gb    NC_002737    529631    569288
   genomes/strep1.gb    NC_002737    778642    820599
   genomes/strep2.gb    NC_003028    100000    145000
   genomes/lacto1.gb    NC_004567    234567    289012
   EOF

Step 2: Mark Prophage Features
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   mark_prophage_features.py -i prophages.txt -o marked_genomes/

Step 3: Create Groups File
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cat > groups.txt << EOF
   strep1.gb    Streptococcus
   strep2.gb    Streptococcus
   lacto1.gb    Lactobacillus
   EOF

Step 4: Generate Training Sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   make_training_sets.py \
       -d marked_genomes/ \
       -g groups.txt \
       --use_taxonomy \
       --phmms pVOGs.hmm \
       --threads 4 \
       --retrain

This creates training sets for each group in ``PhiSpyModules/data/``.

Step 5: Use Training Sets
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # List available training sets
   PhiSpy.py --list short
   
   # Use your new training set
   PhiSpy.py new_genome.gb -o results -t data/trainSet_Streptococcus.txt

Advanced Options
----------------

Using Taxonomy
^^^^^^^^^^^^^^

The ``--use_taxonomy`` flag automatically groups genomes by taxonomy:

.. code-block:: bash

   make_training_sets.py \
       -d marked_genomes/ \
       --use_taxonomy \
       --retrain

PhiSpy reads taxonomy from GenBank files and creates training sets for:

- Genus level (if available)
- Family level (if genus not available)
- "Bacteria" (if no taxonomy information)

HMM-Enhanced Training
^^^^^^^^^^^^^^^^^^^^^

Include HMM signals in training:

.. code-block:: bash

   # Download pVOGs
   wget http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz
   tar -zxvf AllvogHMMprofiles.tar.gz
   cat AllvogHMMprofiles/* > pVOGs.hmm
   
   # Train with HMM
   make_training_sets.py \
       -d marked_genomes/ \
       -g groups.txt \
       --phmms pVOGs.hmm \
       --threads 8 \
       --retrain

Absolute Retraining
^^^^^^^^^^^^^^^^^^^

To use ONLY your data (ignore PhiSpy defaults):

.. code-block:: bash

   make_training_sets.py \
       -d marked_genomes/ \
       -g groups.txt \
       --absolute_retrain

This is useful for:

- Highly specialized organisms
- Quality control with known prophages
- Specific research questions

Best Practices
--------------

Selecting Genomes
^^^^^^^^^^^^^^^^^

- **Use diverse prophages**: Include genomes with different prophage types
- **Avoid too many genomes**: ~5-20 genomes is often sufficient
- **Include remnants**: Mark even degraded prophage regions
- **Quality over quantity**: Better to have fewer well-annotated genomes

Marking Prophages
^^^^^^^^^^^^^^^^^

- **Be thorough**: Mark all prophage proteins, including small ones
- **Include boundaries**: Capture integration and excision genes
- **Mark remnants**: Even partial prophages provide useful signal
- **Verify manually**: Check that marked regions make biological sense

K-mer Selection
^^^^^^^^^^^^^^^

- **Default (12-mers, all)**: Works well for most cases
- **Codon k-mers**: Better for very AT-rich or GC-rich genomes
- **Larger k-mers**: More specific but require more data
- **Smaller k-mers**: Less specific but work with less data

Training Set Management
^^^^^^^^^^^^^^^^^^^^^^^

- **Use consistent naming**: Name sets after taxonomic groups
- **Document parameters**: Record k-mer size, type, and source genomes
- **Version control**: Keep track of training set versions
- **Test performance**: Validate on held-out genomes

Storage Location
^^^^^^^^^^^^^^^^

Training sets are stored in::

   PhiSpyModules/data/

Temporary files during training are in::

   PhiSpyModules/data/testSets/

These are preserved to speed up retraining.

Troubleshooting
---------------

Not Enough Phage-Specific K-mers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Problem**: Warning about small k-mer sets.

**Solutions**:

- Add more diverse genomes
- Use smaller k-mer size
- Include more prophage regions
- Try ``--kmers_type simple``

Training Set Not Found
^^^^^^^^^^^^^^^^^^^^^^^

**Problem**: PhiSpy can't find your training set.

**Solutions**:

.. code-block:: bash

   # List available sets
   PhiSpy.py --list long
   
   # Use full path if needed
   PhiSpy.py genome.gb -o results -t /path/to/trainSet_Custom.txt

Poor Performance
^^^^^^^^^^^^^^^^

**Problem**: Custom training set performs worse than default.

**Possible causes**:

- Not enough training genomes
- Poorly marked prophage regions
- Overfitting to specific prophage types
- K-mer parameters not optimal

**Solutions**:

- Add more diverse training genomes
- Verify prophage annotations
- Try default training set for comparison
- Experiment with k-mer parameters

Example: Training for Novel Genus
----------------------------------

Complete workflow for a novel bacterial genus:

.. code-block:: bash

   # 1. Annotate genomes
   for genome in genomes/*.fasta; do
       prokka --outdir prokka_output/$(basename $genome .fasta) $genome
   done
   
   # 2. Run PhiSpy with default settings to get initial predictions
   for gb in prokka_output/*/*.gbk; do
       name=$(basename $gb .gbk)
       PhiSpy.py $gb -o initial_predictions/$name --output_choice 11
   done
   
   # 3. Manually review and create prophage coordinate file
   # (Review prophage_information.tsv files and create prophages.txt)
   
   # 4. Mark features
   mark_prophage_features.py -i prophages.txt -o marked_genomes/
   
   # 5. Create training set
   make_training_sets.py \
       -d marked_genomes/ \
       --use_taxonomy \
       --phmms pVOGs.hmm \
       --threads 8 \
       --retrain
   
   # 6. Re-run with custom training set
   for gb in marked_genomes/*.gb; do
       name=$(basename $gb .gb)
       PhiSpy.py $gb -o final_predictions/$name \
           -t data/trainSet_YourGenus.txt \
           --phmms pVOGs.hmm \
           --output_choice 11
   done

Further Reading
---------------

- See the `Parameters Reference <parameters.html>`_ for details on all options
- Check the `Output Files <output.html>`_ page for result interpretation
- Visit the `GitHub repository <https://github.com/linsalrob/PhiSpy>`_ for examples
