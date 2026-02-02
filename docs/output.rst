Output Files
============

PhiSpy generates multiple output files containing prophage predictions and related data.

Choosing Output Files
----------------------

Use the ``--output_choice`` parameter to control which files are created. Each file has a code, and you add codes together to get multiple files.

.. list-table:: Output File Codes
   :header-rows: 1
   :widths: 10 50

   * - Code
     - File
   * - 1
     - prophage_coordinates.tsv
   * - 2
     - GenBank format output
   * - 4
     - Prophage and bacterial sequences
   * - 8
     - prophage_information.tsv
   * - 16
     - prophage.tsv
   * - 32
     - GFF3 format output of just the prophages
   * - 64
     - prophage.tbl
   * - 128
     - Test data used in the random forest
   * - 256
     - GFF3 format output for annotated genomic contigs

**Examples:**

- ``--output_choice 3`` (default): Files 1 + 2 (coordinates + GenBank)
- ``--output_choice 10``: Files 2 + 8 (GenBank + information)
- ``--output_choice 512``: All files

Output File Descriptions
-------------------------

1. prophage_coordinates.tsv (Code: 1)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tab-separated file with prophage coordinates and attachment (att) sites.

**Basic Columns:**

1. Prophage number
2. Contig name
3. Start location
4. Stop location

**If att sites are detected (additional columns):**

11. Start of attL
12. End of attL
13. Start of attR
14. End of attR
15. Sequence of attL
16. Sequence of attR
17. Explanation of why this att site was chosen

**Example:**

.. code-block:: text

   pp_1    NC_002737    529631    569288    ...
   pp_2    NC_002737    778642    820599    ...

2. GenBank Format Output (Code: 2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A duplicate of the input GenBank file with prophage information inserted, including att sites.

**Features:**

- Maintains all original annotations
- Adds prophage region annotations
- Includes att site information when found
- If input is gzipped, output is also gzipped

**Usage:** View in genome browsers like Artemis or use for downstream analysis.

3. Prophage and Bacterial Sequences (Code: 4)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Separates DNA sequences into prophage and bacterial components.

**GenBank Files:**

- ``*_bacteria.gb``: Bacterial regions
- ``*_prophage.gb``: Prophage regions

**FASTA Files:**

- ``*.fasta``: Complete genome with prophage regions masked with N's
- Prophage regions are replaced with N's (not removed) to:
  
  - Allow easy conversion to separate contigs if needed
  - Maintain genome structure for read mapping
  - Show insertion points clearly

4. prophage_information.tsv (Code: 8)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**KEY FILE for assessing predictions**. Contains all genes in the genome, one per line, with prophage predictions.

**Columns:**

1. Gene ID
2. Function (product from GenBank)
3. Contig
4. Start location
5. Stop location
6. Position (sequential number starting at 1)
7. Rank (from random forest)
8. my_status (gene status from random forest)
9. pp (classification based on function)
10. Final_status (0 = bacterial, >0 = prophage number)

**If att sites detected (additional columns):**

11-16. Same as in prophage_coordinates.tsv

**Interpreting Final_status:**

- **0**: Bacterial gene
- **1, 2, 3, etc.**: Prophage number
- Higher non-zero values indicate stronger confidence

**Usage:** Import into Excel/LibreOffice/Google Sheets for analysis.

5. prophage.tsv (Code: 16)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Simplified version of prophage_coordinates.tsv with only:

1. Prophage number
2. Contig
3. Start
4. Stop

**Example:**

.. code-block:: text

   pp_1    NC_002737    529631    569288
   pp_2    NC_002737    778642    820599

6. GFF3 Format (Code: 32)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Prophage information in GFF3 format for insertion into genome browsers.

**Note:** This is a legacy format. Contains only prophage coordinates. If more complete GFF3 files are needed, please open a GitHub issue.

7. prophage.tbl (Code: 64)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Legacy format with two tab-separated columns:

1. Prophage number
2. Location (format: contig_start_stop)

**Example:**

.. code-block:: text

   1    NC_002737_529631_569288
   2    NC_002737_778642_820599

8. Test Data (Code: 128)
^^^^^^^^^^^^^^^^^^^^^^^^^

Raw data used in the random forest classification.

**Columns:**

- Identifier
- Median ORF length
- Shannon slope
- Adjusted AT skew
- Adjusted GC skew
- Maximum number of ORFs in the same direction
- PHMM matches
- Status

Numbers are averaged across the window size (``--window_size``).

**Usage:** For debugging or understanding how PhiSpy made predictions.

9. GFF3 Annotated Contigs (Code: 256)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Full genome annotation in GFF3 format, including prophages.

**Usage:** Best choice for loading into Artemis. Handles multiple contigs correctly.

Assessing Predictions
----------------------

It's critical to assess PhiSpy output to ensure predictions make biological sense.

Spreadsheet Analysis
^^^^^^^^^^^^^^^^^^^^

1. Generate the prophage_information.tsv file::

      PhiSpy.py genome.gb -o results --output_choice 11

2. Open prophage_information.tsv in a spreadsheet program

3. **Recommended workflow:**

   a. Freeze the first row (column headers)
   b. Sort by *my_status* column and color rows red where value > 0
   c. Sort by *Final_status* column and color prophage rows green (value > 0)
   d. Sort by *position* column to restore genome order

4. **Review the results:**

   - Green rows = predicted prophages
   - Red rows = potential phage genes not included in prophages
   - White rows = bacterial genes
   - Look for patterns and evaluate if excluded regions should be prophages

5. **Adjust parameters** if needed based on your assessment

Visual Inspection
^^^^^^^^^^^^^^^^^

Use genome browsers to visualize predictions:

**Artemis:**

.. code-block:: bash

   PhiSpy.py genome.gb -o results --output_choice 258 --color

Then load the output GenBank file in Artemis to see:

- Prophage regions highlighted
- CDS colors based on function
- Gene organization

**Contig Boundaries:**

Pay attention to contig boundaries in draft genomes. Prophages often span contig breaks because:

- They contain repeated sequences
- Assembly algorithms struggle with repeats
- This is a known issue (see `GitHub issue #33 <https://github.com/linsalrob/PhiSpy/issues/33>`_)

Common Patterns
^^^^^^^^^^^^^^^

**True Prophages:**

- Cluster of phage-related genes
- Clear att sites
- Integration at tRNA genes (common)
- Distinct GC content or codon usage
- Size: typically 20-100 kb

**False Positives:**

- Single or few scattered genes
- Genomic islands with mobile element genes
- Ribosomal RNA operons (if ``--phage_genes 0``)
- Very small regions (<10 kb)

**Borderline Cases:**

- Prophage remnants (degraded prophages)
- Satellite prophages
- Gene transfer agents
- Defective prophages

Adjusting Parameters
^^^^^^^^^^^^^^^^^^^^

Based on your assessment:

**Too many predictions:**

.. code-block:: bash

   # Increase stringency
   PhiSpy.py genome.gb -o results --phage_genes 5

**Too few predictions:**

.. code-block:: bash

   # Decrease stringency
   PhiSpy.py genome.gb -o results --phage_genes 1

**Fragmented predictions:**

.. code-block:: bash

   # Allow more non-phage genes between regions
   PhiSpy.py genome.gb -o results --nonprophage_genegaps 15

**Better specificity:**

.. code-block:: bash

   # Use specific metrics and HMM
   PhiSpy.py genome.gb -o results \
       --metrics shannon_slope gc_skew \
       --phmms pVOGs.hmm \
       --phage_genes 3

Example Output
--------------

For *Streptococcus pyogenes* M1 GAS::

   PhiSpy.py -o output_directory -t data/trainSet_160490.61.txt tests/Streptococcus_pyogenes_M1_GAS.gb.gz

**Expected prophage.tsv output:**

.. list-table:: Prophages in S. pyogenes M1 GAS
   :header-rows: 1

   * - Prophage
     - Contig
     - Start
     - Stop
   * - pp_1
     - NC_002737
     - 529631
     - 569288
   * - pp_2
     - NC_002737
     - 778642
     - 820599
   * - pp_3
     - NC_002737
     - 1192630
     - 1222549
   * - pp_4
     - NC_002737
     - 1775862
     - 1782822

Best Practices
--------------

1. **Always generate prophage_information.tsv** for manual review
2. **Use genome browsers** for visual confirmation
3. **Check for att sites** - their presence increases confidence
4. **Consider genome context** - location near tRNAs is common
5. **Validate key predictions** with BLAST or HMM searches
6. **Document your assessment** for reproducibility
7. **Adjust parameters iteratively** based on results
