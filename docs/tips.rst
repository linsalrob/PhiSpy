Tips and Tricks
===============

This page contains helpful tips, troubleshooting advice, and solutions to common problems.

PATH Issues
-----------

Command Not Found
^^^^^^^^^^^^^^^^^

If you get this error::

   $ PhiSpy.py -v
   -bash: PhiSpy.py: command not found

**Solution 1: Use full path**

.. code-block:: bash

   ~/.local/bin/PhiSpy.py -v

**Solution 2: Add to PATH** (recommended)

.. code-block:: bash

   echo "export PATH=\$HOME/.local/bin:\$PATH" >> ~/.bashrc
   source ~/.bashrc
   PhiSpy.py -v

Installation Tips
-----------------

Quick Installation
^^^^^^^^^^^^^^^^^^

The simplest installation (requires sudo):

.. code-block:: bash

   sudo apt install -y python3-pip
   python3 -m pip install --user phispy

Note: ``python3-pip`` automatically installs ``build-essential`` and ``python3-dev``.

Verify Installation
^^^^^^^^^^^^^^^^^^^

Check that PhiSpy is installed correctly:

.. code-block:: bash

   PhiSpy.py --version
   PhiSpy.py --list short

Running PhiSpy
--------------

Handling Large Genomes
^^^^^^^^^^^^^^^^^^^^^^

For large or draft genomes:

1. **Filter small contigs**::

      PhiSpy.py genome.gb -o results --min_contig_size 10000

2. **Process contigs separately** if memory is an issue

3. **Use appropriate resources**: 
   
   - Memory: ~2-4 GB per genome typically
   - CPU: Use ``--threads`` for HMM searches

Working with Draft Genomes
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Draft genomes often have prophages spanning contig breaks.

**Best practices:**

- Pay attention to prophages at contig ends
- Consider using assembly tools to close gaps
- Review contig boundaries manually
- Expect some fragmented predictions

Gzip Files
^^^^^^^^^^

PhiSpy handles gzip automatically::

   # These all work
   PhiSpy.py genome.gb -o results
   PhiSpy.py genome.gb.gz -o results
   
   # Output matches input format
   # If input is .gz, output files are also .gz

Multiple Genomes
^^^^^^^^^^^^^^^^

Process multiple genomes efficiently:

**Method 1: Use file prefixes**

.. code-block:: bash

   for genome in *.gb; do
       name=$(basename "$genome" .gb)
       PhiSpy.py "$genome" -o results -p "${name}_"
   done

**Method 2: Separate output directories**

.. code-block:: bash

   for genome in *.gb; do
       name=$(basename "$genome" .gb)
       PhiSpy.py "$genome" -o "results_${name}"
   done

Parameter Tuning
----------------

Finding the Right Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^

Start with defaults and adjust based on results:

.. code-block:: bash

   # 1. Run with defaults
   PhiSpy.py genome.gb -o results_default
   
   # 2. Review output
   less results_default/prophage_coordinates.tsv
   
   # 3. Adjust if needed
   # Too many predictions:
   PhiSpy.py genome.gb -o results_strict --phage_genes 5
   
   # Too few predictions:
   PhiSpy.py genome.gb -o results_relaxed --phage_genes 0

Systematic Testing
^^^^^^^^^^^^^^^^^^

Test multiple parameter combinations:

.. code-block:: bash

   for pg in 0 1 3 5; do
       PhiSpy.py genome.gb -o results_pg${pg} --phage_genes $pg
   done
   
   # Compare results
   wc -l results_pg*/prophage_coordinates.tsv

Using Jupyter Notebook
^^^^^^^^^^^^^^^^^^^^^^

The `interactive notebook <https://github.com/linsalrob/PhiSpy/blob/master/jupyter_notebooks/PhiSpy.ipynb>`_ is excellent for parameter exploration:

1. Clone the repository
2. Install Jupyter: ``pip install jupyter``
3. Run: ``jupyter notebook jupyter_notebooks/PhiSpy.ipynb``
4. Adjust parameters interactively

Performance
-----------

Speed Up Processing
^^^^^^^^^^^^^^^^^^^

1. **Use multiple threads**::

      PhiSpy.py genome.gb -o results --threads 8

2. **Skip HMM search on reruns**::

      PhiSpy.py genome.gb -o results --phmms pVOGs.hmm --skip_search

3. **Use appropriate training sets** (faster than generic)

4. **Filter small contigs** (reduces processing time)

Parallel Processing
^^^^^^^^^^^^^^^^^^^

Process multiple genomes in parallel:

.. code-block:: bash

   # Using GNU parallel
   parallel -j 4 'PhiSpy.py {} -o results/{/.} --threads 2' ::: genomes/*.gb
   
   # Or simple background jobs
   for genome in genomes/*.gb; do
       name=$(basename "$genome" .gb)
       PhiSpy.py "$genome" -o "results_${name}" --threads 2 &
   done
   wait

Debugging
---------

Verbose Output
^^^^^^^^^^^^^^

For detailed logging:

.. code-block:: bash

   PhiSpy.py genome.gb -o results --log phispy_detailed.log

Review the log file for:

- Which metrics were used
- Training set loaded
- Number of genes/contigs processed
- Decisions made during prediction

Keep Temporary Files
^^^^^^^^^^^^^^^^^^^^

Preserve intermediate files for inspection:

.. code-block:: bash

   PhiSpy.py genome.gb -o results --keep

This keeps all temporary files for debugging.

Evaluation Mode
^^^^^^^^^^^^^^^

Re-run evaluation without reprocessing:

.. code-block:: bash

   PhiSpy.py genome.gb -o results --evaluate

Useful when testing different thresholds.

Common Issues
-------------

No Prophages Found
^^^^^^^^^^^^^^^^^^

**Possible causes:**

1. Genome truly has no prophages
2. Parameters too strict
3. Poor annotation quality
4. Contigs too small

**Solutions:**

.. code-block:: bash

   # Try relaxed parameters
   PhiSpy.py genome.gb -o results --phage_genes 0
   
   # Check for annotation
   grep -c "CDS" genome.gb
   
   # Verify contig sizes
   PhiSpy.py genome.gb -o results --min_contig_size 1000

Too Many False Positives
^^^^^^^^^^^^^^^^^^^^^^^^^

**Symptoms:**

- Many small predictions
- Ribosomal operons called as prophages
- Pathogenicity islands included

**Solutions:**

.. code-block:: bash

   # Increase stringency
   PhiSpy.py genome.gb -o results --phage_genes 5
   
   # Use specific metrics
   PhiSpy.py genome.gb -o results --metrics shannon_slope gc_skew
   
   # Use HMM database
   PhiSpy.py genome.gb -o results --phmms pVOGs.hmm

Fragmented Predictions
^^^^^^^^^^^^^^^^^^^^^^

**Symptom:** One prophage split into multiple predictions.

**Solution:**

.. code-block:: bash

   # Allow more non-phage genes between regions
   PhiSpy.py genome.gb -o results --nonprophage_genegaps 20

Memory Issues
^^^^^^^^^^^^^

**For very large genomes:**

1. Filter contigs::

      PhiSpy.py genome.gb -o results --min_contig_size 10000

2. Reduce threads::

      PhiSpy.py genome.gb -o results --threads 1

3. Process on high-memory machine

4. Split genome into smaller files

File Format Issues
------------------

GenBank Format Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PhiSpy expects:

- Valid GenBank format
- CDS features with locations
- DNA sequence included

**Verify format:**

.. code-block:: bash

   # Check for required elements
   grep "LOCUS" genome.gb
   grep "CDS" genome.gb
   grep "ORIGIN" genome.gb

Converting Formats
^^^^^^^^^^^^^^^^^^

From GFF + FASTA to GenBank:

.. code-block:: bash

   # Using BioPython
   python -c "from BCBio import GFF; from Bio import SeqIO; \
   SeqIO.write(GFF.parse('genome.gff', 'genome.fasta'), 'genome.gb', 'genbank')"

From EMBL to GenBank:

.. code-block:: bash

   # Using BioPython
   python -c "from Bio import SeqIO; \
   SeqIO.convert('genome.embl', 'embl', 'genome.gb', 'genbank')"

Quality Control
---------------

Validating Predictions
^^^^^^^^^^^^^^^^^^^^^^

Always validate key predictions:

1. **BLAST search** prophage genes
2. **Check for att sites** (strong evidence)
3. **Verify integration sites** (often at tRNAs)
4. **Compare to known prophages** in related species
5. **Check GC content** (often different from host)

Comparing Versions
^^^^^^^^^^^^^^^^^^

When testing parameters::

   # Generate comparable outputs
   PhiSpy.py genome.gb -o v1 --phage_genes 1
   PhiSpy.py genome.gb -o v2 --phage_genes 3
   
   # Compare
   diff v1/prophage_coordinates.tsv v2/prophage_coordinates.tsv

Batch Analysis
--------------

Create Summary Statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # Count prophages per genome
   for dir in results_*/; do
       count=$(wc -l < "$dir/prophage_coordinates.tsv")
       echo "$(basename $dir): $count prophages"
   done

Merge Results
^^^^^^^^^^^^^

.. code-block:: bash

   # Combine all prophage coordinates
   echo -e "Genome\tProphage\tContig\tStart\tStop" > all_prophages.tsv
   for dir in results_*/; do
       genome=$(basename $dir)
       awk -v g="$genome" '{print g"\t"$0}' "$dir/prophage_coordinates.tsv" >> all_prophages.tsv
   done

Best Practices Summary
----------------------

1. **Start simple**: Use default parameters first
2. **Review carefully**: Always manually inspect results
3. **Use HMM databases**: Improves accuracy significantly
4. **Choose appropriate training sets**: Use organism-specific when available
5. **Document parameters**: Record what you used for reproducibility
6. **Validate predictions**: Use independent evidence
7. **Iterate**: Adjust parameters based on results
8. **Keep logs**: They're invaluable for troubleshooting

Getting Help
------------

If you encounter issues:

1. **Check this documentation** thoroughly
2. **Review the log file** for error messages
3. **Search GitHub issues**: https://github.com/linsalrob/PhiSpy/issues
4. **Open a new issue** if your problem is novel:
   
   - Include PhiSpy version
   - Describe the problem clearly
   - Provide example data if possible
   - Include command used and error message

5. **Join the community**: Participate in discussions on GitHub

Useful Resources
----------------

- **GitHub Repository**: https://github.com/linsalrob/PhiSpy
- **Original Paper**: https://doi.org/10.1093/nar/gks406
- **RAST Annotation**: http://rast.nmpdr.org/
- **PROKKA Annotation**: https://github.com/tseemann/prokka
- **pVOG Database**: http://dmk-brain.ecn.uiowa.edu/pVOGs
- **VOGdb Database**: http://vogdb.org/
- **Artemis Genome Browser**: https://sanger-pathogens.github.io/Artemis/
