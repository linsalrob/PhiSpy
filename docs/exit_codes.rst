Exit Codes
==========

PhiSpy uses specific exit codes to indicate different error conditions. This helps with scripting and debugging.

Exit Code Reference
-------------------

.. list-table:: PhiSpy Exit Codes
   :header-rows: 1
   :widths: 10 30 60

   * - Code
     - Meaning
     - Suggested Solution
   * - 0
     - Success
     - Normal completion
   * - 2
     - No input file provided
     - Provide a GenBank file as input
   * - 3
     - No output directory provided
     - Specify output directory with ``-o`` flag
   * - 10
     - No training sets available
     - Check PhiSpy installation
   * - 11
     - Specific training set not available
     - Verify training set name with ``--list``
   * - 13
     - No kmers file found
     - Check PhiSpy installation
   * - 20
     - IO Error
     - Check input file exists and is readable
   * - 25
     - Non-nucleotide base found
     - Check for non-standard bases in sequence
   * - 26
     - An ORF with no bases
     - Remove very short ORFs from annotation
   * - 30
     - No contigs
     - Reduce ``--min_contig_size`` parameter
   * - 40
     - No ORFs in GenBank file
     - Annotate genome with RAST or PROKKA
   * - 41
     - Less than 100 ORFs
     - Genome needs more complete annotation

Detailed Descriptions
---------------------

Exit Code 0: Success
^^^^^^^^^^^^^^^^^^^^

**Meaning:** PhiSpy completed successfully.

**Action:** Review output files in the specified output directory.

Exit Code 2: No Input File
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** No GenBank file was provided.

**Example error:**

.. code-block:: bash

   $ PhiSpy.py -o results
   Error: No input file provided

**Solution:**

.. code-block:: bash

   PhiSpy.py genome.gb -o results

Exit Code 3: No Output Directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** Output directory was not specified.

**Example error:**

.. code-block:: bash

   $ PhiSpy.py genome.gb
   Error: No output directory provided

**Solution:**

.. code-block:: bash

   PhiSpy.py genome.gb -o results

Exit Code 10: No Training Sets Available
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** PhiSpy cannot find any training sets. This indicates a problem with the installation.

**Troubleshooting:**

1. **Verify installation:**

   .. code-block:: bash

      PhiSpy.py --list short

2. **Reinstall PhiSpy:**

   .. code-block:: bash

      pip install --force-reinstall PhiSpy

3. **Check installation directory:**

   .. code-block:: bash

      python -c "import PhiSpyModules; print(PhiSpyModules.__file__)"
      # Should show path to PhiSpyModules
      
      ls -la $(python -c "import PhiSpyModules, os; print(os.path.dirname(PhiSpyModules.__file__))")/data/

Exit Code 11: Training Set Not Available
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** The specified training set was not found.

**Example error:**

.. code-block:: bash

   $ PhiSpy.py genome.gb -o results -t data/trainSet_Nonexistent.txt
   Error: Training set not found

**Solution:**

1. **List available training sets:**

   .. code-block:: bash

      PhiSpy.py --list short

2. **Use a valid training set:**

   .. code-block:: bash

      PhiSpy.py genome.gb -o results -t data/trainSet_Streptococcus.txt

3. **Use full path if needed:**

   .. code-block:: bash

      PhiSpy.py genome.gb -o results -t /full/path/to/trainSet_Custom.txt

Exit Code 13: No K-mers File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** K-mers file required for Shannon score calculation is missing.

**Troubleshooting:**

Similar to exit code 10 - indicates installation problem. Reinstall PhiSpy.

Exit Code 20: IO Error
^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** Error reading the input file.

**Common causes:**

- File doesn't exist
- File is corrupted
- Insufficient permissions
- Wrong file format

**Solutions:**

1. **Verify file exists:**

   .. code-block:: bash

      ls -lh genome.gb

2. **Check file format:**

   .. code-block:: bash

      head -n 20 genome.gb
      # Should show GenBank format with LOCUS, DEFINITION, etc.

3. **Check permissions:**

   .. code-block:: bash

      chmod +r genome.gb

4. **Verify it's a valid GenBank file:**

   .. code-block:: bash

      grep "LOCUS" genome.gb

Exit Code 25: Non-Nucleotide Base
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** The sequence contains a non-standard nucleotide.

**Common causes:**

- Ambiguous bases beyond IUPAC codes
- Protein sequences in nucleotide field
- Corrupted file

**Solutions:**

1. **Identify the problematic base:**

   Check log file for details about which base and location.

2. **Clean the sequence:**

   Replace ambiguous codes with N or remove them.

3. **Re-annotate the genome:**

   If the original FASTA had issues, re-annotate with corrected sequence.

Exit Code 26: ORF with No Bases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** An ORF annotation has no associated sequence.

**Common causes:**

- Very short ORF (< 3 bases)
- Annotation error
- Corrupted GenBank file

**Solutions:**

1. **Remove short ORFs:**

   Edit the GenBank file to remove extremely short CDS features.

2. **Re-annotate:**

   Use PROKKA or RAST to generate new annotations.

3. **Fix annotation:**

   Manually correct the problematic CDS feature coordinates.

Exit Code 30: No Contigs
^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** No contigs passed the minimum size filter.

**Example error:**

.. code-block:: bash

   $ PhiSpy.py small_genome.gb -o results
   Error: No contigs meet minimum size requirement

**Solutions:**

1. **Reduce minimum contig size:**

   .. code-block:: bash

      PhiSpy.py genome.gb -o results --min_contig_size 1000

2. **Check contig sizes:**

   .. code-block:: bash

      grep "LOCUS" genome.gb

3. **Note:** Very small contigs (< 5000 bp) are unlikely to contain complete prophages.

Exit Code 40: No ORFs
^^^^^^^^^^^^^^^^^^^^^

**Meaning:** The GenBank file contains no CDS features.

**Example error:**

.. code-block:: bash

   $ PhiSpy.py unannotated_genome.gb -o results
   Error: No ORFs found in GenBank file

**Solutions:**

1. **Annotate the genome:**

   **Using PROKKA:**

   .. code-block:: bash

      prokka --outdir prokka_output genome.fasta
      PhiSpy.py prokka_output/genome.gbk -o results

   **Using RAST:**

   - Upload genome to http://rast.nmpdr.org/
   - Download annotated GenBank file
   - Run PhiSpy on annotated file

2. **Verify file has CDS features:**

   .. code-block:: bash

      grep -c "CDS" genome.gb

Exit Code 41: Insufficient ORFs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Meaning:** The genome has fewer than 100 ORFs. This is insufficient for reliable prophage prediction.

**Common causes:**

- Incomplete annotation
- Very small genome
- Plasmid instead of chromosome
- Partial genome sequence

**Solutions:**

1. **Re-annotate with better tools:**

   .. code-block:: bash

      prokka --force --outdir prokka_output genome.fasta

2. **Check if this is a complete genome:**

   - Is the sequence complete?
   - Is it a plasmid?
   - Is annotation adequate?

3. **For small genomes/plasmids:**

   PhiSpy is designed for bacterial chromosomes (typically > 1 Mbp with thousands of genes). It may not work well on very small genomes.

Handling Exit Codes in Scripts
-------------------------------

Bash Script Example
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   #!/bin/bash
   
   PhiSpy.py genome.gb -o results
   exit_code=$?
   
   case $exit_code in
       0)
           echo "Success! Processing results..."
           ;;
       2)
           echo "Error: No input file provided"
           exit 1
           ;;
       3)
           echo "Error: No output directory provided"
           exit 1
           ;;
       40)
           echo "Error: No ORFs found. Running annotation..."
           prokka --outdir prokka_out genome.fasta
           PhiSpy.py prokka_out/*.gbk -o results
           ;;
       *)
           echo "PhiSpy failed with exit code $exit_code"
           exit $exit_code
           ;;
   esac

Python Script Example
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   import subprocess
   import sys
   
   result = subprocess.run(
       ['PhiSpy.py', 'genome.gb', '-o', 'results'],
       capture_output=True
   )
   
   if result.returncode == 0:
       print("PhiSpy completed successfully")
   elif result.returncode == 40:
       print("No ORFs found. Genome needs annotation.")
       # Run annotation
       subprocess.run(['prokka', '--outdir', 'prokka_out', 'genome.fasta'])
       # Retry PhiSpy
       subprocess.run(['PhiSpy.py', 'prokka_out/genome.gbk', '-o', 'results'])
   else:
       print(f"PhiSpy failed with exit code {result.returncode}")
       print(result.stderr.decode())
       sys.exit(result.returncode)

Workflow Integration
--------------------

Nextflow Example
^^^^^^^^^^^^^^^^

.. code-block:: groovy

   process phispy {
       errorStrategy { task.exitStatus == 30 ? 'ignore' : 'terminate' }
       
       input:
       path genome
       
       output:
       path "results/*", optional: true
       
       script:
       """
       PhiSpy.py ${genome} -o results --min_contig_size 1000
       """
   }

Snakemake Example
^^^^^^^^^^^^^^^^^

.. code-block:: python

   rule phispy:
       input:
           "genomes/{genome}.gb"
       output:
           "results/{genome}/prophage_coordinates.tsv"
       params:
           outdir="results/{genome}"
       log:
           "logs/{genome}.log"
       shell:
           """
           PhiSpy.py {input} -o {params.outdir} 2> {log} || \
           {{ if [ $? -eq 30 ]; then echo "No contigs" > {log}; exit 0; fi; exit 1; }}
           """

Debugging Tips
--------------

When encountering an exit code error:

1. **Check the log file:**

   .. code-block:: bash

      PhiSpy.py genome.gb -o results --log detailed.log
      cat detailed.log

2. **Run with verbose output:**

   .. code-block:: bash

      PhiSpy.py genome.gb -o results 2>&1 | tee phispy_output.txt

3. **Validate input file:**

   .. code-block:: bash

      # For GenBank files
      python -c "from Bio import SeqIO; list(SeqIO.parse('genome.gb', 'genbank'))"

4. **Check installation:**

   .. code-block:: bash

      PhiSpy.py --version
      PhiSpy.py --list short

5. **Report issues:**

   If you believe the exit code is erroneous, report it at:
   https://github.com/linsalrob/PhiSpy/issues
