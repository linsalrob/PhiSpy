Installation
============

PhiSpy can be installed using several methods. We recommend using either conda (easiest) or pip.

Conda (Recommended)
-------------------

The easiest way to install PhiSpy for all users is to use ``bioconda``:

.. code-block:: bash

   conda install -c bioconda phispy

This method handles all dependencies automatically and works across different platforms.

PIP
---

To install PhiSpy using pip, you need a C++ compiler and the Python header files. 

On Ubuntu/Debian systems:

.. code-block:: bash

   sudo apt install -y build-essential python3-dev python3-pip
   python3 -m pip install --user PhiSpy

On macOS (requires Xcode Command Line Tools):

.. code-block:: bash

   xcode-select --install
   python3 -m pip install --user PhiSpy

The ``pip`` installation will place ``PhiSpy.py`` in ``~/.local/bin`` which should be in your ``$PATH``. If it's not, see the `Tips and Tricks <tips.html#path-issues>`_ section.

Advanced Installation
---------------------

For advanced users, you can clone the git repository and install from source:

.. code-block:: bash

   git clone https://github.com/linsalrob/PhiSpy.git
   cd PhiSpy
   python3 -m pip install --user .

If you want to track installed files for easy uninstallation:

.. code-block:: bash

   git clone https://github.com/linsalrob/PhiSpy.git
   cd PhiSpy
   python3 setup.py install --user --record installed_files.txt

To uninstall later:

.. code-block:: bash

   cat installed_files.txt | xargs rm -f

Global Installation (Requires Root)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have root access and want to install globally:

.. code-block:: bash

   git clone https://github.com/linsalrob/PhiSpy.git
   cd PhiSpy
   python3 -m pip install .

Software Requirements
---------------------

PhiSpy requires the following programs to be installed. Most of these are installed automatically when using conda or pip:

1. **Python** - version 3.8 or later
2. **Biopython** - version 1.74 or later
3. **numpy** - version 1.16.0 or later
4. **scikit-learn** - version 0.21.3 or later
5. **bcbio-gff** - version 0.6.6 or later
6. **gcc** - GNU project C and C++ compiler - version 4.4.1 or later
7. **Python.h** header file - included in ``python3-dev`` package

Optional Dependencies
---------------------

For HMM-based prophage detection:

- **HMMER3** - for using pVOGs or other HMM databases

Verifying Installation
----------------------

To verify that PhiSpy is installed correctly:

.. code-block:: bash

   PhiSpy.py --version

You should see the version number printed. If you get a "command not found" error, see the `Tips and Tricks <tips.html#path-issues>`_ section.
