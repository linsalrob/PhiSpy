PhiSpy Documentation
====================

.. image:: https://img.shields.io/badge/Bioinformatics-EdwardsLab-03A9F4
   :target: https://edwards.flinders.edu.au/
   :alt: Edwards Lab

.. image:: https://www.zenodo.org/badge/60999054.svg
   :target: https://www.zenodo.org/badge/latestdoi/60999054
   :alt: DOI

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

.. image:: https://img.shields.io/pypi/pyversions/phispy.svg?style=flat-square&label=PyPi%20Versions
   :target: https://pypi.org/project/PhiSpy/
   :alt: PyPI

.. image:: https://img.shields.io/conda/dn/bioconda/phispy.svg?style=flat-square&label=BioConda%20install
   :target: https://anaconda.org/bioconda/phispy
   :alt: BioConda Install

What is PhiSpy?
---------------

PhiSpy identifies prophages in Bacterial (and probably Archaeal) genomes. Given an annotated genome, it will use several approaches to identify the most likely prophage regions.

PhiSpy was originally written by:

* Sajia Akhter (sajia@stanford.edu)
* Edwards Bioinformatics Lab (http://edwards.sdsu.edu/research/)

Improvements, bug fixes, and other changes were made by:

* Katelyn McNair - Edwards Bioinformatics Lab
* Przemysław Decewicz - DEMB at the University of Warsaw

Citation
--------

If you use PhiSpy in your work, please cite:

Sajia Akhter, Ramy K. Aziz, Robert A. Edwards; PhiSpy: a novel algorithm for finding prophages in bacterial genomes that combines similarity- and composition-based strategies. *Nucleic Acids Research* 2012; 40 (16): e126. doi: `10.1093/nar/gks406 <https://doi.org/10.1093/nar/gks406>`_

For the newest versions, please also cite:

McNair, K., Decewicz, P., Akhter, S., Aziz, R.K., Daniel, S., Edwards, R.A. 2019. PhiSpy. https://github.com/linsalrob/PhiSpy/ doi://10.5281/zenodo.3475717

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   installation
   quickstart
   usage
   output
   training
   tips

.. toctree::
   :maxdepth: 2
   :caption: Reference

   parameters
   exit_codes

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
