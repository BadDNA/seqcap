.. _chap:introduction:

****************
Introduction
****************

This `repository <http://github.com/baddna/seqcap>`_ holds computer code (and, temporarily, data) used as part of [McCormack:20??]_.

In essence, the software included within the repository are meant to:

#. identify ultra-conserved regions of genomic DNA (AKA ultra-conserved elements or UCEs) in [vertebrate] organisms
#. design *in silico* (or *in vitro*) "probes" for the enrichment of UCEs
#. align *in silico* probes to genomic DNA of numerous sources and extract those reads + some flanking DNA
#. cluster *in silico* probes + flanking sequence from individual species into "loci" representing entire UCEs.  We can easily translate this step to an *in vitro* approach by clustering reads we derive from target enrichment where the enriching agents are the probes designed above - then we just need cluster these probes + reads into the corresponding UCEs from which they were derived.
#. align, across species, UCE "loci" from above
#. generate gene trees from these inter-species alignments (using PhyML)
#. generate the species trees from individual gene trees (using STAR, STEAC, or MP-EST)

Because of the complex nature of the above steps, there are a number of independent and interdependent programs within this repository, that we hopefully better explain in :ref:`workflow`. The fruits of these programs labor are found in the `Downloads <http://github.com/BadDNA/seqcap/archives/master>`_ as described in the :ref:`data` section. Since many people are mostly interested in these data, we discuss these data first.

.. _dependencies:

Dependencies
============

Hardware Dependencies
*********************

For many of the programs within, you should at least have access to a multi-core computer with a sufficient amount of RAM (>8GB).  A number of the programs allow you to to parallelize jobs using the multiprocessing module of Python, and this provides and almost linear increase in processing speed.

Additionally, several programs (particularly those within phylo/*) require access to a cluster-computing system.  We used a large cluster running the `LSF queuing system <http://www.platform.com/workload-management/high-performance-computing>`_.  

Several of these programs cannot, realistically, be run on very small systems unless you have a very long time to wait.


Software Dependencies
*********************

There are numerous software dependencies that you must satisfy if you wish to run all of the code contained within this repository.  These include:

* A working `Kent Source tree < http://genome-source.cse.ucsc.edu/kent.git>`_.  In particular, you need the following binaries
    * faFrag
    * faToTwoBit
    * ??
* `Lastz <http://www.bx.psu.edu/~rsharris/lastz/>`_
* `Muscle <http://www.drive5.com/muscle>`_
* `PhyML <http://atgc.lirmm.fr/phyml/>`_
* `MP-EST <http://github.com/brantfaircloth/mpest>`_ - we use a modified version of the `original code <http://mpest.googlecode.com>`_
* Python > 2.6.x
    * numpy
    * bx-python
    * biopython (> 1.54)
    * oursql
    * dendropy
    * ??

.. _notes_on_the_code:

Notes on the Code
=========================

The code is available at `http://github.com/baddna/seqcap <http://github.com/baddna/seqcap>`_.

We are in the process of cleaning and standardizing the code, while also increasing the available documentation for the code - both in the source files and here, in the documentation.  Additionally, we have made incremental improvements to a number of routines used herein, and we will be merging these changes into this repository after the initial commit and tagging of files (see :ref:`tagging`).

You will notice, should you scrutinize the code, that some programs write to an `sqlite <http://www.sqlite.org/>`_ database while others write to a `mysql <http://www.mysql.com/>`_ database. The reasons for this additional level of complexity are several-fold. Generally speaking, we started using sqlite as the initial database for holding data generated as part of this project, but we moved to mysql when we decided that we needed a greater level of concurrency (sqlite does not support concurrent **writes**.

As we incrementally improve individual files, we will make the switch to mysql-only (or `ORM <http://en.wikipedia.org/wiki/Object-relational_mapping>`_, which should be more generic) support.

In the development version of the code (currently a private repository - see :ref:`additional_notes`) we have also replaced slower external dependencies (e.g. `BLAST <http://blast.ncbi.nlm.nih.gov/Blast.cgi>`_) with similar, yet speedier alternatives (e.g. `lastz <http://www.bx.psu.edu/~rsharris/lastz/>`_)

.. _additional_notes:

Additional Notes
==================

#. we have moved the code within this repository here from a private repository that I (BCF) maintain for development. You should generally be happy about this, because it has allowed us to do a fair amount of housekeeping. It also allows us to work on some things that we may not want everyone to know about (yet!).  However, if you believe a program is missing that may be in this private repository, please let me know, and I'll attempt to move it over.  **Eventually, we plan to move all code to this repository, and remove the private repository**.  At that point, we will develop from within this repository.  The downside of this approach is that we will lose some history information of particular pieces of code.

#. we have an updated workflow for a number of the steps detailed below (particularly the initial steps of UCE location and probe design). When the time comes, we will tag pertinent files in the current repo (see :ref:`tagging`), and then move in the new bits.

#. some of the methods/code within are likely confusing, particularly if you are trying to piece together what we did without actually reading the code. For the most part, we'll try to give you some guidance, but you'll also need to read the code. It may be helpful to enlist someone with knowledge of Python to aid this process.

.. _tagging:

Tagging
================

Because code is a moving target, we have tagged particular commits that contain files of a certain vintage and/or purpose.  For instance, the `0.1 <http://github.com/baddna/seqcap/tree/v0.1>`_ tag holds the initial version of the code that we use to search MAF files and run parallel blast jobs.  Subsequently, we updated these files to work more efficiently (e.g. the `multiprocessing <http://docs.python.org/library/multiprocessing.html>`_ module in place of `pp <http://www.parallelpython.com/>`_).