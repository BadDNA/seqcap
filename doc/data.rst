.. _data:

***********
Data
***********

We make a number of data files available (temporarily) on this site.  The data files are as follows.

uce_probe.sqlite.bz2
======================

This is a bzip archive of an sqlite database containing the data we used during our probe design process and the probes that results from this process.  The tables are:

* ``cons`` - this table contains matches that we identified in alignments
  (`MAF <http://genome.ucsc.edu/FAQ/FAQformat#format5>`_) between chicken and
  lizard using the initial version (0.1) of the `genome/summary.py`code in
  this repository (see :ref:`tagging` below).

* ``blast`` - this table contains results from a blast of the matches from cons onto the zebra finch genome.  For these matches we used the 0.1 tagged `genome/summaryBlast.py` code

* ``gallus_refseq`` - this is a table of the refseq genes in galGal3.  it would likely be best to download a recent version of these data from NCBI or UCSC rather than using what is here for anything.
    
* ``probes`` - this is a table of the probes that we designed from the UCEs we
  located in chicken, lizard, and zebra finch.  We designed these probes using
  `design/sure_select_tiler.py`

uce_bed.bz2
==============

This is a `BED <http://genome.ucsc.edu/FAQ/FAQformat#format1>`_ format file of the UCEs we located that are shared between chicken, lizard, and zebra finch. The locations within this file are relative to chicken (galGal3).

uce_probes.bed.bz2
===================

This is a `BED <http://genome.ucsc.edu/FAQ/FAQformat#format1>`_ format file of
the locations of the probes we designed.  The locations within this file are
relative to chicken (galGal3).

db1.sql.bz2
==========================================

This is a bzipped dump file of the `mysql <http://www.mysql.com/>`_ database 
in which we stored our alignments of each probe sequence to other genomes.
Initially and for error checking purposes, we conducted the align, clean, and
insert steps by hand (as seen in :ref:`Simulation:steps`).  This process has 
now largely been automated in :ref:`future:run_mutiple_lastz.py`.  The tables are:

.. _group_data_tables:

group data tables
-------------------------------

* `all_orgs` - this table contains data from all of the matches to individual
  species listed in :ref:`individual_species`
* `cons` - this table is identical to the `cons` table from uce_probe.sqlite.  it
  is repeated here for ease of use
* `group_XAB_5` - this contains the list of probes found in all members of 

    * monDom5
    * loxAfr3
    * choHof1
    * hg19
    * mm9
    
* `group_Bats_5` - this contains the list of probes found in all members of

    * eriEur1
    * equCab2
    * bosTau4
    * canFam2
    * pteVam1
    
* `group_Bats_7` - this contains the list of probes found in all members of

    * eriEur1
    * equCab2
    * bosTau4
    * canFam2
    * pteVam1
    * myoLuc1
    * ailMel1

* `group_Bats_8` - this contains the list of probes found in all members of

    * eriEur1
    * pteVam1
    * myoLuc1
    * bosTau4
    * vicPac1
    * equCab2
    * canFam2
    * ailMel1

* `group_Bats_12` - this contains the list of probes found in all members of

    * eriEur1
    * pteVam1
    * myoLuc1
    * bosTau4
    * vicPac1
    * susScr9
    * oviAri1
    * turTru1
    * equCab2
    * canFam2
    * felCat3
    * ailMel1

* `group_elephants_7` - this contains the list of probes found in all members of
    * loxAfr3
    * canFam2
    * echTel1
    * choHof1
    * dasNov2
    * hg19
    * monDom5

* `group_size_19` - this contains the list of probes found in all members of

    * anoCar2
    * bosTau4
    * calJac3
    * canFam2
    * cavPor3
    * chinese
    * equCab2
    * gorGor3
    * hg19
    * korean
    * loxAfr3
    * mm9
    * monDom5
    * oryCun2
    * panTro2
    * ponAbe2
    * rheMac2
    * taeGut1
    * venter

* `group_size_23` - this contains the list of probes found in all memebrs of

    * bosTau4
    * canFam2
    * cavPor3
    * chinese
    * dasNov2
    * echTel1
    * equCab2
    * gorGor3
    * hg19
    * korean
    * loxAfr3
    * mm9
    * monDom5
    * ochPri2
    * oryCun2
    * panTro2
    * ponAbe2
    * pteVam1
    * rn4
    * speTri1
    * tarSyr1
    * tupBel1
    * venter

* `group_size_25` - this contains the list of probes found in all members of

    * anoCar2
    * bosTau4
    * calJac3
    * canFam2
    * cavPor3
    * chinese
    * dipOrd1
    * equCab2
    * gorGor3
    * hg19
    * korean
    * loxAfr3
    * mm9
    * monDom5
    * ornAna1
    * oryCun2
    * panTro2
    * ponAbe2
    * pteVam1
    * rheMac2
    * rn4
    * taeGut1
    * tarSyr1
    * venter
    * vicPac1

* `group_size_29` - this contains the list of probes found in all members of

    * anoCar2
    * bosTau4
    * calJac3
    * canFam2
    * cavPor3
    * chinese
    * choHof1
    * dipOrd1
    * echTel1
    * equCab2
    * eriEur1
    * gorGor3
    * hg19
    * korean
    * loxAfr3
    * mm9
    * monDom5
    * ornAna1
    * oryCun2
    * panTro2
    * ponAbe2
    * pteVam1
    * rheMac2
    * rn4
    * taeGut1
    * tarSyr1
    * tupBel1
    * venter
    * vicPac1

* `probe_distribution` - this contains a binary "matrix" indicating
  presence/absence (1/0) of probe matches by species
* `probes` - this is a table providing the ids of the probes we designed
* `species` - this table provides information on the genome build of each
  organism to which we aligned probes
* `sureselect` - this table is identical to the `probes` table from
  uce_probe.sqlite

.. _individual_species:

probe matches to individual species
------------------------------------ 

These tables provide the individual matches of probes to difference organisms, as found by using lastz to align our probe sequences to the individual genomes.  We provide the build version, name, etc. of each genome sequence in the *species* table of the database dump described in :ref:`group_data_tables`.

* `ailMel1`
* `anoCar2`
* `bosTau4`
* `calJac3`
* `canFam2`
* `cavPor3`
* `chinese`
* `choHof1`
* `danRer6`
* `dasNov2`
* `dipOrd1`
* `echTel1`
* `equCab2`
* `eriEur1`
* `felCat3`
* `gasAcu1`
* `gorGor3`
* `hg19`
* `korean`
* `loxAfr3`
* `macEug1`
* `micMur1`
* `mm9`
* `monDom5`
* `myoLuc1`
* `ochPri2`
* `ornAna1`
* `oryCun2`
* `otoGar1`
* `oviAri1`
* `panTro2`
* `ponAbe2`
* `proCap1`
* `pteVam1`
* `rheMac2`
* `rn4`
* `sorAra1`
* `speTri1`
* `susScr9`
* `taeGut1`
* `tarSyr1`
* `tetNig1`
* `tupBel1`
* `turTru1`
* `venter`
* `vicPac1`
* `xenTro2`