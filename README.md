# seqcap - initial thoughts

This repository holds computer code used as part of McCormack et al. XXXX.

[CITE]

There are a very large number of independent programs and interdependent
programs within.  You will notice, likely at first glance, that some of the
code seems all over the place (style-wise).  That's because, to some degree, it
is.  We began this project during 2008, and it has stretched to the present.

There are several programs that were use-once-and-forget, and others that
became indispensable or are newer/prettier/better/etc. I would say, too, that
you can see some evolution in the code itself.  About 3/4 of the way into the
project, I (BCF) also started to better follow
[PEP8](http://www.python.org/dev/peps/pep-0008/) which you'll also likely
notice.  Nick Crawford (NGC) was better about following PEP8 than I.

You will also notice, should you scrutinize the code, that some programs write
to an sqlite database while other write to a mysql database.  The reasons for
this additional level of complexity are several-fold.  Generally speaking, we
started using sqlite as the initial database for holding data generated as
part of this project, but we moved to mysql when demands for concurrency
required that we use a database supporting concurrent writes (sqlite does
not).

Three additional notes:  

1. we have moved the code within this repository here from a private repository
   that I (BCF) maintain for the development portions of this project.  You
   should generally be happy about this, beacuse it has allowed use to do
   a fair amount of housekeeping.  If you believe a program is missing that
   may be in this private repository, please let me know, and I'll attempt to
   move it over.

2. we have an updated workflow for a number of the steps detailed below
   (particularly the initial steps of UCE location and probe design.  when the
   time comes, we will tag pertinent files in the current repo, and then move
   in the new bits.

3. some of the methods/code within are likely confusing to others, particularly
   if you are trying to piece together what we did without actually reading the
   code.  For the most part, we'll try to give you some guidance, but
   you'll also need to read the code.  It may be helpful to enlist someone
   with knowledge of [Python](http://www.python.org) to aid this process.


## code

The code is available at http://github.com/BadDNA/seqcap/.  This file is the
top-level README.


## workflow

We provide a general overview of the overall workflow in WORKFLOW.md.  We also
provide task-specific workflows in:

* x
* y
* z


## data

We make a number of data files available (at present) on this site.  In
essence, these data are the fruits of our labor.  the data files are as
follows.

### uce_probe.sqlite.bz2

This is an sqlite database containing the data that are central to our probe
design process.  the tables are:

* `cons` - this table contains matches that we identified in alignments
  ([MAF](http://genome.ucsc.edu/FAQ/FAQformat#format5)) between chicken and
  lizard using the initial version (0.1) of the `genome/summary.py`code in this
  repository (see "tags" below).
* `blast` - this table contains results from a blast of the matches from cons
  onto the zebra finch genome.  For these matches we used the 0.1 tagged
  `genome/summaryBlast.py` code
* `gallus_refseq` - this is a table of the refseq genes in galGal3.  it would
  likely be best to download a recent version of these data from NCBI or UCSC
  rather than using what is here for anything.
* `probes` - this is a table of the probes that we designed from the UCEs we
  located in chicken, lizard, and zebra finch.  We designed these probes using
  `design/sure_select_tiler.py

### uce_bed.bz2

This is a [BED](http://genome.ucsc.edu/FAQ/FAQformat#format1) format file of
the UCEs we located that are shared between chicken, lizard, and zebra finch.
The locations within this file are relative to chicken (galGal3).

### uce_probes.bed.bz2

This is a [BED](http://genome.ucsc.edu/FAQ/FAQformat#format1) format file of
the locations of the probes we designed.  The locations within this file are
relative to chicken (galGal3).

### probe_matches_to_other_genomes.sql.bz2

This is a bzipped dump file of the [mysql](http://www.mysql.com/) database in
which we stored our alignments of each probe sequence to other genomes.
Initially and for error checking purposes, we conducted the align, clean, and
insert steps by hand (as seen in Simulation/STEPS.md).  This process has now
largely been automated in `Future/run_mutiple_lastz.py`.  the tables are:

#### organizational/metadata tables

* `all_orgs` - this table contains data from all of the matches to individual
  species listed below
* `cons` - this table is identical to the `cons` table from uce_probe.sqlite.  it
  is repeated here for ease of use
* `group_XAB_5` - this contains the list of probes found in all members of
  [monDom5, loxAfr3, choHof1, hg19, mm9]
* `group_Bats_5` - this contains the list of probes found in all members of
  [eriEur1, equCab2, bosTau4, canFam2, pteVam1]
* `group_Bats_7` - this contains the list of probes found in all members of
  [eriEur1, equCab2, bosTau4, canFam2, pteVam1, myoLuc1, ailMel1]
* `group_Bats_8` - this contains the list of probes found in all members of
  [eriEur1, pteVam1, myoLuc1, bosTau4, vicPac1, equCab2, canFam2, ailMel1]
* `group_Bats_12` - this contains the list of probes found in all members of
  [eriEur1, pteVam1, myoLuc1, bosTau4, vicPac1, susScr9, oviAri1, turTru1,
  equCab2, canFam2, felCat3, ailMel1]
* `group_elephants_7` - this contains the list of probes found in all members of
  [loxAfr3,  canFam2,  echTel1,  choHof1,  dasNov2,  hg19, monDom5]
* `group_size_19` - this contains the list of probes found in all members of
  [anoCar2, bosTau4, calJac3, canFam2, cavPor3, chinese, equCab2, gorGor3,
  hg19, korean, loxAfr3, mm9, monDom5, oryCun2, panTro2, ponAbe2, rheMac2,
  taeGut1, venter]
* `group_size_23` - this contains the list of probes found in all memebrs of
  [bosTau4, canFam2, cavPor3, chinese, dasNov2, echTel1, equCab2, gorGor3,
  hg19, korean, loxAfr3, mm9, monDom5, ochPri2, oryCun2, panTro2, ponAbe2,
  pteVam1, rn4, speTri1, tarSyr1, tupBel1, venter]
* `group_size_25` - this contains the list of probes found in all members of
  [anoCar2, bosTau4, calJac3, canFam2, cavPor3, chinese, dipOrd1, equCab2,
  gorGor3, hg19, korean, loxAfr3, mm9, monDom5, ornAna1, oryCun2, panTro2,
  ponAbe2, pteVam1, rheMac2, rn4, taeGut1, tarSyr1, venter, vicPac1]
* `group_size_29` - this contains the list of probes found in all members of
  [anoCar2, bosTau4, calJac3, canFam2, cavPor3, chinese, choHof1, dipOrd1,
  echTel1, equCab2, eriEur1, gorGor3, hg19, korean, loxAfr3, mm9, monDom5,
  ornAna1, oryCun2, panTro2, ponAbe2, pteVam1, rheMac2, rn4, taeGut1, tarSyr1,
  tupBel1, venter, vicPac1]
* `probe_distribution` - this contains a binary "matrix" indicating
  presence/absence (1/0) of probe matches by species
* `probes` - this is a table providing the ids of the probes we designed
* `species` - this table provides information on the genome build of each
  organism to which we aligned probes
* `sureselect` - this table is identical to the `probes` table from
  uce_probe.sqlite

#### lastz matches of probes to individual species (build version, name, etc in
`species` above)

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


## acknowledgments

We thank the [UCSC genome browser](http://genome.ucsc.edu), in particular, for
being an awesome resource that enables much of the work within.  We also thank
all of the organizations that have made genomic sequences available for the
many organisms we've used as part of this and other work.
