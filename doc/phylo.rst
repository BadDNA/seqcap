Phylo
======

Script Descriptions
*******************

-  **bootstrap_locus_alignments.py**

   *Description:* This script creates bootstrap replicates of
   individual locus alignments using a two stage bootstrap method (Seo
   *et al.* 2005; Seo 2008). First the genes (e.g., UCEs) are
   resampled with replacement. Second, aligned positions within each
   gene are resampled with replacement.

   The script takes as input a directory of individual alignments in
   either phylip or nexus format. The output is series of
   multi-alignment files equal to the number of bootstrap replicates.
   Each of these multi-alignment files contain individual bootstrapped
   alignments totaling the number of individual alignments in the
   original input directory.

   *Requirements:* Python version 2.6 with modules pylab and biopython
   installed.

   *Command line:*

   ::

       >>> python bootstrap_locus_alignments.py --input-dir /some/dir/ \
       --output-dir /some/out/dir --input-file-format nexus \
       --output-file-format phylip --bootstrap-reps 500 --multi_out

   Detailed information about the command line arguments may be
   obtained by typing:
   
   ::
   
        >>> python bootstrap_locus_alignments.py --help

-  **bootstrap_concat_alignments.py**

   *Description:* This script takes a single alignment and generates
   alignments bootstrapped and aligned position. It can handle large
   alignments, but is a bit slow.

   *Requirements:* Python version 2.6 with modules pylab and biopython
   installed.

   *Command line:*

   ::

       >>> python bootstrap_concat_alignments --input-file /some/file.align \
       --output-dir /some/out/dir --input-file-format nexus \  
       --output-file-format phylip --bootstrap-reps 500

   Detailed information about the command line arguments may be
   obtained by typing:
   
   ::
   
        >>> python bootstrap_locus_alignments.py --help

-  **execute_PhyML.py**

   *Description:* This is a template for executing multiple PhyML jobs
   on an LSF cluster. Just modify the path and the command line as
   necessary and then run the code from within the python interpreter
   (e.g., iPython). Because every computer cluster is different, it is
   impossible to provide code that can be run everywhere.

-  **phybase.py**

   *Description:* Script to run Phybase on a set of trees. Easier than
   loading R everytime.

   *Requirements:* Python version 2.6 with modules rpy2, dendropy, and
   biopython installed.

   *Command line:*

   ::

       >>> python phybase.py --input-file /some/file.trees \
       --output-dir /some/out/dir --outgroup fish \
       --taxa 'mouse human chimp etc'

   Detailed information about the command line arguments may be
   obtained by typing: 
   
   ::
   
        >>> python phybase.py --help

-  **phylo.py**

   *Description:* A home brewed phylogenetic module. Required by some
   of the preceding scripts.

   *Requirements:* Python version 2.6 with module Numpy.

-  **get_mpest_format.py**

   *Description:* Given directory or file input, read the contents of
   the either and root the trees within (using --root and --outgroup=
   options). Additionally, build a control file for our slightly
   customized version of mpest.

   *Requirements:* `ETE2 <http://ete.cgenomics.org/>`_ module and
   numpy

   *Command line:*

   ::

       >>> python get_mpest_format.py --input=file_of_.trees \
       --output=ouput.tree --root --outgroup=ASpeciesName --build-control

   Detailed information about the command line arguments may be
   obtained by typing ``>>> get_mpest_format.py --help``

-  **run_mpest.py**

   *Description:* Run an input file through our slightly customized
   version of mpest either using a single core or multiple cores, and
   feeding mpest a random (integer) seed drawn from a uniform
   distribution.

   *Requirements:* Python version 2.6

   *Command line:*

   ::

       >>> python run_mpest.py --control-file=903_loci_5_species_bats.control \
       --iterations=1000 --cores=7 --output=903_loci_5_species_bats_mpest.tree

-  **add_taxa_names_to_mpest.py**

   *Description:* Adds human-readable labels to MP-EST output.

   *Requirements:* Python version 2.6

   *Command line:*

   ::

       >>> python /n/home06/ngcrawford/data1/seqcap/Phylo/add_taxa_names_to_mpest.py \
       --input=/path/dir/containing/trees/and/control/files


Pipeline Examples
*****************


**STEAC or STAR Trees:**

    1.) Execute PhyML on each locus alignment using the template in
    ``execute_PhyML.py``

    2.) Run ``phybase.py`` on the gene trees to make the species tree
    (STEC or STAR).

**Concatenated Alignment Boostrapping:**

    1.) Generate bootstrapped alignments with
    ``bootstrap_concat_alignments.py``

    2.) Execute PhyML on each replicate using the template in
    ``execute_PhyML.py``

    3.) Concatenate the trees into a single file with the linux 'cat'
    command.

    ::

        >>> cat *.tree > concat.trees

4.) Generate a consensus tree on the species trees with PAUP or
similar.

**Per Locus Bootstrapping:**

    1.) Generate bootstrapped alignments with
    ``bootstrap_locus_alignments.py``

    2.) Execute PhyML on each replicate using the template in
    ``execute_PhyML.py``

    3.) Run ``phybase.py`` on each set of gene trees to make the
    species trees.

    4.) Concatenate the species trees into a single file with the linux
    'cat' command.

    ::

        >>> cat star*.tree > all_star_species.trees

5.) Generate a consensus tree on the species trees with PAUP or
similar.

**MP-EST Trees:**

1. Single MP-EST tree from a collection of PhyML trees.

    a.) Run get `get_mpest_format.py` on tree file.
    
    b.) Run, the appropriately named, `run_mpest.py`
    
    c.) Run `add_taxa_names_to_mpest.py` to added readable taxa names
    
    d.) Concatenate output:
     
    ::
    
        >>> cat *.mpest.named.trees > all_mpest_named_species.trees
    
    e.) Generate consensus tree in `PAUP <http://paup.csit.fsu.edu/>`_ and visualize with PAUP or `FigTree <http://tree.bio.ed.ac.uk/software/figtree/>`_

2. Bootstrap MP-EST trees from a collection of PhyML trees.

    a.) Generate bootstrapped alignments with `bootstrap_locus_alignments.py`
    
    b.) Execute PhyML on each replicate using the template in `execute_PhyML.py`
    
    c.) Make control files using `get_mpest_format.py`
    
    d.) Run MP-EST. I used the following PYTHON script on an `LSF <http://en.wikipedia.org/wiki/Platform_LSF>`_ cluster to automate the process:
    
    
    ::

            import os, glob, shlex, subprocess
            paths = glob.glob('/path/to/bootreps/*.control')
            for count, p in enumerate(paths):
                path, fname = os.path.split(p)
                fout = fname.split('.')[0] + '.mpest.tree'
                fout = os.path.join(path,fout)
                command = "bsub -o sterr.out -q normal_serial -R 'select[mem>3000]' \
                    'python /path/to/run_mpest.py \
                    --control-file=%s \
                    --iterations=100 \
                    --cores=1 \
                    --output=%s'" % (p, fout)
                command = shlex.split(command)
                subprocess.Popen(command)
    
    e.) Run `add_taxa_names_to_mpest.py` to added readable taxa names
    
    d.) Concatenate output e.g, `>>> cat *.mpest.named.trees > all_mpest_named_species.trees`
    
    f.) Generate consensus tree in `PAUP <http://paup.csit.fsu.edu/>`_ and visualize with PAUP or `FigTree <http://tree.bio.ed.ac.uk/software/figtree/>`_

Citations
*********

Seo et al. Incorporating gene-specific variation when inferring and
evaluating optimal evolutionary tree topologies from multilocus
sequence data. Proc. Natl. Acad. Sci. U.S.A. (2005) vol. 102 (12)
pp. 4436-41

Seo. Calculating bootstrap probabilities of phylogeny using
multilocus sequence data. Mol. Biol. Evol. (2008) vol. 25 (5) pp.
960-71