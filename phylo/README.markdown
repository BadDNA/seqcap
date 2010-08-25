Scripts:
========

* **`bootstrap_locus_alignments.py`**

    *Description:* This script creates bootstrap replicates of individual locus alignments using a two stage bootstrap method (Seo *et al.* 2005; Seo 2008). First the genes (e.g., UCEs) are resampled with replacement..  Second, aligned positions within each gene are resampled with replacement.
    
    The script takes as input a directory of individual alignments in either phylip or nexus format.  The output is series of multi-alignment files equal to the number of bootstrap replicates.  Each of these multi-alignment files contain individual bootstrapped alignments totaling the number of individual alignments in the original input directory.
    
    *Requirements:* Python version 2.6 with modules pylab and biopython installed. 
    
    *Command line:* 
    
        >>> python bootstrap_locus_alignments.py --input-dir /some/dir/ --output-dir /some/out/dir --input-file-format nexus --output-file-format phylip --bootstrap-reps 500 --multi_out
        
    Detailed information about the command line arguments may be obtained by typing `>>> python bootstrap_locus_alignments.py --help`

* **`bootstrap_concat_alignments.py`**

    *Description:* This script takes a single alignment and generates alignments bootstrapped and aligned position. It can handle large alignments, but is a bit slow.
    
    *Requirements:* Python version 2.6 with modules pylab and biopython installed. 
    
    *Command line:*
    
        >>> python bootstrap_concat_alignments --input-file /some/file.align --output-dir /some/out/dir --input-file-format nexus --output-file-format phylip --bootstrap-reps 500
    
    Detailed information about the command line arguments may be obtained by typing `>>> python bootstrap_locus_alignments.py --help`

* **`execute_PhyML.py`** *

    *Description:* This is a template for executing multiple PhyML jobs on an LSF cluster.  Just modify the path and the command line as necessary and then run the code from within the python interpreter (e.g., iPython).  Because every computer cluster is different, it is impossible to provide code that can be run everywhere.
    
* **`phybase.py`** *

    *Description:* Script to run Phybase on a set of trees. Easier than loading R everytime.
    
    *Requirements:* Python version 2.6 with modules rpy2, dendropy, and biopython installed. 

    *Command line:*
    
        >>> python phybase.py --input-file /some/file.trees --output-dir /some/out/dir --outgroup fish --taxa 'mouse human chimp etc'

    Detailed information about the command line arguments may be obtained by typing `>>> python phybase.py --help`


* **`phylo.py`**

    *Description:* A home brewed phylogenetic module. Required by some of the preceding scripts.
    
    *Requirements:* Python version 2.6 with module Numpy.

Methods:
=======

**Species Tree:**

1.) Execute PhyML on each locus alignment using the template in `execute_PhyML.py`

2.) Run `phybase.py` on the gene trees to make the species tree.


**Concatenated Alignment Boostrapping:**

1.) Generate bootstrapped alignments with `bootstrap_concat_alignments.py`

2.) Execute PhyML on each replicate using the template in `execute_PhyML.py`

3.) Concatenate the trees into a single file with the linux 'cat' command. 

        >>> cat *.tree > concat.trees

4.) Generate a consensus tree on the species trees with PAUP or similar.


**Per Locus Bootstrapping:**

1.) Generate bootstrapped alignments with `bootstrap_locus_alignments.py`

2.) Execute PhyML on each replicate using the template in `execute_PhyML.py`

3.) Run `phybase.py` on each set of gene trees to make the species trees.

4.) Concatenate the species trees into a single file with the linux 'cat' command. 

        >>> cat star*.tree > all_star_species.trees

5.) Generate a consensus tree on the species trees with PAUP or similar. 


Citations:
==========

    Seo et al. Incorporating gene-specific variation when inferring and evaluating optimal evolutionary tree topologies from multilocus sequence data. Proc. Natl. Acad. Sci. U.S.A. (2005) vol. 102 (12) pp. 4436-41
    
    Seo. Calculating bootstrap probabilities of phylogeny using multilocus sequence data. Mol. Biol. Evol. (2008) vol. 25 (5) pp. 960-71
    
