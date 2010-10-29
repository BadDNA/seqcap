#!/usr/bin/env python
# encoding: utf-8
"""
bootstrap.py

Created by Nick Crawford on 2010-05-16.
Copyright (c) 2010

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com
"""

import os
import sys
import glob
import phylo
import optparse
import subprocess
from pylab import *
from Bio.Seq import Seq 
from Bio import AlignIO, SeqIO
from copy import deepcopy, copy
from Bio.Align.Generic import Alignment 
from Bio.Alphabet import IUPAC, Gapped

def interface():
    """ create commandline interface for script"""
    p = optparse.OptionParser()
    
    p.add_option('--input-dir','-i', type='string',
        help='Path to input directory.')
    p.add_option('--output-dir','-o', type='string',
        help='Path to output directory.')
    p.add_option('--input-file-format','-f', type='string',
        help="Either 'phylip' or 'nexus'.")
    p.add_option('--output-file-format','-F', type='string',
        help="Either 'phylip' or 'nexus'.")
    p.add_option('--bootstrap-reps','-b', type='int',
        help="Number of bootstrap replicates to generate.")
    p.add_option('--multi_out','-m', action="store_true",
        help="Create individual files for each bootstrap replicate.")
    
    options, arg = p.parse_args()
    
    # check options for errors, etc.
    if options.input_dir == None:
        print "Input directory required."
        print "Type 'bootstrap.py -h' for details" 
        sys.exit()
    
    if options.output_dir == None:
        print "No output directory supplied."
        print "Output will be created in input directory."
        options.output_dir = os.path.expanduser(options.output_dir)
        
    if os.path.isdir(options.output_dir) == False:
        print 'Output directory does not exist.'
        os.mkdir(options.output_dir)
        print 'Created new directory at %s' % options.output_dir
    
    if options.output_file_format == None:
        print "Output file format required."
        print "Type 'bootstrap.py -h' for details"
        sys.exit()
    
    if options.input_file_format == None:
        print "Input file format required."
        print "Type 'bootstrap.py -h' for details"
        sys.exit()
    
    if options.bootstrap_reps == None:
        print 'An integer representing the number of'
        print 'bootstrap replicates is required.'
        print '(defaults to 10)'
        print '\n'
        print "Type 'bootstrap.py -h' for details" 
        sys.exit()
        
    return options # not returning arguments

def biopy2strarray(biopython_alignment):
    """convert biopython alignment into an numpy 
    character array with list of taxon IDs"""
    
    str_alignment = []
    ids = []
    for count, sequence in enumerate(biopython_alignment):
        bases = sequence.seq
	id_ = sequence.id.split('_')[-1]                    # quick fix for weird taxa names
        bases = bases.tostring()
        bases_str_array = fromstring(bases,dtype='|S1') # convert char string to numpy array
        str_alignment.append(bases_str_array)
        ids.append(id_)
    str_alignment = array(str_alignment)
    return [str_alignment, ids]

def strarray2biopy(align):
    """ take a 2d character array with an associated ID list 
    and convert it into a biopython DNA alignment."""
        
    seqs = align[0]
    ids = align[1]

    alphabet = Gapped(IUPAC.unambiguous_dna) 
    alignment = Alignment(alphabet) 
       
    for count, array_seq in enumerate(seqs):
        bases = ''
            
        for base in array_seq:
            bases += base
                
        alignment.add_sequence(ids[count],bases)
        
    return alignment

    
def make_align_list(in_dir,file_type):
    """ From the supplied directory parse the appropriate file types
    into as list of Biopython alignments
    """
        
    in_dir = os.path.expanduser(in_dir) # properly deal with ~ in file path 
    print in_dir
        
    # open the correct type of file
    if file_type.lower() == 'nexus':
        alignments = glob.glob( os.path.join(in_dir,'*.nex'))
    if file_type.lower() == 'phylip':
        alignments = glob.glob( os.path.join(in_dir,'*.phylip'))
  
    concat_alignments = []
    for count, alignment_file in enumerate(alignments):
        biopy_alignment = AlignIO.read(open(alignment_file), file_type)
        # print biopy_alignment
        numeric_alignment = biopy2strarray(biopy_alignment)
        concat_alignments.append(numeric_alignment)
        
    return concat_alignments
    
def executePhyML(fname):
    phyml_cmd = 'PhyML_3.0_linux64 --input %s --multiple 1100 --quiet' % fname
                
    bsub_cmd = ['bsub','-q','normal_serial','-o','/n/home06/ngcrawford/sterr.out','-R',"'select[mem>3000]'", phyml_cmd]
    print subprocess.list2cmdline(bsub_cmd)
    p = subprocess.Popen(bsub_cmd)
    

def main():
    options = interface()
    array_alignments = make_align_list(options.input_dir, options.input_file_format)
    
    # DO BOOTSTRAPPING AND WRITE REPLICATES TO OUTFILES
    # options.bootstrap_reps += 1
    boot_alignments = phylo.bootstrap(array_alignments, options.bootstrap_reps)
    bootstrapped_datasets = []
    for rep_count, boot_rep in enumerate(boot_alignments):
        boot_bases = []
        
        # setup outfile names (make into function?)
        if options.output_file_format == 'nexus':
            fname = 'bootrep_%s.nex' % rep_count
        if options.output_file_format == 'phylip':
            fname = 'bootrep_%s.phylip' % rep_count
        final_path = os.path.join(options.output_dir,fname)
        fout = open(final_path,'a')
    
        for count, align in enumerate(boot_rep):
            seqs = copy(align[0])                               # lots of copying to be 'safe'
            ids = copy(align[1])
            bases_by_col = np.column_stack(seqs)                # flip rows and columns
            bs_bases = phylo.bootstrap(bases_by_col, 1)         # bootstrap the bases within the bootstrapped alignments
            bs_bases = np.column_stack(bs_bases[0])             # [0] corrects weirdnesss due to extra set of brackets
            bs_bases = bs_bases.copy()                          # copy modified replicate
            pair = [bs_bases, ids]                              
            biopy_align = [strarray2biopy(pair)]
            AlignIO.write(biopy_align, fout, options.output_file_format)
            fout.write('\n')
        
        fout.close()
        # print 'created bootstrap replicate %s' % rep_count
        # # calculate trees using PhyML
        # executePhyML(final_path)
        # print 'executed PhyML on replicate %s' % rep_count
        # print '\n'


if __name__ == '__main__':
    main()

