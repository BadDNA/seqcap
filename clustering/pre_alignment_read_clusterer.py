#!/usr/bin/env python
# encoding: utf-8

"""
clusterer.py

Created by Brant Faircloth on 01 May 2010 15:12 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import pdb
import os
import re
import sys
import glob
import tempfile
import optparse
import subprocess
import seqcap.lib.lastz
import seqcap.lib.muscle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--input', dest = 'input', action='store', 
type='string', default = None, help='The path to the input directory')

    p.add_option('--probes-2bit', dest = 'probes', action='store', 
type='string', default = None, help='The path to probes 2bit file', \
metavar='FILE')    

    p.add_option('--output', dest = 'output', action='store', 
type='string', default = None, help='The pathname to the output file.')

    p.add_option('--matchcount', dest = 'matchcount', action='store', \
type='int', default = 80, help='The number of matched bases in the \
alignment.')

    p.add_option('--identity', dest = 'identity', action='store', \
type='float', default = 60, help='The fraction of aligned bases \
(excluding columns containing gaps or non-ACGT characters) that are \
matches, expressed as a percentage')

    p.add_option('--fish', dest = 'fish', action='store_true', default=False, 
help='If pre-clustering fish reads')

    (options,arg) = p.parse_args()
    if not options.input:
        p.print_help()
        sys.exit(2)
    if not os.path.isdir(options.input):
        print "You must provide a valid path to the input directory."
        p.print_help()
        sys.exit(2)
    return options, arg 

def main():
    options, args = interface()
    if not options.fish:
        regex = re.compile('(.*)_chr.*_probe_[0-9]')
    # alter things a bit for fish since their locus names are slightly
    # different
    else:
        regex = re.compile('(.*)_[chr|ultracontig|scaffold].*_[0-9]+_[0-9]+.*')
    # all output is going into a single file, so go ahead and open it up
    #pdb.set_trace()
    outp_path = os.path.expanduser(options.output)
    outp = open(outp_path, 'w')
    for infile in glob.glob(os.path.join(os.path.expanduser(options.input), '*.fa')):
        #pdb.set_trace()
        # get the species name from the fasta file
        species = os.path.basename(infile).split('.')[0]
        # make a tempfile for the lastz output
        fd, lastz_temp_file = tempfile.mkstemp(suffix='.lastz')
        os.close(fd)
        # lastz this file against the probe 2bit file
        lastz = seqcap.lib.lastz.Align(options.probes, infile, options.matchcount, \
            options.identity, lastz_temp_file)
        lastz_stdout, lastz_stderr = lastz.run()
        # parse the lastz result to place a read with it's best matching
        # **UCE** (by parsing the match header) versus it's best matching
        # probe
        #pdb.set_trace()
        if not lastz_stderr:
            match = {}
            for line in open(lastz_temp_file, 'rU'):
                line_list = line.strip('\n').split('\t')
                match_score = int(line_list[0])
                match_name  = line_list[1]
                seq_id      = line_list[6].strip('>').split(' ')[0]
                if seq_id in match:
                    match[seq_id][match_score] = match_name
                else:
                    match[seq_id] = {match_score:match_name}
            # iterate over the match_dict to get the best match for a read
            best_match = {}
            cons_match = {}
            #pdb.set_trace()
            for m in match:
                temp_best = match[m][max(match[m].keys())]
                temp2_best = regex.search(temp_best)
                best = temp2_best.groups()[0]
                if not best:
                    print 'no best match'
                    pdb.set_trace()
                best_match[m] = best
            #pdb.set_trace()
            for k,v in best_match.iteritems():
                if v in cons_match:
                    cons_match[v].append(k)
                else:
                    cons_match[v] = [k]
            multi_probe_holder = {}
            sequences_to_write = []
            #pdb.set_trace()
            for record in SeqIO.parse(open(infile, 'rU'), 'fasta'):
                try:
                    locus_name = best_match[record.id]
                    if len(cons_match[locus_name]) == 1:
                        #pdb.set_trace()
                        # reformat the record.id
                        record.id = '{1}'.format(locus_name, species)
                        record.name = record.id
                        record.description = \
                            '|{1}|{0} Locus {1} singleton sequence'.format(species, 
                            locus_name)
                        # we just write the record
                        sequences_to_write.append(record)
                    else:
                        if best_match[record.id] in multi_probe_holder:
                            multi_probe_holder[best_match[record.id]].append(record)
                        else:
                            multi_probe_holder[best_match[record.id]] = [record]
                except KeyError:
                    print 'Sequence {0} skipped'.format(record.id)
            # now, for each multi-record, we need to align the records, snip
            # the alignment, generate a consensus, and write that out to 
            # wherever the singleton sequences are
            for locus_name, record_set in multi_probe_holder.iteritems():
                fd, fasta_temp_file = tempfile.mkstemp(suffix='.fasta')
                os.close(fd)
                SeqIO.write(record_set, open(fasta_temp_file, 'w'), 'fasta')
                muscle = seqcap.lib.muscle.Align(fasta_temp_file)
                muscle.run_alignment()
                locus_consensus_seq_record = SeqRecord(muscle.alignment_consensus)
                locus_consensus_seq_record.id = '{0}_{1}'.format(locus_name, species)
                locus_consensus_seq_record.name = locus_consensus_seq_record.id
                locus_consensus_seq_record.description = \
                    '|{1}|{0} Locus {1} consensus sequence'.format(species, 
                    locus_name)
                sequences_to_write.append(locus_consensus_seq_record)
            # write the records
            SeqIO.write(sequences_to_write, outp, 'fasta')
        else:
            print 'Lastz: ', lastz_stderr
        os.remove(lastz_temp_file)
    outp.close()
                

if __name__ == '__main__':
    main()
