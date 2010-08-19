#!/usr/bin/env python
# encoding: utf-8

"""
muscle.py

Created by Brant Faircloth on 02 May 2010 12:10 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import pdb
import sys
import os
import re
import tempfile
import subprocess
from numpy import sum, sqrt, argmax
from numpy.matlib import repmat, repeat
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment
from Bio.Align.Applications import MuscleCommandline



class Align(object):
    """docstring for Align"""
    def __init__(self, input):
        self.input = input
        self.alignment = None
        self.trimmed_alignment = None
        self.perfect_trimmed_alignment = None
    
    def _clean(self, outtemp):
        # cleanup temp file
        os.remove(outtemp)
        # cleanup input file
        os.remove(self.input)
    
    def _find_ends(self, forward=True):
        """determine the first (or last) position where all reads in an alignment 
        start/stop matching"""
        if forward:
            theRange = xrange(self.alignment.get_alignment_length())
        else:
            theRange = reversed(xrange(self.alignment.get_alignment_length()))
        for col in theRange:
            if '-' in self.alignment.get_column(col):
                pass
            else:
                break
        return col
    
    def _base_checker(self, bases, sequence, loc):
        """ensure that any trimming that occurs does not start beyong the
        end of the sequence being trimmed"""
        # deal with the case where we just want to measure out from the
        # middle of a particular sequence
        if len(loc) == 1:
            loc = (loc, loc)
        if not bases > len(sequence.seq[:loc[0]]) and \
            not bases > len(sequence.seq[loc[1]:]):
            return True
    
    def _record_formatter(self, temp):
        """return a string formatted as a biopython sequence record"""
        temp_record = SeqRecord(temp)
        temp_record.id = sequence.id
        temp_record.name = sequence.name
        temp_record.description = sequence.description
        return temp_record
    
    def _alignment_summary(self, alignment):
        """return summary data for an alignment object using the AlignInfo
        class from BioPython"""
        summary = AlignInfo.SummaryInfo(alignment)
        consensus = summary.dumb_consensus()
        return summary, consensus
    
    def _read(self, format):
        """read an alignment from the CLI - largely for testing purposes"""
        self.alignment = AlignIO.read(open(self.input,'rU'), format)
    
    def get_probe_location(self):
        '''Pull the probe sequence from an alignment object and determine its position
        within the read'''
        # probe at bottom => reverse order
        for record in self.alignment[::-1]:
            if record.id == 'probe':
                start = re.search('^-*', str(record.seq))
                end   = re.search('-*$', str(record.seq))
                # should be first record
                break
        # ooh, this seems so very backwards
        self.ploc = (start.end(), end.start(),)
    
    def run_alignment(self, clean = True, consensus = True):
        """Align, as originally written gets bogged down. Add communicate method 
        and move away from pipes for holding information (this has always been 
        problematic for me with multiprocessing).  Move to tempfile-based
        output."""
        # create results file
        fd, outtemp = tempfile.mkstemp(suffix='.align')
        os.close(fd)
        # run MUSCLE on the temp file
        cline = MuscleCommandline(input=self.input, out=outtemp)
        stdout, stderr = subprocess.Popen(str(cline),
                                 stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 shell=True).communicate(None)
        self.alignment = AlignIO.read(open(outtemp,'rU'), "fasta", alphabet = Gapped(IUPAC.unambiguous_dna, "-"))
        # build a dumb consensus
        if consensus:
            self.alignment_summary, self.alignment_consensus = \
                self._alignment_summary(self.alignment)
        # cleanup temp files
        if clean:
            self._clean(outtemp)
    
    def running_average(self, window_size, threshold):
        # iterate across the columns of the alignment and determine presence
        # or absence of base-identity in the column
        differences = []
        for column in xrange(self.alignment.get_alignment_length()):
            column_values = self.alignment.get_column(column)
            # get the count of different bases in a column (converting
            # it to a set gets only the unique values)
            if len(set(list(column_values))) > 1:
                differences.append(0)
            else:
                differences.append(1)
        # compute the running average from the start => end of the sequence
        forward_average = []
        for start in xrange(len(differences)):
            end = start + window_size
            if end < len(differences):
                forward_average.append(sum(differences[start:end])/float(len(differences[start:end])))
        # compute the running average from the end => start of the sequence
        # we do this, because, otherwise, this end would be neglected.
        reverse_average = []
        for end in reversed(xrange(-len(differences), 0)):
            start = end - window_size
            if start > -len(differences):
                reverse_average.append(sum(differences[start:end])/float(len(differences[start:end])))
        # find where each running average first reaches some threshold 
        # identity over the run span chosen.
        for start_clip, avg in enumerate(forward_average):
            if round(avg, 1) >= float(threshold):
                break
        for temp_end_clip, avg in enumerate(reverse_average):
            if round(avg, 1) >= float(threshold):
                end_clip = len(differences) - temp_end_clip
                break
        return start_clip, end_clip
    
    def trim_alignment(self, method = 'edges', remove_probe = None, bases = None, consensus = True, window_size = 20, threshold = 0.5):
        """Trim the alignment"""
        if method == 'edges':
            # find edges of the alignment
            start   = self._find_ends(forward=True)
            end     = self._find_ends(forward=False)
        elif method == 'running':
            start, end = self.running_average(window_size, threshold)
        # create a new alignment object to hold our alignment
        self.trimmed_alignment = Alignment(Gapped(IUPAC.ambiguous_dna, "-"))
        for sequence in self.alignment:
            # ignore the probe sequence we added
            if (method == 'edges' or method == 'running') and not remove_probe:
                # it is totally retarded that biopython only gives us the option to
                # pass the Alignment object a name and str(sequence).  Given this 
                # level of retardation, we'll fudge and use their private method
                self.trimmed_alignment._records.append(sequence[start:end])  
            elif method == 'static' and not remove_probe and bases:
                # get middle of alignment and trim out from that - there's a
                # weakness here in that we are not actually locating the probe
                # region, we're just locating the middle of the alignment
                mid_point = len(sequence)/2
                if self._base_checker(bases, sequence, mid_point):
                    self.trimmed_alignment._records.append(
                        sequence[mid_point-bases:mid_point+bases]
                        )
                else:
                    self.trimmed_alignment = None
            elif method == 'static' and not remove_probe and bases and self.ploc:
                # get middle of alignment and trim out from that - there's a
                # weakness here in that we are not actually locating the probe
                # region, we're just locating the middle of the alignment
                if self._base_checker(bases, sequence, self.ploc):
                    self.trimmed_alignment._records.append(
                        sequence[self.ploc[0]-bases:self.ploc[1]+bases]
                        )
                else:
                    self.trimmed_alignment = None
            elif remove_probe and self.ploc:
                # we have to drop to sequence level to add sequence slices
                # where we basically slice around the probes location
                temp = sequence.seq[:self.ploc[0]] + sequence.seq[self.ploc[1]:]
                self.trimmed_alignment._records.append( \
                    self._record_formatter(temp)
                    )
            elif method == 'static' and remove_probe and bases and self.ploc:
                if self._base_checker(bases, sequence, self.ploc):
                    temp = sequence.seq[self.ploc[0]-bases:self.ploc[0]] + \
                        sequence.seq[self.ploc[1]:self.ploc[1]+bases]
                    self.trimmed_alignment._records.append( \
                        self._record_formatter(temp)
                        )
                else:
                    self.trimmed_alignment = None
        # build a dumb consensus
        if consensus:
            self.trimmed_alignment_summary, self.trimmed_alignment_consensus = \
                self._alignment_summary(self.trimmed_alignment)
    
    def trim_ambiguous_bases(self):
        """snip ambiguous bases from a trimmed_alignment"""
        ambiguous_bases = []
        # do this by finaing all ambiguous bases and then snipping the largest
        # chunk with no ambiguous bases from the entire alignment
        for column in xrange(0, self.trimmed_alignment.get_alignment_length()):
            if 'N' in self.trimmed_alignment.get_column(column):
                ambiguous_bases.append(column)
        maximum = 0
        maximum_pos = None
        #pdb.set_trace()
        if ambiguous_bases:
            # prepend and append the start and end of the sequence so consider
            # those chunks outside the stop and start of ambiguous base runs.
            ambiguous_bases.insert(0,0)
            ambiguous_bases.append(self.trimmed_alignment.get_alignment_length() - 1)
            # create a new alignment object to hold our alignment
            self.perfect_trimmed_alignment = \
                Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
            for pos in xrange(len(ambiguous_bases)):
                if pos + 1 < len(ambiguous_bases):
                    difference = ambiguous_bases[pos + 1] - \
                        ambiguous_bases[pos]
                    if difference > maximum:
                        maximum = difference
                        maximum_pos = (pos, pos+1)
                else:
                    pass
            # make sure we catch cases where there is not best block
            if maximum_pos:
                for sequence in self.trimmed_alignment:
                    self.perfect_trimmed_alignment._records.append(
                        sequence[ambiguous_bases[maximum_pos[0]] + 1
                            :ambiguous_bases[maximum_pos[1]]]
                            )
            else:
                self.perfect_trimmed_alignment = None
        else:
            self.perfect_trimmed_alignment = self.trimmed_alignment
            


if __name__ == '__main__':
    #test_alignment = Align('/Users/bcf/git/brant/seqcap/Test/align/new/chrZ_8059.nex')
    test_alignment = Align('/Users/bcf/git/brant/seqcap/Test/align/amb_trim_concat/concat.nex')
    test_alignment._read('nexus')
    pdb.set_trace()
    test_alignment.trim_alignment(method='running')
    test_alignment.trim_ambiguous_bases()
    pdb.set_trace()