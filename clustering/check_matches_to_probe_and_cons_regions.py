#!/usr/bin/env python
# encoding: utf-8

"""
targetMapperChecker.py

Created by Brant Faircloth on 1 May 2010 11:00 PDT (-0700).
Copyright (c) 2010 Brant Faircloth. All rights reserved.

PURPOSE:  

The purpose of this program is to check the results of a lastz alignment
between (simulated) sequence reads and the 2bit file of probes designed from
UCE regions.  This is slightly different from targetMapperChecker.py in that
this program maps reads to probes from UCEs rather than UCEs themselves.  I
actually think this is more accurate, because we **know** that each read
should contain a singificant portion of a particular probe, but we do not know
how much of each conserved region is should contain (this depends on the size
of the UCE).

The program takes output from ../Simulation/runLastz.py or
../Simulation/easyLastz.py

USAGE:
    
    python targetMapperCheckerII.py \
        --lastz=testAnoCar2ToConservedRegions.lastz \
        --reads=../Simulation/johns_sequence/group_size_29/anoCar2.simSequence.500_flank.fa

"""

import pdb
import os
import re
import sys
#import time
import optparse
from Bio import SeqIO


def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"
    
    p = optparse.OptionParser(usage)
    
    p.add_option('--lastz', dest = 'lastz', action='store',
type='string', default = None, help='The path to the lastz file.',
metavar='FILE')
    
    p.add_option('--reads', dest = 'reads', action='store',
type='string', default = None, help='The path to the lastz file.',
metavar='FILE')
    
    (options,arg) = p.parse_args()
    if not options.lastz or not options.reads:
        p.print_help()
        sys.exit(2)
    if not os.path.isfile(options.lastz) or not os.path.isfile(options.reads):
        print 'You must provide a valid path to the configuration file.'
        p.print_help()
        sys.exit(2)
    return options, arg

def buildTagDictionary(conserved):
    '''docstring for buildTagDictionary'''
    ts = set()
    for record in SeqIO.parse(open(conserved,'rU'), 'fasta'):
        ts.add(record.id)
    return ts

def main():
    options, args = interface()
    tagSet = buildTagDictionary(options.reads)
    probe_error     = 0
    locus_error     = 0
    errorList       = []
    probe_correct   = 0
    locus_correct   = 0
    probe_count     = 0
    locus_count     = 0
    readSet         = set()
    readMatch       = {}
    readDict        = {}
    for line in open(options.lastz, 'rU'):
        bits = line.strip('\n').split('\t')
        bitsList = bits[6].strip('>').split(' ')
        readUniqueName, readProbeName = bitsList[0], bitsList[-1].split('|')[1]
        # we can get multiple matches for a read.  Determine the best match
        # by building a dictionary for each read containing all the matches
        # indexed by their score.
        if readUniqueName in readMatch:
            readDict[readUniqueName] = readProbeName
            readMatch[readUniqueName][int(bits[0])] = bits[1]
        else:
            readDict[readUniqueName] = readProbeName
            readMatch[readUniqueName] = {int(bits[0]):bits[1]}
        readSet.add(readUniqueName)
    regex = re.compile('(.*)_chr.*_probe_[0-9]')
    for read in readMatch:
        #pdb.set_trace()
        m = readMatch[read]
        bestMatch = m[max(m.keys())]
        temp = regex.search(bestMatch)
        bestMatchLocus = temp.groups()[0]
        realName  = readDict[read]
        temp = regex.search(realName)
        realNameLocus = temp.groups()[0]
        if bestMatch != realName:
            probe_error     += 1
            errorList.append('{0}:{1}'.format(bestMatch, realName))
        else:
            probe_correct   += 1
        probe_count += 1
        if bestMatchLocus != realNameLocus:
            locus_error     += 1
        else:
            locus_correct   += 1
        locus_count +=1
    # we use count rather than len(readSet) because some reads may be matched
    # to their source UCE > 1 time, which makes them right, but duplicates.
    # As duplicates, they wouldn't show up in the set.
    print "Correct (Probe)\t\t= {0:2} %".format(probe_correct/float(probe_count) * 100)
    print "Incorrect (Probe)\t= {0:2} %".format(probe_error/float(probe_count) * 100)
    print "Correct (locus)\t\t= {0:2} %".format(locus_correct/float(locus_count) * 100)
    print "Incorrect (locus)\t= {0:2} %".format(locus_error/float(locus_count) * 100)
    #for e in errorList:
    #    print '\t\tUCE:{0}'.format(e)
    #pdb.set_trace()
    missing = tagSet.difference(readSet)
    print "Missed\t\t\t= {0:.2} %".format(len(missing)/float(len(tagSet)) * 100)
    if len(missing):
        print "Missing from lastz groupings:"
        for missed in missing:
            print '\t\t{0}'.format(missed)

if __name__ == '__main__':
    main()