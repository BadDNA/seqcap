#!/usr/bin/env python
# encoding: utf-8

"""
targetMapperChecker.py

Created by Brant Faircloth on 24 March 2010 21:55 PDT (-0700).
Copyright (c) 2010 Brant Faircloth. All rights reserved.

The purpose of this program is to check the results of a lastz alignment
between (simulated) sequences reads and the 2bit file of UCE regions from
which we designed probes (we generated the 2bit file using
buildLociFromProbes.py).

The program takes output from ../Simulation/runLastz.py or
../Simulation/easyLastz.py

USAGE:
    
    python ../../targetMapperChecker.py \
    --lastz=amniotesMappedToSelectedConservedLociPlusFlank.lastz \
    --reads=../amniotes.ProbeSimSeq.100.entire.fa > results.txt

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
    td = {}
    ts = set()
    regex = re.compile('(.*)_chr.*_probe_[0-9]')
    for record in SeqIO.parse(open(conserved,'rU'), 'fasta'):
        temp = regex.search(record.description.split('|')[-1])
        seq = temp.groups()[0]
        td[record.description] = seq
        ts.add(seq)
    invTd = dict((v,k) for k, v in td.iteritems())
    return td, invTd, ts

def main():
    options, args = interface()
    tagDictionary, invTagDictionary, tagSet = \
        buildTagDictionary(options.reads)
    error       = 0
    errorList   = []
    correct     = 0
    count       = 0
    readSet     = set()
    readMatch   = {}
    for line in open(options.lastz, 'rU'):
        bits = line.strip('\n').split('\t')
        # we can get multiple matches for a read.  Determine the best match
        # by building a dictionary for each read containing all the matches
        # indexed by their score.
        if bits[1] in readMatch:
            readMatch[bits[1]][int(bits[0])] = bits[6]
        else:
            readMatch[bits[1]] = {int(bits[0]):bits[6]}
        # alternate method ~ same speed
        # readMatch[bits[1]] = dict(readMatch.get(bits[1],{}), \
        #    **{int(bits[0]):bits[6]})
    for read in readMatch:
        m = readMatch[read]
        bestMatch = m[max(m.keys())].lstrip('>')
        if read not in invTagDictionary.keys():
            error       += 1
            #errorList.append('{0}:{1}'.format(tagDictionary[read], bestMatch))
        elif invTagDictionary[read] == bestMatch:
            correct     += 1
        elif invTagDictionary[read] != bestMatch:
            error       += 1
            #errorList.append('{0}:{1}'.format(tagDictionary[read], bestMatch))
        readSet.add(read)
        count += 1
    # we use count rather than len(readSet) because some reads may be matched
    # to their source UCE > 1 time, which makes them right, but duplicates.
    # As duplicates, they wouldn't show up in the set.
    print "Correct\t\t= {0:2} %".format(correct/float(count) * 100)
    print "Incorrect\t= {0:2} %".format(error/float(count) * 100)
    for e in errorList:
        #pdb.set_trace()
        inv = invTagDictionary[e.split(':')[0]]
        print '\t\t{0} UCE:{1}'.format(inv, e)
    #pdb.set_trace()
    missing = tagSet.difference(readSet)
    print "Missed\t\t= {0:.2} %".format(len(missing)/float(len(tagSet)) * 100)
    if len(missing):
        print "Missing from lastz groupings:"
        for missed in missing:
            print '\t\t{0}'.format(missed)

if __name__ == '__main__':
    main()