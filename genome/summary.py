#!/usr/bin/env python
# encoding: utf-8
"""
summary.py

Created by Brant Faircloth on 2008-07-05.
Copyright (c) 2008 Brant Faircloth. All rights reserved.

This program scans MAF files for conserved elements and stores 
those results in an sqlite database

"""

import os
import pp
import re
import pdb
import time
import numpy
import sqlite3
import sequence
import bx.align.maf

def spScreen(a, minAlignLength):
    '''screen alignments to ensure minSpecies and minAlignLength'''
    for spp in a.components:
        if len(a.components[0].text) > minAlignLength:
            return a

def alignMetadata(counter, candAlign, cons, refPosition, altPosition, metadataKey):
    '''get metdata for alignment based on species in metadataKey'''
    for seq in candAlign.components:
        name = seq.src
        metadata = {}
        if name[0:7] == metadataKey:
            metadata['spp']         = name[0:7]
            metadata['chromo']      = name[10:]
            metadata['aln_start']   = seq.forward_strand_start
            metadata['aln_len']     = seq.size
            metadata['aln_end']     = seq.forward_strand_end
            metadata['cons']        = cons
            metadata['cons_len']    = len(cons)
            # add values to metadata, making up for 0 indexing
            metadata['cons_start']  = metadata['aln_start'] + 1 + refPosition[0]
            metadata['cons_end']    = metadata['aln_start'] + refPosition[1]
            metadata['alt_start']   = candAlign.components[1].start + 1 + refPosition[0]
            metadata['alt_end']     = candAlign.components[1].start + refPosition[1]
            metadata['alt_chromo']  = candAlign.components[1].src[7:]
            metadata['map']         = (('chr%s:%s-%s') % (metadata['chromo'], metadata['cons_start'], metadata['cons_end']))
            metadata['seq']         = (('chr%s_%s') % (metadata['chromo'], counter))
        break
    #pdb.set_trace()
    return metadata

def createCons(candAlign):
    '''stack sequence and return dumb (but smart!) consensus with
    metadata'''
    for seq in range(len(candAlign.components)):
        if seq == 0:
            zString = candAlign.components[seq].text
            zString = numpy.array(list(zString))
            seqArray = zString
        else:
            nzString = candAlign.components[seq].text
            nzString = numpy.array(list(nzString))
            seqArray = numpy.vstack((seqArray,nzString))
    #pdb.set_trace()
    seqStack = sequence.stack(seqArray)
    consensus = seqStack.consensus()
    return consensus

def filterCons(unfilteredConsensus, minConsensusLength):
    '''filter out alignments with short, gappy, mismatching shit (most of them)'''
    # find masked|unmasked block > minConsensusLength
    searchString = (('[acgtACGT]{%i,}') % (minConsensusLength))
    pattern = re.compile(searchString)
    masked = pattern.search(unfilteredConsensus)
    if masked:
        return masked.group()
    else:
        return False

def positioner(candAlign, cons):
    '''return correct positions of the conserved area relative to the reference seq
    by degapping while also dealing with repeat-masked sequence in the conserved area'''
    # strip gap character from reference seq
    pattern = re.compile('-+')
    cleanCandAlign = pattern.sub('', candAlign.text)
    # deal with upper/lowercase issues btw reference <--> alt and
    # repeat-masked bases
    caseUnawareCons = []
    for letter in cons:
        if letter.isupper():
            bracket = (('[%s%s]') % (letter, letter.lower()))
            caseUnawareCons.append(bracket)
        else:
            bracket = (('[%s%s]') % (letter, letter.upper()))
            caseUnawareCons.append(bracket)
    caseUnawareCons = ''.join(caseUnawareCons)
    # find position of conserved sequence relative to gapless
    # candAlign
    pattern = re.compile(caseUnawareCons)
    position = pattern.search(cleanCandAlign)
    return position.span()
    
def createDbase(name):
    '''create a database to hold the results'''
    c = sqlite3.connect(name)
    try:
        # if previous tables exist, drop them
        # TODO: fix createDbase() to drop tables safely
        c.execute('''drop table cons''')
        # commit the change
        c.commit()
    except:
        pass
    # create the primers results table
    c.execute('''create table cons (
    id integer primary key autoincrement,
    seq text unique,
    spp text,
    chromo text,
    aln_start int,
    aln_end int,
    aln_len int,
    cons text,
    cons_start int,
    cons_end int,
    cons_len int,
    map text,
    alt_chromo text,
    alt_start,
    alt_end
    )''')
    c.close()
    
def store(dBase, META):
    '''store the results in sqlite'''
    #pdb.set_trace()
    c = sqlite3.connect(dBase)
    for metadata in META:
        c.execute("insert into cons (seq, spp, chromo, aln_start, aln_end, aln_len, cons, cons_start, cons_end, cons_len, map, alt_chromo, alt_start, alt_end) values (:seq, :spp, :chromo, :aln_start, :aln_end, :aln_len, :cons, :cons_start, :cons_end, :cons_len, :map, :alt_chromo, :alt_start, :alt_end)", META[metadata])
        c.commit()
    c.close()

#TODO:  Move contents of main loop to another function such that we 
# can use pp to run multiple jobs.

def processor(input, minConsensusLength=60, minAlignLength=60, metadataKey='galGal3'):
    file = open(input,'rU')
    parser = bx.align.maf.Reader(file)
    a = parser.next()
    # select only those alignments of > minSpecies
    counter = 0
    META = {}
    while a:
        counter += 1
        candAlign = spScreen(a, minAlignLength)
        if candAlign:
            # create sequence stack and stack -> dumb consensus
            unfilteredConsensus = createCons(candAlign)
            # filter out consensi with < 1 contiguous block of minConsensus
            cons = filterCons(unfilteredConsensus, minConsensusLength)
            if cons:
                #print '%s: ****Valid consensus****' % counter
                #print cons
                # find 'real' positions in reference sequence (galGal3 here)
                # by degapping
                refPosition = positioner(candAlign.components[0], cons)
                # find 'real' positions in alternate sequence (anoCar1 here)
                # by degapping
                altPosition = positioner(candAlign.components[1], cons)
                # get sequence metadata
                metadata = alignMetadata(counter, candAlign, cons, refPosition, altPosition, metadataKey)
                # store start, totalLength, end, consensus somewhere
                META[counter] = metadata
                #try:
                #    store(c, metadata)
                #except:
                #    print 'error storing data'
                    #pdb.set_trace()
        a = parser.next()
    # close the MAF reader
    parser.close()
    # close the file
    file.close()
    return META

def main():
    #################
    directory = '/Users/bcf/Git/seqcap/'
    cpu = 1
    dBase = 'probe.sqlite'
    #################
    startTime = time.time()
    # create a dbase.  there is no return here.
    createDbase(dBase)
    files = []
    # get a list *.maf files from directory 
    [files.append(os.path.join(directory, f)) for f in os.listdir(directory) 
        if os.path.isfile(os.path.join(directory, f)) and f.split('.')[1]=='maf']
    # create the jobServer
    jobServer = pp.Server(ncpus=cpu)
    print "Starting pp with %s workers" % jobServer.get_ncpus()
    # create the worker process stuff for nproc
    jobs = [(file, jobServer.submit(processor, args=(file,), depfuncs=(spScreen, alignMetadata, createCons, filterCons, positioner,), modules=("re", "numpy", "sequence", "bx.align.maf",), callback=store, callbackargs=(dBase,))) for file in files]
    # run the processes
    for file, job in jobs:
        print (("File processed = %s (%s)") % (file, job()))
    jobServer.print_stats()

    endTime = (time.time() - startTime)/60.
    print 'Time for execution = %f min.' % (endTime)
    
if __name__ == '__main__':
    main()
    
'''
Single Job Code:

job = jobServer.submit(processor, args=(os.path.join(directory,files[0]),), depfuncs=(spScreen, alignMetadata, createCons, filterCons, positioner,), modules=("re", "numpy", "sequence", "bx.align.maf",), callback=store, callbackargs=(dBase,))

'''