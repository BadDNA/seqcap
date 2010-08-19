#!/usr/bin/env python
# encoding: utf-8
"""
summaryBlast.py

Created by Brant Faircloth on 2008-08-03.
Copyright (c) 2008 Brant Faircloth. All rights reserved.
"""
import pp
import re
import os
import pdb
import time
import sqlite3
import Bio.Blast.NCBIXML
import Bio.Blast.NCBIStandalone

def sqlQuery(sqlDb):
    '''get sequences from sqlite and return'''
    c = sqlite3.connect(sqlDb)
    query = c.execute('''select seq, cons from cons''')
    conservedSeq = query.fetchall()
    return conservedSeq

def tempFile(seqId, strand):
    '''write sequences and their ids out to a temporary fasta file'''
    tempName = '.%s.fa' % seqId
    temp__handle =  open(tempName, 'w')
    temp__handle.write(('>%s\n%s') % (seqId, strand))
    temp__handle.close()
    return tempName

def blastQuery(tempName, blastExe, blastDb):
    '''query blast for matching sequences'''
    result__handle, error__handle = Bio.Blast.NCBIStandalone.blastall(blastExe, "blastn", blastDb, tempName)
    blastResults = Bio.Blast.NCBIXML.parse(result__handle)
    return blastResults.next()

def parseBlast(seqId, strand, blastResults, eThresh):
    alnScore = {}
    counter = 0
    for alignment in blastResults.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < eThresh:
                print '****Alignment****'
                print 'sequence:', seqId
                print 'e value:', hsp.expect
                print hsp.query
                print hsp.match
                print hsp.sbjct
                # we want sequences that match in their entirety
                if len(hsp.query) == hsp.identities == len(strand):
                    # exact match btw. subject and query
                    perfect     = True
                    types       = None
                    positions   = None

                else:
                    pattern     = re.compile('\s+')
                    match       = pattern.search(hsp.match)
                    if not match and hsp.identities != len(strand):
                        perfect     = True
                        types       = 'partial'
                        positions   = None
                    else:
                        iterator    = pattern.finditer(hsp.match)
                        # gappy or wobbly
                        perfect     = False
                        positions, types   = [],[]
                        for match in iterator:
                            location = match.span()[0]
                            positions.append(str(location))
                            if hsp.query[location]=='-' or hsp.sbjct[location]=='-':
                                types.append('g')
                            else:
                                types.append('w')
                        types = ','.join(types)
                        positions = ','.join(positions)
                alnScore[counter] = (seqId, counter, hsp.expect, perfect, types, positions)
                counter += 1
    return alnScore

def addTable(sqlDb):
    '''add a table to the database'''
    c = sqlite3.connect(sqlDb)
    try:
        # if previous tables exist, drop them
        # TODO: fix createDbase() to drop tables safely
        c.execute('''drop table blast''')
        # commit the change
        c.commit()
    except:
        pass
    # create the primers results table
    c.execute('''create table blast (
    id integer primary key autoincrement,
    seq text,
    matches integer,
    match_cnt integer,
    e_value float,
    perfect text,
    types text,
    positions text
    )''')
    c.close()

def store(sqlDb, alnScore):
    '''store the results in sqlite'''
    #pdb.set_trace()
    if alnScore:
        c = sqlite3.connect(sqlDb)
        for aln in alnScore:
            row = alnScore[aln]
            row += (len(alnScore),)
            print row
            print '\n\n'
            try:
                c.execute('insert into blast (seq, match_cnt, e_value, perfect, types, positions, matches) values (?,?,?,?,?,?,?)', row)
            except:
                #rollback the pending transaction
                c.rollback()
                # wait for the dbase lock
                time.wait(0.2)
                c.execute('insert into blast (seq, match_cnt, e_value, perfect, types, positions) values (?,?,?,?,?,?)', row)
            c.commit()
        c.close()
    else:
        pass


def processor(seq):
    blastDb     = '/Users/bcf/lib/blast/data/taeGut3'
    blastExe    = '/Users/bcf/bin/blast-2.2.18/bin/blastall'
    # set e-value low to drop spurios alignments
    eThresh     = 1e-15
    seqId, strand  = seq
    # write a tempfile for the sequence
    tempName    = tempFile(seqId, strand)
    # blast the query against the subject
    blastResults = blastQuery(tempName, blastExe, blastDb)
    # search the results for mismatches
    alnScore = parseBlast(seqId, strand, blastResults, eThresh)
    # erase the tempfile
    os.remove(tempName)
    return alnScore

def main():
    startTime       = time.time()
    sqlDb           = 'probe.sqlite'
    jobServer       = pp.Server(ncpus=6)
    print "Starting pp with %s workers" % jobServer.get_ncpus()
    conservedSeq    = sqlQuery(sqlDb)
    addTable(sqlDb)
    jobs            = [(seq, jobServer.submit(processor, args=(seq,), depfuncs=(tempFile, blastQuery, parseBlast,), modules=("Bio.Blast.NCBIXML", "Bio.Blast.NCBIStandalone","re","os","time"), callback=store, callbackargs=(sqlDb,))) for seq in conservedSeq]
    #
    for seq, job in jobs:
        print (("File processed = %s (%s)") % (seq, job()))
    jobServer.print_stats()
    endTime         = (time.time() - startTime)/60.
    print 'Time for execution = %f min.' % (endTime)

if __name__ == '__main__':
    main()