#!/usr/bin/env python
# encoding: utf-8

"""
Lastz.py

Created by Brant Faircloth on 01 May 2010 19:01 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  This is a helper class to simplify calls to lastz.

USAGE:  

    import Lastz
    
    lastz = Lastz.Align(target, query, matchcount, identity, output)
    lastz_stderr, lastz_stdout = lastz.run()
    results_file = lastz.output

"""

import os
import tempfile
import subprocess

class Align():
    '''docstring for lastz'''
    def __init__(self, target, query, matchcount, identity, out=False):
        # if not an output file, create a temp file to hold output
        if not out:
            fd, self.output = tempfile.mkstemp(suffix='.lastz')
            os.close(fd)
        else:
            self.output = out
        self.cli = 'lastz {0}[multiple,nameparse=full] {1}[nameparse=full]\
            --strand=both \
            --seed=12of19 \
            --transition \
            --nogfextend \
            --nochain \
            --gap=400,30 \
            --xdrop=910 \
            --ydrop=8370 \
            --hspthresh=3000 \
            --gappedthresh=3000 \
            --noentropy \
            --coverage={2} \
            --identity={3} \
            --output={4} \
            --format=general-:score,name1,strand1,zstart1,end1,length1,name2,\
strand2,zstart2,end2,length2,diff,cigar,identity,\
continuity'.format(target, query, matchcount, identity, self.output)
    
    def run(self):
        lastz_stdout, lastz_stderr = subprocess.Popen(self.cli, shell=True, \
            stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        return lastz_stdout, lastz_stderr
        

if __name__ == '__main__':
    pass