#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Nick Crawford on 2010-05-18.
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
import subprocess


def main():
    for count in range(0,1000):
        
        phybase_cmd = "python phybase.py -i /n/home06/ngcrawford/data1/seqcap/Manuscripts/Tetrapods/Data/917Loci_19Species/alignments_from_loci/917Loci_19Species_alignments_bootreps/bootrep_%s.phylip_phyml_tree.txt -o /n/home06/ngcrawford/data1/seqcap/Manuscripts/Tetrapods/Data/917Loci_19Species/alignments_from_loci/917Loci_19Species_alignments_bootreps -g anoCar2 -t 'anoCar2 taeGut1 monDom5 mm9 cavPor3 loxAfr3 oryCun2 bosTau4 calJac3 rheMac2 ponAbe2 gorGor3 chinese hg19 korean venter panTro2 canFam2 equCab2'" % count

        bsub_cmd = ['bsub',
         '-q',
         'short_serial',
         '-R',
         'select[mem>3000]',
         '-o',
         '/n/home06/ngcrawford/sterr.out',
         phybase_cmd]
         
        p = subprocess.Popen(bsub_cmd)


if __name__ == '__main__':
    main()
