#!/usr/bin/env python
# encoding: utf-8
"""
Phylo.py

Created by Nicholas Crawford on 2010-03-29.
Copyright (c) 2010 Boston Univeristy. All rights reserved.

Phylogenetic module.
"""

import sys
import os
import numpy as np
from copy import deepcopy 

def hasAllTaxa(tree, expected_taxa):
    """ checks if tree has all expected taxa. Returns 'False' if taxa are missing"""
    has_all_taxa = True
    for taxon in expected_taxa:
        if tree.count(taxon) != 1:
            print 'tree missing %s' % taxon
            has_all_taxa = False
    return has_all_taxa
    
    
def bootstrap(sample, replicates=10):
    """bootstrap(initial_sample, replicates=10)

    Create boostrapped replicates of an a numpy array.  Results
    are returned in a python list.
    
    Parameters
    ----------
    sample : array_like
        An array, any object exposing the array interface, an
        object whose __array__ method returns an array, or any
        (nested) sequence. 
    
    replicates : int
        The number of bootstrap replicates to produce.
    
    Examples
    --------
    >>> bootstrap(array([1,2,3,4,5]),5)
    [array([1, 5, 1, 2, 2]),
     array([1, 5, 1, 4, 2]),
     array([2, 2, 2, 4, 4]),
     array([2, 2, 4, 4, 2]),
     array([2, 4, 3, 3, 3])]
    
    
    """
    replicates = int(replicates)
    sample_size = len(sample)
    
    bootstrap_replicates = []
    for boot_rep_num in range(0,replicates):
        choices = np.random.random_integers(0, sample_size-1, sample_size)  # generate index array of random choices
        
        if type(sample) == list:
            boot_rep = []
            for choice in choices:
                element = deepcopy(sample[choice])
                boot_rep.append(sample[choice])
        else:
            boot_rep = sample[choices]
        
        bootstrap_replicates.append(boot_rep)
    return bootstrap_replicates

def branch_lengths_2_decimals(str_newick_tree):
    """replaces branch lengths in scientific notation with decimals"""
    colon_s = 0
    comma_back_paren_s = 0
    num = ''
    new_tree = ''
    for count, char in enumerate(str_newick_tree):
        if char == ':': 
            colon_s = count
            continue

        if char in (')',','): 
            comma_back_paren_s = 1
            num = '%f' % float(num)
            new_tree += ":" + num 
            colon_s = 0
            num = ''

        if colon_s != 0:
            num = num + char

        if colon_s == 0:
            new_tree += char
    new_tree += ";"
    return new_tree

def main():
    pass


if __name__ == '__main__':
    main()
