#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 12:37:34 2023

@author: darwin
"""

"""
Python script to filter KofamScan results based on evalue and score,
retaining only the best hit for each gene. Input should be in the non-tabular format returned by KofamScan.

KofamScan is a tool for functional gene annotation using Hidden Markov Models (HMM) profiles.
This script ensures that only the most significant results for each gene are retained,
optimizing the interpretation and use of KofamScan-generated data.

Example of use:
    filter-kofamscan-for-best-hits.py kofamscan-output.txt
"""

import sys

infilename = sys.argv[1]

evalue_threshold = 1e-5 # all results above this threshold will be discarded
factor_threshold = 1 # all results with a score below the threshold multiplied by this factor will be discarded


print('\n\nWorking with {}...\n'.format(infilename))

outfilename = infilename.replace('.txt', '_filtered.txt')

infile = open(infilename)
filew = open(outfilename, 'w')

new_header = 'Gene_name\tKO\tthrshld\tscore\tE-value\tKO_def\n'
filew.write(new_header)

infile.readline()
infile.readline()

hits = {}

for line in infile:
    if line[0:2] == '* ':
        line = line.strip('* ')
    l = line.split(' ')
    l = [x for x in l if x != '']
    l_info = l[:5]
    
    gene_name = l_info[0]
    KO = l_info[1]
    threshold = l_info[2]
    score = l_info[3]
    evalue = l_info[4]
    
    KO_def = ' '.join(l[5:]).strip('\n')
    
    add = False
    
    
    
    if threshold != '-':
        if float(score) > factor_threshold * float(threshold):
            if float(evalue) < evalue_threshold:
                add = True
    else:
        if float(score) > 100:
            if float(evalue) < evalue_threshold:
                add = True
    if add == True:
        try:
            hits[gene_name].append([KO, threshold, score, evalue, KO_def])
        except KeyError:
            hits[gene_name] = [[KO, threshold, score, evalue, KO_def]]
            
best_hits = {}

for gene in hits:
    best_hit = []
    high_score = -100
    
    for hit in hits[gene]:
        if float(hit[2]) > high_score:
            best_hit = hit
            high_score = float(hit[2])
    best_hits[gene] = best_hit
    
for hit in best_hits:
    string = '\t'.join(best_hits[hit])
    filew.write('{}\t{}\n'.format(hit, string))


infile.close()
filew.close()

print('\n***{} has been created.***\n'.format(outfilename))