#!/usr/bin/env python
'''
overlap2mutprobs.py was written to determine the potential significance of
multiple de novo mutations found in any given gene across or within studies.
It takes the overlap files produced by multi_dnm.py (or similar programs), a
file containing the mutation probabilities, and the number of individuals in
all studies to create an output file that contains the probability of observing
the number of mutations given the expected number based on gene mutability and
sample size.

Required inputs:
    1) Overlap file created from multi_dnm.py
    2) Per gene mutation probability file
    3) The number of individuals in the study

Output contains: gene, mutation list, # loss of function (LoF), # missense,
prob(LoF), prob(mis), prob(LoF+mis), 2 times the comparison probability,
the expected number of mutations given the number of individuals,
the poisson probability of the observed number of mutations, and the comparison
probability.

Version 1.2 handles probabilities of mutation adjusted for depth and
divergence (mut_prob_fs_adjdepdiv.txt)
'''

__version__ = 1.2
__author__ = 'Kaitlin E. Samocha <ksamocha@fas.harvard.edu>'
__date__ = 'February 12th, 2013'

import sys
import os
import argparse
from rpy2.robjects.packages import importr

def process_overlaps(file1):
    'Open the overlap file and extract data of interest'
    multi_hits = {}

    with open(file1, 'r') as multi_f:
        for line in multi_f:
            line = line.split()
            gene = line[0]
            mutations = line[1].split('/')

            nlof = 0
            nmis = 0
            func_mut = []
            for mut in mutations:
                if mut == 'missense':
                    nmis += 1
                    func_mut.append(mut)
                elif mut in ('nonsense', 'splice', 'frameshift'):
                    nlof += 1
                    func_mut.append(mut)

            if (nmis+nlof) < 2:
                continue

            muts = ', '.join(func_mut)

            multi_hits[gene] = [muts, nlof, nmis]

    return multi_hits


def get_mut_probs(file2):
    'Open the file and save the mutation probabilities'
    probs = {}

    with open(file2, 'r') as probs_f:
        for line in probs_f:
            if line.startswith('transcript'):
                continue
            line = line.split()
            gene = line[1]
            
            if line[6] == 'NA':
                pmis = 0.0
            else:
                pmis = 10**float(line[6])

            if line[7] == 'NA':
                pnon = 0.0
            else:
                pnon = 10**float(line[7])

            if line[9] == 'NA':
                psplice = 0.0
            else:
                psplice = 10**float(line[9])

            if line[10] == 'NA':
                pfs = 0.0
            else:
                pfs = 10**float(line[10])

            plof = pnon + psplice + pfs
            probs[gene] = [plof, pmis]

    return probs


def determine_significance(results, mutprob, num_i, type):
    'Determines the probability of observing the number of mutations'
    if mutprob != '.':
        doubleprob = 2*mutprob
        expnum = num_i*doubleprob

        if type == 'lofmis':
            total_num = (results[2] + results[3]) - 1
        elif type == 'lof':
            total_num = results[2] - 1
        stats = importr('stats')
        pois_prob = stats.ppois(total_num, expnum, lower=False)[0]
        results.extend([doubleprob, expnum, pois_prob])
    else:
        results.extend(['.', '.', '.'])

    return results


def main(args):
    'Control the flow of the script'
    
    thingsToReturn = ''#added, string to collect output
    multi_hits = process_overlaps(args[0])#changing args from argparse to a list
    probs = get_mut_probs(args[1])#same as above

    labels = ['#gene', 'mutations', '#LoF', '#mis', 'prob(LoF)', 'prob(mis)',
              'prob(LoF+mis)', '2*prob', 'exp#[{0}]'.format(args[2]), 'ppois',
              'compared_to']#changed arg format as above
    #print('\t'.join(labels))# labels don't need to be printed

    for gene in multi_hits.keys():
        results = multi_hits[gene]
        results.insert(0, gene)

        try:
            (plof, pmis) = probs[gene]
            plofmis = plof + pmis
        except KeyError:
            plof = '.'
            pmis = '.'
            plofmis = '.'

        results.extend([plof, plofmis, pmis])#changed order from plof, pmis, plofmis to this

        #if results[3] > 0:
        results = determine_significance(results, plofmis, args[2], 'lofmis')
        results.append('LoF+mis')
        thingsToReturn += ('\t'.join(map(str, results)))+'\n'#return string
            #print('\t'.join(map(str, results)))# doesn't need to print

        #if results[2] > 1:
        if results[-1] == 'LoF+mis':
            results = results[:-4]
        l_results = determine_significance(results, plof, args[2], 'lof')
        l_results.append('LoF')
        thingsToReturn += ('\t'.join(map(str, l_results)))+'\n'#returned string
            #print('\t'.join(map(str, l_results)))# doesn't need to print
    return thingsToReturn

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
        Pipeline to take the genes with multiple mutations and determine
        their probabilities for mutation and the probability of observing
        the number of mutations seen
    ''')
    parser.add_argument('multi', action='store', type=str,
                        help='File that contains multiply hit genes')
    parser.add_argument('probs', action='store', type=str,
                        help='Mutation probability file')
    parser.add_argument('num', action='store', type=int,
                        help='Number of individuals in the study')

    args = parser.parse_args()

    if not os.path.exists(args.multi):
        sys.exit('{0}: No such file or directory'.format(args.multi))
    if not os.path.exists(args.probs):
        sys.exit('{0}: No such file or directory'.format(args.probs))

    main(args)
