from Bio import SeqIO
import json
import os
import argparse
from math import e
import random
import csv

def arguments():

    parser = argparse.ArgumentParser()
    
    parser.add_argument('--jsons', '-j', required = True, help = "Path to JSON directory")
    parser.add_argument('--global', action = 'store_true', help = 'Consider all known alleles')
    parser.add_argument('--test', '-t', required = True, help = "Name of MIST test")
    parser.add_argument('--replicates', default = 1000, help = "Number of replicates to perform")
    parser.add_argument('--out', required = True, help = "Output file path")

    return parser.parse_args()

def load_json(path):

    with open(path, 'r') as f:
        data = json.load(f)
    return data

def gather_sequences(testname, jsondir):
    """Walks through all jsons in jsondir and creates a dict of
    known alleles for each gene.
    
    Alleles are collected using the SubjAln field.
    """

    gene_alleles = {}

    jsons = (load_json(os.path.join(jsondir, j)) for j in os.listdir(jsondir) if '.json' in j)

    for strain in jsons:

        genes = strain["Results"][0]["TestResults"][testname]

        for gene in genes:
            
            if genes[gene]["BlastResults"] != None:
                current = genes[gene]["BlastResults"]["SubjAln"]

                try:
            
                    gene_alleles[gene].add(current)
            
                except KeyError:
                
                    gene_alleles[gene] = set([current])

    return gene_alleles


def bootstrap(gene_alleles, replicates):
    """Generates fragments of observed alleles to be synthetic
    contig truncations.
    
    This is used to determine how many observed alleles
    the fragment matches and the function returns the reciprocal.
    """
    
    def evalue(m, n):

        # constants for blastn
        K = 0.460 # Kappa
        L = 1.280 # Lambda
        S = m     # true if 100% identity

        E = K * m * n * (e**(-L * S))

        return E

    def truncate_and_match(alleles):
       
        n = 1641481 # NCTC11168 length
        tuple_alleles = tuple(alleles)

        good = False
        
        while not good:
            
            current = random.choice(tuple_alleles)
            pivot = random.randrange(len(current))

            if random.randrange(2):

                m = current[pivot:]
            
            else:
                
                m = current[:pivot]
            
            if evalue(len(m), n) <= 10.0:
                good = True

        return 1.0 / sum(m in a for a in alleles)
    
    return {gene : [truncate_and_match(gene_alleles[gene]) for x in range(replicates)]
            for gene in gene_alleles}

def write_out(values, outpath):

    with open(outpath, 'w') as f:
        
        out = csv.writer(f)
        
        header = values.keys()
        out.writerow(header)

        for row in range(len(values[header[0]])):

            out.writerow( [values[x][row] for x in values] )
        
        out.writerow(['' for x in headers])

def main():

    args = arguments()

    sequences = gather_sequences(args.test, args.jsons)
    
    values = bootstrap(sequences, args.replicates)

    write_out(values, args.out)

if __name__ == '__main__':
    main()
