from Bio import SeqIO
from math import log
import argparse
import json
import numpy
import os
def shannon(seqs):

    N = float(len(seqs))
    S = set(seqs)

    p_i = lambda i: len(filter(lambda x: x == i, seqs)) / N

    return -sum(map(lambda n: p_i(n) * log(p_i(n), 2), S))

def slice_generator(width):

    start = 0
    stop = start + width

    while True:
        yield start, stop
        
        start += width
        stop = start + width

def shannon_kmers(msadir, k):

    msas = [os.path.join(msadir, m) for m in os.listdir(msadir) if '.fasta' in m]
    shannons = { os.path.basename(msa).split('.')[0]: walk_msa(msa, k) for msa in msas }

    return shannons

def walk_msa(msa, k):
    """Walks through the MSA of alleles for a particular gene.

    For each kmer, find the Shannon entropy for that region
    of the alignment.
    """

    shannons = []
    with open(msa, 'r') as f:

        # inelegant hack to get MSA length
        length = len(list(SeqIO.parse(f, 'fasta'))[0]) 
       
        start = 0
        while start < length:
            start, stop = walk_msa(msa, k)
            current = []
            
            for record in SeqIO.parse(f, 'fasta'):
                current.append(record.seq[start, stop])
            
            shannons.append(shannon(current))
    
    return shannons, length

def contig_locs(json_dir, test_name):
    
    def load_json(jsonpath):
        with open(jsonpath, 'r') as f:
            data = json.load(f)
        return data
    
    def combine_dictionaries(dict_list):

        d = {}

        for i in dict_list:
            for elem in i:
                try:
                    d[elem].append(i[elem])
                except KeyError:
                    d[elem] = [i[elem]]

        return d                            

    def parse_data(data, test_name):
        
        gene_locations = {}

        genes = [data]["Results"][0]["TestResults"][test_name]

        for gene in genes:

            if genes[gene]["BlastResults"] != None:

                if genes[gene]["IsContigTruncation"]:

                    qlen = genes[gene]["BlastResults"]["QueryLength"]
                    alnlen = genes[gene]["BlastAlignmentLength"]
                    start_index = genes[gene]["StartIndex"]

                
                    if start_index <  0:

                        loc = -start_index

                    else:
                        loc = alnlen

                    gene_locations[gene] = loc 

    jsons = [os.path.join(json_dir, j) for j in os.listdir(json_dir) if '.json' in j]
    
    locs = [parse_data(load_json(j), test_name) for j in jsons]

    return combine_dictionaries(locs)

def combine_contigs_shannon(shannons, contigs):
    """Combine shannon and contig sequences into 2D - tuples 
    for easy creation of Genetic Algorithm individuals later.
    """

    combined = {}

    for gene in shannons:

        bins = numpy.linspace(0, len(gene), len(gene))
        contig_kmer = numpy.histogram(contigs[gene], bins = bins)[0]
        
        z = zip(shannon[gene], contig_kmer)

        combined[gene] = z

    return combined

def arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument('--jsondir', required = True)
    parser.add_argument('--msadir', required = True)
    parser.add_argument('--test', required = True)
    parser.add_argument('-k', type = int, required = True)

    return parser.parse_args()

def main():

    args = arguments()
    
    shannons = shannon_kmers(args.msadir, args.k)
    contigs = contig_locs(args.jsondir, args.test)

    combined = combine_contigs_shannon(shannons, contigs)

if __name__ == '__main__':
    main()
