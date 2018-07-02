#!/usr/bin/env python

'''
Given a gene sequence file from 1KP and a corresponding table mapping the 4-letter codes
to different taxonomic groups, return the number of parsimony informative characters
for a specific clade.
'''

import os,sys
from Bio import SeqIO

def parsimony_informative(seqs):
    pars_inf_sites = 0
    num_sites = 0
    num_var_sites = 0
    for pos in range(len(seqs[0])):
        column = "".join([seqs[i].seq[pos] for i in range(len(seqs))]).upper()
        a = column.count('A')
        c = column.count('C')
        t = column.count("T")
        g = column.count("G")
        if sum([x > 1 for x in [a,c,t,g]]) > 1:
            pars_inf_sites += 1
        if sum([x > 0 for x in [a,c,t,g]]) > 1:
            num_var_sites  += 1
        if sum([a,t,c,g]) > 0:
            num_sites += 1

    return len(seqs), num_sites,num_var_sites#, num_variable_sites #pars_inf_sites
    
def main():
    mapping_filename = "/Users/mjohnson/onedrive/Projects/AngiospermHybSeq/Analysis/fine_annotations_MJ.txt"
    order_col = 1
    family_col = 2
    genus_col = 3


    order_dict = {}
    family_dict = {}
    genus_dict = {}

    id_to_taxon = open(mapping_filename).readlines()

    header = id_to_taxon.pop(0)
    for line in id_to_taxon:
        line = line.split("\t")
    #    print line
        try:
            order_dict[line[order_col]].append(line[0])
        except KeyError:
            order_dict[line[order_col]] = [line[0]]
        try:
            family_dict[line[family_col]].append(line[0])
        except KeyError:
            family_dict[line[family_col]] = [line[0]]
        try:
            genus_dict[line[genus_col]].append(line[0])
        except KeyError:
            genus_dict[line[genus_col]] = [line[0]]

    

    #print order_dict

    order = sys.argv[1]
    family = sys.argv[2]
    genus = sys.argv[3]

    order_species_list = order_dict[order]
    family_species_list = family_dict[family]
    genus_species_list = genus_dict[genus]

    #gene_list = [x for x in os.listdir("genes") if os.path.isdir(os.path.join("genes",x))]
    gene_list = [x.rstrip() for x in open("genes_for_probes.txt")]

    for gene in gene_list:
        seqs = [seq for seq in SeqIO.parse(os.path.join("genes",gene,"FNA2AA-upp-masked.fasta"),'fasta')]
        order_seqs = [seq for seq in seqs if seq.id in order_species_list]
        family_seqs = [seq for seq in seqs if seq.id in family_species_list]
        genus_seqs = [seq for seq in seqs if seq.id in genus_species_list]
        gene_len = "NA"
        if len(order_seqs) > 3:
            taxa, gene_len, order_parsinf = parsimony_informative(order_seqs)
        else:
            order_parsinf = "NA"
        if len(family_seqs) > 3:
            taxa, gene_len, family_parsinf = parsimony_informative(family_seqs)
        else:
            family_parsinf = "NA"
        if len(genus_seqs) > 3:
            taxa, gene_len, genus_parsinf = parsimony_informative(genus_seqs)
        else:
            genus_parsinf = "NA"
        
        #sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(gene,gene_len,order_parsinf,family_parsinf,genus_parsinf))    
        #sys.stderr.write("{}\t{}\t{}\t{}\n".format(gene,len(order_seqs),len(family_seqs),len(genus_seqs)))
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(gene,gene_len,len(genus_seqs),genus_parsinf))
        
        
if __name__ == "__main__":main()
