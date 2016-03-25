'''
Created on Mar 24, 2016

@author: Mingxun Wang

This file contains functions to transform a GFF file into a proteome
'''


def gff_to_proteins(gff_dict):
    for chromosome in gff_dict:
        print "Chromosome - " + chromosome

        all_genes_in_chromosome = gff_dict[chromosome]


        for gene in all_genes_in_chromosome:
            print "gene - " + str(gene[2])
