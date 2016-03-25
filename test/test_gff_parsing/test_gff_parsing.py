#!/usr/bin/python


import sys
import getopt
import os
import DNA_AA_utilities
import GFF_utilities

INPUT_DATA_LOCATION = "../../data/test_human.gff"

def main():
    gff_dict = GFF_utilities.parseGFF(INPUT_DATA_LOCATION)

    for chromosome in gff_dict:
        print "Chromosome - " + chromosome

        all_genes_in_chromosome = gff_dict[chromosome]


        for gene in all_genes_in_chromosome:
            print "gene - " + str(gene[2])




if __name__ == "__main__":
    main()
