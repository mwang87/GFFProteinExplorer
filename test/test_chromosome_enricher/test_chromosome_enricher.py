#!/usr/bin/python


import sys
import json
import os
import GFF_utilities

INPUT_GFF_LOCATION = "../../data/test_human.gff"
INPUT_CHROMOSOME_FOLDER = "../../data/chromosomes"

def main():
    gff_dict = GFF_utilities.parseGFF(INPUT_GFF_LOCATION)
    gff_dict2 = GFF_utilities.transform_gff_dict(gff_dict)

    chromosome_map = GFF_utilities.parse_chromosomes(INPUT_CHROMOSOME_FOLDER)
    gff_dict_with_chromosomes = GFF_utilities.enrich_gff_dict_with_chromosome_data(gff_dict2, chromosome_map)




if __name__ == "__main__":
    main()
