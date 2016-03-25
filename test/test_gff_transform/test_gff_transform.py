#!/usr/bin/python


import sys
import json
import os
import GFF_utilities

INPUT_DATA_LOCATION = "../../data/test_human.gff"

def main():
    gff_dict = GFF_utilities.parseGFF(INPUT_DATA_LOCATION)
    gff_dict2 = GFF_utilities.transform_gff_dict(gff_dict)





if __name__ == "__main__":
    main()
