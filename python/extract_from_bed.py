#!/usr/bin/python3

################################################################################
# repeat_summary.py
# Author: Linnea Smeds
# Date: 27-August-2024
# Description: Takes a bed file and a list of names or regions, and only returns
# the lines that matches the list. Which column to match is given as input
# (starting from 1).
################################################################################
##### Import libraries
import argparse
import re
from datetime import datetime


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script matches bed entries \
with a list of desired elements")
parser.add_argument('-b', '--bedfile', help = 'Bed input file name', type = str, required = True)
parser.add_argument('-l', '--list', help = 'List of wanted elements', type = str, required = True)
parser.add_argument('-c', '--column', help = 'Column to match', type = int, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
args = parser.parse_args()

################################################################################
##### Main code

#COLUMN
col=args.column-1

# SAVE NAMES IN DICT
l_dict={}
with open(args.list, 'r') as infile:
    for line in infile:
        tabs = line.rstrip().split("\t")
        if tabs[0] not in l_dict:
            l_dict[tabs[0]]=1

# GO THROUGH THE BED
with open(args.output, 'w') as outfile:
    with open(args.bedfile, 'r') as infile:
        for line in infile:
            s=re.match("^#", line)
            if s:
                print(line.rstrip(), file=outfile)
            else:
                tabs = line.rstrip().split("\t")
                if tabs[col] in l_dict:
                    print(line.rstrip(), file=outfile)
