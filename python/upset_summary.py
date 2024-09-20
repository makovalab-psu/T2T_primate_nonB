#!/usr/bin/python3

################################################################################
# upset_summary.py
# Author: Linnea Smeds
# Date: 13-September-2024
# Description: Goes through bedfiles of each non-B motif type (listed below)
# adding each base to a joint matrix with the types as columns and positions 
# as rows.
################################################################################
##### Import libraries
import argparse
import re
from datetime import datetime


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "Finds overlaps of non-B motifs \
from multiple bedfiles")
parser.add_argument('-b', '--bedprefix', help = 'Prefix of the bed files', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
args = parser.parse_args()

################################################################################
##### Main code

nonB=["APR", "DR", "GQ", "IR", "MR", "STR", "Z"]

# prepare dict with positions 
p_dict={}

# Go through the bed files 
for i in range(len(nonB)): 
    file=args.bedprefix+nonB[i]+".bed"
    print("looking at the file: "+file)
    with open(file, 'r') as infile:
        for line in infile:
            tabs = line.rstrip().split("\t")
            for j in range(int(tabs[1]),int(tabs[2])):
                pos=tabs[0]+":"+str(j)
                if pos in p_dict:
                    p_dict[pos][i]="1"
                else:
                    p_dict[pos]=["0"]*len(nonB)
                    p_dict[pos][i]="1"
            

# Go through dict and count the different types
c_dict={}
for k in p_dict:
    type="-".join(p_dict[k])
    if type in c_dict:
        c_dict[type]=c_dict[type]+1
    else:
        c_dict[type]=1

# Delete original dict and print summary
p_dict={}
with open(args.output, 'w') as outfile:
    for k in c_dict:
        name=""
        new_list=k.split("-")
        for i in range(len(new_list)):
            if new_list[i] == "1":
                if name == "":
                    name=nonB[i]
                else:
                    name=name+"-"+nonB[i]

        print(name+"\t"+str(c_dict[k]), file=outfile)