#!/usr/bin/python3

################################################################################
# print_lines.py
# Author: Linnea Smeds
# Date: 27-September-2024
# Description: Takes a file with or without header, and a comma separated list 
# of numbers, then it prints the lines corresponding to these numbers in the 
# input file. 
################################################################################
##### Import libraries
import argparse
import re


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script prints desired lines \
    in an input file, given as a comma separated list of numbers. ")
parser.add_argument('-i', '--input', help = 'A text file with or without header', type = str, required = True)
parser.add_argument('-s', '--skipheader', help = 'Set this if input file as header \
    that should be excluded', action='store_true')
parser.add_argument('-r', '--rownumbers', help = 'List of numbers separated by ","', type = str, required = True)
args = parser.parse_args()

################################################################################
##### Main code

#Save a list with the columns numbers 
col=list(map(int, args.rownumbers.split(",")))

# Go through list, starting on first line or not depending on skipheader flag 
cnt=1
if args.skipheader:
    cnt=0

with open(args.input, 'r') as infile:
    for line in infile:
        if cnt in col:
            print(line.rstrip())
        cnt=cnt+1