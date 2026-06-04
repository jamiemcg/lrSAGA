#!/usr/bin/env python

import sys
import argparse
import os

# Convert GFA assembly graph to fasta file, optionally rename sequence identifiers

parser = argparse.ArgumentParser(description = "Convert GFA assembly graph to fasta file")

parser.add_argument("-i", "--input", type=str, help="Input GFA file", required=True)
parser.add_argument("-o", "--output", type=str, help="Output GFA file", required=True)
parser.add_argument("--rename", help="Rename sequence identifiers", action="store_true")

args = parser.parse_args()

if not os.path.isfile(args.input):
    print("ERROR, input file does not exist:", args.input)
    sys.exit(1)    

if os.path.exists(args.output):
    print("ERROR, output file already exists:", args.output)
    sys.exit(1)

output_file = open(args.output, "w")

if not args.rename:
    with open(args.input, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            if line[0] == "S":
                output_file.write(">" + line[1] + "\n")
                output_file.write(line[2] + "\n")

else:
    n = 1

    # Store original names for tracking
    output_file_names = open(args.output + "_names.tsv", "w")
    output_file_names.write("NEW_NAME\tORIGINAL_NAME\n")

    with open(args.input, "r") as f:
        for line in f:
            line = line.strip().split("\t")
            if line[0] == "S":
                output_file.write(">NODE_" + str(n) + "\n")
                output_file.write(line[2] + "\n")

                output_file_names.write("NODE_" + str(n) + "\t" + line[1] + "\n")
                
                n += 1
    
    output_file_names.close()

output_file.close()
