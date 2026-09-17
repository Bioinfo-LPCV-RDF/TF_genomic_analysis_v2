#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
TFFM Sequence Scanner
=====================

DESCRIPTION:
This script scans DNA sequences using a Transcription Factor Flexible Model (TFFM)
to identify potential binding sites. It loads a pre-trained TFFM model from an XML file
and scans a FASTA file containing positive sequences to find matches.

RELATED TO compute_ROCS

USAGE:
The script is designed to be called from the command line with the following arguments:
```bash
python tffm_scanner.py --output <output_file> --fasta_pos <fasta_file> --tffm <tffm_xml>
```

ARGUMENTS: 
--output, -o        : Path to the output file for hit positions. Required.
--fasta_pos, -pos   : Path to the FASTA file with positive sequences. Required.
--tffm, -t          : Path to the TFFM XML model file. Required.

OUTPUT:
- <output_file>     : Text file containing positions of all identified binding sites

NOTES:
- The script uses the TFFM's scan_sequences() method with only_best=True to find
  the best hit per sequence.
- Each hit position is written as a separate line in the output file.
- The TFFM model should be compatible with the sequences being scanned (same organism,
  appropriate model order, etc.).
- For large FASTA files, the scanning process may take significant time and memory.
"""

# IMPORTS
from os import system, mkdir
from os.path import isdir, basename, dirname, join
import sys
sys.path.append("/home/312.3-StrucDev/312.3.1-Commun/lib/TFFM")
import tffm_module
from constants import TFFM_KIND
import argparse

# LOAD ARGUMENTS
parser = argparse.ArgumentParser()

parser.add_argument("--output", "-o", type=str)
parser.add_argument("--fasta_pos", "-pos", type=str)
parser.add_argument("--tffm", "-t", type=str)

args = parser.parse_args()
output = args.output
# DEBUG PRINTS
# print("TFFM: {}".format(args.tffm))
# print("FASTA Pos: {}".format(args.fasta_pos))

# MAIN 
tffm = args.tffm
fasta_pos = args.fasta_pos

hit_pos_list = []
tffm_first_order = tffm_module.tffm_from_xml(tffm, TFFM_KIND.FIRST_ORDER)

try:
    for hit in tffm_first_order.scan_sequences(fasta_pos, only_best=True):
        if hit is None:
            print("No hit found for this sequence.")
        else:
            # print("Hit Position: {}".format(hit)) #DEBUG
            hit_pos_list.append(str(hit))
except Exception as e:
    print("Error during sequence scanning: {}".format(e))

# DEBUG
# print("Hit positions list: {}".format(hit_pos_list))

try:
    with open(output, "w") as f1:
        for hit_pos in hit_pos_list:
            # print("Writing Hit Position: {}".format(hit_pos)) #DEBUG
            f1.write(hit_pos + "\n")
except Exception as e:
    print("Error writing to output file: {}".format(e))
