#!/usr/bin/env python
#
#
# 
#
import argparse
import pysam


def main():
    parser = argparse.ArgumentParser(description="Convert alignment of FASTA to variants for further analysis")
    parser.add_argument("alignment_file",dest="alignment_file", help="Alignment file")



