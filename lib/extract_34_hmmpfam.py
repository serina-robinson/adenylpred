#!/usr/bin/env python

import subprocess
from sys import argv
import os

from Bio import SearchIO

import subprocess
import os
import argparse

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

from lib.run_hmmpfam import *
from lib.parse_hmmpfam import *

def define_arguments():
    parser = argparse.ArgumentParser(description = "Extract 34 aa using hmmalign.")
    parser.add_argument("-i", type = str, required = True,
                        help = "Input fasta directory.")
    parser.add_argument("-o", type = str, required = True,
                        help = "Output fasta directory.")
    return parser

def extract_34_aa(hmm_dir, fasta_dir, out_dir):
    run_hmmpfam(hmm_dir, fasta_dir, 'temp_hmmsearch2.txt')

    out_file = open(out_dir, 'w')

    for result in SearchIO.parse('temp_hmmsearch2.txt', 'hmmer2-text'):
        for i, hit in enumerate(result.hits):
            if hit.id == 'AMP-binding':
                hsp = result.hsps[i]
                seq = hsp.query.seq
                seq = remove_insertions(seq)
                seq34 = extract_34(seq)
                header = make_header(result.id, result.description)
                out_file.write(">%s\n%s\n" % (header, seq34))
    subprocess.check_call(["rm", "temp_hmmsearch2.txt"])

    
