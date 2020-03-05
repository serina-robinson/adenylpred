#!/usr/bin/env python

import subprocess
import os
import argparse
from sys import argv
from Bio import SearchIO

def define_arguments():
    parser = argparse.ArgumentParser(description = "Extract 34 aa using hmmalign.")
    parser.add_argument("-i", type = str, required = True,
                        help = "Input fasta directory.")
    parser.add_argument("-o", type = str, required = True,
                        help = "Output fasta directory.")
    return parser

def run_hmmpfam(hmm_dir, fasta_dir, out_dir):
    program = 'hmmpfam2'
    options = '-T 50'
    arguments = [hmm_dir, fasta_dir]
    out_file = open(out_dir, 'w')

    command = [program] + [options] + arguments
    subprocess.call(command, stdout = out_file)
    out_file.close()
