#!/usr/bin/env python

#Libraries to import
import sys
import os
import string
import subprocess

from io import StringIO   
from Bio import SearchIO

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

from antismash.common import fasta

#Global settings
global cpu
##Number of CPU cores used for multiprocessing
cpu = 4


def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    proc = subprocess.Popen(commands, stdin=stdin_redir,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    try:
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError as e:
        raise

def run_hmmsearch(query_hmmfile, fasta_file):
    "Run hmmsearch"
    command = ["hmmsearch", "--cpu", "3","--domT", "50", "--domtblout", "temp.out", query_hmmfile, fasta_file]
    try:
        out, err, retcode = execute(command)
        out = open("temp.out","r").read()
    except OSError:
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmsearch3-domtab'))
    return results

def read_fasta_file(fastafile):
    "Get sequences from fasta file and store in dictionary"
    infile = open(fastafile).read()
    entries = infile.split(">")[1:]
    fastadict = {} #Accession as key, sequence as value
    for entry in entries:
        accession = entry.partition("\n")[0].replace("\r","").partition(" ")[0]
        #print accession
        sequence = entry.partition("\n")[2].replace("\n","").replace("\r","")
        #print sequence
        fastadict[accession] = sequence
    return fastadict

def run_HMMer(hmmfile, fastafile):
    "Run HMMer to find start and end sites of AMP-binding domains"
    runresults = run_hmmsearch(hmmfile, fastafile)
    results_by_id = {}
    for runresult in runresults:
        #Store results in dictionary by accession
        for hsp in runresult.hsps:
            if not hsp.hit_id in results_by_id:
                results_by_id[hsp.hit_id] = [hsp]
            else:
                results_by_id[hsp.hit_id].append(hsp)
    return results_by_id

def write_fasta(fastadict, outfile, hmmer_results, size_threshold, standalone=True):
    #For each HMM, print a FASTA file with all domain hits
    domaindict = {}
    size_threshold = int(size_threshold)
    out_file = open(outfile, "w")
    domnr = 1
    for cds in hmmer_results.keys():
        #Get sequence from fastadicts
        cds_sequence = fastadict[cds]

        #For each hit, write out domain name + sequence in FASTA format
        for hit in hmmer_results[cds]:
            domain_name = "%s|dom%s" % (cds, domnr)
            domain_sequence = cds_sequence[hit.hit_start:hit.hit_end]
            if len(domain_sequence) > size_threshold: ##Don't get truncated domains
                out_file.write(">%s\n%s\n" % (domain_name, domain_sequence))
                if not standalone:
                    domaindict[domain_name] = int(domain_name.partition("|")[0].replace("genome",""))
            domnr += 1
        domnr = 1
    out_file.close()
    return domaindict

def xtract_doms(fastafile, hmmfile, outfile, size_threshold, standalone=True):
    """Main function to extract domains from FASTA sequence"""
    #Read FASTA file

    fastadict = read_fasta_file(fastafile)
    
    #Run HMMer
    hmmer_results = run_HMMer(hmmfile, fastafile)
    #Write FASTA file with domains
    
    domaindict = write_fasta(fastadict, outfile, hmmer_results, size_threshold, standalone)
    return domaindict
    

