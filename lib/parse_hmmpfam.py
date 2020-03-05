#!/usr/bin/env python

from Bio import SearchIO
from sys import argv

def extract_34(seq):
    new_seq = []
    indices = [163, 166, 167, 183, 187, 188, 189, 190, 191, 192, 193, 196, 231, 232, 258, 259, 260, 261, 262, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296]

    for index in indices:
        new_seq.append(seq[index])

    seq34 = ''.join(new_seq)
    return seq34

def remove_insertions(seq):
    new_seq = []
    for character in seq:
        if not character.islower():
            new_seq.append(character)

    new_seq = ''.join(new_seq)
    return new_seq
    

def make_header(ID, description):
    description_elements = description.split()
    header = [ID] + description_elements
    return '_'.join(header)

if __name__ == "__main__":

    out_dir = argv[2]
    out_file = open(out_dir, 'w')
    
    for result in SearchIO.parse(argv[1], 'hmmer2-text'):
        print(result.id)
        print(result.description)
        for i, hit in enumerate(result.hits):
            if hit.id == 'AMP-binding':
                hsp = result.hsps[i]
                seq = hsp.query.seq
                seq = remove_insertions(seq)
                print(seq)
                seq34 = extract_34(seq)
                header = make_header(result.id, result.description)
                out_file.write(">%s\n%s\n" % (header, seq34))

    out_file.close()
