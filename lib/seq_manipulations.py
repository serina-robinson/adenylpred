#!/usr/bin/env python
import Bio.SeqIO
import sys
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from lib import fasta

def gbk_to_faa(gbk_dir, out_file):
    out_file = open(out_file, "w")
    for seq_record in Bio.SeqIO.parse(gbk_dir, 'genbank'):
	    for seq_feature in seq_record.features:
		    if seq_feature.type == 'CDS':
			     x = seq_feature.qualifiers
			     assert len(x['protein_id']) == 1
			     assert len(x['product']) == 1
			     assert len(x['translation']) == 1
			     out_file.write('>%s %s\n%s\n' % (
			         x['protein_id'][0],
				     x['product'][0],
				     x['translation'][0]))
    out_file.close()
    return(out_file.name)

def nuc_to_aa(input_dir, out_file):
    seqs = fasta.read_fasta(input_dir)
    translated = []
    for seq in seqs.values():
        aa = Seq(str(seq), IUPAC.unambiguous_dna)
        translated.append(aa.translate())
    fasta.write_fasta(seqs.keys(), translated, out_file)
    return(out_file)