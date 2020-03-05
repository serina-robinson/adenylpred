#!/usr/bin/env python
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.


import Bio.SeqIO
import sys


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

if __name__ == "__main__":

    # Parse arguments
    parser = define_arguments()
    args = parser.parse_args()
    fasta_dir = args.input
    silent = args.silent

    # Suppress stdout 
    if silent:
        text_trap = io.StringIO()
        sys.stdout = text_trap

    # Checks if GenBank file and converts to FASTA
    if args.genbank_input:
        out_file = "%s/data/gbk_to_fasta.faa" % parent_folder
        faa_converted = gbk_to_faa(fasta_dir, out_file)
        fasta_dir = faa_converted

    # Checks if nucleotide sequence and converts to amino acid
    if args.nucleotide:
        try:
            out_file = "%s/data/nuc_to_aa.faa" % parent_folder
            seqs = fasta.read_fasta(fasta_dir)
            translated = []
            for seq in seqs.values():
                aa = Seq(str(seq), IUPAC.unambiguous_dna)
                translated.append(aa.translate())
            fasta.write_fasta(seqs.keys(), translated, out_file)
            fasta_dir = out_file
        except:
            print("Error: please check your file is a valid FASTA or GenBank file with the appropriate file suffix (.fasta, .faa, or .gbk). \n If your input is a nucleotide sequence, please check the -n option is set to 1")
            sys.exit(1)

    # Extracts AMP-binding hits from a multi-FASTA file
    if not silent:
        print("##### Extracting AMP-binding domains... #####")
    

    if not args.do_not_xtract_domains:
        out_file = "%s/data/AMP_binding_extracted.fasta" % parent_folder
        try:
            xtract_doms(fasta_dir, HMM_FILE, out_file, 50, True)
        except:
            print("Error: please check your file is a valid FASTA or GenBank file")
            sys.exit(1)
    else:
        out_file = fasta_dir

    # Reads in file
    try:
        create_domain_fa = fasta.read_fasta(out_file)
    except:
        print("Error: please check your file is a valid FASTA or GenBank file")
        sys.exit(1)

    # Create antiSMASH-like domain objects from AMP-binding hits
    domain_list = []
    
    for i, domain in enumerate(create_domain_fa):
        domain_list.append(AntismashDomain(FeatureLocation(1, 1, 1), tool="test")) # arbitrary feature location for testing
        domain_list[i].domain_id = list(create_domain_fa.keys())[i]
        domain_list[i].translation = list(create_domain_fa.values())[i]
    
    # Extract the active site residues
    res, nms, sqs = [], [], []
    new_path = "%s/data/34_aa_xtracted.fasta" % parent_folder
    
    if not silent:
        print("##### Extracting active site residues... #####")
    for i, domain in enumerate(tqdm.tqdm(domain_list, disable = silent)):
        res.append(get_34_aa_signature(domain))
        nms.append(domain.domain_id)
        sqs.append(res[i])

    fasta.write_fasta(nms, sqs, new_path)

    # Make predictions based on 34 active site residues
    if not silent:
        print("##### \n Making predictions... #####")
    results = make_prediction(new_path, silent = silent)

    if silent:
        sys.stdout = sys.__stdout__
    print("Query_name\tPrediction\tProbability_score")
    for x in range(len(results)):
        print("%s\t%s\t%.3f\n" % \
        (results[x].name, results[x].prediction, results[x].probability))
  
    # Write predictions to file
    if args.output is not None:
        with open(args.output, 'w') as output_file:
             output_file.write("Query_name\tPrediction\tProbability_score")
             output_file.write("\n")
             for result in results:
                result.write_result(output_file)
        output_file.close()