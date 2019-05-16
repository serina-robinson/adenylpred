#!/usr/bin/env python
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
AdenylPred: script to predict substrate given an amino acid sequence
"""
from sys import argv, path
from statistics import stdev
from typing import Any, Dict, List, Tuple, Optional, List, Set
from sklearn.ensemble import RandomForestClassifier
from joblib import dump, load
from collections import defaultdict
from helperlibs.wrappers.io import TemporaryDirectory
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import random
import argparse
import os
import copy
import pickle
import numpy
import tqdm
import sys
import warnings
import io

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

from lib.get_seq_properties import *
from lib.make_test_set import *
from lib.test_classifier import *
from lib.hmmsearch_domains import *
from lib.get_34_aa_signature import *
from lib.gbk_to_faa import *

from antismash.common.secmet import AntismashDomain, FeatureLocation
from antismash.common import path, subprocessing

PROPERTIES_15 = "%s/data/15_aa_properties.txt" % parent_folder
PROPERTIES_4 = "%s/data/aa_properties.txt" % parent_folder
REF_SEQUENCE = "P0C062_A1"
A34_POSITIONS_FILENAME = "%s/data/A34positions.txt" % parent_folder
APOSITION_FILENAME = "%s/data/Apositions.txt" % parent_folder
KNOWN_CODES = "%s/data/knowncodes.fasta" % parent_folder
ILLEGAL_CHARS = "!@#$%^&*(){}:\"<>?/.,';][`~1234567890*-+-=_\\|"
ADOMAINS_FILENAME = "%s/data/A_domains_muscle.fasta" % parent_folder
START_POSITION = 66
HMM_FILE = "%s/data/AMP-binding.hmm" % parent_folder


def define_arguments() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description = "Prediction tool for adenylate-forming enzyme substrate specificity")
    parser.add_argument("-i", "--input", required = True, type = str, help = "Input file (FASTA or GenBank format).")
    parser.add_argument("-o", "--output", required = False, type = str, help = "Output file directory. Default is stdout")
    parser.add_argument("-s", "--silent", required = False, action = "store_true", help = "Silences all progress updates to stdout")
    parser.add_argument("-n", "--nucleotide", required = False, action = "store_true", help = "Nucleotide sequence")
    parser.add_argument("-g", "--genbank_input", required = False, action = "store_true", help = "Input is in GenBank format")
    return parser

class PredictRFResult():
    """ Holds all relevant results from Random Forest substrate prediction """
    def __init__(self, seqname, prediction, probability):
        self.name = seqname
        self.prediction = prediction
        self.probability = probability

    def write_result(self, out_file):
        out_file.write("%s\t%s\t%.3f\n" % \
                       (self.name, self.prediction, self.probability))

def get_feature_matrix(fasta_dir: str, property_dict: Dict[str, PhysicochemicalProps]) -> Tuple[List[List[float]], List[str]]:
    """ Converts 34 amino acids to physicochemical properties

    Arguments: 
        fasta_dir: directory of input fasta file
        property_dict: dict of {aa: PhysicochemicalProps, ->}, with aa str
    """
    fasta_dict = fasta.read_fasta(fasta_dir)

    features = []
    response = []
    seq_IDs = []
    
    for ID in fasta_dict:
        sequence = fasta_dict[ID]
        specificity = extract_specificity(ID)
        property_vector = extract_features(property_dict, sequence)
        features.append(property_vector)
        response.append(specificity)
        seq_IDs.append(ID)

    feature_list = []
    feature_list.append(seq_IDs)
    feature_list.append(features)
    return(feature_list)

def extract_features(property_dict: Dict[str, PhysicochemicalProps], sequence: str) -> List[float]:
    """Return list of float, with each float a physicochemical property

    Arguments:
        property_dict: dict of {aa: PhysicochemicalProps, ->}, with aa str
        sequence: str, amino acid sequence of length 34
    """

    property_vector = []

    for aa in sequence:
        for aa_property in property_dict[aa].get_all():
            property_vector.append(aa_property)

    return property_vector

def make_prediction(fasta_dir, silent):
    """ Makes a substrate prediction given an input of 34 active site residues

        Arguments:
            fasta_dir: directory of input fasta file 
    """

    filename = "%s/data/models/finalized_rf_model.sav" % parent_folder
    rf_mod = pickle.load(open(filename, 'rb'))

    properties_4 = parse_4_properties(PROPERTIES_4)
    feature_list = get_feature_matrix(fasta_dir, properties_4)

    seqname = feature_list[0]
    features = feature_list[1]
    prediction = rf_mod.predict(features)
    probability = numpy.amax(rf_mod.predict_proba(features), axis = 1)

    results = []
    for i,res in enumerate(tqdm.tqdm(prediction, disable = silent)):
        result = PredictRFResult(seqname[i], prediction[i], probability[i])
        results.append(result)

    return results

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
    out_file = "%s/data/AMP_binding_extracted.fasta" % parent_folder

    try:
        xtract_doms(fasta_dir, HMM_FILE, out_file, 50, True)
    except:
         print("Error: please check your file is a valid FASTA or GenBank file")
         sys.exit(1)

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