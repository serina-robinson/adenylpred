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

import random
import argparse
import os
import copy
import pickle
import numpy
import tqdm
import sys
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

from lib.get_seq_properties import *
from lib.make_test_set import *
from lib.test_classifier import *
from lib.hmmsearch_domains import *
from lib.get_34_aa_signature import *

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
    parser = argparse.ArgumentParser(description = "Run random forest prediction.")
    parser.add_argument("-i", "--input", required = True, type = str, help = "Input file (FASTA format).")
    parser.add_argument("-o", "--output", required = False, type = str, help = "Output file directory. Default is stdout")
    parser.add_argument("-x", "--xtract_A_domains", required = False, default = 1, type = int, help = "[0|1 extract AMP-binding domains?].")
    return parser

class PredictRFResult():
    def __init__(self, seqname, prediction, probability):
    	self.name = seqname
    	self.prediction = prediction
    	self.probability = probability

    def write_result(self, out_file):
        out_file.write("%s\t%s\t%.3f\n" % \
                       (self.name, self.prediction, self.probability))

def get_feature_matrix(fasta_dir: str, properties: Dict[str, PhysicochemicalProps]) -> Tuple[List[List[float]], List[str]]:
    fasta_dict = fasta.read_fasta(fasta_dir)

    features = []
    response = []
    seq_IDs = []
    
    for ID in fasta_dict:
        sequence = fasta_dict[ID]
        specificity = extract_specificity(ID)
        property_vector = extract_features(properties, sequence)
        features.append(property_vector)
        response.append(specificity)
        seq_IDs.append(ID)

    feature_list = []
    feature_list.append(seq_IDs)
    feature_list.append(features)
    return(feature_list)

def extract_features(property_dict: Dict[str, PhysicochemicalProps], sequence: str) -> List[float]:
    """Return list of float, with each float a physicochemical property

    Input:
    property_dict: dict of {aa: PhysicochemicalProps, ->}, with aa str
    sequence: str, amino acid sequence of length 34
    """

    property_vector = []

    for aa in sequence:
        for aa_property in property_dict[aa].get_all():
            property_vector.append(aa_property)

    return property_vector

def make_prediction(fasta_dir):
    filename = "%s/data/models/finalized_rf_model.sav" % parent_folder
    rf_mod = pickle.load(open(filename, 'rb'))

    properties_4 = parse_4_properties(PROPERTIES_4)
    feature_list = get_feature_matrix(fasta_dir, properties_4)

    seqname = feature_list[0]
    features = feature_list[1]
    prediction = rf_mod.predict(features)
    probability = numpy.amax(rf_mod.predict_proba(features), axis = 1)

    results = []
    for i,res in enumerate(tqdm.tqdm(prediction)):
        result = PredictRFResult(seqname[i], prediction[i], probability[i])
        results.append(result)

    return(results)

if __name__ == "__main__":
    parser = define_arguments()
    args = parser.parse_args()

    fasta_dir = args.input

    # Optional: extract AMP-binding hits from a FASTA file
    """
    Arguments:
        fastafile: FASTA file 
        hmmfile: AMP-binding.hmm
        outfile: path for output FASTA file
        size_threshold: threshold for the minimum size so you don't get truncated domains (recommendation is 50 - 100)
    """

    if args.xtract_A_domains == 1:
        print("##### Extracting AMP-binding domains... #####")
        out_file = "%s/data/AMP_binding_extracted.fasta" % parent_folder
        xtract_doms(fasta_dir, HMM_FILE, out_file, 50, True)

    # Create antiSMASH-like domain objects from AMP-binding hits
    domain_list = []
    create_domain_fa = fasta.read_fasta(fasta_dir)
    for i, domain in enumerate(create_domain_fa):
        domain_list.append(AntismashDomain(FeatureLocation(1, 1, 1), tool="test")) # arbitrary feature location for testing
        domain_list[i].domain_id = list(create_domain_fa.keys())[i]
        domain_list[i].translation = list(create_domain_fa.values())[i]
    
    # Extract the active site residues
    res = []
    new_path = "%s/data/34_aa_xtracted.fasta" % parent_folder
    nms = []
    sqs = []
    print("##### Extracting active site residues... #####")
    for i, domain in enumerate(tqdm.tqdm(domain_list)):
        res.append(get_34_aa_signature(domain))
        nms.append(domain.domain_id)
        sqs.append(res[i])

    fasta.write_fasta(nms, sqs, new_path)

    print("##### \n Making predictions... #####")
    results = make_prediction(new_path)
    print("Query_name\tPrediction\tProbability")
    for x in range(len(results)):
        print("%s\t%s\t%.3f\n" % \
        (results[x].name, results[x].prediction, results[x].probability))
    print("########################")
    
    if args.output is not None:
        with open(args.output, 'w') as output_file:
             output_file.write("Query_name\tPrediction\tProbability")
             output_file.write("\n")
             for result in results:
                result.write_result(output_file)
        output_file.close()
    


    

