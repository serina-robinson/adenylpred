#!/usr/bin/env python
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
AdenylPred: script to predict substrate given a protein sequence
"""

from sys import argv, path
from typing import Any, Dict, List, Tuple, Optional, List, Set
from sklearn.ensemble import RandomForestClassifier
from collections import defaultdict

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
import tempfile

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

from lib.get_seq_properties import *
from lib.seq_manipulations import *
from lib.extract_34_hmmpfam import *

PROPERTIES_15 = "%s/data/15_aa_properties.txt" % parent_folder
PROPERTIES_4 = "%s/data/aa_properties.txt" % parent_folder
HMM_FILE = "%s/data/AMP-binding_hmm2.hmm" % parent_folder

def define_arguments() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description = "Prediction tool for adenylate-forming enzyme substrate specificity")
    parser.add_argument("-i", "--input", required = True, type = str, help = "Input file (FASTA or GenBank format).")
    parser.add_argument("-o", "--output", required = True, type = str, help = "Output file directory. Default is stdout")
    parser.add_argument("-v", "--verbose", required = False, action = "store_true", help = "Verbose. Prints progress to stdout")
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

def make_prediction(fasta_dir, verbose):
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
    for i,res in enumerate(tqdm.tqdm(prediction, disable = not verbose)): # this seems weird
        result = PredictRFResult(seqname[i], prediction[i], probability[i])
        results.append(result)
    return results

def print_results(results, output_handle, verbose):
    if verbose:
        print("Query_name\tPrediction\tProbability_score")
        for result in results:
            result.write_result(sys.stdout)
    with open(output_handle, 'w') as output_file:
            output_file.write("Query_name\tPrediction\tProbability_score")
            output_file.write("\n")
            for result in results:
                result.write_result(output_file)
    output_file.close()

def clean_up():
    folder_path = "%s/data/" % parent_folder
    for file_name in os.listdir(folder_path):
        if file_name.endswith('_tmp.faa'):
            os.remove(folder_path + file_name)

def main(args):
    """Executes the main functions of AdenylPred: file format conversion, extraction of active site residues, and prediction"""

    verbose = args.verbose

    # Checks if GenBank file and converts to FASTA
    if args.genbank_input:
        try:
            out_file = '%s/data/gbk_to_fasta_tmp.faa' % parent_folder
            fasta_dir = gbk_to_faa(args.input, out_file)
        except:
           print('Error: please check your file is a valid FASTA or GenBank file. \n If your input is a nucleotide sequence, please check the -n option is set to 1')
           return

    # Checks if nucleotide sequence and converts to amino acid FASTA
    if args.nucleotide: # elif
        try:
            out_file = '%s/data/nuc_to_aa_tmp.faa' % parent_folder
            fasta_dir = nuc_to_aa(args.input, out_file)
        except:
            print('Error: please check your file is a valid FASTA or GenBank file. \n If your input is a nucleotide sequence, please check the -n option is set to 1')
            return

    # If neither nucleotide or GenBank then use input directly
    if not args.genbank_input and not args.nucleotide:  #check this # else
            fasta_dir = args.input

    # Align and extract 34 residues from alignment
    extract_file = '%s/data/34_aa_xtracted_tmp.faa' % parent_folder
    try:
        extract_34_aa(HMM_FILE, fasta_dir, extract_file)
    except:
        print('Could not extract active site residues. Please check sequence length(s) of input file.')
        return

    # Make predictions using 34 active site residues as features 
    if verbose:
        print('##### \n Making predictions... #####')
    try:
        results = make_prediction(extract_file, verbose = verbose)
    except:
         print('Error: could not make a prediction. Please check your input file. \n If your input is a nucleotide sequence, please check the -n option is set to 1')
         return

    # Clean up intermediate files
    clean_up()
    
    # Write predictions to file (and stdout if verbose)
    print_results(results, output_handle = args.output, verbose = verbose) 

if __name__ == '__main__':

    # Parse arguments
    parser = define_arguments()
    args = parser.parse_args()

    # Check if file or dir here 
    #for:
    main(args)
        