#!/usr/bin/env python

"""
Script to extract sequence features from a fasta file

Make sure to have the antismash5 environment activated
"""

import sys
import os

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

from statistics import stdev, mean
from typing import Any, Dict, List, Tuple, Optional
from sys import argv
import copy
import random
from pprint import pprint

from lib import fasta

class DataSet():
    chuck_dict = {"ser": 'pac',
                  "oh-asn": 'pac',
                  "o-iv": 'pac',
                  "noh-val": 'pac',
                  "noh-ala": 'pac',
                  "glu": 'pac',
                  "gln": 'pac',
                  "arg": 'pac',
                  "arg|leu": 'pac',
                  "dab": 'pac',
                  "oh-orn": 'pac',
                  "hser": 'pac',
                  "cysa": 'pac',
                  "asp|bme-asp": 'pac',
                  "asp": 'pac',
                  "asn|gln": 'pac',
                  "asn": 'pac',
                  "arg|lys": 'pac',
                  "oh-asp": 'pac',
                  "orn": 'pac',
                  "oh-orn|orn": 'pac',

                  "cys": 'cys',

                  "tyr": 'bulky',
                  "trp": 'bulky',
                  "phe|trp": 'bulky',
                  "phe": 'bulky',
                  "htyr": 'bulky',
                  "his": 'bulky',
                  "gln|tyr": 'bulky',
                  "lys": 'bulky',
                  "glyco-cya": 'bulky',
                  "nme-asp": 'bulky',
                  "fo-oh-orn": 'bulky',
                  "dhh-tyr": 'bulky',
                  "asp|nme-asp": 'bulky',
                  "aad": 'bulky',
                  "boh-tyr|tyr": 'bulky',
                  "boh-tyr": 'bulky',
                  "bme-phe": 'bulky',
                  "bhceala": 'bulky',
                  "6cl-trp|trp": 'bulky',
                  "3oh-gln": 'bulky',
                  "3me-tyr|ddh-arg": 'bulky',
                  "3me-glu|glu": 'bulky',
                  "3me-asp": 'bulky',
                  "33p-ala": 'bulky',
                  "met": 'bulky',
                  "homm-tyr": 'bulky',
                  "3me-tyr": 'bulky',
                  "kyn": 'bulky',
                  #fix in dataset
                  "boh|tyr": 'bulky',
                  "phe|tyr": 'bulky',
                  "phgly": 'bulky',
                  "hpg": 'bulky',
                  "hcn": 'bulky',
                  "dhpg": 'bulky',
                  "cl2-hpg|hpg": 'bulky',
                  #fix in dataset
                  "2cl-hpg|hpg": 'bulky',

                  "val": 'smh',
                  "norcma": 'smh',
                  "nme-gly": 'smh',
                  "ile|val": 'smh',
                  "ile": 'smh',
                  "goh-val": 'smh',
                  "ala|val": 'smh',
                  "aile|ile": 'smh',
                  "abu": 'smh',
                  "3oh-leu": 'smh',
                  "hiv": 'smh',
                  "valol": 'smh',
                  "ncp-ala": 'smh',
                  "leu": 'smh',
                  "akic": 'smh',
                  "ala|leu|val": 'smh',
                  "4oh-5oh-orn": 'smh',
                  "3oh-3me-pro": 'smh',
                  #fix in dataset
                  "iso": 'smh',

                  "thr": 'lwat',
                  "piz": 'lwat',
                  "me-dpr": 'lwat',
                  "di-oh-bz": 'lwat',
                  "di-oh-bz|sal": 'lwat',
                  "dhdab": 'lwat',
                  "dhabu|dhtyr|thr": 'lwat',
                  "dhabu": 'lwat',
                  #check in source
                  "dh-thr": 'lwat',
                  "bmt": 'lwat',
                  "athr|thr": 'lwat',
                  "athr": 'lwat',
                  "alaol": 'lwat',
                  "ala|gly": 'lwat',
                  "ala": 'lwat',
                  "akiv": 'lwat',
                  "hse": 'lwat',

                  "pro": 'ringy',
                  "pip": 'ringy',
                  "leu|val": 'ringy',
                  "gly|nme-gly": 'ringy',
                  "gac": 'ringy',
                  "ac-phe": 'ringy',
                  "abu|ival": 'ringy',
                  "6cl-hic": 'ringy',
                  "4pe-pro": 'ringy',
                  "4p-pro": 'ringy',
                  "4oh-pro": 'ringy',
                  "4me-pro|pro": 'ringy',
                  "3me-pro": 'ringy',
                  "4me-pro": 'ringy',

                  "gly": 'tiny',
                  "gca": 'tiny',
                  "d-ala": 'tiny',

                  "d-lea": 'reject',
                  "bala": 'reject',
                  "blys": 'reject',
                  "ahnpb": 'reject',
                  "btrp": 'reject',
                  "athr|ser": 'reject',
                  "aeo": 'reject',
                  "ame-ser": 'reject',
                  "ser|thr": 'reject',
                  "pac": 'reject',
                  "ddh-arg|3me-tyr": 'reject',
                  "me-bala": 'reject'}
    
    def __init__(self, fasta, properties, groups = False, cutoff = 9):
        if groups:
            cutoff = 0
            
        self.features, self.responses, self.sequences = get_feature_matrix(\
            fasta, properties, cutoff)

        if groups:
            self.change_responses()

        self.group_data()
        print(self.classes.keys())

    def change_responses(self):
        new_responses = []
        new_features = []
        for i, response in enumerate(self.responses):
            if not self.chuck_dict[response] == 'reject':
                new_responses.append(self.chuck_dict[response])
                new_features.append(self.features[i])

        self.features = new_features
        self.responses = new_responses

    def group_data(self):
        response_dict = {}
        for i, feature_set in enumerate(self.features):
            response = self.responses[i]
            if not response in response_dict:
                response_dict[response] = []
            response_dict[response].append((feature_set, self.sequences[i]))

        self.classes = response_dict

    def decide_train_size(self, train_fraction: float, size: int) -> int:
        """Return number of samples to be drawn for training set

        Input:
        train_fraction: float, number between 0 and 1, indicates percentage of
            data points which will be included in training data

        size: int, number of data points in training data
        """
        if (size * train_fraction) % 1 == 0:
            train_nr = int(size * train_fraction)

        else:
            train_nr = int(size * train_fraction) + 1

        return train_nr

    def split_randomly(self, choices, class_1_size):
        random.shuffle(choices)
        class_1 = choices[:class_1_size]
        class_2 = choices[class_1_size:]

        return class_1, class_2

    def stratify_data(self, fraction):
        self.test_features = []
        self.test_response = []
        self.test_seqs = []

        self.training_features = []
        self.training_response = []
        self.training_seqs = []
        

        for response in self.classes:
            choices = copy.copy(self.classes[response])
            train_nr = self.decide_train_size(fraction, len(choices))
            training, test = self.split_randomly(choices, train_nr)

            for feature in test:
                self.test_features.append(feature[0])
                self.test_seqs.append(feature[1])
                self.test_response.append(response)

            for feature in training:
                self.training_features.append(feature[0])
                self.training_seqs.append(feature[1])
                self.training_response.append(response)

    def write_test_training(self, fraction, suffix):
        out_dir = "training_and_test_%.2f_%s.txt" % (fraction, suffix)
        out_file = open(out_dir, 'w')
        
        for seq in self.test_seqs:
            out_file.write("%s\t%s\n" % (seq, 'test'))
            
        for seq in self.training_seqs:
            out_file.write("%s\t%s\n" % (seq, 'training'))

        out_file.close()

class PhysicochemicalProps():
    def __init__(self, polarity: float, structure: float, size: float, codon: float, charge: float):
        self.polarity = polarity
        self.structure = structure
        self.size = size
        self.codon = codon
        self.charge = charge

    def get_all(self):
        return [self.polarity, self.structure, self.size, self.charge]

    def get_all_2(self):
        return [self.polarity, self.structure, self.size, self.codon, self.charge]


class PhysicochemicalPropsBig():
    def __init__(self, prop_1, prop_2, prop_3, prop_4, prop_5, prop_6, prop_7,
                 prop_8, prop_9, prop_10, prop_11, prop_12, prop_13, prop_14,
                 prop_15):
        self.prop_1 = prop_1
        self.prop_2 = prop_2
        self.prop_3 = prop_3
        self.prop_4 = prop_4
        self.prop_5 = prop_5
        self.prop_6 = prop_6
        self.prop_7 = prop_7
        self.prop_8 = prop_8
        self.prop_9 = prop_9
        self.prop_10 = prop_10
        self.prop_11 = prop_11
        self.prop_12 = prop_12
        self.prop_13 = prop_13
        self.prop_14 = prop_14
        self.prop_15 = prop_15

    def get_all(self):
        return [self.prop_1, self.prop_2, self.prop_3, self.prop_4,
                self.prop_5, self.prop_6, self.prop_7, self.prop_8,
                self.prop_9, self.prop_10, self.prop_11, self.prop_12,
                self.prop_13, self.prop_14, self.prop_15]
    
def parse_15_properties(property_dir: str) -> List[float]:
    property_file = open(property_dir, 'r')
    property_file.readline()


    matrix = []
    aa_list = []

    for line in property_file:
        line = line.strip()
        aa, abbr, prop_1, prop_2, prop_3, prop_4, prop_5, prop_6, prop_7, \
            prop_8, prop_9, prop_10, prop_11, prop_12, prop_13, prop_14, \
            prop_15 = line.split()

        aa = abbr.strip()
        aa_list.append(aa)
        prop_1 = float(prop_1.strip())
        prop_2 = float(prop_2.strip())
        prop_3 = float(prop_3.strip())
        prop_4 = float(prop_4.strip())
        prop_5 = float(prop_5.strip())
        prop_6 = float(prop_6.strip())
        prop_7 = float(prop_7.strip())
        prop_8 = float(prop_8.strip())
        prop_9 = float(prop_9.strip())
        prop_10 = float(prop_10.strip())
        prop_11 = float(prop_11.strip())
        prop_12 = float(prop_12.strip())
        prop_13 = float(prop_13.strip())
        prop_14 = float(prop_14.strip())
        prop_15 = float(prop_15.strip())

        matrix.append([prop_1, prop_2, prop_3, prop_4, prop_5, prop_6,
                       prop_7, prop_8, prop_9, prop_10, prop_11, prop_12,
                       prop_13, prop_14, prop_15])

    property_file.close()

    meanstddev = get_mean_stddev(matrix)
    normalise(matrix, meanstddev)

    pprint(matrix)

    property_dict = {}

    for i, aa in enumerate(aa_list):
        prop_1, prop_2, prop_3, prop_4, prop_5, prop_6, prop_7, prop_8, \
                prop_9, prop_10, prop_11, prop_12, prop_13, prop_14, prop_15 \
                = matrix[i]
        properties = PhysicochemicalPropsBig(prop_1, prop_2, prop_3, prop_4,
                                             prop_5, prop_6, prop_7, prop_8,
                                             prop_9, prop_10, prop_11, prop_12,
                                             prop_13, prop_14, prop_15)
        property_dict[aa] = properties

    property_dict['-'] = PhysicochemicalPropsBig(0.0, 0.0, 0.0, 0.0, 0.0,
                                                 0.0, 0.0, 0.0, 0.0, 0.0,
                                                 0.0, 0.0, 0.0, 0.0, 0.0)


    return property_dict

def get_mean_stddev(matrix):
    meanstddev = []

    for i in range(len(matrix[0])):
        col_vals = [row[i] for row in matrix]
        stddev = stdev(col_vals)
        mean_val = mean(col_vals)
        meanstddev.append((mean_val, stddev))

    return meanstddev

def normalise(matrix, meanstddev):
    for row in matrix:
        for i in range(len(row)):
            mean_val = meanstddev[i][0]
            stddev = meanstddev[i][1]
            row[i] = (row[i] - mean_val) / (stddev)


def make_chuck_dict():

    chuck_dict = {"ser": 'pac',
                  "oh-asn": 'pac',
                  "o-iv": 'pac',
                  "noh-val": 'pac',
                  "noh-ala": 'pac',
                  "glu": 'pac',
                  "gln": 'pac',
                  "arg": 'pac',
                  "arg|leu": 'pac',
                  "dab": 'pac',
                  "oh-orn": 'pac',
                  "hser": 'pac',
                  "cysa": 'pac',
                  "asp|bme-asp": 'pac',
                  "asp": 'pac',
                  "asn|gln": 'pac',
                  "asn": 'pac',
                  "arg|lys": 'pac',
                  "oh-asp": 'pac',
                  "orn": 'pac',
                  "oh-orn|orn": 'pac',

                  "cys": 'cys',

                  "tyr": 'bulky',
                  "trp": 'bulky',
                  "phe|trp": 'bulky',
                  "phe": 'bulky',
                  "htyr": 'bulky',
                  "his": 'bulky',
                  "gln|tyr": 'bulky',
                  "lys": 'bulky',
                  "glyco-cya": 'bulky',
                  "nme-asp": 'bulky',
                  "fo-oh-orn": 'bulky',
                  "dhh-tyr": 'bulky',
                  "asp|nme-asp": 'bulky',
                  "aad": 'bulky',
                  "boh-tyr|tyr": 'bulky',
                  "boh-tyr": 'bulky',
                  "bme-phe": 'bulky',
                  "bhceala": 'bulky',
                  "6cl-trp|trp": 'bulky',
                  "3oh-gln": 'bulky',
                  "3me-tyr|ddh-arg": 'bulky',
                  "3me-glu|glu": 'bulky',
                  "3me-asp": 'bulky',
                  "33p-ala": 'bulky',
                  "met": 'bulky',
                  "homm-tyr": 'bulky',
                  "3me-tyr": 'bulky',
                  "kyn": 'bulky',
                  #fix in dataset
                  "boh|tyr": 'bulky',
                  "phe|tyr": 'bulky',
                  "phgly": 'bulky',
                  "hpg": 'bulky',
                  "hcn": 'bulky',
                  "dhpg": 'bulky',
                  "cl2-hpg|hpg": 'bulky',
                  #fix in dataset
                  "2cl-hpg|hpg": 'bulky',

                  "val": 'smh',
                  "norcma": 'smh',
                  "nme-gly": 'smh',
                  "ile|val": 'smh',
                  "ile": 'smh',
                  "goh-val": 'smh',
                  "ala|val": 'smh',
                  "aile|ile": 'smh',
                  "abu": 'smh',
                  "3oh-leu": 'smh',
                  "hiv": 'smh',
                  "valol": 'smh',
                  "ncp-ala": 'smh',
                  "leu": 'smh',
                  "akic": 'smh',
                  "ala|leu|val": 'smh',
                  "4oh-5oh-orn": 'smh',
                  "3oh-3me-pro": 'smh',
                  #fix in dataset
                  "iso": 'smh',

                  "thr": 'lwat',
                  "piz": 'lwat',
                  "me-dpr": 'lwat',
                  "di-oh-bz": 'lwat',
                  "di-oh-bz|sal": 'lwat',
                  "dhdab": 'lwat',
                  "dhabu|dhtyr|thr": 'lwat',
                  "dhabu": 'lwat',
                  #check in source
                  "dh-thr": 'lwat',
                  "bmt": 'lwat',
                  "athr|thr": 'lwat',
                  "athr": 'lwat',
                  "alaol": 'lwat',
                  "ala|gly": 'lwat',
                  "ala": 'lwat',
                  "akiv": 'lwat',
                  "hse": 'lwat',

                  "pro": 'ringy',
                  "pip": 'ringy',
                  "leu|val": 'ringy',
                  "gly|nme-gly": 'ringy',
                  "gac": 'ringy',
                  "ac-phe": 'ringy',
                  "abu|ival": 'ringy',
                  "6cl-hic": 'ringy',
                  "4pe-pro": 'ringy',
                  "4p-pro": 'ringy',
                  "4oh-pro": 'ringy',
                  "4me-pro|pro": 'ringy',
                  "3me-pro": 'ringy',
                  "4me-pro": 'ringy',

                  "gly": 'tiny',
                  "gca": 'tiny',
                  "d-ala": 'tiny',

                  "d-lea": 'reject',
                  "bala": 'reject',
                  "blys": 'reject',
                  "ahnpb": 'reject',
                  "btrp": 'reject',
                  "athr|ser": 'reject',
                  "aeo": 'reject',
                  "ame-ser": 'reject',
                  "ser|thr": 'reject',
                  "pac": 'reject',
                  "ddh-arg|3me-tyr": 'reject',
                  "me-bala": 'reject'}

    return chuck_dict
     
            

def parse_4_properties(property_dir: str) -> Dict[str, PhysicochemicalProps]:
    """Return dict of {aa: PhysocochemicalProps, ->}

    Input:
    property_dir: str, directory of file containing data on scaled amino
        acid properties

    Output:
    property_dict: dict of {aa: PhysicochemicalProps, ->}, with aa str
    """

    property_file = open(property_dir, 'r')
    property_file.readline()
    
    property_dict = {}
    for line in property_file:
        line = line.strip()
        aa, polarity, structure, size, codon, charge = line.split()
        aa = aa.strip()
        polarity = float(polarity.strip())
        structure = float(structure.strip())
        size = float(size.strip())
        codon = float(codon.strip())
        charge = float(charge.strip())

        property_dict[aa] = PhysicochemicalProps(polarity, structure, size, codon, charge)
    property_dict['-'] = PhysicochemicalProps(0.0, 0.0, 0.0, 0.0, 0.0)

    property_file.close()
    return property_dict

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

def extract_specificity(fasta_ID: str) -> str:
    specificity = "unknown" #fasta_ID.split('_')[1]
    specificity = specificity.lower()

    if '|' in specificity:
        specificities = specificity.split('|')
        specificities.sort()
        specificity = '|'.join(specificities)

    return specificity

def count_specificities(fasta_dict: Dict[str, str]) -> Dict[str, int]:
    specificity_dict = {}
    for ID in fasta_dict:
        specificity = extract_specificity(ID)
        if not specificity in specificity_dict:
            specificity_dict[specificity] = 1
        else:
            specificity_dict[specificity] += 1

    return specificity_dict


def select_sequences(fasta_dict: Dict[str, str], cutoff: int = 9) -> Dict[str, bool]:
    
    specificity_dict = count_specificities(fasta_dict)
    print(cutoff)

    selected_sequences = {}
    excluded = 0
    omitted = []

    for ID in fasta_dict:
        specificity = extract_specificity(ID)
        support = specificity_dict[specificity]
        if support > cutoff:
            selected_sequences[ID] = True
        else:
            selected_sequences[ID] = False
            excluded += 1
            omitted.append(specificity)

    print("%d data points omitted due to lack of support." % excluded)
    print(set(omitted))

    return selected_sequences

def get_feature_matrix(fasta_dir: str, properties: Dict[str, PhysicochemicalProps], cutoff: int = 9) -> Tuple[List[List[float]], List[str]]:
    fasta_dict = fasta.read_fasta(fasta_dir)

    features = []
    response = []
    seq_IDs = []

    selected_sequences = select_sequences(fasta_dict, cutoff = cutoff)
    
    for ID in fasta_dict:
        if selected_sequences[ID]:
            sequence = fasta_dict[ID]
            specificity = "unknown" #extract_specificity(ID
            property_vector = extract_features(properties, sequence)
            features.append(property_vector)
            response.append(specificity)
            seq_IDs.append(ID)

    return features, response, seq_IDs

def write_property_matrix(property_dict, out_dir):
    aas = list(property_dict.keys())
    aas.sort()

    out_file = open(out_dir, 'w')
    for aa in aas:
        properties = property_dict[aa]
        prop_1 = properties.prop_1
        prop_2 = properties.prop_2
        prop_3 = properties.prop_3
        prop_4 = properties.prop_4
        prop_5 = properties.prop_5
        prop_6 = properties.prop_6
        prop_7 = properties.prop_7
        prop_8 = properties.prop_8
        prop_9 = properties.prop_9
        prop_10 = properties.prop_10
        prop_11 = properties.prop_11
        prop_12 = properties.prop_12
        prop_13 = properties.prop_13
        prop_14 = properties.prop_14
        prop_15 = properties.prop_15

        out_file.write("%s\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n" %\
                       (aa, prop_1, prop_2, prop_3, prop_4, prop_5, prop_6, prop_7, prop_8,
                        prop_9, prop_10, prop_11, prop_12, prop_13, prop_14, prop_15))

    out_file.close()

def write_to_file(property_dict, fasta_dir):
    fasta_dict = fasta.read_fasta(fasta_dir)
    out_dir = "%s_features.txt" % fasta_dir.split('.')[0]
    out_file = open(out_dir, 'w')
    
    for ID in fasta_dict:
        sequence = fasta_dict[ID]
        # specificity = extract_specificity(ID)
        specificity = "unknown"
        property_vector = extract_features(property_dict, sequence)
        
        out_file.write("%s\t%s" % (ID, specificity))
        for prop in property_vector:
            out_file.write("\t%.10f" % prop)
        out_file.write('\n')

    out_file.close()

def make_training_dict(test_training_dir):
    test_training_file = open(test_training_dir, 'r')
    training_dict = {}
    for line in test_training_file:
        line = line.strip()
        ID, decision = line.split()
        training_dict[ID] = decision

    test_training_file.close()

    return training_dict
        
        
def split_test_training_from_file(test_training_dir, properties, fasta_dir):
    fasta_dict = fasta.read_fasta(fasta_dir)
    training_dict = make_training_dict(test_training_dir)

    features_test = []
    response_test = []
    features_training = []
    response_training = []

    
    for ID in fasta_dict:
        if ID in training_dict:
            
            sequence = fasta_dict[ID]
            specificity = extract_specificity(ID)
            property_vector = extract_features(properties, sequence)

            if training_dict[ID] == 'training':
                features_training.append(property_vector)
                response_training.append(specificity)
            else:
                features_test.append(property_vector)
                response_test.append(specificity)

    return features_test, response_test, features_training, response_training
            
    

if __name__ == "__main__":
    #PHYSICOCHEMICAL = parse_4_properties("aa_properties.txt")
    PHYSICOCHEMICAL = parse_15_properties("15_aa_properties.txt")

    
 #   write_to_file(PHYSICOCHEMICAL, argv[1])
    #features, response = get_feature_matrix(argv[1], PHYSICOCHEMICAL)
   # print(features[0])
    #print(len(features[0]))
    #print(response[0])
    write_property_matrix(PHYSICOCHEMICAL, argv[1])
    
            
    
    
    
