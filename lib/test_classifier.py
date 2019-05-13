#!/usr/bin/env python

"""
Script for testing a classifier
"""

from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
import os
from sys import path

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

from joblib import load
from typing import Any, Dict, List, Tuple, Optional
from lib.make_test_set import inventorise_data_per_class

def determine_accuracy(predictions, responses):
    matches = 0
    
    for i, prediction in enumerate(predictions):

        if prediction == responses[i]:
            matches += 1

    return float(matches)/len(responses)

def test(classifier, features, responses):
    predictions = classifier.predict(features)
    accuracy = determine_accuracy(predictions, responses)

    return accuracy

def test_per_class(classifier, features, responses):
    response_dict = inventorise_data_per_class(features, responses)
    accuracy_dict = {}

    for response in response_dict:
        feature_list = response_dict[response]
        response_list = [response] * len(feature_list)
        print(response, len(feature_list))
        accuracy = test(classifier, feature_list, response_list)
        accuracy_dict[response] = accuracy

    return accuracy_dict

def print_accuracy_dict(accuracy_dict):
    aas = list(accuracy_dict.keys())
    aas.sort()

    for aa in aas:
        print("%s: %.3f" % (aa, accuracy_dict[aa]))
