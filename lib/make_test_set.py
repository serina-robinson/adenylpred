#!/usr/bin/env python

"""
Functions to split up a data set into test and training data
"""

import random

def select_with_probability(probability: float) -> bool:
    """Return True if probability <= random roll, False if it doesn't

    Input:
    probability: float, value between 0 and 1

    Output:
    bool: True if data point is assigned to test set, False if it is assigned
        to the training set
    """

    if random.random <= probability:
        return True
    else:
        return False

def decide_train_size(train_fraction: float, size: int) -> int:
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

def split_randomly(choices, class_1_size):
    random.shuffle(choices)
    class_1 = choices[:class_1_size]
    class_2 = choices[class_1_size:]

    return class_1, class_2
    

def inventorise_data_per_class(features, responses):
    response_dict = {}
    for i, feature_set in enumerate(features):
        response = responses[i]
        if not response in response_dict:
            response_dict[response] = []
        response_dict[response].append(feature_set)

    return response_dict

def inventorise_data_per_class_seqs(features, responses, seq_IDs):
    response_dict = {}
    for i, feature_set in enumerate(features):
        response = responses[i]
        if not response in response_dict:
            response_dict[response] = []
        response_dict[response].append((feature_set, seq_IDs[i]))

    return response_dict

def write_test_training(test_seqs, training_seqs, out_file):
    for seq in test_seqs:
        out_file.write("%s\t%s\n" % (seq, 'test'))
    for seq in training_seqs:
        out_file.write("%s\t%s\n" % (seq, 'training'))

def stratify_data_test(response_dict, fraction):
    clean_test_features = []
    clean_test_response = []
    clean_test_seqs = []

    training_features = []
    training_response = []
    training_seqs = []
    

    for response in response_dict:
        choices = response_dict[response]
 #       test, training = split_randomly(choices, test_nr)
        train_nr = decide_train_size(fraction, len(choices))
        training, test = split_randomly(choices, train_nr)

        for feature in test:
            clean_test_features.append(feature[0])
            clean_test_seqs.append(feature[1])
            clean_test_response.append(response)

        for feature in training:
            training_features.append(feature[0])
            training_seqs.append(feature[1])
            training_response.append(response)

    return clean_test_features, clean_test_response, clean_test_seqs, \
           training_features, training_response, training_seqs
