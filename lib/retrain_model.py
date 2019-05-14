#!/usr/bin/env python

"""
Script to train, cross-validate and test a random forest classifier
"""

from sys import argv, path
from statistics import stdev
from typing import Any, Dict, List, Tuple, Optional

import random
import argparse
import os
import copy
import pickle

from sklearn.ensemble import RandomForestClassifier
from joblib import dump, load

from train_RF import train_RF, store_model


current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

path.append("%s/Data_preparation/" % parent_folder)
path.append("%s/Testing/" % parent_folder)

from get_seq_properties import *
from make_test_set import *
from test_classifier import *

PROPERTIES_15 = "%s/Data_preparation/15_aa_properties.txt" % parent_folder
PROPERTIES_4 = "%s/Data_preparation/aa_properties.txt" % parent_folder

def define_arguments() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description = "Run random forest.")
    
    parser.add_argument("--fasta", required = True, type = str,
                        help = "Fasta directory.")
    parser.add_argument("--trees", required = True, nargs = '+', type = int,
                        help = "Forest sizes to try.")
    parser.add_argument("--depth", required = True, nargs = '+',
                        help = "Maximum tree depths to try.")
    parser.add_argument("--features", required = True, nargs = '+',
                        help = "Feature numbers to try.")
    parser.add_argument("--iterations", required = True, type = int,
                        help = "Number of times to run cross-validation.")
    parser.add_argument("--train_size", required = True, nargs = '+',
                        type = float, help = "Test sizes to try.")
    parser.add_argument("--groups", type = string_to_bool, default = False,
                        help = "Cross-validate group specificities only.")
    parser.add_argument("--cutoffs", nargs = '+', type = int,
                        help = "Minimum support of data points for a class.")
    parser.add_argument("--clean", required = True, type = string_to_bool,
                        default = True, help = "True if you want to start from\
scratch, False if you want to include previous results")

    return parser

def parse_depths(depths):
    depth_list = []
    for depth in depths:
        if depth == 'None':
            depth_list.append(None)
        else:
            depth_list.append(int(depth))

    return depth_list

def parse_max_features(max_features):
    feature_list = []

    for max_feature in max_features:
        if max_feature == 'log2' or max_feature == 'auto':
            feature_list.append(max_feature)
        else:
            feature_list.append(float(max_feature))

    return feature_list

def run_rf(training_features, training_responses, trees, depths, max_features,
             iterations, train_size, cutoff):

    forest = train_RF(training_features, training_responses,
                                      forest_size, max_feature, depth)
    filename = '../serina_scratch/finalized_rf_model.sav' # save finalized model 
    pickle.dump(forest, open(filename, 'wb'))

def crossval(training_features, training_responses, trees, depths, max_features,
             iterations, train_size, cutoff):

    results = []

    
    for forest_size in trees:
        for depth in depths:
            for max_feature in max_features:
                result = CrossvalResult(forest_size, depth, max_feature,
                                        iterations, train_size, cutoff)
                for i in range(iterations):
                    
                    forest = train_RF(training_features, training_responses,
                                      forest_size, max_feature, depth)
                    filename = '../serina_scratch/finalized_rf_model.sav' # save finalized model 
                    pickle.dump(forest, open(filename, 'wb'))

                    accuracy = forest.oob_score_
                    result.add_accuracy(accuracy)

                result.calc_properties()
                results.append(result)

    results.sort(key = lambda result: result.average_accuracy, reverse = True)

    return results
                    
def write_results(results, out_dir):
    out_file = open(out_dir, 'a+')
    out_file.write("\tTraining size\tCutoff\tTrees\tDepth\tMax_features\tIterations\tOOB_score\tStd_dev")
    out_file.write("\n")
    
    for result in results:
        result.write_result(out_file)

    out_file.close()

def parse_results(results_dir):
    results_file = open(results_dir, 'r')
    results = []

    results_file.readline()
    
    for line in results_file:
        line = line.strip()
        features = line.split('\t')

        train_size = float(features[0])
        cutoff = int(features[1])
        trees = int(features[2])
        depth = features[3]
        max_features = features[4]
        iterations = int(features[5])
        accuracy = float(features[6])
        stddev = float(features[7])

        result = CrossvalResult(trees, depth, max_features, iterations,
                                train_size, cutoff)
        result.average_accuracy = accuracy
        result.stddev = stddev

        results.append(result)

    results_file.close()

    return results
    
class CrossvalResult():
    def __init__(self, trees, depth, max_features, iterations, train_size,
                 cutoff):
        self.trees = trees
        self.depth = depth
        self.max_features = max_features
        self.iterations = iterations
        self.train_size = train_size
        self.cutoff = cutoff


        self.depth_ID = self.get_depth_ID()
        self.features_ID = self.get_features_ID()
        
        self.ID = "%d_%s_%s_%d_%.2f_%d" % (self.trees, self.depth_ID,
                                           self.features_ID, self.iterations,
                                           self.train_size, self.cutoff)

        self.accuracies = []

    def __hash__(self):
        return self.ID

    def write_result(self, out_file):
        out_file.write("%.2f\t%d\t%d\t%s\t%s\t%d\t%.3f\t%.3f\n" % \
                       (self.train_size, self.cutoff, self.trees,
                        self.depth_ID, self.features_ID, self.iterations,
                        self.average_accuracy, self.stddev))

    def get_features_ID(self):
        features_ID = str(self.max_features)
        return features_ID

    def set_suffix(self, suffix):
        self.ID += "_%s" % suffix

    def get_depth_ID(self):
        if not self.depth:
            depth_ID = 'max'
        else:
            depth_ID = str(self.depth)

        return depth_ID

    def add_accuracy(self, accuracy):
        self.accuracies.append(accuracy)

    def set_stdev(self):
        try:
            self.stddev = stdev(self.accuracies)
        except Exception:
            self.stddev = 0.0

    def set_average_accuracy(self):
        assert self.iterations == len(self.accuracies)
        self.average_accuracy = sum(self.accuracies)/self.iterations

    def calc_properties(self):
        self.set_stdev()
        self.set_average_accuracy()

    def get_properties(self):
        return self.ID, self.average_accuracy, self.stddev



def main(fasta, trees, depths, max_features, iterations, train_sizes, groups,
         cutoffs, clean):
    depths = parse_depths(depths)
    max_features = parse_max_features(max_features)
    
    properties_15 = parse_15_properties(PROPERTIES_15)
    properties_4 = parse_4_properties(PROPERTIES_4)

    if groups:

        dir_15 = 'results_15_groups.txt'
        dir_4 = 'results_4_groups.txt'

        if clean:
            
            all_results_15 = []
            all_results_4 = []
            
        else:
            assert os.path.isfile(dir_15)
            
            all_results_15 = parse_results(dir_15)
            all_results_4 = parse_results(dir_4)

        if os.path.isfile(dir_15):
            os.remove(dir_15)

        if os.path.isfile(dir_4):
            os.remove(dir_4)
            
        
        for train_size in train_sizes:
            data_15 = DataSet(fasta, properties_15, True)
            data_4 = DataSet(fasta, properties_4, True)

            data_15.stratify_data(train_size)
            data_4.stratify_data(train_size)

            data_15.write_test_training(train_size, "15_groups")
            data_4.write_test_training(train_size, "4_groups")

            results_15 = crossval(data_15.training_features,
                                  data_15.training_response,
                                  trees, depths, max_features,
                                  iterations, train_size, 0)

            results_4 = crossval(data_4.training_features,
                                 data_4.training_response,
                                 trees, depths, max_features,
                                 iterations, train_size, 0)

            all_results_15 += results_15
            all_results_4 += results_4

        all_results_15.sort(key = lambda result: result.average_accuracy, reverse = True)
        all_results_4.sort(key = lambda result: result.average_accuracy, reverse = True)

        
        
            
            

        write_results(all_results_15, dir_15)
        write_results(all_results_4, dir_4)
    else:

        dir_15 = 'results_15.txt'
        dir_4 = 'results_4.txt'

        if clean:
            
            all_results_15 = []
            all_results_4 = []
            
        else:
            assert os.path.isfile(dir_15)
            
            all_results_15 = parse_results(dir_15)
            all_results_4 = parse_results(dir_4)

        if os.path.isfile(dir_15):
            os.remove(dir_15)

        if os.path.isfile(dir_4):
            os.remove(dir_4)

        for train_size in train_sizes:
            for cutoff in cutoffs:
                data_15 = DataSet(fasta, properties_15, False, cutoff)
                data_4 = DataSet(fasta, properties_4, False, cutoff)

                data_15.stratify_data(train_size)
                data_4.stratify_data(train_size)

                data_15.write_test_training(train_size, "15_cutoff%d" % cutoff)
                data_4.write_test_training(train_size, "4_cutoff%d" % cutoff)

                results_15 = crossval(data_15.training_features,
                                      data_15.training_response,
                                      trees, depths, max_features,
                                      iterations, train_size, cutoff)

                results_4 = crossval(data_4.training_features,
                                     data_4.training_response,
                                     trees, depths, max_features,
                                     iterations, train_size, cutoff)

                all_results_15 += results_15
                all_results_4 += results_4


        all_results_15.sort(key = lambda result: result.average_accuracy, reverse = True)
        all_results_4.sort(key = lambda result: result.average_accuracy, reverse = True)

        write_results(all_results_15, dir_15)
        write_results(all_results_4, dir_4)

def string_to_bool(value):
    return value.lower() == 'true'

if __name__ == "__main__":
    parser = define_arguments()
    args = parser.parse_args()

    fasta = args.fasta
    trees = args.trees
    depths = args.depth
    max_features = args.features
    iterations = args.iterations
    train_sizes = args.train_size
    groups = args.groups
    cutoffs = args.cutoffs
    clean = args.clean

    print(type(groups))
    print(groups)

    main(fasta, trees, depths, max_features, iterations, train_sizes, groups,
         cutoffs, clean)

    




        
