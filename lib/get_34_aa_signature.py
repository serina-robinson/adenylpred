#!/usr/bin/env python
from sys import argv, path
from statistics import stdev
from typing import Any, Dict, List, Tuple, Optional, List, Set

import random
import argparse
import os
import copy
import pickle
import numpy
import tqdm
import sys
import warnings

current_folder = os.path.dirname(os.path.abspath(__file__))
parent_folder = os.path.dirname(current_folder)

from sklearn.ensemble import RandomForestClassifier
from joblib import dump, load
from collections import defaultdict
from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import path, subprocessing
from antismash.modules.nrps_pks.nrps_predictor import build_position_list, read_positions, extract, verify_good_sequence

from lib.get_seq_properties import *
from lib.make_test_set import *
from lib.test_classifier import *
from lib.hmmsearch_domains import *

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

def build_position_list(positions: List[int], reference_seq: str) -> List[int]:
    """ Adjusts a list of positions to account for gaps in the reference sequence

        Arguments:
            positions: a list of ints that represent positions of interest in
                       the reference sequence
            reference_seq: the (aligned) reference sequence

        Returns:
            a new list of positions, each >= the original position
    """
    poslist = []
    position = 0
    for i, ref in enumerate(reference_seq):
        if ref != "-":
            if position in positions:
                poslist.append(i)
            position += 1
    return poslist

def get_34_aa_signature(domain: AntismashDomain) -> str:
    """ Extract 10 / 34 AA NRPS signatures from A domains """
    assert " " not in domain.get_name()
    # assert verify_good_sequence(domain.translation) # TODO uncomment this
    # Run muscle and collect sequence positions from file
    alignments = subprocessing.run_muscle_single(domain.get_name(), domain.translation, ADOMAINS_FILENAME)

    domain_alignment = alignments[domain.get_name()]
    reference_alignment = alignments[REF_SEQUENCE]

    positions = read_positions(APOSITION_FILENAME, START_POSITION)
    # Count residues in ref sequence and put positions in list
    poslist = build_position_list(positions, reference_alignment)

    # Extract positions from query sequence
    query_sig_seq = extract(domain_alignment, poslist)
    # Add fixed lysine 517
    query_sig_seq += "K"

    # repeat with 34 AA codes
    angpositions = read_positions(A34_POSITIONS_FILENAME, START_POSITION)
    poslist = build_position_list(angpositions, reference_alignment)

    return extract(domain_alignment, poslist)