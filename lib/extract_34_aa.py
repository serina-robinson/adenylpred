""" Provides a collection of functions and classes to run the external Java program
    NRPSPredictor2 and interpret the results
"""

from collections import defaultdict
import os
import sys
from typing import Any, Dict, List, Set

from helperlibs.wrappers.io import TemporaryDirectory

sys.path.append(os.path.abspath('/Users/robi0916/Documents/Wageningen_UR/github/as5'))
from antismash.common import fasta, path, subprocessing
from antismash.common.secmet import AntismashDomain
from antismash.config import ConfigType
from antismash.common.secmet import AntismashDomain, FeatureLocation
from antismash.modules.nrps_pks.nrps_predictor import get_34_aa_signature

REF_SEQUENCE = "P0C062_A1"
A34_POSITIONS_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "A34positions.txt")
APOSITION_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "Apositions.txt")
KNOWN_CODES = path.get_full_path(__file__, "knowncodes.fasta")
ILLEGAL_CHARS = "!@#$%^&*(){}:\"<>?/.,';][`~1234567890*-+-=_\\|"
ADOMAINS_FILENAME = path.get_full_path(__file__, "external", "NRPSPredictor2", "A_domains_muscle.fasta")
START_POSITION = 66

# create_domain_fa = fasta.read_fasta('../flat/669_training_set_to_extract_34aa_20192803.faa')
domain_list = []
for i, domain in enumerate(create_domain_fa):
    domain_list.append(AntismashDomain(FeatureLocation(1, 1, 1), tool="test")) # arbitrary feature location for testing
    domain_list[i].domain_id = list(create_domain_fa.keys())[i]
    domain_list[i].translation = list(create_domain_fa.values())[i]

res = []
for i, domain in enumerate(domain_list):
	res.append(get_34_aa_signature(domain))
	print('>' + domain.domain_id)
	print(res[i]) 
