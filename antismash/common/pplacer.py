# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" pplacer parsing and clade assignment. """

import logging
import os
from typing import Any, Dict, List, Optional, Tuple

from antismash.common import fasta, module_results, pfamdb, subprocessing
from antismash.common.secmet import Record, CDSFeature
from antismash.common.secmet.features.feature import FeatureLocation #PFAMDomain

from ete3 import Tree

class CladeAssignmentResults(module_results.ModuleResults):
    """ Results for clade assignments """

    def __init__(self, record_id: str, clade_assignment: str, tool: str, snn_score: float) -> None:
        super().__init__(record_id)
        self.clade_assignment = str(clade_assignment)
        self.tool = str(tool)
        self.snn_score = float(snn_score)


## Note:  Key of returned dictionary may eventually be changed to class Record
## Note2: this was tested only for TransAT KS domains without an SNN score. This
##        will need to be adapted to fit SANDPUMA
def parse_pplacer(leaf2clade: Dict[str, str], pplacer_tree: str, masscutoff: float, mode: str, calc_snn: bool) -> Dict[str, CladeAssignmentResults]:
    """ Parses pplacer trees (guppy sing) and returns clade assignments and SNN scores

        Arguments:
            leaf2clade: dictionary of leaves in tree to clade
            pplacer_tree: pplacer tree in newick format
            masscutoff: threshold of mass to assign a clade ranging from 0 (low) to 1 (high), default=0.9
            mode: clade assignment mode
            calc_snn: True=calculate SNN score, False=don't
        Returns:
            dictionary of query_name -> clade_assignment
    """
    
    monoclade = {} ## Dictionary for final clade assignments
    tool = 'pplacer_'+mode
    with open(pplacer_tree) as f:
        for ln in f.read().splitlines():
            t = Tree(ln)
            tree_hits = {}
            ## Identify the list of pplacer query placements
            for leaf in t.get_leaves():
                ln = leaf.name.split("_")
                if re.match("^#\d+$", ln[-2]) is not None: ## matches pplacer query formatting
                    n = re.sub(r"^#(\d+)$", "\g<1>", ln[-2])
                    tree_hits[n] = leaf
                    leaf2clade[leaf.name] = 'query_seq'
                    ## Different methods
                    if mode == 'top_down':
                        to_check = []
                        running_mass = 0
                        for placement in sorted(tree_hits): ## NOTE: unsure about behavior when n>=10
                            running_mass += float(re.sub(r"^.+#\d+_M=(\d+?\.?\d*)$", "\g<1>", tree_hits[placement].name))
                            to_check.append(placement)
                            if(running_mass >= masscutoff): ## threshold reached
                                break
                        clade_assignments = {} ## Clade to count
                        for q in to_check:
                            pref, clade_assignment = get_clade(q, tree_hits, t, leaf2clade)
                            if clade_assignment in clade_assignments:
                                clade_assignments[clade_assignment] += 1
                            else:
                                clade_assignments[clade_assignment] = 1;
                        if len(clade_assignments) == 1: ## A single assignment made
                            for a in clade_assignments:
                                monoclade[pref] = CladeAssignmentResults(pref, a, tool, None)
                        else: ## None or more than one assignment made
                            monoclade[pref] = CladeAssignmentResults(pref, 'clade_not_conserved', tool, None)
                    elif mode == 'additive':
                        totalmass = {}
                        pref = ''
                        for placement in tree_hits:
                            mass = float(re.sub(r"^.+#\d+_M=(\d+?\.?\d*)$", "\g<1>", tree_hits[placement].name))
                            pref, clade_assignment = get_clade(placement, tree_hits, t, leaf2clade)
                            if clade_assignment in totalmass:
                                totalmass[clade_assignment] += mass
                            else:
                                totalmass[clade_assignment] = mass
                        best_clade = max(totalmass, key=totalmass.get)
                        
                        if totalmass[best_clade] >= masscutoff: ## Best clade over threshold
                            monoclade[pref] = CladeAssignmentResults(pref, best_clade, tool, None)
                        else: ## None or more than one assignment made
                            monoclade[pref] = CladeAssignmentResults(pref, 'clade_not_conserved', tool, None)
                    else:
                        raise RuntimeError("Unrecognized mode in parse_pplacer() (in common/pplacer.py)")
    return monoclade

def pplacer_clade_assignment(leaf2clade_tbl: str, reference_tree: str, reference_alignment: str, reference_package: str, alignment: Dict[str, str], masscutoff: float) -> Dict[str, str]: ## TODO: update with return type
    """ Pipleline for placing query sequences onto a precalculated phylogeny

        Example pplacer reference package creation (requires taxtastic, available on bioconda):
        taxit create --aln-fasta example.afa --tree-stats RAxML_info.example.nwk --tree-file RAxML_bestTree.example.nwk -P example.refpkg -l example

        Arguments:
            leaf2clade_tbl: path to the leaf to clade map
                            tab separated, columns= Leafname, Cladename, Cladedesc
                            assumes leaf names match reference tree names
            reference_tree: path to precomputed tree (newick format)
            reference_alignment: path to the reference alignment
            reference_package: path to the pplacer reference package, see above comment for generation
            alignment: fasta Dictionary of seq names (str) to seqs (str)
            masscutoff: pplacer parameter that dials confidence from 0 (low) to 1 (high). Default=0.9

        Returns:                                                                                                                             dictionary of query_name -> clade_assignment
    """
    config = get_config
    
    ## Runs pplacer and guppy to summarize to a single tree
    pplacer_tree = subprocessing.run_pplacer(reference_tree, reference_alignment, reference_package, alignment)
    
    leaf2clade = {}
    with open(leaf2clade_tbl) as c:
        for ln in c.read().splitlines():
            leafname, clade, ann = ln.split("\t")
            leaf2clade[leafname] = clade

    ## Currently there are two methods for clade assignment based on the mass cutoff
    ## mode=additive computes total masses for each clade assignment and returns the
    ## best assignment above the cutoff
    ## Alternatively, mode=top_down starts with the clade of the placement with the
    ## highest mass, and then adds only if the next highest placement is of the same
    ## clade; returns best clade over cutoff
    ## additive is the default, could potentially build this in as a switch later on
    clade_assignment = parse_pplacer(leaf2clade, pplacer_tree, masscutoff, 'additive', False)

    return clade_assignment
