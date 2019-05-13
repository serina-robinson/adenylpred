AdenylPred
===========
Substrate and function prediction tool for class I adenylate-forming enzymes.

Class I adenylate-forming enzymes including acyl-CoA synthetases, NRPS A-domains, firefly luciferases, fatty acyl-AMP ligases, and β-lactone synthetases. They are critical for 

As the number of sequences for adenylate-forming enzymes in databases far outnumber our capacity for experimental characterization, there is a need for adenylation enzyme substrate prediction tools. As a result, AdenylPred was developed using a random forest machine learning approach to predict substrate specificity from amino acid sequence.

This tool is the product of a collaboration between Dr. Larry Wackett’s lab at the University of Minnesota and Dr. Marnix Medema’s lab at Wageningen University and Research, the Netherlands. It is supported a National Science Foundation Graduate Research Fellowship under Grant No. 00039202. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

Installation
------------
Prerequisite software and packages:
* python (version 3.7.3 tested, any version >= 3.5.0 should work)
* python-virtualenv (not needed, but highly recommended)
* [muscle](http://www.drive5.com/muscle/downloads.htm) (version 3.8.1551 tested)
* [hmmer3](http://hmmer.org/) (version 3.2.1 tested)

First, we recommend you create a python virtualenv for installing AdenylPred dependencies

```
virtualenv --python $(which python3) adenylpred_env
source activate adenylpred_env
```

Next obtain a copy of AdenylPred source code. You will need to clone the AdenylPred git repository:

```
git clone https://github.com/serina-robinson/adenylpred/adenylpred.git
cd adenylpred
```

Or, from github.com, click the "clone or download" button and "Download ZIP"

Next, all python dependencies are specified in the `setup.py` file

License
-------
AdenylPred is an open source tool available under the GNU Affero General Public
License version 3.0 or greater. See the [`LICENSE.txt`](LICENSE.txt) file for
details.
