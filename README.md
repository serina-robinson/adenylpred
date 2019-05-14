AdenylPred
===========
Welcome to AdenylPred! AdenylPred is a substrate and function prediction tool for class I adenylate-forming enzymes.

Class I adenylate-forming enzymes have diverse functions including the acyl-CoA synthetases, NRPS A-domains, firefly luciferases, fatty acyl-AMP ligases, and β-lactone synthetases. These enzymes play critical roles in primary and secondary metabolism in all branches of the tree of life.

The number of adenylate-forming enzymes in sequence databases (>700,000) far outnumber our capacity to experimentally characterize them. Since these enzymes activate a variety of fatty, aryl, and amino acid precursors in biosynthetic pathways, prediction of substrate can help inform the chemical structure of downstream metabolites. To meet this challenge, AdenylPred was developed using a random forest machine learning approach to predict substrate specificity from amino acid sequence.

![](https://github.com/serina-robinson/adenylpred/blob/master/data/ml_workflow.png)

Installation
------------
For running small queries (<50 sequences), an AdenylPred web app is available at [z.umn.edu/adenylpred](z.umn.edu/adenylpred). For computationally-intensive queries we request use of the command-line version here.

Prerequisite software and packages:
* python (version 3.7.3 tested, any version >= 3.5.0 should work)
* python-virtualenv (not needed, but highly recommended)
* [muscle](http://www.drive5.com/muscle/downloads.htm) (version 3.8.1551 tested) or if you have conda installed you can use `conda install -c bioconda muscle `
* [hmmer3](http://hmmer.org/) (version 3.2.1 tested) or if you have conda installed you can use `conda install -c biocore hmmer`

Create a python virtualenv for installing AdenylPred dependencies

```
python3 -m venv /path/to/new/virtual/env
``` 
OR if you are using conda:
```
conda create --name /path/to/new/virtual/env
```
Then activate your environment using
```
source activate /path/to/new/virtual/env
```

Next navigate to a location on your local drive where you would like AdenylPred instaled. Obtain a copy of AdenylPred source code. You will need to clone the AdenylPred git repository:

```
git clone https://github.com/serina-robinson/adenylpred.git
cd adenylpred
```

Or, from github.com, click the "clone or download" button and "Download ZIP"

All python dependencies are specified in the `requirements.txt` file. To load all the requirements simply run:
```
pip install -r requirements.txt
```

You can then run adenylpred as follows:

```
usage: python3 adenylpred.py [-h] -i INPUT [-o OUTPUT] [-x XTRACT_A_DOMAINS]

  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file (FASTA format).
  -o OUTPUT, --output OUTPUT
                        Output file directory. Default is stdout
  -x XTRACT_A_DOMAINS, --xtract_A_domains XTRACT_A_DOMAINS
                        [1 = extract AMP-binding hits using AMP-binding.hmm
                         0 = AMP-binding domains already extracted]
```

Example Usages
--------------

Predict substrates of A domains from a complete biosynthetic gene cluster FASTA file
```
python3 adenylpred.py -i examples/lipstatin.fasta -o lipstatin_predictions.txt
```

Predict substrates for example NRPS sequences where the AMP-binding domains have already been extracted
```
python3 adenylpred.py -i examples/test5.fasta -x 0
```

Acknowledgements
-------
This tool is the product of a collaboration between Dr. Larry Wackett’s lab at the University of Minnesota and Dr. Marnix Medema’s lab at Wageningen University and Research, the Netherlands. It is supported a National Science Foundation Graduate Research Fellowship under Grant No. 00039202. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

Key contributors
-------
* Serina Robinson (robi0916@umn.edu)
* Barbara Terlouw
* Marnix Medema
* Larry Wackett

License
-------
AdenylPred is an open source tool available under the GNU Affero General Public
License version 3.0 or greater. See the [`LICENSE.txt`](LICENSE.txt) file for
details.
