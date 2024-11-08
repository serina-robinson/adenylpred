AdenylPred
------------

Welcome to the AdenylPred github repo! AdenylPred is a substrate and function prediction tool for class I adenylate-forming enzymes. Class I adenylate-forming enzymes have diverse functions including the acyl-CoA synthetases, NRPS A-domains, firefly luciferases, fatty acyl-AMP ligases, and β-lactone synthetases. These enzymes play critical roles in primary and secondary metabolism in all branches of the tree of life.

The number of adenylate-forming enzymes in sequence databases (>700,000) far outnumber our capacity to experimentally characterize them. Since these enzymes activate a variety of fatty, aryl, and amino acid precursors in biosynthetic pathways, prediction of substrate can help inform the chemical structure of downstream metabolites. To meet this challenge, AdenylPred was developed using a random forest machine learning approach to predict substrate specificity from amino acid sequence. For more information, see our paper in JBC:

Robinson SL, Terlouw BR, Smith MD, Pidot SJ, Stinear TP, Medema MH, Wackett LP. Global analysis of adenylate-forming enzymes reveals β-lactone biosynthesis pathway in pathogenic Nocardia. J Biol Chem. 2020 Oct 30;295(44):14826-14839. doi: [10.1074/jbc.RA120.013528](https://pubmed.ncbi.nlm.nih.gov/32826316/). PMID: 32826316.

![](https://github.com/serina-robinson/adenylpred/blob/master/data/ml_workflow.png)

Installation
------------
For running small queries (<50 sequences), an AdenylPred web app is available at [z.umn.edu/adenylpred](https://srobinson.shinyapps.io/AdenylPred/). For computationally-intensive queries we request use of the command-line version here.

Prerequisite software and packages:
* python (version 3.7.3 tested, any version >= 3.5.0 should work)
* python-virtualenv (not needed, but highly recommended)
* [muscle](http://www.drive5.com/muscle/downloads.htm) (version 3.8.1551 tested) or if you have conda installed you can use `conda install -c bioconda muscle `
* [hmmer2](http://hmmer.org/) (HMMER 2.3.2 tested) or if you have conda installed you can use `conda install -c bioconda hmmer2`

Create a python virtual environment for installing AdenylPred dependencies. If conda is not installed, see [venv documentation](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) for more information. If conda is install you can use the following commands: 

```
conda create --name adenylpred_env
```
Then activate your environment using:
```
source activate adenylpred_env
```

Next navigate to a location on your local drive where you would like AdenylPred installed. To obtain a copy of the AdenylPred source code you will need to clone the AdenylPred git repository:

```
git clone https://github.com/serina-robinson/adenylpred.git
cd adenylpred
```

Or, from github.com, click the "clone or download" button and "Download ZIP"

All python dependencies are specified in the `requirements.txt` file. To load all the requirements simply run:
```
conda install -c bioconda hmmer2
conda install -c bioconda muscle
pip3 install -r requirements.txt
```

You should then be able to run adenylpred as follows:

```
usage: python3 adenylpred.py [-h] -i INPUT [-o OUTPUT] [-v] [-s] [-n] [-g]

Prediction tool for adenylate-forming enzyme substrate specificity

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input file (FASTA or GenBank format).
  -o OUTPUT, --output OUTPUT
                        Output file directory. Default is stdout
  -v, --verbose         Verbose. Prints progress to stdout
  -n, --nucleotide      Nucleotide sequence
  -g, --genbank_input   Input is in GenBank format
```

Common issues include not having hmmer2 or muscle installed, so please check the dependencies before using AdenylPred.

Example Usages
--------------

Predict substrates of A domains from a complete biosynthetic gene cluster FASTA file
```
python3 adenylpred.py -v -i examples/lipstatin.fasta -o examples/lipstatin_predictions.txt
```

Predict substrates for AMP-binding enzymes in a nucleotide FASTA file
```
python3 adenylpred.py -v -n -i examples/lipstatin_nucleotide.fasta -o examples/lipstatin_nucleotide_predictions.txt
```

Predict substrates from a GenBank file
```
python3 adenylpred.py -v -g -i examples/daptomycin.gbk -o examples/daptomycin_predictions.txt
```

Note
-------
This command-line version of AdenylPred does not have all the features of [z.umn.edu/adenylpred](https://srobinson.shinyapps.io/AdenylPred/). Predictions and probability scores may differ from the web application. Please contact Serina Robinson robinsonserinalee[at]gmail.com with questions.

Acknowledgements
-------
This tool is the product of a collaboration between Dr. Larry Wackett’s lab at the University of Minnesota and Dr. Marnix Medema’s lab at Wageningen University and Research, the Netherlands. It is supported a National Science Foundation Graduate Research Fellowship under Grant No. 00039202. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

Cite us
-------
If you use either the web or command-line versions of AdenylPred, or any of the code associated, please cite us at: 

Robinson SL, Terlouw BR, Smith MD, Pidot SJ, Stinear TP, Medema MH, Wackett LP. Global analysis of adenylate-forming enzymes reveals β-lactone biosynthesis pathway in pathogenic Nocardia. J Biol Chem. 2020 Oct 30;295(44):14826-14839. doi: [10.1074/jbc.RA120.013528](http://dx.doi.org/10.1074/jbc.RA120.013528). PMID: 32826316

Key contributors
-------
* Serina Robinson robinsonserinalee[at]gmail.com
* Barbara Terlouw
* Marnix Medema
* Larry Wackett

License
-------
AdenylPred is an open source tool available under the GNU Affero General Public
License version 3.0 or greater. See the [`LICENSE.txt`](LICENSE.txt) file for
details.
