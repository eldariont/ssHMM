## ssHMM - Sequence-structure hidden Markov model

A motif finder for sequence-structure binding preferences of RNA-binding proteins.

RNA-binding proteins (RBPs) play a vital role in the post-transcriptional control of RNAs. They are known to recognize RNA molecules by their nucleotide sequence as well as their three-dimensional structure. ssHMM is an RNA motif finder that combines a hidden Markov model (HMM) with Gibbs sampling to learn the joint sequence and structure binding preferences of RBPs from high-throughput RNA-binding experiments, such as CLIP-Seq. The model can be visualized as an intuitive graph illustrating the interplay between RNA sequence and structure.

#### Scripts

ssHMM consists of 3 main scripts:
- **preprocess_dataset**: Prepares a CLIP-Seq dataset in BED format for ssHMM. It filters the BED file, fetches the genomic sequences, and predicts RNA secondary structures.
- **train_seqstructhmm**: Trains ssHMM on a given CLIP-Seq dataset and produces an intuitive visualization of the recovered motif.
- **batch_seqstructhmm**: Trains ssHMM on several different CLIP-Seq datasets and recovers one motif for each dataset.

#### Installation

##### Prerequisites:
- GHMM (http://ghmm.org/)
- GraphViz (http://www.graphviz.org/)

##### Installation steps:

1.  Install prerequisites for GHMM as described on http://ghmm.sourceforge.net/installation.html. The commands for Ubuntu are:
  
  ```bash
  sudo apt-get update
  sudo apt-get install build-essential automake autoconf libtool
  sudo apt-get install python-dev
  sudo apt-get install libxml++2.6-dev
  sudo apt-get install swig
  ```
2.  Download and unpack GHMM from https://sourceforge.net/projects/ghmm/
3.  Install GHMM as described on http://ghmm.sourceforge.net/installation.html. The commands for Ubuntu are:
  
  ```bash
  cd ghmm
  sh autogen.sh
  sudo ./configure
  sudo make
  sudo make install
  sudo ldconfig
  ```
4.  Install GraphViz. On Ubuntu:
  
  ```bash
  sudo apt-get install graphviz
  sudo apt-get install libgraphviz-dev
  ```
6.  Install pip if not already installed. On Ubuntu:
  
  ```bash
  sudo apt-get install python-pip
  ```
5.  Install PyGraphViz:
  
  ```bash
  sudo PKG_CONFIG_ALLOW_SYSTEM_LIBS=OHYESPLEASE pip install pygraphviz
  ```
5.  Download ssHMM from this page

7.  Install ssHMM:
  
  ```bash
  sudo python setup.py install`
  ```
  This will install ssHMM and the following python package dependencies: numpy, graphviz, pygraphviz, weblogo, forgi. If setuptools fails to install any of the dependencies, try to install it separately (e.g. with `sudo pip install numpy`).
  
