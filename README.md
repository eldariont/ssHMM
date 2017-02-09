## ssHMM - Sequence-structure hidden Markov model

ssHMM is an RNA motif finder. It recovers sequence-structure motifs from RNA-binding protein data, such as CLIP-Seq data.

### Background

RNA-binding proteins (RBPs) play a vital role in the post-transcriptional control of RNAs. They are known to recognize RNA molecules by their nucleotide sequence as well as their three-dimensional structure. ssHMM combines a hidden Markov model (HMM) with Gibbs sampling to learn the joint sequence and structure binding preferences of RBPs from high-throughput RNA-binding experiments, such as CLIP-Seq. The model can be visualized as an intuitive graph illustrating the interplay between RNA sequence and structure.

### Scope

ssHMM was developed for the analysis of data from RNA-binding assays. Its aim is to help biologists to derive a binding motif for one or a number of RNA-binding proteins. ssHMM was written in Python and is a pure command-line tool.

### Documentation

Check out our documentation: http://sshmm.readthedocs.io

### Citation

Heller, D., Krestel, R., Ohler, U., Vingron, M., & Marsico, A. (2016). ssHMM: [Extracting intuitive sequence-structure motifs from high-throughput RNA-binding protein data](http://dx.doi.org/10.1101/076034). bioRxiv, 076034.
