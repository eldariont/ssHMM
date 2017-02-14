from setuptools import setup
from setuptools.extension import Extension

long_description = """RNA-binding proteins (RBPs) play a vital role in the post-transcriptional control of RNAs. They are known to recognize RNA molecules by their nucleotide sequence as well as their three-dimensional structure. ssHMM is an RNA motif finder that combines a hidden Markov model (HMM) with Gibbs sampling to learn the joint sequence and structure binding preferences of RBPs from high-throughput RNA-binding experiments, such as CLIP-Seq. The model can be visualized as an intuitive graph illustrating the interplay between RNA sequence and structure."""

setup(name='sshmm',
      version='1.0.0b4',
      description='A sequence-structure hidden Markov model for the analysis of RNA-binding protein data.',
      long_description=long_description,
      url='https://github.molgen.mpg.de/heller/ssHMM',
      author='David Heller, Annalisa Marsico',
      author_email='heller_d@molgen.mpg.de',
      license='GPLv3',
      classifiers=[
      'Development Status :: 4 - Beta',
      'Environment :: Console',
      'Intended Audience :: Science/Research',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      'Programming Language :: Python :: 2.7'
      ],
      keywords='sshmm CLIP-Seq RBP RNA RNA-binding protein motif HMM sequence structure',
      packages=['sshmm'],
      package_data={'sshmm': ['img/*.png'],},
      zip_safe=False,
      install_requires=['numpy', 'graphviz', 'pygraphviz', 'weblogo', 'forgi'],
      scripts=['bin/preprocess_dataset', 'bin/train_seqstructhmm', 'bin/batch_seqstructhmm'],
      ext_modules = [Extension("cEstimate", ["sshmm/cEstimate.c"])])
