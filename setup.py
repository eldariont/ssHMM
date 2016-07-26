from setuptools import setup
from setuptools.extension import Extension

setup(name='sshmm',
      version='0.1',
      description='A sequence-structure hidden Markov model for the analysis of RNA-binding protein data.',
      url='https://github.molgen.mpg.de/heller/RBPbinding_Paper',
      author='David Heller, Annalisa Marsico',
      author_email='heller_d@molgen.mpg.de',
      license='MIT',
      packages=['sshmm'],
      package_data={'sshmm': ['img/*.png'],},
      zip_safe=False,
      install_requires=['numpy', 'graphviz', 'pygraphviz', 'weblogo', 'forgi'],
      scripts=['bin/preprocess_dataset', 'bin/train_seqstructhmm', 'bin/batch_seqstructhmm'],
      ext_modules = [Extension("cEstimate", ["sshmm/cEstimate.c"])])
