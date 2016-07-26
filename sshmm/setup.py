from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("cEstimate.pyx", gdb_debug=True)
)