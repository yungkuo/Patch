"""
Compile the extension with:

    python setup.py build_ext --inplace

"""

from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy as np


ext_modules = [Extension("bsearch_c", ["bsearch_c.pyx"])]

setup(
    name = 'Burst search',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],
    ext_modules = ext_modules
)
