from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy
ext_modules = [Extension("_genomicfeatures", ["_genomicfeatures.pyx"], include_dirs=[numpy.get_include()])]

setup(
  name = 'genomicfeatures',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
