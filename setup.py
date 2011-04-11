from setuptools import setup
from setuptools.extension import Extension
import os

try:
    from Cython.Distutils import build_ext
except:
    print "You don't seem to have Cython installed. Please get a"
    print "copy from www.cython.org and install it,"
    print "or if you have easy_install, try"
    print 
    print "  easy_install cython"
    print
    sys.exit(1)

from Cython.Distutils import build_ext

import numpy

def scandir(dir, files=[]):
    """
    scan the directory for extension files, converting
    them to extension names in dotted notation
    """
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if os.path.isfile(path) and path.endswith(".pyx"):
            files.append(path.replace(os.path.sep, ".")[:-4])
        elif os.path.isdir(path):
            scandir(path, files)
    return files

def make_extension(extName):
    """
    generate an Extension obj from dotted name
    """
    extPath = extName.replace(".", os.path.sep)+".pyx"
    return Extension(
        extName,
        [extPath],
        include_dirs = [numpy.get_include(), "."],   # adding the '.' to include_dirs is CRUCIAL!!
        )

# See http://wiki.cython.org/PackageHierarchy
#ext_modules = [
#                Extension("genomicfeatures._AbstractBaseClasses", ["src/_AbstractBaseClasses.pyx"], include_dirs=[numpy.get_include(),'.']),
#                Extension("genomicfeatures._genomicfeatures", ["src/_genomicfeatures.pyx"], include_dirs=[numpy.get_include(),'.']),
#                Extension("genomicfeatures._Window", ["src/_Window.pyx"], include_dirs=[numpy.get_include(),'.']),
#                Extension("genomicfeatures._Scores", ["src/_Scores.pyx"], include_dirs=[numpy.get_include(),'.']),
#
#              ]

ext_names = scandir('genomicfeatures')
extensions = [make_extension(i) for i in ext_names]

setup(
  name = 'genomicfeatures',
  cmdclass = {'build_ext': build_ext},
  ext_modules = extensions,
  packages = ['genomicfeatures',],
)
