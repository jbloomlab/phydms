"""Setup script for ``phydms``.

Written by Jesse Bloom.
"""

import sys
import os
import re
import glob
try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    raise ImportError("You must install setuptools")

if not ((sys.version_info[0] == 2 and sys.version_info[1] == 7) or
        (sys.version_info[0] == 3 and sys.version_info[1] >= 4)):
    raise RuntimeError('phydms requires Python 2.7, or Python 3.4 or higher.\n'
            'You are using Python {0}.{1}'.format(
            sys.version_info[0], sys.version_info[1]))

# get metadata, which is specified in another file
metadata = {}
lines = open('phydmslib/_metadata.py').readlines()
for dataname in ['version', 'author', 'author_email', 'url']:
    for line in lines:
        entries = line.split('=')
        assert len(entries) == 2, "Failed to parse metadata:\n{0}".format(line)
        if entries[0].strip() == '__{0}__'.format(dataname):
            if dataname in metadata:
                raise ValueError("Duplicate metadata for {0}".format(dataname))
            else:
                metadata[dataname] = entries[1].strip()[1 : -1]
    assert dataname in metadata, "Failed to find metadata {0}".format(dataname)

with open('README.rst') as f:
    readme = f.read()

class lazy_cythonize(list):
    """Lazy evaluation of cythonize so it isn't needed until installed.

    Following this:
    http://stackoverflow.com/questions/11010151/distributing-a-shared-library-and-some-c-code-with-a-cython-extension-module
    """
    def __init__(self, callback):
        self._list = None
        self.callback = callback
    def c_list(self):
        if self._list is None:
            self._list = self.callback()
        return self._list
    def __iter__(self):
        for e in self.c_list(): yield e
    def __getitem__(self, ii):
        return self.c_list()[ii]
    def __len__(self):
        return len(self.c_list())

def extensions():
    """Returns list of `cython` extensions for `lazy_cythonize`."""
    import numpy
    from Cython.Build import cythonize
    ext = [
            Extension('phydmslib.numutils', ['phydmslib/numutils.pyx'],
                    include_dirs=[numpy.get_include()],
                    extra_compile_args=['-Wno-unused-function']),
          ]      
    return cythonize(ext)

# main setup command
setup(
    name = 'phydms',
    version = metadata['version'],
    author = metadata['author'],
    author_email = metadata['author_email'],
    url = metadata['url'],
    download_url = 'https://github.com/jbloomlab/phydms/tarball/{0}'.format(
            metadata['version']), # assumes tagged version is on GitHub
    description = 'Phylogenetic analyses informed by deep mutational scanning data.',
    long_description = readme,
    license = 'GPLv3',
    setup_requires = [
        'cython>=0.21',
        'numpy>=0.18',
        ],
    install_requires = [
        'biopython>=1.67',
        'cython>=0.21',
        'numpy>=1.11',
        'scipy>=0.18',
        'matplotlib>=1.5.1',
        'natsort>=5.0.1',
        'sympy>=1.0',
        'six>=1.10',
        'pandas>=0.18',
        'pyvolve>=0.8.4',
        'statsmodels>=0.8',
        'weblogo>=3.4, <3.6',
        'PyPDF2>=1.26',
        ],
    packages = ['phydmslib'],
    package_dir = {'phydmslib':'phydmslib'},
    package_data = {'phydmslib':['_weblogo_template.eps']},
    ext_modules = lazy_cythonize(extensions),
    scripts = [
            'scripts/phydms',
            'scripts/phydms_comprehensive',
            'scripts/phydms_prepalignment',
            'scripts/phydms_logoplot',
            'scripts/phydms_divpressure',
            ],
)
