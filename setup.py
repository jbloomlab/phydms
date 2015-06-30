"""Setup script for ``phydms``.

Written by Jesse Bloom.
"""

import sys
import os
import re
import fnmatch
try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    raise ImportError("You must install setuptools")

if (sys.version_info[0], sys.version_info[1]) != (2, 7):
    raise RuntimeError('phydms is currently only compatible with Python 2.7.\nYou are using Python %d.%d' % (sys.version_info[0], sys.version_info[1]))


# get metadata, which is specified in another file
metadata = {}
lines = open('phydmslib/_metadata.py').readlines()
for dataname in ['version', 'author', 'author_email', 'url']:
    for line in lines:
        entries = line.split('=')
        assert len(entries) == 2, "Failed to parse metadata:\n%s" % line
        if entries[0].strip() == '__%s__' % dataname:
            if dataname in metadata:
                raise ValueError("Duplicate metadata for %s" % dataname)
            else:
                metadata[dataname] = entries[1].strip()[1 : -1]
    assert dataname in metadata, "Failed to find metadata for %s" % dataname


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
    """List of cython extensions for *lazy_cythonize*"""
    from Cython.Build import cythonize
    include_dirs = []
    bpp_sources = []
    for bpplib in ['bpp-core', 'bpp-seq', 'bpp-phyl']:
        include_dirs.append('phydmslib/Bpp/%s/src/' % bpplib)
        for (path, dirs, files) in os.walk('phydmslib/Bpp/%s/src' % bpplib):
            for fname in fnmatch.filter(files, '*.cpp'):
                bpp_sources.append('%s/%s' % (path, fname))
    ext = [\
        Extension(\
            'phydmslib.pybpp',\
            sources=['phydmslib/pybpp.pyx', 'phydmslib/BppExtensions/BppTreeLikelihood.cpp', 'phydmslib/BppExtensions/ExperimentallyInformedCodonModel.cpp'] + bpp_sources,\
            language='c++',\
            include_dirs=include_dirs,\
            extra_compile_args=['-O2'],\
            ),\
        ]
    #commenting out command for linking to existing Bpp libraries:        
    #ext = [Extension('phydmslib.pybpp', sources=['phydmslib/pybpp.pyx', 'phydmslib/BppExtensions/BppTreeLikelihood.cpp', 'phydmslib/BppExtensions/ExperimentallyInformedCodonModel.cpp'], language='c++', extra_compile_args=['-I%s/.local/include/' % os.path.expanduser("~")], extra_link_args=['-L%s/.local/lib/' % os.path.expanduser('~'), '-lbpp-core', '-lbpp-seq', '-lbpp-phyl']),]\
    return cythonize(ext)


with open('README.rst') as f:
    readme = f.read()


# main setup command
setup(
    name = 'phydms', 
    version = metadata['version'],
    author = metadata['author'],
    author_email = metadata['author_email'],
    url = metadata['url'],
    download_url = 'https://github.com/jbloom/phydms/tarball/%s' % metadata['version'], # assumes appropriate tagged version is on GitHub
    description = 'Phylogenetic analyses informed by deep mutational scanning data.',
    long_description = readme,
    license = 'GPLv3',
    install_requires = [\
        'biopython>=1.65',\
        'cython>=0.21',\
        'dms_tools>=1.1.2',\
        ],
    platforms = 'Linux',
    packages = ['phydmslib'],
    package_dir = {'phydmslib':'phydmslib'},
    ext_modules = lazy_cythonize(extensions),
    scripts = [
            'scripts/phydms',
            'scripts/phydms_comprehensive',
            ],
)
