"""Setup script for ``phydms``.

Written by Jesse Bloom.
"""

import sys
import os
import re
import glob
import fnmatch
try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    raise ImportError("You must install setuptools")

if not ((sys.version_info[0] == 2 and sys.version_info[1] == 7) or
        (sys.version_info[0] == 3 and sys.version_info[1] >= 4)):
    raise RuntimeError('phydms is currently only tested with Python 2.7, or Python 3.4 or higher.\nYou are using Python %d.%d' % (sys.version_info[0], sys.version_info[1]))

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

with open('README.rst') as f:
    readme = f.read()

# main setup command
setup(
    name = 'phydms', 
    version = metadata['version'],
    author = metadata['author'],
    author_email = metadata['author_email'],
    url = metadata['url'],
    download_url = 'https://github.com/jbloomlab/phydms/tarball/%s' % metadata['version'], # assumes appropriate tagged version is on GitHub
    description = 'Phylogenetic analyses informed by deep mutational scanning data.',
    long_description = readme,
    license = 'GPLv3',
    install_requires = [
        'biopython>=1.67',
        'scipy>=0.16',
        'matplotlib>=1.5.1',
        'natsort>=5.0.1',
        ],
    packages = ['phydmslib'],
    package_dir = {'phydmslib':'phydmslib'},
    scripts = [
            'scripts/phydms',
#            'scripts/phydms_comprehensive',
            'scripts/phydms_prepalignment',
#            'scripts/phydms_renumber',
#            'scripts/phydms_plotselection',
#            'scripts/phydms_analyzeselection',
#            'scripts/phydms_testdivpressure',
            ],
)
