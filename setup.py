'''pygeo
'''

# pygeo - a distribution of tools for managing geophysical data
# Copyright (C) 2015 Brendan Smithyman

# This file is part of pygeo.

# pygeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# pygeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with pygeo.  If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------

from distutils.core import setup
from distutils.extension import Extension
from setuptools import find_packages
from Cython.Build import cythonize
import numpy as np

CLASSIFIERS = [
'Development Status :: 4 - Beta',
'Intended Audience :: Developers',
'Intended Audience :: Science/Research',
'License :: OSI Approved :: License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
'Programming Language :: Python',
'Topic :: Scientific/Engineering',
'Topic :: Scientific/Engineering :: Mathematics',
'Topic :: Scientific/Engineering :: Physics',
'Operating System :: Microsoft :: Windows',
'Operating System :: POSIX',
'Operating System :: Unix',
'Operating System :: MacOS',
'Natural Language :: English',
]

import os, os.path

with open('README.md') as fp:
    LONG_DESCRIPTION = ''.join(fp.readlines())

extensions = [
    Extension('analysis', ['analysis.pyx', 'filter.c']),
    Extension('autopick', ['autopick.pyx', 'filter.c']),
    Extension('dipfilt', ['dipfilt.pyx']),
    Extension('segyread', ['segyread/segyread.pyx', 'segyread/fpconvert.c']),
]
    

sourcefiles = ['segyread/segyread.pyx', 'segyread/fpconvert.c']
extensions.append(Extension('segyread', sourcefiles))

setup(
    name = 'Pygeo',
    version = '0.1',
    packages = find_packages(),
    install_requires = ['numpy',
                        'scipy',
                        'Cython',
                        'pyopencl',
                        ],
    author = 'Brendan Smithyman',
    author_email = 'brendan@bitsmithy.net',
    description = 'pygeo',
    long_description = LONG_DESCRIPTION,
    license = 'LGPL',
    keywords = 'exploration seismology geophysics',
    download_url = 'http://github.com/bsmithyman/pygeo',
    classifiers = CLASSIFIERS,
    platforms = ['Linux', 'Solaris', 'Mac OS-X', 'Unix'],
    use_2to3 = False,
    ext_modules = cythonize(extensions),
)
    
