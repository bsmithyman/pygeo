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

import os
from setuptools import find_packages, setup
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize

CLASSIFIERS = [
'Development Status :: 4 - Beta',
'Intended Audience :: Developers',
'Intended Audience :: Science/Research',
'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
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

with open('README.md') as fp:
    LONG_DESCRIPTION = ''.join(fp.readlines())

NAME = 'pygeo'

genPath = lambda fns: ['%s/__init__.py'%(NAME,)] + ['%s/%s'%(NAME, fn) for fn in fns]
genName = lambda subname: '%s.%s'%(NAME, subname)

extensions = [
    Extension(genName('analysis'),  genPath(['analysis.pyx', 'filter.c'])),
    Extension(genName('autopick'),  genPath(['autopick.pyx', 'filter.c'])),
    Extension(genName('coord'),     genPath(['coord.py'])),
    Extension(genName('dipfilt'),   genPath(['dipfilt.pyx'])),
    Extension(genName('fastio'),    genPath(['fastio.py'])),
    Extension(genName('fullpy'),    genPath(['fullpy.py'])),
    Extension(genName('rsfread'),   genPath(['rsfread.py'])),
    Extension(genName('segyarray'), genPath(['segyarray.py'])),
    Extension(genName('segyread'),  genPath(['segyread.pyx', 'fpconvert.c'])),
    Extension(genName('testing'),   genPath(['testing.py'])),
]

class NumpyBuild(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

setup(
    name = NAME,
    version = '0.1.0',
    packages = find_packages(),
    install_requires = ['numpy>=1.7',
                        'scipy>=0.13',
                        'Cython',
                        ],
    author = 'Brendan Smithyman',
    author_email = 'brendan@bitsmithy.net',
    description = 'pygeo',
    long_description = LONG_DESCRIPTION,
    license = 'LGPL',
    keywords = 'exploration seismology geophysics',
    download_url = 'https://github.com/bsmithyman/pygeo',
    url = 'https://github.com/bsmithyman/pygeo',
    classifiers = CLASSIFIERS,
    platforms = ['Linux', 'Solaris', 'Mac OS-X', 'Unix'],
    use_2to3 = False,
    cmdclass={'build_ext': NumpyBuild},
    setup_requires=['numpy'],
    extra_compile_args=['-w','-fopenmp'],
    extra_link_args=['-fopenmp'],
    ext_modules = cythonize(extensions),
)

