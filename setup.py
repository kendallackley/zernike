"""Setup script for the Zernike package."""
import sys
from setuptools import setup
from codecs import open
import os
from re import compile as re_compile

PACKAGES = ['zernike'
            ]

REQUIRES = ['astromatic-wrapper']
# REQUIRES = ['os',
#             'sys',v
#             'glob',
#             'subprocess',
#             'shutil',
#             'math',
#             'numpy',
#             'astropy',
#             'hotpants',
#             'scipy',
#             'time',
#             'logging',
#             'configobj',
#             'spalipy',
            # 'multiprocessing'
            # ]

# Get the version string
__version__ = None
with open('zernike/version.py') as f:
    exec(f.read())  # Should set __version__


def read(filename):
    kwds = {"encoding": "utf-8"} if sys.version_info[0] >= 3 else {}
    with open(filename, **kwds) as fp:
        contents = fp.read()
    return contents


setup(name='zern',
      version=__version__,
      description='Zernike Shapelet analysis',
      url='http://github.com/kendallackley/zernike',
      author='Kendall Ackley',
      author_email='kendall.ackley@monash.edu',
      python_requires = '>=3.6.0',
      license="MIT",
      classifiers=[
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Topic :: Scientific/Engineering :: Astronomy",
          "Topic :: Scientific/Engineering :: Physics"
      ],
      install_requires=REQUIRES,
      packages=PACKAGES,
      )
