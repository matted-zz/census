#!/usr/bin/env python
from setuptools import setup

setup(name="census",
      version="0.9.1",
      description="Sequencing library complexity estimation.",
      author="Matt Edwards",
      author_email="matted@mit.edu",
      license="MIT",
      url="https://github.com/matted/census",
      packages=[],
      scripts=["calculate_libsize.py", "bam_to_histo.py"],
      zip_safe=True,
      install_requires=["scipy", "numpy", "pysam"],
     )
