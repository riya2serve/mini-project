#!/usr/bin/env python

"""
Call `pip install -e .` to install package locally for testing.
"""

from setuptools import setup

# build command
setup(
    name="mini-project",
    version="0.0.1",
    author="Riya Rampalli",
    author_email="rr3491@columbia.edu",
    license="GPLv3",
    description="A package for simulating a VCF file",
    classifiers=["Programming Language :: Python :: 3"],
    entry_points={
        "console_scripts": ["mini-project = mini-project.__main__:main"]
    },
)
