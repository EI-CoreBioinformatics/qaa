# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys


here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
    description = long_description = description.read()

name = "qaa"
version = "0.4"

if sys.version_info.major != 3:
    raise EnvironmentError("""{} is a python module that requires python3, and is not compatible with python2.""".format(name))


setup(
	name=name,
	version=version,
	description=description,
	long_description=long_description,
	url="https://github.com/EI-CoreBioinformatics/qaa",
	author="Christian Schudoma",
	author_email="christian.schudoma@earlham.ac.uk",
	license="MIT",
	classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific Engineering :: Bio/Informatics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        'Programming Language :: Python :: 3.4',
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6"
    ],
	zip_safe=False,
	keywords="genome assembly quality control",
	packages=find_packages(exclude=["test"]),
	scripts=["bin/qaa_sub"],
	install_requires=[
    	'snakemake>=4.4.0',
    	'drmaa'
    ],
	entry_points={
		"console_scripts": [
			"qaa = qaa.__main__:main"
		]
	},
	package_data={
		"qaa.zzz": ["*.smk.py"]
	},
	include_package_data=True,
	data_files=[
		("etc", ["etc/qaa_config.yaml", "etc/hpc_config.json"]), 
		("etc/util", ["etc/util/busco_init_dir"])
	]
)
