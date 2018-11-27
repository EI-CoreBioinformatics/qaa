# coding: utf-8

"""Setup file for PyPI"""

import sys
import os
import shutil

from setuptools import setup, find_packages, Command
from subprocess import check_call

name = "qaa"
version = "0.3"
release = version + ".0"


if sys.version_info.major != 3:
	raise EnvironmentError("""qaa is a python module that requires python3,
    and is not compatible with python2.""")


class InstallAll(Command):
	description="Build qaa, documentation and set python paths"
	user_options=[('prefix=', 'd', "directory to install the files to")]

	def initialize_options(self):
		self.prefix = None

	def finalize_options(self):
		pass

	def run(self):

		pythonpath = None
		env = os.environ.copy()
		if self.prefix:
			pythonpath = self.prefix + "/lib/python" + str(sys.version_info.major) + "." + str(
				sys.version_info.minor) + "/site-packages"
			print()
			print("Custom install prefix provided: " + self.prefix)
			print("PYTHONPATH=" + pythonpath)
			os.makedirs(pythonpath, exist_ok=True)
			env["PYTHONPATH"] = pythonpath

		if os.path.exists('.git'):
			print()
			print("Updating git submodules")
			check_call(['git', 'submodule', 'init'])
			check_call(['git', 'submodule', 'update'])

			print()
			print("Building eicore")
			mqc_cmd = [sys.executable, "setup.py", "install"]
			if self.prefix:
				mqc_cmd.append("--prefix=" + self.prefix)

			check_call(mqc_cmd, cwd=os.path.join("deps", "eicore"), env=env)

		print()
		print("Building qaa")
		eim_cmd = [sys.executable, "setup.py", "install"]
		if self.prefix:
			eim_cmd.append("--prefix=" + self.prefix)
		check_call(eim_cmd, env=env)

		print()
		"""
		print("Making documentation")

		mkfile = os.path.join("doc", "Makefile")
		bakmkfile = os.path.join("doc", "Makefile.in")
		with open(bakmkfile, 'r') as mkin:
			with open(mkfile, 'w') as mkout:
				for line in mkin:
					if "<PREFIX>" in line:
						newline = line.replace("<PREFIX>", self.prefix if self.prefix else "/usr/local")
						print(newline, end="", file=mkout)
					else:
						print(line, end="", file=mkout)

		check_call(["make", "html", "man"], cwd=os.path.join("doc"))
		check_call(["make", "install"], cwd=os.path.join("doc"))
		"""


# External python modules that can be gathered from
install_requires = [
    'snakemake',
    'biopython',
    'drmaa'
    ]

setup(
	name=name,
	version=release,
	description="QAA - Quality Assessment for Assemblies",
	long_description='''QAA runs a couple of 3rd party tools and generates reports about the quality of an assembly''',
	author="Christian Schudoma",
	author_email="christian.schudoma@earlham.ac.uk/cschu1981@gmail.com",
	license="MIT",
	zip_safe=False,
	keywords="genome assembly quality control",
	packages=find_packages(exclude=["test"]),
	entry_points={"console_scripts": ["qaa = qaa.__main__:main"]},
	test_suite="nose.collector",

	install_requires=install_requires,
	tests_require = [
		'nose',
    ],
	package_data={
		"qaa.zzz": ["*.smk.py"]
	},
	include_package_data=True,
	data_files=[("etc", ["etc/qaa_config.yaml", "etc/hpc_config.json"]), ("etc/util", ["etc/util/busco_init_dir"])],
	scripts=["bin/qaa_sub"],
	cmdclass={
		'install_all': InstallAll
	}
)
