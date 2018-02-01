import pkg_resources

__title__ = "qaa"
__author__ = 'Christian Schudoma (cschu)'
__license__ = 'MIT'
__copyright__ = 'Copyright 2018 Earlham Institute'
__version__ = pkg_resources.require("qaa")[0].version

import os

ETC_DIR = os.path.join(os.path.dirname(__file__), "..", "etc")
DEFAULT_HPC_CONFIG_FILE = os.path.join(ETC_DIR, "hpc_config.json")
DEFAULT_CONFIG_FILE = os.path.join(ETC_DIR, "qaa_config.yaml")

import csv
from collections import namedtuple

import yaml

QAA_SAMPLESHEET_COLUMNS = "id assembly bamfile r1 r2 busco_id transcripts proteins".split(" ")
QAA_Sample = namedtuple("QAA_Sample", QAA_SAMPLESHEET_COLUMNS)

def readSamplesheet(_in):
	import csv
	for row in csv.reader(_in, delimiter=","):
		if row and row[0]:
			while len(row) != len(QAA_SAMPLESHEET_COLUMNS):
				row.append("")
			yield (row[0], QAA_Sample(*row))


TIME_CMD = " /usr/bin/time -v"

def loadPreCmd(command, is_dependency=True):
        '''
        Used to prefix a shell command that utilises some external software with another command used to load that software
        '''
        if command:
                cc = command.strip()
                if cc != "":
                        if is_dependency:
                                return "set +u && {} &&".format(cc)
                        else:
                                return " {} ".format(cc)

        return ""

from eicore.external_process.snakemake_helper import *

class QAA_Runner(object):
	def __init__(self, qaa_args):
		print("Configuring execution environment ... ", end="", flush=True)
		self.output_dir = qaa_args.output_dir
		self.logs_dir = os.path.join(self.output_dir, "hpc_logs")
		self.exe_env = ExecutionEnvironment(qaa_args, NOW, job_suffix=qaa_args.input + "_" + self.output_dir, log_dir=self.logs_dir)
		print("done.")
		print(str(self.exe_env))
		print()

		if qaa_args.config:
			print("Custom configuration file specified, overriding defaults ... ", end="", flush=True)
			self.config = yaml.load(open(qaa_args.config))
			self.config_file = qaa_args.config
			print("done.")
			print()
		else:
			print("Loading default configuration ... ", end="", flush=True)
			self.config = yaml.load(open(DEFAULT_CONFIG_FILE))
			self.config_file = DEFAULT_CONFIG_FILE
			print("done.")
			print()

		if os.path.exists(self.output_dir):
			print("Output directory already exists, attempting to resume ...")
		else:
			print("Output directory does not exist. Creating ... ", end="", flush=True)
			os.makedirs(self.output_dir)
			print("done.")

		if not os.path.exists(self.logs_dir) and self.exe_env.use_scheduler:
			print("HPC log dir does not exist. Creating ... ", end="", flush=True)
			os.makedirs(self.logs_dir)
			print("done.")

		print()
			
		self.config["samplesheet"] = qaa_args.input
		self.config["out_dir"] = self.output_dir
		self.config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
		self.config["cwd"] = os.getcwd()

		self.config["blobtools_run_bwa"] = not qaa_args.blobtools_no_bwa
		requested_modes = qaa_args.qaa_mode.split(",")
		modes = ("geno", "genome", "tran", "transcriptome", "prot", "proteome")
		invalid_modes = list(filter(lambda s:s not in modes, requested_modes))
		if invalid_modes:
			valid_modes = list(set(requested_modes).difference(invalid_modes))
			if valid_modes:
				print("--qaa-mode: Found invalid modes {}. Proceeding with {}.".format(invalid_modes, valid_modes))
				requested_modes = valid_modes
			else:
				print("--qaa-mode: No valid mode provided. Valid modes are {}. Exiting.".format(",".join(modes)))
				sys.exit(1)

		self.config["survey_assembly"] = qaa_args.survey_assembly
		self.config["run_genome_module"] = qaa_args.survey_assembly or ("geno" in requested_modes or "genome" in requested_modes)
		self.config["run_transcriptome_module"] = not qaa_args.survey_assembly and ("tran" in requested_modes or "transcriptome" in requested_modes)
		self.config["run_proteome_module"] = not qaa_args.survey_assembly and ("prot" in requested_modes or "proteome" in requested_modes)

		self.new_config_file = os.path.join(self.output_dir, "qaa.conf.xml")
		with open(self.new_config_file, "w") as conf_out:
			yaml.dump(self.config, conf_out, default_flow_style=False)

		self.unlock = qaa_args.unlock

	def run(self):		
		return run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "qaa.smk.py"), self.output_dir, self.new_config_file, self.exe_env, dryrun=False, unlock=self.unlock)
	

class QAA_ArgumentsAdapter(object):
	def __init__(self):
		pass
