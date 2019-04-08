import pkg_resources

__title__ = "qaa"
__author__ = "Christian Schudoma (cschu)"
__email__ = "christian.schudoma@earlham.ac.uk"
__license__ = "MIT"
__copyright__ = "Copyright 2018 Earlham Institute"
__version__ = pkg_resources.require("qaa")[0].version

import os
from os.path import join, dirname, basename
import csv
from collections import namedtuple
import yaml
import sys
from copy import copy
import pathlib

from .snakemake_helper import *
from .workflow_runner import WorkflowRunner
from .qaa_environment import QAA_Environment

from qaa.reporting.busco_report import compileBUSCOReport
from qaa.reporting.quast_report import compileQUASTReport
from qaa.reporting.blobtools_report import compileBlobReport
from qaa.reporting.sixteen_s_reporter import SixteenSReporter, HEADER as sixteenS_header

ETC_DIR = join(dirname(__file__), "etc")
DEFAULT_HPC_CONFIG_FILE = join(ETC_DIR, "hpc_config.json")
DEFAULT_CONFIG_FILE = join(ETC_DIR, "qaa_config.yaml")
QAA_ID = str(__version__) + "XXX"


class QAA_Runner(WorkflowRunner):

	@staticmethod
	def create_samplesheet_stream(args, annotation="prokka"):
		stream, asm_path = list(), ""
		for row in csv.reader(open(args.input), delimiter=","):
			if args.runmode == "asm":
				asm_path = join(args.output_dir, "assembly", row[0], row[0] + ".assembly.fasta")
			elif args.runmode == "survey":
				asm_path = join(args.output_dir, "qc", "tadpole", row[0], row[0] + "_tadpole_contigs.fasta")
			elif args.runmode == "asm,ann":
				asm_path = join(args.output_dir, "annotation", "prokka", row[0], row[0] + ".fna")
			
			if len(row) == 3: # ann-formatted sheet
				r1, r2 = "", ""
			else:
				r1, r2 = row[2:4]
		
			new_row = [row[0], asm_path, "", r1, r2, args.busco_db]

			annotation_path = join(args.output_dir, "annotation", annotation)
			
			if args.runmode in ("ann", "asm,ann"):
				new_row.extend([
					join(annotation_path, row[0], row[0] + ".ffn"), 
					join(annotation_path, row[0], row[0] + ".faa")
				])
			stream.append(new_row)
		return stream


	def __calc_qualimap_memory_setting(self):
		# from eimethyl: Extract qualimap mem setting from hpc_config, convert to gigabytes and subtract a bit to account for memory above
		# java heap
		import json
		print("HPCCONFIG", self.hpc_config_file, file=sys.stderr)
		data = json.load(open(self.hpc_config_file))
		try:
			qmem = str((int(int(data["qaa_qualimap"]["memory"]) / 1000)) - 2) + "G"
			self.config["qualimap_mem"] = qmem
		except:
			raise ValueError("Could not find qaa_qualimamp:memory-record in {}.".format(self.hpc_config_file))


	def __init__(self, args):
		self.args = copy(args)

		args.alt_hpc_config_warning = "Please provide a valid HPC configuration file with --hpc_config."
		args.alt_config_warning = "Please provide a valid configuration file with --config."

		super().__init__(args)
		self.__calc_qualimap_memory_setting()

		self.report_dir = join(args.output_dir, "reports")
		if not os.path.exists(self.report_dir):			
			pathlib.Path(self.report_dir).mkdir(parents=True, exist_ok=True)

		self.runmode = args.runmode if hasattr(args, "runmode") else "asm"

		# todo: change how this is done...
		if hasattr(args, "make_input_stream"):
			self.config["samplesheet"] = QAA_Runner.create_samplesheet_stream(args)
			self.config["has_stream"] = True
			self.config["full_samplesheet"] = args.input
		else:
			self.config["has_stream"] = False
			self.config["full_samplesheet"] = self.config["samplesheet"] = args.input

		self.config["out_dir"] = args.output_dir
		self.config["etc"] = ETC_DIR
		self.config["cwd"] = os.getcwd()
		self.config["create_bam"] = args.align_reads != "no"
		self.config["align_reads"] = args.align_reads if args.align_reads != "no" else False

		requested_modes = args.qaa_mode.split(",")
		modes = ["geno", "genome", "tran", "transcriptome", "prot", "proteome", "all"]
		invalid_modes = list(filter(lambda s:s not in modes, requested_modes))
		if invalid_modes:
			valid_modes = list(set(requested_modes).difference(invalid_modes))
			if valid_modes:
				print("--qaa-mode: Found invalid modes {}. Proceeding with {}.".format(invalid_modes, valid_modes))
				if "all" in valid_modes:
					requested_modes = modes[:-1]
				else:
					requested_modes = valid_modes
			else:
				print("--qaa-mode: No valid mode provided. Valid modes are {}. Exiting.".format(",".join(modes)))
				sys.exit(1)

		# these are coming in via api: #api
		self.config["normalized"] = args.normalized if hasattr(args, "normalized") else False # api
		try:
			self.config["project_prefix"] = args.project_prefix # api
		except:
			pass

		self.config["run_multiqc"] = args.run_multiqc if hasattr(args, "run_multiqc") else True # api
		self.config["multiqc_dir"] = args.multiqc_dir if hasattr(args, "multiqc_dir") else "."
		self.config["multiqc_config"] = args.multiqc_config if hasattr(args, "multiqc_config") else "."
		self.config["survey_assembly"] = self.runmode == "survey" # api
		self.config["run_blobtools"] = args.run_blobtools if hasattr(args, "run_blobtools") else True # api
		self.config["run_qualimap"] = args.run_qualimap if hasattr(args, "run_qualimap") else True
		self.config["run_busco"] = args.run_busco if hasattr(args, "run_busco") else True # api
		self.config["run_genome_module"] = self.runmode == "survey" or {"geno", "genome"}.intersection(requested_modes)
		self.config["run_transcriptome_module"] = self.runmode != "survey" and {"tran", "transcriptome"}.intersection(requested_modes)
		self.config["run_proteome_module"] = self.runmode != "survey" and {"prot", "proteome"}.intersection(requested_modes)
		self.config["quast_mincontiglen"] = args.quast_mincontiglen if hasattr(args, "quast_mincontiglen") else 0
		self.config["annotation"] = args.annotation if hasattr(args, "annotation") else None

		self.new_config_file = join(args.output_dir, "qaa.conf.yaml")
		with open(self.new_config_file, "w") as conf_out:
			yaml.dump(self.config, conf_out, default_flow_style=False)

		self.qaa_env = QAA_Environment(self.config)

		self.unlock = args.unlock
		self.output_dir = args.output_dir

	def __clean_blobtools_trash(self):
		import glob	
		print("QAA_CLEANUP:", os.getcwd())	
		for f in glob.glob(join(os.getcwd(), "*.bam.cov")):
			try:
				if os.path.isfile(f):
					os.unlink(f)
			except Exception as e:
				print(e)

	def run(self):
		run_result = run_snakemake(join(dirname(__file__), "zzz", "qaa.smk.py"), self.output_dir, self.new_config_file, self.exe_env, dryrun=False, unlock=self.unlock)
		print("QAA_RUN_RESULT=", run_result)
		self.report()
		self.__clean_blobtools_trash()
		return run_result

	def report(self):
		def report_func(qa_dir, report_out, rfunc):
			with open(report_out, "w") as rep_out:
				rfunc(qa_dir, out=rep_out)

		print("QAA_CONFIG_!!!")
		print(self.config)

		if self.runmode == "survey":
			report_func(self.qaa_env.quast_dir, join(self.report_dir, "quast_survey_report.tsv"), compileQUASTReport)
			if self.config["run_busco"]:
				report_func(self.qaa_env.busco_geno_dir, join(self.report_dir, "busco_survey_report.tsv"), compileBUSCOReport)
			if self.config["run_blobtools"]:
				report_func(self.qaa_env.blob_dir, join(self.report_dir, "blobtools_survey_report.tsv"), compileBlobReport)
		else:
			if "asm" in self.runmode:
				if self.config["run_genome_module"]:
					report_func(self.qaa_env.quast_dir, join(self.report_dir, "quast_report.tsv"), compileQUASTReport)
				if self.config["run_busco"]:
					report_func(self.qaa_env.busco_geno_dir, join(self.report_dir, "busco_genome_report.tsv"), compileBUSCOReport)
				if self.config["run_blobtools"]:
					report_func(self.qaa_env.blob_dir, join(self.report_dir, "blobtools_report.tsv"), compileBlobReport)
			if "ann" in self.runmode:
				if self.config["run_transcriptome_module"] or self.config["run_proteome_module"]:
					report_func(self.qaa_env.busco_dir, join(self.report_dir, "busco_report.tsv"), compileBUSCOReport)
				with open(join(self.report_dir, "16S_report.tsv"), "wt") as ostream:
					SixteenSReporter(self.qaa_env.prokka_dir, header=sixteenS_header, out=ostream).generateReport()
