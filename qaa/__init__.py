import pkg_resources
import os
from os.path import join, dirname, basename
import csv
from collections import namedtuple
import yaml

from eicore.external_process.snakemake_helper import *
from qaa.reporting.busco_report import compileBUSCOReport
from qaa.reporting.quast_report import compileQUASTReport
from qaa.reporting.blobtools_report import compileBlobReport

__title__ = "qaa"
__author__ = 'Christian Schudoma (cschu)'
__license__ = 'MIT'
__copyright__ = 'Copyright 2018 Earlham Institute'
__version__ = pkg_resources.require("qaa")[0].version


ETC_DIR = os.path.join(os.path.dirname(__file__), "..", "etc")
DEFAULT_HPC_CONFIG_FILE = os.path.join(ETC_DIR, "hpc_config.json")
DEFAULT_CONFIG_FILE = os.path.join(ETC_DIR, "qaa_config.yaml")

QAA_SAMPLESHEET_COLUMNS = "id assembly bamfile r1 r2 busco_id transcripts proteins".split(" ")
QAA_Sample = namedtuple("QAA_Sample", QAA_SAMPLESHEET_COLUMNS)

TIME_CMD = " /usr/bin/time -v"
QAA_ID = "XXX"


def readSamplesheet(_in):
	import csv
	# for row in csv.reader(_in, delimiter=","):
	for row in _in:
		if row and row[0]:
			while len(row) != len(QAA_SAMPLESHEET_COLUMNS):
				row.append("")
			yield (row[0], QAA_Sample(*row))

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


class QAA_Runner(object):
	def __init__(self, args, **kwargs):
		def _create_input_stream(args):
			stream, asm_path = list(), ""
			for row in csv.reader(open(args.input), delimiter=","):
				if args.runmode == "asm":
					asm_path = join(args.output_dir, "assembly", row[0], row[0] + ".assembly.fasta")
				elif args.runmode == "survey":
					asm_path = join(args.output_dir, "qc", "tadpole", row[0], row[0] + "_tadpole_contigs.fasta")
				new_row = [row[0], asm_path, "", row[2], row[3], args.busco_db]
				if args.runmode == "ann":
					ann = "prokka" if args.annotation in ("prokka", "both") else "ratt"
					new_row.extend([join(args.output_dir, "annotation", ann, row[0], row[0] + ".ffn"), join(args.output_dir, "annotation", ann, row[0], row[0] + ".faa")])
				stream.append(new_row)
			return stream


		print("Assimilating kwargs ... ", end="", flush=True)
		from copy import copy
		args = copy(args)
		for k in kwargs:
			setattr(args, k, kwargs[k])
		print(" done.")
		print()

		print("Configuring execution environment ... ", end="", flush=True)
		self.output_dir = args.output_dir
		self.report_dir = os.path.join(self.output_dir, "reports")
		if not os.path.exists(self.report_dir):
			os.makedirs(self.report_dir)
		self.logs_dir = os.path.join(self.output_dir, "hpc_logs")
		self.exe_env = ExecutionEnvironment(args, NOW, job_suffix=args.input + "_" + self.output_dir, log_dir=self.logs_dir)
		print("done.")
		print(str(self.exe_env))
		print()
		self.runmode = args.runmode

		if args.config:
			print("Custom configuration file specified, overriding defaults ... ", end="", flush=True)
			self.config = yaml.load(open(args.config))
			self.config_file = args.config
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
		
		if "make_input_stream" in args:
			self.config["samplesheet"] = _create_input_stream(args)
			self.config["has_stream"] = True
		else:
			self.config["has_stream"] = False
			self.config["samplesheet"] = args.input

		self.config["out_dir"] = self.output_dir
		self.config["etc"] = os.path.join(os.path.dirname(__file__), "..", "etc")
		self.config["cwd"] = os.getcwd()

		self.config["blobtools_run_bwa"] = not args.blobtools_no_bwa
		requested_modes = args.qaa_mode.split(",")
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

		self.config["project_prefix"] = args.project_prefix if "project_prefix" in args else "dummy_project"
		self.config["no_multiqc"] = args.no_multiqc if "no_multiqc" in args else False
		self.config["multiqc_dir"] = args.multiqc_dir if "multiqc_dir" in args else "."
		self.config["survey_assembly"] = args.survey_assembly
		self.config["no_blobtools"] = args.no_blobtools if "no_blobtools" in args else False
		self.config["no_busco"] = args.no_busco if "no_busco" in args else False
		self.config["run_genome_module"] = args.survey_assembly or ("geno" in requested_modes or "genome" in requested_modes)
		self.config["run_transcriptome_module"] = not args.survey_assembly and ("tran" in requested_modes or "transcriptome" in requested_modes)
		self.config["run_proteome_module"] = not args.survey_assembly and ("prot" in requested_modes or "proteome" in requested_modes)
		self.config["quast_mincontiglen"] = args.quast_mincontiglen if "quast_mincontiglen" in args else 0
		self.new_config_file = os.path.join(self.output_dir, "qaa.conf.xml")
		with open(self.new_config_file, "w") as conf_out:
			yaml.dump(self.config, conf_out, default_flow_style=False)

		self.unlock = args.unlock

	def _clean_blobtools_trash(self):
		import glob		
		for f in glob.glob(os.path.join(os.getcwd(), "*.blob_bwa.bam.cov")):
			try:
				if os.path.isfile(f):
					os.unlink(f)
			except Exception as e:
				print(e)
				

	def run(self):
		# return run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "qaa.smk.py"), self.output_dir, self.new_config_file, self.exe_env, dryrun=False, unlock=self.unlock)
		"""
		# was worth a try...
		print("RRRRRRRUN...")
		from snakemake import snakemake
		cluster_cfg = self.exe_env.hpc_config if self.exe_env.use_scheduler else None
		cluster = (self.exe_env.sub_cmd + self.exe_env.res_cmd if not self.exe_env.use_drmaa else None) if self.exe_env.use_scheduler else None
		drmaa = self.exe_env.res_cmd if self.exe_env.use_drmaa else None
		snakefile = os.path.join(os.path.dirname(__file__), "zzz", "qaa.smk.py")
		
		return snakemake(snakefile, cores=self.exe_env.max_cores, local_cores=self.exe_env.max_cores, nodes=self.exe_env.max_nodes, config=self.config, workdir=".", cluster_config=cluster_cfg, cluster=cluster, drmaa=drmaa, unlock=self.unlock, printshellcmds=True, printreason=True, stats=os.path.join(self.output_dir, os.path.basename(snakefile) + "-" + NOW + ".stats"), jobname="qaa.{rulename}.{jobid}", force_incomplete=True, latency_wait=60 if self.exe_env.use_scheduler or self.exe_env.use_drmaa else 1, printdag=False, dryrun=False, forceall=False, verbose=True)
		"""

		run_result = run_snakemake(os.path.join(os.path.dirname(__file__), "zzz", "qaa.smk.py"), self.output_dir, self.new_config_file, self.exe_env, dryrun=False, unlock=self.unlock)

		self.report()

		self._clean_blobtools_trash()
		

		return run_result

	def report(self):
		def report_func(qa_dir, report_out, rfunc):
			with open(report_out, "w") as rep_out:
				rfunc(qa_dir, out=rep_out)
			
		# if self.runmode == "survey":
		if self.config["survey_assembly"]:
			report_func(os.path.join(self.output_dir, "qa", "survey", "quast"), os.path.join(self.report_dir, "quast_survey_report.tsv"), compileQUASTReport)
			if not self.config["no_busco"]:
				report_func(os.path.join(self.output_dir, "qa", "survey", "busco", "geno"), os.path.join(self.report_dir, "busco_survey_report.tsv"), compileBUSCOReport)
			if not self.config["no_blobtools"]:
				report_func(os.path.join(self.output_dir, "qa", "survey", "blobtools", "blob"), os.path.join(self.report_dir, "blobtools_survey_report.tsv"), compileBlobReport)
		# elif self.runmode == "asm":
		if self.config["run_genome_module"]:
			report_func(os.path.join(self.output_dir, "qa", "asm", "quast"), os.path.join(self.report_dir, "quast_report.tsv"), compileQUASTReport)
			if not self.config["no_busco"]:
				report_func(os.path.join(self.output_dir, "qa", "asm", "busco", "geno"), os.path.join(self.report_dir, "busco_genome_report.tsv"), compileBUSCOReport)
			if not self.config["no_blobtools"]:
				report_func(os.path.join(self.output_dir, "qa", "asm", "blobtools", "blob"), os.path.join(self.report_dir, "blobtools_report.tsv"), compileBlobReport)
				

		# elif self.runmode == "ann":
		if self.config["run_transcriptome_module"] or self.config["run_proteome_module"]:
			report_func(os.path.join(self.output_dir, "qa", "asm", "busco"), os.path.join(self.report_dir, "busco_report.tsv"), compileBUSCOReport)
		
	

class QAA_ArgumentsAdapter(object):
	def __init__(self):
		pass
