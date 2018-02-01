#!/usr/bin/env python3

import os
from os.path import join, dirname, basename
import shutil
from argparse import ArgumentParser, RawDescriptionHelpFormatter as RDHF
from collections import defaultdict, namedtuple

import yaml
from snakemake.utils import min_version

from eicore import NOW
from eicore.external_process.snakemake_helper import *
from . import DEFAULT_CONFIG_FILE, DEFAULT_HPC_CONFIG_FILE, __version__

min_version("4.0")

# https://www.ascii-art.net/content.php?id=disney
kaa = [ 
"                                     ..:::::;'|",
"             __   ___             / \:::::;'  ;",
"           ,::::`'::,`.          :   ___     /",
"          :_ `,`.::::)|          | ,'SSt`.  /",
"          |(` :\)):::`;          |:::::::| :",
"          : \ ,`'`:::::`.        |:::::::| |",
"           \ \  ,' `:::::`.      :\::::::; |",
"            \ `.  ,' ` ,--.)     : `----'  |",
"             :  `-.._,'__.'      :   ____  |",
"             |     |              :,'::::\ |",
"             :  _  |              :::::::::|",
"             ;     :              |:::::::||",
"            :      |              |:::::::;|",
"            :  __  :              |\:::::; |",
"            |      |              | _____  |",
"            |      :              |':::::\ |",
"            : .--  |\             |::::::::|",
"            :      :.\            |:::::::;|",
"             :      \:\           |::::::/ |",
"              \ .-'  \'`.         ;`----' /;",
"               \      \|::-...__,'._     //",
"                `.  ,' `:::/ |::::::`. ,'/",
"                  `.     `-._;:::::::.','",
"                    `-..__,   ``--'' ,'",
"                          ``---....-'",
]

ASM_Sample = namedtuple("ASM_Sample", "id assembly bamfile r1 r2 busco_id".split(" "))

def readSamplesheet(_in):
	import csv
	for row in csv.reader(_in, delimiter=","):
		yield (row[0], ASM_Sample(*row))



def main():
	print("Starting QAA V" + __version__)
	print()
	print("\n".join(kaa))
	print()


	parser = ArgumentParser("The Earlham Institute Quality Assessment for Assemblies (QAA)",
	                        description="""...""",
				formatter_class=RDHF)

	parser.add_argument("input", help="""Path to assembly samplesheet.""")
	parser.add_argument("--output-dir", "-o", type=str, default=".", help="QAA will output data to this directory.")
	parser.add_argument("--blobtools-no-bwa", action="store_true", help="Use this switch to avoid aligning your reads against the reference, if you are providing your own bam files in the samplesheet. [off/False]")
	parser.add_argument("--qaa-mode", type=str, default="geno", help="Comma-separated list of run modes (geno[me], tran[scriptome], prot[eome]). If tran and/or prot are chosen, the samplesheet must include paths to transcriptome/proteome data. Omitting the genome mode allows to run QAA purely on transcriptome/proteome data. In this case, only BUSCO analyses will be performed. QAA'ing a transcriptome assembly can be done by supplying the respective fasta file in the assembly column of the samplesheet. [geno]")
	parser.add_argument("--survey-assembly", action="store_true", help="If True, then this qaa will only run quast and blobtools. This mode is for assessing assemble-ability of a set of reads as well as a closer look at taxonomic composition of the data. genome-mode is implied and prot/tran mode requests are ignored. [off/False]")

	parser.add_argument("--config", help="Configuration file for the pipeline. This file specifies details for accessing services and commands to be executed prior to running each pipeline tool.  Default config file is: " + DEFAULT_CONFIG_FILE)
	

	make_exeenv_arg_group(parser) #, hpc_config=DEFAULT_HPC_CONFIG_FILE)	# Add in cluster and DRMAA options

	args = parser.parse_args()
	print(args)
	
	print("Configuring execution environment ... ", end="", flush=True)
	logs_dir = join(args.output_dir, "hpc_logs")
	exe_env = ExecutionEnvironment(args, NOW, job_suffix=args.input + "_" + args.output_dir, log_dir=logs_dir)
	print("done.")
	print(str(exe_env))
	print()

	if args.config:
		print("Custom configuration file specified, overriding defaults ... ", end="", flush=True)
		config = yaml.load(open(args.config))
		config_file = args.config
		print("done.")
		print()
	else:
		print("Loading default configuration ... ", end="", flush=True)
		config = yaml.load(open(DEFAULT_CONFIG_FILE))
		config_file = DEFAULT_CONFIG_FILE
		print("done.")
		print()

	if os.path.exists(args.output_dir):
		print("Output directory already exists, attempting to resume...")
	else:
		print("Output director does not exist. Creating ... ", end="", flush=True)
		os.makedirs(args.output_dir)
		print("done.")
	
	if not os.path.exists(logs_dir) and exe_env.use_scheduler:
		print("HPC log dir does not exist. Creating ... ", end="", flush=True)
		os.makedirs(logs_dir)
		print("done.")
	
	print()
	
	print(config)
	config["samplesheet"] = args.input
	config["out_dir"] = args.output_dir
	config["etc"] = join(dirname(__file__), "..", "etc")
	config["cwd"] = os.getcwd()
	config["blobtools_run_bwa"] = not args.blobtools_no_bwa

	"""if args.contig_minlen:
        	config["use_asm_lengthfilter"] = True
	        config["asm_lengthfilter_contig_minlen"] = args.contig_minlen
	else:
        	config["use_asm_lengthfilter"] = False
	        config["asm_lengthfilter_contig_minlen"] = 0
	"""

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

	config["survey_assembly"] = args.survey_assembly
	config["run_genome_module"] = args.survey_assembly or ("geno" in requested_modes or "genome" in requested_modes)
	config["run_transcriptome_module"] = not args.survey_assembly and ("tran" in requested_modes or "transcriptome" in requested_modes)
	config["run_proteome_module"] = not args.survey_assembly and ("prot" in requested_modes or "proteome" in requested_modes)

	new_config_file = join(args.output_dir, "qaa.conf.xml")
	with open(new_config_file, "w") as conf_out:
		yaml.dump(config, conf_out, default_flow_style=False)
	
	run_result = run_snakemake(join(dirname(__file__), "zzz", "qaa.smk.py"), args.output_dir, new_config_file, exe_env, dryrun=False, unlock=args.unlock)
	
		

if __name__ == "__main__":
	main()
