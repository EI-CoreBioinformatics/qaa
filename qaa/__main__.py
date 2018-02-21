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
from . import DEFAULT_CONFIG_FILE, DEFAULT_HPC_CONFIG_FILE, __version__, QAA_Runner

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
	parser.add_argument("--no-blobtools", action="store_true", help="Avoids failing blobtools on low quality assemblies. BGRRL helper function [False]")
	parser.add_argument("--quast-mincontiglen", type=int, default=1000, help="Minimum contig length [bp] for quast to consider [1000]")
	parser.add_argument("--config", help="Configuration file for the pipeline. This file specifies details for accessing services and commands to be executed prior to running each pipeline tool.  Default config file is: " + DEFAULT_CONFIG_FILE)
	

	make_exeenv_arg_group(parser) #, hpc_config=DEFAULT_HPC_CONFIG_FILE)	# Add in cluster and DRMAA options

	args = parser.parse_args()
	print(args)

	qaa_runner = QAA_Runner(args)
	run_result = qaa_runner.run()

	qaa_runner.report(os.path.join(args.output_dir, "reports"), os.path.join(args.output_dir, "qa", "asm", "busco"))

if __name__ == "__main__":
	main()
