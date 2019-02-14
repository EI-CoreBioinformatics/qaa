#!/usr/bin/env python3

import os
from os.path import join, dirname, basename
import shutil
from argparse import ArgumentParser, RawDescriptionHelpFormatter as RDHF
from collections import defaultdict, namedtuple

import yaml
from snakemake.utils import min_version

from .snakemake_helper import *
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



def main():
	print("Starting QAA V" + __version__)
	print()
	print("\n".join(kaa))
	print()


	parser = ArgumentParser(
		"The Earlham Institute Quality Assessment for Assemblies (QAA)",
		description="""...""",
		formatter_class=RDHF
	)

	parser.add_argument(
		"input", 
		help="""Path to assembly samplesheet."""
	)
	
	parser.add_argument(
		"--output-dir", "-o", 
		type=str, 
		default=".", 
		help="""QAA will output data to this directory."""
	)

	parser.add_argument(
		"--align-reads", 
		choices=("bwa", "bowtie2", "no"), 
		default="bowtie2", 
		help="""blobtools and qualimap require bam-input. You can either have qaa handle it with bwa or provide your own bam file. [bowtie2]"""
	)

	parser.add_argument(
		"--qaa-mode", 
		type=str, 
		default="all", 
		help="""Comma-separated list of run modes (geno[me], tran[scriptome], prot[eome], all). 
				If tran and/or prot are chosen, the samplesheet must include paths to transcriptome/proteome data. 
				Omitting the genome mode allows to run QAA purely on transcriptome/proteome data. In this case, only BUSCO analyses will be performed. 
				QAA'ing a transcriptome assembly can be done by supplying the respective fasta file in the assembly column of the samplesheet. [all]"""
	)

	parser.add_argument(
		"--multiqc-dir", 
		type=str, 
		default="."
	)

	parser.add_argument(
		"--quast-mincontiglen", 
		type=int, 
		default=1000, 
		help="""Minimum contig length [bp] for quast to consider [1000]"""
	)

	parser.add_argument(
		"--config", 
		help="""Configuration file for the pipeline. This file specifies details for accessing services and commands to be executed prior to running each pipeline tool. 
				Default config file is: """ + DEFAULT_CONFIG_FILE
	)
	
	make_exeenv_arg_group(parser) 


	args = parser.parse_args()

	qaa_runner = QAA_Runner(args)
	run_result = qaa_runner.run()
	qaa_runner.report() 

if __name__ == "__main__":
	main()
