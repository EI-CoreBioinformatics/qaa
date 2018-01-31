import sys
import csv
import os
import glob
from os.path import join, basename, dirname

from qaa import readSamplesheet, loadPreCmd, TIME_CMD

# setup I/O
INPUTFILES = dict(readSamplesheet(open(config["samplesheet"])))
with open("qaa-inputfiles.txt", "w") as input_out:
	print(*INPUTFILES.values(), sep="\n", file=input_out)
OUTPUTDIR = config["out_dir"]
QA_DIR = join(OUTPUTDIR, "qa")
LOG_DIR = join(QA_DIR, "log")

QUAST_DIR = join(QA_DIR, "quast")
BUSCO_DIR = join(QA_DIR, "busco")
BLOB_DIR = join(QA_DIR, "blobtools", "blob")
BLOB_BWA_DIR = join(QA_DIR, "blobtools", "bwa")
BLOB_BLAST_DIR = join(QA_DIR, "blobtools", "blast")


BBNORM_DIR = join(OUTPUTDIR, "qc", "bbnorm")
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")

# BUSCO environment
BUSCO_INIT_DIR = join(config["etc"], "util", "busco_init_dir")
BUSCO_DATA_PATH = config["resources"]["busco_databases"]
BUSCO_DATA = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl/data/busco/bacteria_odb9"
CWD = os.getcwd()

def getBUSCOData(sample):
	return join(BUSCO_DATA_PATH, INPUTFILES[sample].busco_id)
def getBAM(wildcards):
	return INPUTFILES[wildcards.sample].bamfile
def getReads(wildcards):
	return INPUTFILES[wildcards.sample].r1, INPUTFILES[wildcards.sample].r2
def getAssembly(wildcards):
	return INPUTFILES[wildcards.sample].assembly


# ASM_Sample = namedtuple("ASM_Sample", "id assembly bamfile r1 r2 busco_id".split(" "))
TARGETS = list()
for sample in INPUTFILES:
	TARGETS.extend([join(QUAST_DIR, sample, "quast.log"),
			join(BUSCO_DIR, sample, "short_summary_{}.txt".format(sample)),
			# join(BLOB_BWA_DIR, sample, sample + ".blob_bwa.bam")
			join(BLOB_DIR, sample, sample + ".blobDB.table.txt")
			])

print("CONFIG")
print(config)

with open("qaa-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

localrules: all

rule all:
	input: TARGETS


rule qa_busco_geno:
	input:
		assembly = getAssembly
		
	output:
		join(BUSCO_DIR, "{sample}", "short_summary_{sample}.txt")
	log:	
		join(config["cwd"], LOG_DIR, "{sample}_busco_geno.log")
	params:
		outdir = lambda wildcards: join(config["cwd"], BUSCO_DIR, "run_" + wildcards.sample),
		final_outdir = lambda wildcards: join(config["cwd"], BUSCO_DIR, wildcards.sample),
		tmp = lambda wildcards: join(config["cwd"], BUSCO_DIR,  "tmp", wildcards.sample),
		load = loadPreCmd(config["load"]["busco"]),
		busco_data = lambda wildcards: getBUSCOData(wildcards.sample)
	threads:
		8
	shell:
		BUSCO_INIT_DIR + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
		" {params.load}" + TIME_CMD + \
		" run_BUSCO.py -i {input.assembly} -c {threads} -m geno" + \
		" --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample} &> {log} && cd " + CWD + \
		" && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
		" && rm -rf {params.outdir}" + \
		" &> {log}"

rule qa_quast:
	input:
		assembly = getAssembly
	output:
		join(QUAST_DIR, "{sample}", "quast.log")
	params:
		load = loadPreCmd(config["load"]["quast"]),
		outdir = lambda wildcards: join(QUAST_DIR, wildcards.sample),
		cmd = loadPreCmd(config["cmd"]["quast"], is_dependency=False)
	log:
		join(LOG_DIR, "{sample}.asm_quast_assembly.log")
	threads:
		2
	shell:
		"{params.load}" + TIME_CMD + \
		" {params.cmd} -o {params.outdir} -t {threads} -L -s {input.assembly} --min-contig 1000 &> {log}"

if config["blobtools_run_bwa"]:
	rule qa_blob_bwa_mem:
		input:
			reads = getReads,
			assembly = getAssembly
		output:
			bam = join(BLOB_BWA_DIR, "{sample}", "{sample}.blob_bwa.bam")
		log:
			join(LOG_DIR, "{sample}.blob_bwa.log")
		params:
			outdir = lambda wildcards: join(BLOB_BWA_DIR, wildcards.sample),
			index = lambda wildcards: join(BLOB_BWA_DIR, wildcards.sample, wildcards.sample + ".assembly.fasta"),
			outbam = lambda wildcards: join(BLOB_BWA_DIR, wildcards.sample, wildcards.sample + ".blob_bwa.bam"),
			load = loadPreCmd(config["load"]["blob_bwa"])
		threads:
			8
		shell:
			"{params.load}" + TIME_CMD + \
			" bwa index -p {params.index} {input.assembly} &&" + \
			" bwa mem -t {threads} {params.index} {input.reads[0]} {input.reads[1]} | samtools view -buSh - > {output.bam} 2> {log}"

rule qa_blob_blast:
	input:
		assembly = getAssembly
	output:
		tsv = join(BLOB_BLAST_DIR, "{sample}", "{sample}.blob_blast.tsv")
	log:
		join(LOG_DIR, "{sample}.blob_blast.log")
	params:
		outdir = lambda wildcards: join(BLOB_BLAST_DIR, wildcards.sample),
		load = loadPreCmd(config["load"]["blob_blast"]),
		blastdb = config["resources"]["blob_blastdb"]
	threads:
		8
	shell:
		"{params.load}" + TIME_CMD + \
		" blastn -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 " + \
		" -num_threads {threads} -query {input.assembly} -db {params.blastdb} -out {output.tsv} &> {log}"

rule qa_blobtools:
	input:
		bwa = join(BLOB_BWA_DIR, "{sample}", "{sample}.blob_bwa.bam") if config["blobtools_run_bwa"] else getBAM,
		blast = join(BLOB_BLAST_DIR, "{sample}", "{sample}.blob_blast.tsv"),
		assembly = getAssembly
	output:
		blobtable = join(BLOB_DIR, "{sample}", "{sample}.blobDB.table.txt")
	log:
		join(LOG_DIR, "{sample}.blobtools.log")
	params:
		prefix = lambda wildcards: join(BLOB_DIR, wildcards.sample, wildcards.sample),
		load = loadPreCmd(config["load"]["blobtools"])
	threads:
		1
	shell:
		"{params.load}" + TIME_CMD + \
		" blobtools create -t {input.blast} -b {input.bwa} -i {input.assembly} -o {params.prefix} -x bestsumorder &&" + \
		" blobtools view -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/ -x bestsumorder -r species &&" + \
		" blobtools blobplot -r species -l 1000 -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/ &> {log}"

