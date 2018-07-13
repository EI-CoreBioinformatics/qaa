import sys
import csv
import os
import glob
from os.path import join, basename, dirname

from qaa import readSamplesheet, loadPreCmd, TIME_CMD

# setup I/O
if type(config["samplesheet"]) is str:
	INPUTFILES = dict(readSamplesheet(csv.reader(open(config["samplesheet"]), delimiter=",")))
else:
	INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("qaa-inputfiles.txt", "w") as input_out:
	print(*INPUTFILES.values(), sep="\n", file=input_out)
OUTPUTDIR = config["out_dir"]

QA_DIR = join(OUTPUTDIR, "qa", "survey" if config["survey_assembly"] else "asm")
QC_DIR = join(OUTPUTDIR, "qc")
LOG_DIR = join(QA_DIR, "log")

QUAST_DIR = join(QA_DIR, "quast")
BUSCO_DIR = join(QA_DIR, "busco")
BUSCO_GENO_DIR = join(BUSCO_DIR, "geno")
BUSCO_TRAN_DIR = join(BUSCO_DIR, "tran")
BUSCO_PROT_DIR = join(BUSCO_DIR, "prot")
BLOB_DIR = join(QA_DIR, "blobtools", "blob")
BLOB_BWA_DIR = join(QA_DIR, "blobtools", "bwa")
BLOB_BLAST_DIR = join(QA_DIR, "blobtools", "blast")


BBNORM_DIR = join(OUTPUTDIR, "qc", "bbnorm")
ASSEMBLY_DIR = join(OUTPUTDIR, "assembly")

# BUSCO environment
BUSCO_INIT_DIR = join(config["etc"], "util", "busco_init_dir")
BUSCO_DATA_PATH = config["resources"]["busco_databases"]
# BUSCO_DATA = "/tgac/workarea/group-pb/schudomc_bact/bact-grrl/data/busco/bacteria_odb9"
CWD = os.getcwd()

def getBUSCOData(sample):
	return join(BUSCO_DATA_PATH, INPUTFILES[sample].busco_id)
def getBAM(wildcards):
	return INPUTFILES[wildcards.sample].bamfile
def getReads(wildcards):
	return INPUTFILES[wildcards.sample].r1, INPUTFILES[wildcards.sample].r2
def getAssembly(wildcards):
	return INPUTFILES[wildcards.sample].assembly
def getTranscripts(wildcards):
	return INPUTFILES[wildcards.sample].transcripts
def getProteins(wildcards):
	return INPUTFILES[wildcards.sample].proteins



# ASM_Sample = namedtuple("ASM_Sample", "id assembly bamfile r1 r2 busco_id transcripts proteins".split(" "))
TARGETS = list()
# runbusco = not config["survey_assembly"]
for sample in INPUTFILES:
	if config["run_genome_module"]:
		TARGETS.append(join(QUAST_DIR, sample, "quast.log"))
		if not config["no_blobtools"]:
			TARGETS.append(join(BLOB_DIR, sample, sample + ".blobDB.table.txt"))
		if not config["no_busco"]:
			# TARGETS.append(join(BUSCO_GENO_DIR, sample, "short_summary_{}.txt".format(sample))) # if runbusco else "",)
			TARGETS.append(join(BUSCO_GENO_DIR, sample, sample + "_short_summary.txt"))
	if config["run_transcriptome_module"]:
		TARGETS.append(join(BUSCO_TRAN_DIR, sample, sample + "_short_summary.txt"))
	if config["run_proteome_module"]:
		TARGETS.append(join(BUSCO_PROT_DIR, sample, sample + "_short_summary.txt"))

if not config["no_multiqc"]:
	TARGETS.append(join(config["multiqc_dir"], config["misc"]["project"] + ("_survey" if config["survey_assembly"] else "_asm") + "_multiqc_report.html"))

TARGETS = list(filter(lambda t:t, TARGETS))

print("CONFIG")
print(config)

with open("qaa-targets.txt", "w") as targets_out:
	print(*TARGETS, sep="\n", file=targets_out)

localrules: all, qaa_multiqc

rule all:
	input:
		TARGETS

if not config["no_multiqc"]:
	rule qaa_multiqc:
		input: 
			TARGETS[:-1]
		output:
			join(config["multiqc_dir"], config["misc"]["project"] + ("_survey" if config["survey_assembly"] else "_asm") + "_multiqc_report.html")	
		params:
			load = loadPreCmd(config["load"]["multiqc"]),
			mqc_config = config["resources"]["multiqc_config"],
			datadir = OUTPUTDIR, # QC_DIR,
			outdir = config["multiqc_dir"], # join(OUTPUTDIR, "reports", "multiqc", "qc"),
			prefix = config["misc"]["project"] + ("_survey" if config["survey_assembly"] else "_asm"),
			ignore_qc = join(QC_DIR, "fastqc", "bbduk", "*"),
			ignore_qa = join(QA_DIR, "log"),
			mqc_files = "MQC_LIST.txt",
			fastqcdir = join(QC_DIR, "fastqc", "bbnorm" if config["normalized"] else "bbduk"),
			katdir = join(QC_DIR, "kat"),
			buscodir = join(QA_DIR, "busco", "geno"),
			quastdir = join(QA_DIR, "quast"),
			samplesheet = config["full_samplesheet"],
			mode = "survey" if config["survey_assembly"] else "asm"
		log:
			"readqc_multiqc.log"
		shell:
			"{params.load}" + \
			" find {params.buscodir} -name '*short_summary.txt' > {params.mqc_files}.tmp &&" + \
			" find {params.quastdir} -name 'report.tsv' >> {params.mqc_files}.tmp &&" + \
			" if [[ -d \"{params.katdir}\" && {params.mode} == \"survey\" ]]; then find {params.katdir} -name '*.json' >> {params.mqc_files}.tmp; fi &&" + \
			" if [[ -d \"{params.fastqcdir}\" && {params.mode} == \"survey\" ]]; then find {params.fastqcdir} -name 'fastqc_data.txt' >> {params.mqc_files}.tmp; fi &&" + \
			" grep -F -f <(cut -f 1 -d , {params.samplesheet}) {params.mqc_files}.tmp > {params.mqc_files} &&" + \
			" rm {params.mqc_files}.tmp &&" + \
			" multiqc -f -n {params.prefix}_multiqc_report -i {params.prefix} -z -c {params.mqc_config} -o {params.outdir} --file-list {params.mqc_files} > {log}"





if config["run_proteome_module"]:
	rule qa_busco_prot:
		input:
			busco_input = getProteins
		output:
			join(BUSCO_PROT_DIR, "{sample}", "{sample}_short_summary.txt")
		log:
			join(config["cwd"], LOG_DIR, "{sample}_busco_prot.log")
		params:
			# input = lambda wildcards: join(config["cwd"], wildcards.busco_input),
			outdir = lambda wildcards: join(config["cwd"], BUSCO_PROT_DIR, "run_" + wildcards.sample),
			final_outdir = lambda wildcards: join(config["cwd"], BUSCO_PROT_DIR, wildcards.sample),
			tmp = lambda wildcards: join(config["cwd"], BUSCO_PROT_DIR, "tmp", wildcards.sample),
			load = loadPreCmd(config["load"]["busco"]),
			busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
			busco_mode = "prot"
		threads:
			8
		shell:
			BUSCO_INIT_DIR + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
			" {params.load}" + TIME_CMD + \
			" run_BUSCO.py -i " + join(config["cwd"], "{input.busco_input}") + " -c {threads} -m {params.busco_mode}" + \
			" --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample} &> {log} && cd " + CWD + \
			" && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
			" && rm -rf {params.outdir}" + \
			" && mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt" + \
			"" # " &> {log}"





if config["run_transcriptome_module"]:
	rule qa_busco_tran:
		input:
			busco_input = getTranscripts
		output:
			join(BUSCO_TRAN_DIR, "{sample}", "{sample}_short_summary.txt")
		log:
			join(config["cwd"], LOG_DIR, "{sample}_busco_tran.log")
		params:
			# input = lambda wildcards: join(config["cwd"], wildcards.busco_input),
			outdir = lambda wildcards: join(config["cwd"], BUSCO_TRAN_DIR, "run_" + wildcards.sample),
			final_outdir = lambda wildcards: join(config["cwd"], BUSCO_TRAN_DIR, wildcards.sample),
			tmp = lambda wildcards: join(config["cwd"], BUSCO_TRAN_DIR, "tmp", wildcards.sample),
			load = loadPreCmd(config["load"]["busco"]),
			busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
			busco_mode = "tran"
		threads:
			8
		shell:
			BUSCO_INIT_DIR + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
			" {params.load}" + TIME_CMD + \
			" run_BUSCO.py -i " + join(config["cwd"], "{input.busco_input}") + " -c {threads} -m {params.busco_mode}" + \
			" --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample} &> {log} && cd " + CWD + \
			" && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
			" && rm -rf {params.outdir}" + \
			" && mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt" + \
			"" # " &> {log}"
		
if config["run_genome_module"]:
	if not config["no_busco"]: # True: # not config["survey_assembly"]:
		rule qa_busco_geno:
			input:
				busco_input = getAssembly
			output:
				# join(BUSCO_GENO_DIR, "{sample}", "short_summary_{sample}.txt")
				join(BUSCO_GENO_DIR, "{sample}", "{sample}_short_summary.txt")
			log:	
				join(config["cwd"], LOG_DIR, "{sample}_busco_geno.log")
			params:
				# input = lambda wildcards: join(config["cwd"], wildcards.input.busco_input),
				outdir = lambda wildcards: join(config["cwd"], BUSCO_GENO_DIR, "run_" + wildcards.sample),
				final_outdir = lambda wildcards: join(config["cwd"], BUSCO_GENO_DIR, wildcards.sample),
				tmp = lambda wildcards: join(config["cwd"], BUSCO_GENO_DIR,  "tmp", wildcards.sample),
				load = loadPreCmd(config["load"]["busco"]),
				busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
				busco_mode = "geno"
			threads:
				8
			shell:
				BUSCO_INIT_DIR + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
				" {params.load}" + TIME_CMD + \
				" run_BUSCO.py -i " + join(config["cwd"], "{input.busco_input}") + " -c {threads} -m {params.busco_mode}" + \
				" --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample} &> {log} && cd " + CWD + \
				" && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
				" && rm -rf {params.outdir} && touch {params.final_outdir}/short_summary_{wildcards.sample}.txt" + \
				" && mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt" + \
				"" # &> {log}"

	rule qa_quast:
		input:
			assembly = getAssembly
		output:
			join(QUAST_DIR, "{sample}", "quast.log")
		params:
			load = loadPreCmd(config["load"]["quast"]),
			outdir = lambda wildcards: join(QUAST_DIR, wildcards.sample),
			cmd = loadPreCmd(config["cmd"]["quast"], is_dependency=False),
			contiglen = config["quast_mincontiglen"]
		log:
			join(LOG_DIR, "{sample}.asm_quast_assembly.log")
		threads:
			2
		shell:
			"{params.load} (" + TIME_CMD + \
			" {params.cmd} -o {params.outdir} -t {threads} -L -s {input.assembly} --min-contig {params.contiglen}" + \
			" || touch {params.outdir}/transposed_report.tsv) &> {log}"

	if not config["no_blobtools"] and config["blobtools_run_bwa"]:
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
				" bwa mem -t {threads} {params.index} {input.reads[0]} {input.reads[1]} | samtools view -buSh - |" + \
				" samtools sort -o {output.bam} -" + \
				" 2> {log}"
	if not config["no_blobtools"]:
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
				load = loadPreCmd(config["load"]["blobtools"]),
                                taxlevel = "family"
			threads:
				1
			shell:
				"{params.load}" + TIME_CMD + \
				" blobtools create -t {input.blast} -b {input.bwa} -i {input.assembly} -o {params.prefix} &&" + \
				" blobtools view -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/ -r {params.taxlevel} &&" + \
				" blobtools blobplot -r {params.taxlevel} -l 1000 -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/ &> {log}"

