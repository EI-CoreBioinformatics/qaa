import sys
import csv
import os
from os.path import join, basename, dirname

from qaa.samplesheet import readQAASamplesheet
from qaa.qaa_environment import QAA_Environment

DEBUG = config.get("debugmode", False)

TIME_CMD = config.get("tools", dict()).get("time", "/usr/bin/time -v")

qaa_env = QAA_Environment(config)
runmode = "survey" if config["survey_assembly"] else "asm"

# setup i/o
if type(config["samplesheet"]) is str:
	INPUTFILES = dict(readQAASamplesheet(csv.reader(open(config["samplesheet"]), delimiter=",")))
else:
	INPUTFILES = dict(readQAASamplesheet(config["samplesheet"]))

# generate target list
TARGETS = list()
for sample in INPUTFILES:
	if config["run_genome_module"]:
		TARGETS.append(join(qaa_env.quast_dir, sample, "transposed_report.tsv"))
		if config["run_blobtools"]:
			TARGETS.append(join(qaa_env.blob_dir, sample, sample + ".blobDB.table.txt"))
		if config["run_qualimap"]:
			TARGETS.append(join(qaa_env.qualimap_dir, sample, "qualimapReport.html"))
		if config["run_busco"]:
			TARGETS.append(join(qaa_env.busco_geno_dir, sample, sample + "_short_summary.txt"))
	if config["run_transcriptome_module"]:
		TARGETS.append(join(qaa_env.busco_tran_dir, sample, sample + "_short_summary.txt"))
	if config["run_proteome_module"]:
		TARGETS.append(join(qaa_env.busco_prot_dir, sample, sample + "_short_summary.txt"))

if config["run_multiqc"]:
	# this always needs to be last!
	TARGETS.append(join(config["multiqc_dir"], config["project_prefix"] + "_" + runmode + "_multiqc_report.html"))


if DEBUG:
	with open("qaa-inputfiles.txt", "w") as input_out:
		print(*INPUTFILES.values(), sep="\n", file=input_out)

	print("CONFIG", config, sep="\n")
	with open("qaa-targets.txt", "w") as targets_out:
		print(*TARGETS, sep="\n", file=targets_out)

# helpers
def getBUSCOData(sample):
	return join(qaa_env.busco_data_dir, INPUTFILES[sample].busco_id)
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

def get_cmd_call(cfg, container):
	call = cfg.get("singularity", dict()).get(container, "")
	return "singularity exec {0} ".format(call) if call and cfg.get("singularity", dict()).get("use_singularity", False) else ""

CMD_CALL = get_cmd_call(config, "qaa_container")


### RULES ###

localrules: all, qaa_compile_multiqc_inputs

rule all:
	input:
		TARGETS

if config["run_multiqc"]:
	rule qaa_compile_multiqc_inputs:
		input:
			TARGETS[:-1]
		output:
			join(config["multiqc_dir"], config["project_prefix"] + "_" + runmode + ".multiqc_input.txt")
		params:
			fastqcdir = join(qaa_env.qc_dir, "fastqc", "bbnorm" if config["normalized"] else "bbduk"),
			katdir = join(qaa_env.qc_dir, "kat"),
			buscodir = qaa_env.busco_geno_dir,
			quastdir = qaa_env.quast_dir,
			qualimapdir = qaa_env.qualimap_dir,
			samplesheet = config["full_samplesheet"],
		run:
			import os
			import sys
			import glob
			import csv
			input_files = list()
			try:
				input_files.extend(glob.glob(os.path.join(params.buscodir, "*", "*short_summary.txt")))
			except:
				print("Could not find busco output in " + params.buscodir, file=log)
				pass
			try:
				input_files.extend(glob.glob(os.path.join(params.quastdir, "*", "report.tsv")))
			except:
				print("Could not find quast output in " + params.quastdir, file=log)
				pass
			if runmode == "survey":
				try:
					input_files.extend(glob.glob(os.path.join(params.katdir, "*", "*.json")))					
				except:
					print("Could not find kat output in " + params.katdir, file=log)
					pass
				try:
					input_files.extend(glob.glob(os.path.join(params.fastqcdir, "*", "*", "fastqc_data.txt")))
				except:
					print("Could not find fastqc output in " + params.fastqcdir, file=log)			
					pass
			try:
				input_files.extend(glob.glob(os.path.join(params.qualimapdir, "*", "*.txt")))
			except:
				print("Could not find qualimap output in " + params.qualimapdir, file=log)			
				pass

			valid_samples = set(row[0] for row in csv.reader(open(params.samplesheet), delimiter=",") if row[0])

			# print("THESE ARE VALID SAMPLES", *valid_samples, sep="\n", file=sys.stderr)
			# print("THESE ARE INPUT FILES", *input_files, sep="\n", file=sys.stderr)
			with open(output[0], "w") as out:
				for f in input_files:
					if (os.path.basename(os.path.dirname(f)) in valid_samples or os.path.basename(os.path.dirname(os.path.dirname(f)))) and os.stat(f).st_size > 0:
						print(f, file=out)

	rule qaa_multiqc:
		input:
			rules.qaa_compile_multiqc_inputs.output[0]
		output:                                                                                           	
			join(config["multiqc_dir"], config["project_prefix"] + "_" + runmode + "_multiqc_report.html")
		params:
			mqc_config = config["multiqc_config"],
			datadir = qaa_env.output_dir,
			outdir = config["multiqc_dir"],
			prefix = config["project_prefix"] + "_" + runmode,
			mode = runmode,           			
			cmd = CMD_CALL + "multiqc"
		log:
			join(qaa_env.log_dir, runmode + "_readqc_multiqc.log")
		shell:
			"{params.cmd} -f -n {params.prefix}_multiqc_report -i {params.prefix}" + \
			" -z -c {params.mqc_config} -o {params.outdir}" + \
			" --file-list {input}" + \
			" &> {log}"
	"""
	rule qaa_multiqc_old:
		input:
			TARGETS[:-1]
		output:
			join(config["multiqc_dir"], config["project_prefix"] + "_" + runmode + "_multiqc_report.html")
		params:
			mqc_config = config["multiqc_config"],
			datadir = qaa_env.output_dir,
			outdir = config["multiqc_dir"],
			prefix = config["project_prefix"] + "_" + runmode,
			mqc_files = runmode + "_MQC_LIST.txt",
			fastqcdir = join(qaa_env.qc_dir, "fastqc", "bbnorm" if config["normalized"] else "bbduk"),
			katdir = join(qaa_env.qc_dir, "kat"),
			buscodir = qaa_env.busco_geno_dir,
			quastdir = qaa_env.quast_dir,
			qualimapdir = qaa_env.qualimap_dir,
			samplesheet = config["full_samplesheet"],
			mode = runmode,
			cmd = CMD_CALL + "multiqc"
		log:
			join(qaa_env.log_dir, runmode + "_readqc_multiqc.log")
		shell:
			" find {params.buscodir} -name '*short_summary.txt' > {params.mqc_files}.tmp" + \
			" && find {params.quastdir} -name 'report.tsv' >> {params.mqc_files}.tmp" + \
			" && if [[ -d \"{params.katdir}\" && {params.mode} == \"survey\" ]]; then" + \
			"   find {params.katdir} -name '*.json' >> {params.mqc_files}.tmp; fi" + \
			" && if [[ -d \"{params.fastqcdir}\" && {params.mode} == \"survey\" ]]; then" + \
			"   find {params.fastqcdir} -name 'fastqc_data.txt' >> {params.mqc_files}.tmp; fi" + \
			" && find {params.qualimapdir} -name '*.txt' >> {params.mqc_files}.tmp" + \
			" && grep -F -f <(cut -f 1 -d , {params.samplesheet}) {params.mqc_files}.tmp" + \
			" > {params.mqc_files}" + \
			" && rm {params.mqc_files}.tmp" + \
			" && {params.cmd} -f -n {params.prefix}_multiqc_report -i {params.prefix}" + \
			" -z -c {params.mqc_config} -o {params.outdir}" + \
			" --file-list {params.mqc_files}" + \
			" && rm {params.mqc_files}" + \
			" &> {log}"
	"""

BUSCO_CMD = "" + \
	"mkdir -p {params.outdir} && cd {params.outdir} && cd .. && " + \
	"export AUGUSTUS_CONFIG_PATH={params.configdir} && mkdir -p $(dirname {params.configdir}) && " + \
	"{params.cp_init} /opt/miniconda/config {params.configdir} && " + \
	"{params.cmd} -i {params.inputpath} -c {threads} -m {params.busco_mode} " + \
	"--force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample} && " + \
	"touch {params.outdir}/short_summary_{wildcards.sample}.txt &> {log} && " + \
	"mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/ && " + \
	"rm -rvf {params.outdir} {params.configdir} && " + \
	"(mv {params.final_outdir}/short_summary_{wildcards.sample}.txt " + \
	"{params.final_outdir}/{wildcards.sample}_short_summary.txt || " + \
	"touch {params.final_outdir}/{wildcards.sample}_short_summary.txt)"


if config["run_proteome_module"]:
	rule qaa_busco_prot:
		input:
			busco_input = getProteins
		output:
			join(qaa_env.busco_prot_dir, "{sample}", "{sample}_short_summary.txt")
		log:
			join(qaa_env.log_dir, "{sample}_busco_prot.log")
		params:
			outdir = lambda wildcards: join(qaa_env.busco_prot_dir, "run_" + wildcards.sample),
			final_outdir = lambda wildcards: join(qaa_env.busco_prot_dir, wildcards.sample),
			tmp = lambda wildcards: join(qaa_env.busco_prot_dir, "tmp", wildcards.sample),
			busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
			busco_mode = "prot",
			inputpath = lambda wildcards: os.path.abspath(getProteins(wildcards)),
			cp_init = CMD_CALL + "cp -r",
			cmd = CMD_CALL + "run_BUSCO.py",
			configdir = lambda wildcards: join(qaa_env.busco_prot_dir, "config", wildcards.sample)
		threads:
			8
		shell:
			BUSCO_CMD

if config["run_transcriptome_module"]:
	rule qaa_busco_tran:
		input:
			busco_input = getTranscripts
		output:
			join(qaa_env.busco_tran_dir, "{sample}", "{sample}_short_summary.txt")
		log:
			join(qaa_env.log_dir, "{sample}_busco_tran.log")
		params:
			outdir = lambda wildcards: join(qaa_env.busco_tran_dir, "run_" + wildcards.sample),
			final_outdir = lambda wildcards: join(qaa_env.busco_tran_dir, wildcards.sample),
			tmp = lambda wildcards: join(qaa_env.busco_tran_dir, "tmp", wildcards.sample),
			busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
			busco_mode = "tran",
			inputpath = lambda wildcards: os.path.abspath(getTranscripts(wildcards)),
			cp_init = CMD_CALL + "cp -r",
			cmd = CMD_CALL + "run_BUSCO.py",
			configdir = lambda wildcards: join(qaa_env.busco_tran_dir, "config", wildcards.sample)
		threads:
			8
		shell:
			BUSCO_CMD

if config["run_genome_module"]:
	if config["run_busco"]:
		rule qaa_busco_geno:
			input:
				busco_input = getAssembly
			output:
				join(qaa_env.busco_geno_dir, "{sample}", "{sample}_short_summary.txt")
			log:
				join(qaa_env.log_dir, "{sample}_busco_geno.log")
			params:
				outdir = lambda wildcards: join(qaa_env.busco_geno_dir, "run_" + wildcards.sample),
				final_outdir = lambda wildcards: join(qaa_env.busco_geno_dir, wildcards.sample),
				tmp = lambda wildcards: join(qaa_env.busco_geno_dir, "tmp", wildcards.sample),
				busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
				busco_mode = "geno",
				inputpath = lambda wildcards: os.path.abspath(getAssembly(wildcards)),
				cp_init = CMD_CALL + "cp -r",
				cmd = CMD_CALL + "run_BUSCO.py",
				configdir = lambda wildcards: join(qaa_env.busco_geno_dir, "config", wildcards.sample)
			threads:
				8
			shell:
				BUSCO_CMD

	rule qaa_quast:
		input:
			assembly = getAssembly
		output:
				join(qaa_env.quast_dir, "{sample}", "transposed_report.tsv")
		params:
				outdir = lambda wildcards: join(qaa_env.quast_dir, wildcards.sample),
				cmd = CMD_CALL + "quast.py",
				contiglen = config["quast_mincontiglen"]
		log:
				join(qaa_env.log_dir, "{sample}.asm_quast_assembly.log")
		threads:
				2
		shell:				
				" ({params.cmd} -o {params.outdir} -t {threads} -L -s {input.assembly} --min-contig {params.contiglen}" + \
				" || touch {params.outdir}/transposed_report.tsv {params.outdir}/report.tsv)" + \
				" && cut -f 1,2 {params.outdir}/report.tsv > {params.outdir}/report.tsv.12" + \
				" && mv {params.outdir}/report.tsv {params.outdir}/report.tsv.full" + \
				" && mv {params.outdir}/report.tsv.12 {params.outdir}/report.tsv" + \
				" 2> {log}"

	if config["align_reads"]:
		BAM_THREADS = 16
		if config["align_reads"] == "bowtie2":
			QAA_ALIGN_BUILD_INDEX = CMD_CALL + "bowtie2-build --threads {threads} {input.assembly} {params.ref}"
			QAA_ALIGN = CMD_CALL + \
				"bowtie2 --threads {params.align_threads} -x {params.ref} " + \
				"-1 {input.reads[0]} -2 {input.reads[1]} --rg-id {params.sample} " + \
				"--rg LB:{params.sample} --rg PL:illumina --rg SM:{params.sample} --rg PU:{params.sample}"
		else: # bwa
			QAA_ALIGN_BUILD_INDEX = CMD_CALL + "bwa index -p {params.ref} {input.assembly}"
			QAA_ALIGN = CMD_CALL + \
				"bwa mem -t {threads} " + \
				"-R '@RG\\tID:1\\tLB:{params.sample}\\tPL:illumina\\tSM:{params.sample}\\tPU:{params.sample}' " + \
				"{params.ref} {input.reads[0]} {input.reads[1]}"

		rule qaa_align_reads:
			input:
				reads = getReads,
				assembly = getAssembly
			output:
				bam = join(qaa_env.bam_dir, "{sample}", "{sample}.align_reads.bam")
			log:
				join(qaa_env.log_dir, "{sample}.align_reads.log")
			params:
				outdir = lambda wildcards: join(qaa_env.bam_dir, wildcards.sample),
				ref = lambda wildcards: join(qaa_env.bam_dir, wildcards.sample, wildcards.sample + ".assembly.fasta"),
				outbam = lambda wildcards: join(qaa_env.bam_dir, wildcards.sample, wildcards.sample, ".align_reads.bam"),
				align_threads = BAM_THREADS // 2,
				sort_threads = BAM_THREADS // 2,
				samtools_cmd = CMD_CALL + "samtools",
				markdup_cmd = CMD_CALL + "picard MarkDuplicates",
				sample = lambda wildcards: wildcards.sample
			threads:
				BAM_THREADS
			shell:
				TIME_CMD + " " + QAA_ALIGN_BUILD_INDEX + " &&" + \
				" touch {output.bam}.ref_index.done &&" + \
				" " + TIME_CMD + " " + QAA_ALIGN + " |" + \
				" {params.samtools_cmd} sort -@ {params.sort_threads} -o {output.bam}.tmp.bam - &&" + \
				" touch {output.bam}.sort.done &&" + \
				TIME_CMD + " {params.markdup_cmd}" + \
				" INPUT={output.bam}.tmp.bam OUTPUT={output.bam}" + \
				" METRICS_FILE={output.bam}.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate &&" + \
				" touch {output.bam}.markdup.done &&" + \
				" {params.samtools_cmd} index {output.bam} &&" + \
				" touch {output.bam}.bam_index.done &&" + \
				" rm {output.bam}.tmp.bam" + \
				" 2> {log}"

	if config["run_qualimap"]:
		rule qaa_qualimap:
			input:
				bam = join(qaa_env.bam_dir, "{sample}", "{sample}.align_reads.bam") if config["align_reads"] else getBAM,
			output:
				join(qaa_env.qualimap_dir, "{sample}", "qualimapReport.html")
			params:
				outdir = lambda wildcards: join(qaa_env.qualimap_dir, wildcards.sample),
				cmd = CMD_CALL + "qualimap",
				mem = config["qualimap_mem"]  # this is coming from within qaa/__init__.py
			log: 
				join(qaa_env.log_dir, "{sample}.qualimap.log")
			threads: 2
			message: "Using qualimap to collect stats for: {input.bam}"
			shell: 
			   TIME_CMD + " {params.cmd} --java-mem-size={params.mem} bamqc " + \
			   "-bam {input.bam} -nt {threads} -outdir {params.outdir} &> {log}"

		rule qaa_blastn:
			input:
				assembly = getAssembly
			output:
				tsv = join(qaa_env.blast_dir, "{sample}", "{sample}.blast.tsv")
			log:
				join(qaa_env.log_dir, "{sample}.blast.log")
			params:
				outdir = lambda wildcards: join(qaa_env.blast_dir, wildcards.sample),
				cmd = CMD_CALL + "blastn",
				blastdb = config["resources"]["blob_blastdb"]
			threads:
				8
			shell: 
				TIME_CMD + \
				" {params.cmd} -db {params.blastdb} -task megablast -query {input.assembly}" + \
				" -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle'" + \
				" -culling_limit 5 -num_threads {threads} -evalue 1e-25 -max_hsps 1 -max_target_seqs 999999 -perc_identity 75" + \
				" | sort -k1,1 -k3,3gr | awk -v FS=\"\\t\" -v OFS=\"\\t\" '{{ if (seen[$1] < 10) {{ seen[$1] += 1; print $0; }} }}'" + \
				" > {output.tsv} 2> {log}"

	if config["run_blobtools"]:
		rule qaa_blobtools:
			input:
				bam = join(qaa_env.bam_dir, "{sample}", "{sample}.align_reads.bam") if config["align_reads"] else getBAM,
				blast = join(qaa_env.blast_dir, "{sample}", "{sample}.blast.tsv"),
				assembly = getAssembly
			output:
				blobtable = join(qaa_env.blob_dir, "{sample}", "{sample}.blobDB.table.txt")
			log:
				join(qaa_env.log_dir, "{sample}.blobtools.log")
			params:
				prefix = lambda wildcards: join(qaa_env.blob_dir, wildcards.sample, wildcards.sample),
				taxlevel = "genus",
				cmd = CMD_CALL + "blobtools",
				min_contiglen = config["quast_mincontiglen"] # 1000
			threads:
				1
			shell:
				" {params.cmd} create -i {input.assembly} -b {input.bam} -t {input.blast} -o {params.prefix} &&" + \
				" {params.cmd} view -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/ -r {params.taxlevel} &&" + \
				" {params.cmd} plot -r {params.taxlevel} -l {params.min_contiglen} -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/" + \
				" &> {log}"				
