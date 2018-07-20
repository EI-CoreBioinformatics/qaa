import sys
import csv
import os
from os.path import join, basename, dirname

from qaa import readSamplesheet, loadPreCmd, TIME_CMD
from qaa import QAA_Environment

# setup I/O
if type(config["samplesheet"]) is str:
        INPUTFILES = dict(readSamplesheet(csv.reader(open(config["samplesheet"]), delimiter=",")))
else:
        INPUTFILES = dict(readSamplesheet(config["samplesheet"]))
with open("qaa-inputfiles.txt", "w") as input_out:
        print(*INPUTFILES.values(), sep="\n", file=input_out)


qaa_env = QAA_Environment(config)
"""
class QAA_Environment(object):
    def __init__(self, config):
        self.output_dir = config.get("out_dir", ".")
        self.qa_dir = join(self.outputdir, "qa", "survey" if config["survey_assembly"] else "asm")
        self.qc_dir = join(self.outputdir, "qc")
        self.log_dir = join(self.qa_dir, "log")

        self.quast_dir = join(self.qa_dir, "quast")
        self.busco_dir = join(self.qa_dir, "busco")
        self.busco_geno_dir = join(self.busco_dir, "geno")
        self.busco_tran_dir = join(self.busco_dir, "tran")
        self.busco_prot_dir = join(self.busco_dir, "prot")
        # busco environment
        self.busco_init_dir = join(config["etc"], "util", "busco_init_dir")
        self.busco_data_dir = config["resources"]["busco_databases"]

        self.blob_dir = join(self.qa_dir, "blobtools")

        self.blast_dir = join(self.qa_dir, "blast")
        self.bam_dir = join(self.qa_dir, "bam")

        self.cwd = os.getcwd()
"""

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


runmode = "survey" if config["survey_assembly"] else "asm"


TARGETS = list()
for sample in INPUTFILES:
    if config["run_genome_module"]:
        TARGETS.append(join(qaa_env.quast_dir, sample, "quast.log"))
        if config["run_blobtools"]:
            TARGETS.append(join(qaa_env.blob_dir, sample, sample + ".blobDB.table.txt"))
        if config["run_busco"]:
            TARGETS.append(join(qaa_env.busco_geno_dir, sample, sample + "_short_summary.txt"))
    if config["run_transcriptome_module"]:
        TARGETS.append(join(qaa_env.busco_tran_dir, sample, sample + "_short_summary.txt"))
    if config["run_proteome_module"]:
        TARGETS.append(join(qaa_env.busco_prot_dir, sample, sample + "_short_summary.txt"))

if config["run_multiqc"]:
    TARGETS.append(join(config["multiqc_dir", config["misc"]["project"] + "_" + runmode + "_multiqc_report.html"))

TARGETS = list(filter(lambda t:t, TARGETS)) # why? 2018-07-20 I don't think this is needed anymore.

print("CONFIG", config, sep="\n")
with open("qaa-targets.txt", "w") as targets_out:
    print(*TARGETS, sep="\n", file=targets_out)


localrules: all, qaa_multiqc


rule all:
    input:
        TARGETS


if config["run_multiqc"]:
    rule qaa_multiqc:
        input = TARGETS[:-1] # TARGETS[-1] is multiqc report
    output:
        join(config["multiqc_dir"], config["misc"]["project"] + "_" + runmode + "_multiqc_report.html")
    params:
        load = loadPreCmd(config["load"]["multiqc"]),
        mqc_config = config["resources"]["multiqc_config"],
        datadir = qaa_env.output_dir,
        outdir = config["multiqc_dir"],
        prefix = config["misc"]["project"] + "_" + runmode,
        ignore_qc = join(qaa_env.qc_dir, "fastqc", "bbduk", "*"),
        ignore_qa = join(qaa_env.qa_dir, "log"),
        mqc_files = "MQC_LIST.txt",
        fastqcdir = join(qaa_env.qc_dir, "fastqc", "bbnorm" if config["normalized"] else "bbduk"),
        katdir = join(qaa_env.qc_dir, "kat"),
        buscodir = qaa_env.busco_geno_dir,
        quastdir = qaa_env.quast_dir,
        samplesheet = config["full_samplesheet"],
        mode = runmode
    log:
        "readqc_multiqc.log"
    shell:
        "{params.load}" + \
        " find {params.buscodir} -name '*short_summary.txt' > {params.mqc_files}.tmp &&" + \
        " find {params.quastdir} -name 'report.tsv' >> {params.mqc_files}.tmp &&" + \
        " if [[ -d \"{params.katdir}\" && {params.mode} == \"survey\" ]]; then" + \
        "   find {params.katdir} -name '*.json' >> {params.mqc_files}.tmp; fi &&" + \
        " if [[ -d \"{params.fastqcdir}\" && {params.mode} == \"survey\" ]]; then" + \
        "   find {params.fastqcdir} -name 'fastqc_data.txt' >> {params.mqc_files}.tmp; fi &&" + \
        " grep -F -f <(cut -f 1 -d , {params.samplesheet}) {params.mqc_files}.tmp > {params.mqc_files} &&" + \
        " rm {params.mqc_files}.tmp &&" + \
        " multiqc -f -n {params.prefix}_multiqc_report -i {params.prefix} -z -c {params.mqc_config} -o {params.outdir} " + \
        " --file-list {params.mqc_files} > {log}"

if config["run_proteome_module"]:
    rule qa_busco_prot:
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
            load = loadPreCmd(config["load"]["busco"]),
            busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
            busco_mode = "prot"
        threads:
            8
        shell:
            qaa_env.busco_init_dir + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
            " {params.load}" + TIME_CMD + \
            " run_BUSCO.py -i " + join(qaa_env.cwd, "{input.busco_input}") + " -c {threads} -m {params.busco_mode}" + \
            " --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample} &> {log} && cd " + qaa_env.cwd + \
            " && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
            " && rm -rf {params.outdir}" + \
            " && mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt"


if config["run_transcriptome_module"]:
    rule qa_busco_tran:
        input:
            busco_input = getProteins
        output:
            join(qaa_env.busco_tran_dir, "{sample}", "{sample}_short_summary.txt")
        log:
            join(qaa_env.log_dir, "{sample}_busco_tran.log")
        params:
            outdir = lambda wildcards: join(qaa_env.busco_tran_dir, "run_" + wildcards.sample),
            final_outdir = lambda wildcards: join(qaa_env.busco_tran_dir, wildcards.sample),
            tmp = lambda wildcards: join(qaa_env.busco_tran_dir, "tmp", wildcards.sample),
            load = loadPreCmd(config["load"]["busco"]),
            busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
            busco_mode = "tran"
        threads:
            8
        shell:
            qaa_env.busco_init_dir + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
            " {params.load}" + TIME_CMD + \
            " run_BUSCO.py -i " + join(qaa_env.cwd, "{input.busco_input}") + " -c {threads} -m {params.busco_mode}" + \
            " --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample} &> {log} && cd " + qaa_env.cwd + \
            " && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
            " && rm -rf {params.outdir}" + \
            " && mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt"


if config["run_genome_module"]:
    if config["run_busco"]:
        rule qa_busco_tran:
            input:
                busco_input = getProteins
            output:
                join(qaa_env.busco_tran_dir, "{sample}", "{sample}_short_summary.txt")
            log:
                join(qaa_env.log_dir, "{sample}_busco_tran.log")
            params:
                outdir = lambda wildcards: join(qaa_env.busco_tran_dir, "run_" + wildcards.sample),
                final_outdir = lambda wildcards: join(qaa_env.busco_tran_dir, wildcards.sample),
                tmp = lambda wildcards: join(qaa_env.busco_tran_dir, "tmp", wildcards.sample),
                load = loadPreCmd(config["load"]["busco"]),
                busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
                busco_mode = "tran"
            threads:
                8
            shell:
                qaa_env.busco_init_dir + " {params.outdir} && cd {params.outdir} && cd .. &&" + \
                " {params.load}" + TIME_CMD + \
                " run_BUSCO.py -i " + join(qaa_env.cwd, "{input.busco_input}") + " -c {threads} -m {params.busco_mode}" + \
                " --force -t {params.tmp} -l {params.busco_data} -o {wildcards.sample} &> {log} && cd " + qaa_env.cwd + \
                " && mkdir -p {params.final_outdir} && mv -v {params.outdir}/* {params.final_outdir}/" + \
                " && rm -rf {params.outdir}" + \
                " && mv {params.final_outdir}/short_summary_{wildcards.sample}.txt {params.final_outdir}/{wildcards.sample}_short_summary.txt"

    rule qa_quast:
        input:
            assembly = getAssembly
        output:
                join(qaa_env.quast_dir, "{sample}", "quast.log")
        params:
                load = loadPreCmd(config["load"]["quast"]),
                outdir = lambda wildcards: join(qaa_env.quast_dir, wildcards.sample),
                cmd = loadPreCmd(config["cmd"]["quast"], is_dependency=False),
                contiglen = config["quast_mincontiglen"]
        log:
                join(qaa_env.log_dir, "{sample}.asm_quast_assembly.log")
        threads:
                2
        shell:
                "{params.load} (" + TIME_CMD + \
                " {params.cmd} -o {params.outdir} -t {threads} -L -s {input.assembly} --min-contig {params.contiglen}" + \
                " || touch {params.outdir}/transposed_report.tsv) &> {log}"

    if config["create_bam"]:

        """
        source bowtie-2.3.4 && /usr/bin/time -v bowtie2-build --threads 16 genome.1kb.fasta genome.1kb.fasta"
        source bowtie-2.3.4 && /usr/bin/time -v bowtie2 --threads 32 -x /tgac/workarea/group-pb/ENQ-1378_PIP-1073_Trevor_Wang_JIC.TW.ENQ-1378_01/annotation/Annotation_Aug2017/Analysis/bowtie-2.3.4/index/genome.1kb.fasta -1 1753_LIB21060_LDI18243_CGATGT_L001_R1.fastq.gz -2 1753_LIB21060_LDI18243_CGATGT_L001_R2.fastq.gz -S GrassPea-vs-GrassPea_output.sam && source samtools-1.5 && /usr/bin/time -v samtools sort -@16 -O BAM -o GrassPea-vs-GrassPea_output.bam GrassPea-vs-GrassPea_output.sam && /usr/bin/time -v samtools index GrassPea-vs-GrassPea_output.bam && /usr/bin/time -v samtools rmdup GrassPea-vs-GrassPea_output.bam nodup_GrassPea-vs-GrassPea_output.bam"
        """

        #!TODO: better bam generation
        rule qa_bwa_mem:
            input: 
                reads = getReads,
                assembly = getAssembly
            output:
                bam = join(qaa_env.bam_dir, "{sample}", "{sample}.create_bam.bam")
            log:
                join(qaa_env.log_dir, "{sample}.create_bam.log")
            params:
                outdir = lambda wildcards: join(qaa_env.bam_dir, wildcards.sample),
                ref = lambda wildcards: join(qaa_env.bam_dir, wildcards.sample, wildcards.sample + ".assembly.fasta"),
                outbam = lambda wildcards: join(qaa_env.bam_dir, wildcards.sample, wildcards.sample + ".create_bam.bam"),
                load = loadPreCmd(config["load"]["blob_bwa"])
            threads:
                8
            shell:
                "{params.load}" + TIME_CMD + \
                "bwa index -p {params.ref} {input.assembly} &&" + \
                "bwa mem -t {threads} {params.ref} {input.reads[0]} {input.reads[1]}" + \
                " | samtools view -buSh - | samtools sort -o {output.bam} - &&" + \
                " samtools index {output.bam} 2> {log}"


    if config["run_blobtools"]:
        rule qa_blast:
            input:
                assembly = getAssembly
            output:
                tsv = join(qaa_env.blast_dir, "{sample}", "{sample}.blast.tsv")
            log:
                join(qaa_env.log_dir, "{sample}.blast.log")
            params:
                outdir = lambda wildcards: join(qaa_env.blast_dir, wildcards.sample),
                load = loadPreCmd(config["load"]["blob_blast"],
                blastdb = config["resources"]["blob_blastdb"]
            threads:
                8
            shell: 
                "{params.load}" + TIME_CMD + \
                " blastn -db {params.blastdb} -task megablast -query {input.assembly}" + \
                " -outfmt '6 qseqid staxids score staxid bitscore std sscinames sskingdoms stitle'" + \
                " -culling_limit 5 -num_threads {threads} -evalue 1e-25 -max_target_seqs 10 -max_hsps 1 -perc_identity 75" + \
                " -out {output.tsv} &> {log}"

        rule qa_blobtools:
            input:
                bam = join(qaa_env.bam_dir, "{sample}", "{sample}.create_bam.bam") if config["create_bam"] else getBAM,
                blast = join(qaa_env.blast_dir, "{sample}", "{sample}.blast.tsv"),
                assembly = getAssembly
            output:
                blobtable = join(qaa_env.blob_dir, "{sample}", "{sample}.blobDB.table.txt")
            log:
                join(qaa_env.log_dir, "{sample}.blobtools.log")
            params:
                prefix = lambda wildcards: join(qaa_env.blob_dir, wildcards.sample, wildcards.sample),
                load = loadPreCmd(config["load"]["blobtools"]),
                taxlevel = "family"
            threads:
                1
            shell:
                "{params.load}" + TIME_CMD + \
                " blobtools create -i {input.assembly} -b {input.bam} -t {input.blast} -o {params.prefix} &&" + \
                " blobtools view -i {params.prefix}.blobDB.json -o $(dirname {params.prefix}) -r {params.taxlevel} &&" + \
                " blobtools blobplot -r {params.taxlevel} -l 1000 -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})" + \
                " &> {log}"                

