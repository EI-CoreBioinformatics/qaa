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
            # blobtable = join(qaa_env.blob_dir, "{sample}", "{sample}.blobDB.table.txt")
            TARGETS.append(join(qaa_env.blob_dir, sample, sample + ".blobDB.table.txt"))
            #    join(qaa_env.qualimap_dir, "{sample}", "qualimapReport.html")
            TARGETS.append(join(qaa_env.qualimap_dir, sample, "qualimapReport.html"))
        if config["run_busco"]:
            TARGETS.append(join(qaa_env.busco_geno_dir, sample, sample + "_short_summary.txt"))
    if config["run_transcriptome_module"]:
        TARGETS.append(join(qaa_env.busco_tran_dir, sample, sample + "_short_summary.txt"))
    if config["run_proteome_module"]:
        TARGETS.append(join(qaa_env.busco_prot_dir, sample, sample + "_short_summary.txt"))

#if config["run_multiqc"]:
#    TARGETS.append(join(config["multiqc_dir", config["misc"]["project"] + "_" + runmode + "_multiqc_report.html"))

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
        # TARGETS[-1] is multiqc report
        input:
            TARGETS[:-1]
        output:
            join(config["multiqc_dir"], config["misc"]["project"] + "_" + runmode + "_multiqc_report.html")
        params:
            load = loadPreCmd(config["load"]["multiqc"]),
            mqc_config = config["resources"]["multiqc_config"],
            datadir = qaa_env.output_dir,
            outdir = config["multiqc_dir"],
            prefix = config["misc"]["project"] + "_" + runmode,
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
            busco_input = getTranscripts
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
        rule qa_busco_geno:
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
                load = loadPreCmd(config["load"]["busco"]),
                busco_data = lambda wildcards: getBUSCOData(wildcards.sample),
                busco_mode = "geno"
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
        """# ./PRO1841_E2_PlateA_E2/PRO1841_E2_PlateA_E2.markdup.bam.rgindex
# i=1; for f in $(find . -name '*.markdup.bam'); do  echo $i > $f.rgindex; i=$((i+1)); done
rule picard_add_correct_rg:
    input:
        bam = join(BWADIR, "{sample}", "{sample}.markdup.bam")
    output:
        bam = join(BWADIR, "{sample}", "{sample}.markdup.rg.bam")
    log:
        join(BWADIR, "logs", "{sample}.picard_add_correct_rg.log")
    params:
        sample = basename("{sample}.markdup.bam").strip(".markdup.bam")
    shell:
        "set +u" + \
        " && source jre-8u144" + \
        " && java -jar /tgac/software/testing/picardtools/2.18.4/x86_64/picard.jar AddOrReplaceReadGroups I={input.bam} O={output.bam} RGLB={params.sample} RGPL=illumina RGPU={params.sample} RGSM={params.sample} RGID=$(cat {input.bam}.rgindex) &> {log}"
        """
        """ --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field
  --rg <text>        add <text> ("lab:value") to @RG line of SAM header.
                     Note: @RG line only printed when --rg-id is set.
@RG     ID:1    LB:PRO1841_B1_PlateA_B1 PL:illumina     SM:PRO1841_B1_PlateA_B1 PU:PRO1841_B1_PlateA_B1
-R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
        """
        BAM_THREADS = 16
        QA_ALIGNER = "bowtie2"
        if QA_ALIGNER == "bowtie2":
            QA_ALIGN_LOAD = loadPreCmd(config["load"]["blob_bowtie2"])
            QA_ALIGN_BUILD_INDEX = "bowtie2-build --threads {threads} {input.assembly} {params.ref}"
            QA_ALIGN = "bowtie2 --threads {params.align_threads} -x {params.ref} -1 {input.reads[0]} -2 {input.reads[1]} --rg-id {sample} --rg LB:{sample} --rg PL:illumina --rg SM:{sample} --rg PU:{sample}"
        else: # bwa
            QA_ALIGN_LOAD = loadPreCmd(config["load"]["blob_bwa"])
            QA_ALIGN_BUILD_INDEX = "bwa index -p {params.ref} {input.assembly}"
            QA_ALIGN = "bwa mem -t {threads} -R '@RG\tID:1\tLB:{sample}\tPL:illumina\tSM:{sample}\tPU:{sample}' {params.ref} {input.reads[0]} {input.reads[1]}"

        rule qa_align_reads:
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
                align_threads = BAM_THREADS / 2,
                sort_threads = BAM_THREADS / 2,
                load_align = QA_ALIGN_LOAD, #loadPreCmd(config["load"]["blob_bowtie2"]),
                load_markdup = loadPreCmd(config["load"]["picard"]),
                cmd_markdup = loadPreCmd(config["cmd"]["picard_markdup"], is_dependency=False)
            threads:
                BAM_THREADS
            shell:
                "{params.load_align}" + \
                TIME_CMD + " " + QA_ALIGN_BUILD_INDEX + " &&" + \
                " " + TIME_CMD + " " + QA_ALIGN + " |" + \
                " samtools sort -@ {params.sort_threads} -o {output.bam}.tmp.bam - &&" + \
                " {params.load_markdup}" + \
                TIME_CMD + " {params.cmd_markdup}" + \
                " INPUT={output.bam}.tmp.bam OUTPUT={output.bam} METRICS_FILE={output.bam}.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate &&" + \
                " samtools index {output.bam} &&" + \
                " 2> {log}"

    if config["run_blobtools"]:
        rule qa_qualimap:
            input:
                bam = join(qaa_env.bam_dir, "{sample}", "{sample}.align_reads.bam") if config["create_bam"] else getBAM,
            output:
                join(qaa_env.qualimap_dir, "{sample}", "qualimapReport.html")
            params:
                load = loadPreCmd(config["load"]["qualimap"]),
                outdir = lambda wildcards: join(qaa_env.qualimap_dir, wildcards.sample),
                mem = config["qualimap_mem"]  # this is coming from within qaa/__init__.py
            log: 
                join(qaa_env.log_dir, "{sample}.qualimap.log")
            threads: 2
            message: "Using qualimap to collect stats for: {input.bam}"
            shell: 
               "{params.load}" + \
               TIME_CMD + " qualimap --java-mem-size={params.mem} bamqc " + \
               "-bam {input.bam} -nt {threads} -outdir {params.outdir} &> {log}"




        rule qa_blast:
            input:
                assembly = getAssembly
            output:
                tsv = join(qaa_env.blast_dir, "{sample}", "{sample}.blast.tsv")
            log:
                join(qaa_env.log_dir, "{sample}.blast.log")
            params:
                outdir = lambda wildcards: join(qaa_env.blast_dir, wildcards.sample),
                load = loadPreCmd(config["load"]["blob_blast"]),
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
                bam = join(qaa_env.bam_dir, "{sample}", "{sample}.align_reads.bam") if config["create_bam"] else getBAM,
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
                " blobtools view -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/ -r {params.taxlevel} &&" + \
                " blobtools blobplot -r {params.taxlevel} -l 1000 -i {params.prefix}.blobDB.json -o $(dirname {params.prefix})/" + \
                " &> {log}"                
