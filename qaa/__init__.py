import pkg_resources
import os
from os.path import join, dirname, basename
import csv
from collections import namedtuple
import yaml
import sys
from copy import copy

from .snakemake_helper import *
from qaa.reporting.busco_report import compileBUSCOReport
from qaa.reporting.quast_report import compileQUASTReport
from qaa.reporting.blobtools_report import compileBlobReport
from qaa.reporting.sixteen_s_reporter import SixteenSReporter, HEADER as sixteenS_header

__title__ = "qaa"
__author__ = 'Christian Schudoma (cschu)'
__license__ = 'MIT'
__copyright__ = 'Copyright 2018 Earlham Institute'
__version__ = pkg_resources.require("qaa")[0].version


ETC_DIR = join(dirname(__file__), "..", "etc")
DEFAULT_HPC_CONFIG_FILE = join(ETC_DIR, "hpc_config.json")
DEFAULT_CONFIG_FILE = join(ETC_DIR, "qaa_config.yaml")
QAA_ID = "XXX"

class QAA_ArgumentsAdapter(object):
    def __init__(self, **kwargs):
        for k in kwargs:
            setattr(self, k, kwargs[k])
        pass
    def update(self, **kwargs):
        for k in kwargs:
            setattr(self, k, kwargs[k])
        pass 

class QAA_Environment(object):
    def __init__(self, config):
        self.cwd = config.get("cwd", ".")

        self.output_dir = config.get("out_dir", ".")
        if not self.output_dir.startswith("/"):
            self.output_dir = join(self.cwd, self.output_dir)
        self.qa_dir = join(self.output_dir, "qa", "survey" if config["survey_assembly"] else "asm")
        self.qc_dir = join(self.output_dir, "qc")
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
        self.qualimap_dir = join(self.qa_dir, "qualimap")
 
        self.blast_dir = join(self.qa_dir, "blast")
        self.bam_dir = join(self.qa_dir, "bam")

        self.prokka_dir = join(self.output_dir, "annotation", "prokka")

class QAA_Runner(object):
    def __init__(self, args):
        print("QAA_RUNNER:__init__")
        def __create_input_stream(args):
            print("QAA_RUNNER:__create_input_stream", args, sep="\n")
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


        args = copy(args) # ???
        print("QAAA_ARGS:", args, file=sys.stderr)

        print("Configuring execution environment ... ", end="", flush=True)
        self.output_dir = args.output_dir
        self.report_dir = join(self.output_dir, "reports")
        if not os.path.exists(self.report_dir):
        	os.makedirs(self.report_dir)
        self.logs_dir = join(self.output_dir, "hpc_logs")
        self.exe_env = ExecutionEnvironment(args, NOW, job_suffix=args.input + "_" + self.output_dir, log_dir=self.logs_dir)
        print("done.")
        print(str(self.exe_env))
        print()
        self.runmode = args.runmode if hasattr(args, "runmode") else "asm"

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

        # from eimethyl: Extract qualimap mem setting from hpc_config, convert to gigabytes and subtract a bit to account for memory above
        # java heap
        import json
        print("HPCCONFIG", args.hpc_config, file=sys.stderr)
        data = json.load(open(args.hpc_config))
        qmem = str((int(int(data["qaa_qualimap"]["memory"]) / 1000)) - 2) + "G"
        self.config["qualimap_mem"] = qmem

                	
        if hasattr(args, "make_input_stream"):
            self.config["samplesheet"] = __create_input_stream(args)
            self.config["has_stream"] = True
            self.config["full_samplesheet"] = args.input
        else:
            self.config["has_stream"] = False
            self.config["full_samplesheet"] = self.config["samplesheet"] = args.input

        self.config["out_dir"] = self.output_dir
        self.config["etc"] = join(dirname(__file__), "..", "etc")
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
            self.config["misc"]["project"] = args.project_prefix # api
        except:
            pass

        self.config["run_multiqc"] = args.run_multiqc if hasattr(args, "run_multiqc") else True # api
        self.config["multiqc_dir"] = args.multiqc_dir if hasattr(args, "multiqc_dir") else "."
        self.config["survey_assembly"] = self.runmode == "survey" # api
        self.config["run_blobtools"] = args.run_blobtools if hasattr(args, "run_blobtools") else True # api
        self.config["run_busco"] = args.run_busco if hasattr(args, "run_busco") else True # api
        self.config["run_genome_module"] = self.runmode == "survey" or {"geno", "genome"}.intersection(requested_modes)
        self.config["run_transcriptome_module"] = self.runmode != "survey" and {"tran", "transcriptome"}.intersection(requested_modes)
        self.config["run_proteome_module"] = self.runmode != "survey" and {"prot", "proteome"}.intersection(requested_modes)
        self.config["quast_mincontiglen"] = args.quast_mincontiglen if hasattr(args, "quast_mincontiglen") else 0
        self.config["annotation"] = args.annotation if hasattr(args, "annotation") else None

        self.new_config_file = join(self.output_dir, "qaa.conf.yaml")
        with open(self.new_config_file, "w") as conf_out:
            yaml.dump(self.config, conf_out, default_flow_style=False)

        self.qaa_env = QAA_Environment(self.config)

        self.unlock = args.unlock

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
            # report_func(join(self.output_dir, "qa", "survey", "quast"), join(self.report_dir, "quast_survey_report.tsv"), compileQUASTReport)
            report_func(self.qaa_env.quast_dir, join(self.report_dir, "quast_survey_report.tsv"), compileQUASTReport)
            if self.config["run_busco"]:
                # report_func(join(self.output_dir, "qa", "survey", "busco", "geno"), join(self.report_dir, "busco_survey_report.tsv"), compileBUSCOReport)
                report_func(self.qaa_env.busco_geno_dir, join(self.report_dir, "busco_survey_report.tsv"), compileBUSCOReport)
            if self.config["run_blobtools"]:
                # report_func(join(self.output_dir, "qa", "survey", "blobtools"), join(self.report_dir, "blobtools_survey_report.tsv"), compileBlobReport)
                report_func(self.qaa_env.blob_dir, join(self.report_dir, "blobtools_survey_report.tsv"), compileBlobReport)
        elif self.runmode == "asm":
            if self.config["run_genome_module"]:
                # report_func(join(self.output_dir, "qa", "asm", "quast"), join(self.report_dir, "quast_report.tsv"), compileQUASTReport)
                report_func(self.qaa_env.quast_dir, join(self.report_dir, "quast_report.tsv"), compileQUASTReport)
            if self.config["run_busco"]:
                # report_func(join(self.output_dir, "qa", "asm", "busco", "geno"), join(self.report_dir, "busco_genome_report.tsv"), compileBUSCOReport)
                report_func(self.qaa_env.busco_geno_dir, join(self.report_dir, "busco_genome_report.tsv"), compileBUSCOReport)
            if self.config["run_blobtools"]:
                # report_func(join(self.output_dir, "qa", "asm", "blobtools"), join(self.report_dir, "blobtools_report.tsv"), compileBlobReport)
                report_func(self.qaa_env.blob_dir, join(self.report_dir, "blobtools_report.tsv"), compileBlobReport)
            if self.config["run_transcriptome_module"] or self.config["run_proteome_module"]:
                # report_func(join(self.output_dir, "qa", "asm", "busco"), join(self.report_dir, "busco_report.tsv"), compileBUSCOReport)
                report_func(self.qaa_env.busco_dir, join(self.report_dir, "busco_report.tsv"), compileBUSCOReport)
            if self.config["annotation"] is not None:
                with open(join(self.report_dir, "16S_report.tsv"), "wt") as ostream:
                    # SixteenSReporter(join(self.output_dir, "annotation", "prokka"), header=sixteenS_header, out=ostream).generateReport()
                    SixteenSReporter(self.qaa_env.prokka_dir, header=sixteenS_header, out=ostream).generateReport()
