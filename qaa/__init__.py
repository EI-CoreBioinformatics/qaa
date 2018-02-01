import pkg_resources

__title__ = "qaa"
__author__ = 'Christian Schudoma (cschu)'
__license__ = 'MIT'
__copyright__ = 'Copyright 2018 Earlham Institute'
__version__ = pkg_resources.require("qaa")[0].version

import os

ETC_DIR = os.path.join(os.path.dirname(__file__), "..", "etc")
DEFAULT_HPC_CONFIG_FILE = os.path.join(ETC_DIR, "hpc_config.json")
DEFAULT_CONFIG_FILE = os.path.join(ETC_DIR, "qaa_config.yaml")

import csv
from collections import namedtuple

QAA_SAMPLESHEET_COLUMNS = "id assembly bamfile r1 r2 busco_id transcripts proteins".split(" ")
QAA_Sample = namedtuple("QAA_Sample", QAA_SAMPLESHEET_COLUMNS)

def readSamplesheet(_in):
	import csv
	for row in csv.reader(_in, delimiter=","):
		if row and row[0]:
			while len(row) != len(QAA_SAMPLESHEET_COLUMNS):
				row.append("")
			yield (row[0], QAA_Sample(*row))


TIME_CMD = " /usr/bin/time -v"

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
