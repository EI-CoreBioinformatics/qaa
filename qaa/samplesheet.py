from collections import namedtuple
import csv

QAA_SAMPLESHEET_COLUMNS = "id assembly bamfile r1 r2 busco_id transcripts proteins".split(" ")

QAA_Sample = namedtuple("QAA_Sample", QAA_SAMPLESHEET_COLUMNS)

def readQAASamplesheet(_in):
    for row in _in:
        if row and row[0]:
            row.extend([""] * max(0, len(QAA_SAMPLESHEET_COLUMNS) - len(row)))
            yield (row[0], QAA_Sample(*row))
