import os
import sys
import csv

def compileQUASTReport(quast_dir, out=sys.stdout):
    header = ""
    for cdir, dirs, files in os.walk(quast_dir):
        if "transposed_report.tsv" in files:
            with open(os.path.join(cdir, "transposed_report.tsv")) as qin:                
                for i, row in enumerate(csv.reader(qin, delimiter="\t")):
                    if i == 0 and not header:
                        header = row
                        print(*header, file=out, sep="\t")
                    elif i > 0:
                        print(*row, file=out, sep="\t")
