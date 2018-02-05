import os
import sys
import csv

def compileQUASTReport(quast_dir, out=sys.stdout):
    header = ""
    for cdir, dirs, files in os.walk(quast_dir):
        if "transposed_report.tsv" in files:
            with open(os.path.join(cdir, "transposed_report.tsv")) as qin:
                r = csv.reader(qin, delimiter="\t")
                first = next(r)
                if not header:
                    header = first
                    print(*header, file=out, sep="\t")
                    # print(*header, sep="\t")
                    # yield header
                for row in r:
                    print(*row, file=out, sep="\t")
                    # print(*row, sep="\t")
                    # yield row
