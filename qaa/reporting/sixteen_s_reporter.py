import sys
import os
from os.path import join
import csv
import glob

from abc import ABC, abstractmethod

from collections import Counter, namedtuple

from qaa.reporting.qaa_reporter import QAAReporter

# join -1 1 -2 1 -a 1 <(find . -name '*.16S' -print -exec grep -c ">" {} \; | cut -f 2 -d \/ | awk -v OFS="\t" '{ id=$1; getline; print id,$1; }' | sort -k1,1) <(find . -name '*.16S' -print | xargs grep partial | cut -f 2 -d \/ | awk -v OFS="\t" '{ct[$1]+=1} END { for (k in ct) print k,ct[k] }' | sort -k1,1) | tr " " "\t" | awk -v OFS="\t" 'BEGIN { print "Sample","predicted_16SrRNAs","predicted_16SrRNAs_partial" } { if (NF < 3) print $0,0}'  > ../../reports/16S_report.tsv

HEADER = [
          "Sample",
          "Predicted_16S_rRNAs",
          "Predicted_partial_16S_rRNAs"
         ]

#class QAAReporter(ABC):
#    def __init__(self, indir, header=list(), out=sys.stdout):
#        self.indir = indir
#        self.outstream = out
#        self.header = header
#        super().__init__()
#    @abstractmethod
#    def generateReport(self):
#        pass

class SixteenSReporter(QAAReporter):
    def generateReport(self):
        print(*self.header, sep="\t", file=self.outstream)
        for cdir, dirs, files in os.walk(self.indir):
            gff = glob.glob(join(cdir, "*.prokka.gff"))
            sixteenS, sixteenS_partial = 0, 0
            if gff:
                with open(gff[0]) as gff_in:
                    for row in csv.reader(gff_in, delimiter="\t"):
                        if not row[0].startswith("#"):
                            if row[1].startswith("barrnap") and row[2] == "rRNA":
                                attrib = dict(item.split("=") for item in row[8].split(";"))
                                if attrib.get("product", "").startswith("16S ribosomal RNA"):
                                    if attrib.get("product", "").endswith("(partial)"):
                                        sixteenS_partial += 1
                                    else:
                                        sixteenS += 1
                print(os.path.basename(cdir),
                      sixteenS,
                      sixteenS_partial,
                      sep="\t", file=self.outstream)
            
        pass
