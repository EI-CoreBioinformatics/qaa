import sys
import os
import csv

from collections import Counter, namedtuple

BlobSample = namedtuple("BlobSample", "sample ncontigs dom_org dom_org_ncontigs dom_org_perc dom_org_span subdom_org subdom_org_ncontigs subdom_org_perc subdom_org_span".split(" "))

def compileBlobReport(blob_dir, out=sys.stdout):
    print("Running BLOBREPORT... " + blob_dir)
    print("Sample", "#contigs", "Predominant genus", "#contigs(Predominant genus)", "%(Predominant genus)", "span(Predominant genus)[bp]", "Subdominant genus", "#contigs(Subdominant genus)", "%(Subdominant genus)", "span(Subdominant genus)[bp]", sep="\t", file=out)
    for cdir, dirs, files in os.walk(blob_dir):
        blobtable = list(filter(lambda s:s.endswith(".blobDB.table.txt"), files))
        if blobtable:
            sample = blobtable[0].split(".")[0]
            taxcounter, spancounter, taxmap = Counter(), Counter(), dict()
            with open(join(cdir, blobtable[0])) as tin:
                for row in csv.reader(tin, delimiter="\t"):
                    if not row[0].startswith("#"):
                        genus = row[5].split(" ")[0]
                        taxcounter[genus] += 1
                        spancounter[genus] += int(row[1])
                        taxmap.setdefault(genus, Counter())[row[5]] += 1
            orgs = sorted(taxcounter.items(), key=lambda x:x[1], reverse=True)
            ncontigs = sum(taxcounter.values())
            dom_org, vdom_org = orgs.pop(0), (None, 0)
            if orgs:
                vdom_org = orgs.pop(0)
            blob_data = BlobSample(sample, ncontigs, dom_org[0], dom_org[1], dom_org[1]/ncontigs if ncontigs else None, spancounter[dom_org[0]], vdom_org[0], vdom_org[1], vdom_org[1]/ncontigs if ncontigs else None, spancounter[vdom_org[0]])
            print(*blob_data, sep="\t", file=out)
            # yield blob_data
