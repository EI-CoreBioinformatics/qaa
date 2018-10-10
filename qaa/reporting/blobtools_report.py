import sys
import os
from os.path import join
import csv

from collections import Counter, namedtuple

# BlobSample = namedtuple("BlobSample", "sample ncontigs dom_org dom_org_ncontigs dom_org_perc dom_org_span subdom_org subdom_org_ncontigs subdom_org_perc subdom_org_span".split(" "))

BlobSample = namedtuple("BlobSample", "sample ncontigs" + \
                                      " Pgbp Pgbp_size Pgbp_size_pc Pgbp_ctg Pgbp_ctg_pc" + \
                                      " Pgct Pgct_size Pgct_size_pc Pgct_ctg Pgct_ctg_pc" + \
                                      " Sgbp Sgbp_size Sgbp_size_pc Sgbp_ctg Sgbp_ctg_pc" + \
                                      " Sgct Sgct_size Sgct_size_pc Sgct_ctg Sgct_ctg_pc".split(" "))

HEADER = [
          "Sample", "#contigs", 
          "Predominant genus by size (Pgbp)", "Predominant genus by contig (Pgct)",
          "size (Pgbp) [bp]", "%size (Pgbp)", "#contigs (Pgbp)", "%contigs (Pgbp)",  
          "size (Pgct) [bp]", "%size (Pgct)", "#contigs (Pgct)", "%contigs (Pgct)",  
          "Subdominant genus by size (Sgbp)", "Subdominant genus by contig (Sgct)",
          "size (Sgbp) [bp]", "%size (Sgbp)", "#contigs (Sgbp)", "%contigs (Sgbp)",  
          "size (Sgct) [bp]", "%size (Sgct)", "#contigs (Sgct)", "%contigs (Sgct)",
         ] 


def compileBlobReport(blob_dir, out=sys.stdout):
    print("Running BLOBREPORT... " + blob_dir)
    print(*HEADER, sep="\t", file=out)
    for cdir, dirs, files in os.walk(blob_dir):
        blobtable = list(filter(lambda s:s.endswith(".blobDB.table.txt"), files))
        if blobtable:
            sample = blobtable[0].split(".")[0]
            taxcounter, spancounter, taxmap = Counter(), Counter(), dict()
            # # name  length  GC      N       bam0    family.t.6      family.s.7      family.c.8
            # 1       193272  0.5448  0       172.984 Lactobacillaceae        15529.0 0
            with open(join(cdir, blobtable[0])) as tin:
                for row in csv.reader(tin, delimiter="\t"):
                    if not row[0].startswith("#"):
                        genus = row[5].split(" ")[0]
                        taxcounter[genus] += 1
                        spancounter[genus] += int(row[1])
                        taxmap.setdefault(genus, Counter())[row[5]] += 1
            
            ncontigs, totalsize = sum(taxcounter.values()), sum(spancounter.values())
            orgs_by_size = sorted(spancounter.items(), key=lambda x:x[1], reverse=True)
            orgs_by_nctg = sorted(taxcounter.items(), key=lambda x:x[1], reverse=True)
            pdomorg_by_size = orgs_by_size.pop(0)
            try:
                sdomorg_by_size = orgs_by_size.pop(0)
            except:
                sdomorg_by_size = (None, 0)
            pdomorg_by_nctg = orgs_by_nctg.pop(0)
            try:
                sdomorg_by_nctg = orgs_by_nctg.pop(0)
            except:
                sdomorg_by_nctg = (None, 0)
            blob_data = BlobSample(sample, 
                                   ncontigs,
                                   pdomorg_by_size[0],                      
                                   pdomorg_by_nctg[0],
                                   spancounter[pdomorg_by_size[0]],
                                   spancounter[pdomorg_by_size[0]] / totalsize if totalsize else None,
                                   taxcounter[pdomorg_by_size[0]],                                                                      
                                   taxcounter[pdomorg_by_size[0]] / ncontigs if ncontigs else None,
                                   spancounter[pdomorg_by_nctg[0]],
                                   spancounter[pdomorg_by_nctg[0]] / totalsize if totalsize else None,
                                   taxcounter[pdomorg_by_nctg[0]],                                                                      
                                   taxcounter[pdomorg_by_nctg[0]] / ncontigs if ncontigs else None,
                                   spancounter[sdomorg_by_size[0]],
                                   spancounter[sdomorg_by_size[0]] / totalsize if totalsize else None,
                                   taxcounter[sdomorg_by_size[0]],                                                                      
                                   taxcounter[sdomorg_by_size[0]] / ncontigs if ncontigs else None,
                                   spancounter[sdomorg_by_nctg[0]],
                                   spancounter[sdomorg_by_nctg[0]] / totalsize if totalsize else None,
                                   taxcounter[sdomorg_by_nctg[0]],                                                                      
                                   taxcounter[sdomorg_by_nctg[0]] / ncontigs if ncontigs else None)

            print(*blob_data, sep="\t", file=out)
