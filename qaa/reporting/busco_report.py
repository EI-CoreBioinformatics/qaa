import os
import sys
import re

def compileBUSCO(busco_dir, out=sys.stdout):
    print("Sample", "Mode", "Complete BUSCOs (C)", "Complete and single-copy BUSCOs (S)", "Complete and duplicated BUSCOs (D)", "Fragmented BUSCOs (F)", "Missing BUSCOs (M)", "Total BUSCO groups searched", "C%", "S%", "D%", "F%", "M%", sep="\t", file=out)
    for cdir, dirs, files in os.walk(busco_dir):
        summary = list(filter(lambda s:"short_summary" in s, files))
        if summary:
            with open(os.path.join(cdir, summary[0])) as _in:
                for line in _in:
                    if line.startswith('# BUSCO was run in mode:'):
                        mode = re.search(' (transcriptome|proteins|genome)', line).group().strip()
                        break
                [next(_in), next(_in), next(_in)]
                numbers = [int(line.strip().split()[0]) for line in _in]
                percentages = list(map(lambda x:"{:.1f}".format(int(x)/int(numbers[-1] * 100)), numbers[:-1]))
                print("_".join(summary[0].split("_")[1:3]), mode, *numbers, *percentages, sep="\t", file=out)

