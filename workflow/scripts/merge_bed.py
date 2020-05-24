from functools import reduce
from collections import defaultdict
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("beds", nargs="+")
parser.add_argument("--min_motif_len", type=int, default=3)
parser.add_argument("--min_est_count", type=int, default=30)
args = parser.parse_args()


filler = "."
#header
print("#chrom", "start", "id", "motif", "n", *map(lambda x: os.path.basename(x)[:-len("-genotype.bed")], args.beds), sep="\t")

data = list()

for filename in args.beds:
    with open(filename, "r") as f:
        d = defaultdict(lambda: filler)
        data.append(d)
        for line in f:
            split = line.strip().split("\t")
            chrom, start, stop, motif = split[0:4]
            if len(motif) < args.min_motif_len:
                continue
            if split[4] == "nan":
                continue
            d[(chrom, int(start), ".", motif)] = split[4]

keys = reduce(lambda a, b: a | b, map(lambda x: set(x.keys()), data))
for key in sorted(keys):
    estimations = [d[key] for d in data]
    valid_estimations = [float(x) for x in estimations if x != filler]
    if max(valid_estimations) < args.min_est_count:
        continue
    print(*key, "<DUP>", sum(d[key] != filler for d in data), ".", ".", "RC", *estimations, sep="\t")
