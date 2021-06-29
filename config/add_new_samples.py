#!/usr/bin/env python

import json
import glob
import os
import re

#query="p[FM]5"
query="pOM1"

fastq_dir="/vol/pico/exome_seq/raw/sequencing_runs"
units_tsv="/vol/pico/exome_pools/dna-seq-varlociraptor/config/units_pools.tsv"
samples_tsv="/vol/pico/exome_pools/dna-seq-varlociraptor/config/samples_pools.tsv"

d = {}
sra = ""
adapters = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


print("{:#^60s}".format(" Files "))
print("query:", query)
print("fastq_dir:", fastq_dir)


print("{:#^60s}".format(" Files "))
path = os.path.join(fastq_dir, "**", f"{query}_*_R?_*.fastq.gz")

for file in glob.glob(path, recursive=True):
    print(file)
    e = os.path.basename(file).split(".")[0].split("_")
    sample_name = e[0]
    sample_number = e[1]
    lane_number = e[2]
    read = e[3]
    last_segment = e[4]
    ext = os.path.basename(file).split(".", 1)[1]

    print("sample_name:", sample_name)
    print("sample_number:", sample_number)
    print("lane_number:", lane_number)
    print("read:", read)
    print("last_segment:", last_segment)
    print("ext:", ext)
    print()

    if sample_name not in d:
        d[sample_name] = {}

    if lane_number not in d[sample_name]:
        d[sample_name][lane_number] = {}

    if read == "R1":
        d[sample_name][lane_number]["R1"] = file
    if read == "R2":
        d[sample_name][lane_number]["R2"] = file


print("{:#^60s}".format(" JSON "))
# print(d)
print(json.dumps(d, indent=4, sort_keys=True))


print("{:#^60s}".format(" units.tsv "))
lines = []
for sample_name in sorted(d):
    for lane_number in sorted(d[sample_name]):
        line = "\t".join([sample_name,
                         lane_number,
                         d[sample_name][lane_number]["R1"],
                         d[sample_name][lane_number]["R2"],
                         sra,
                         adapters])
        lines.append(line)


with open(units_tsv, "a") as f:
    lines = [line+"\n" for line in lines]
    f.writelines(lines)
print("Saved (appended) to:", units_tsv)


print("{:#^60s}".format(" samples.tsv "))

lines = []

for sample_name in sorted(d):

    alias = ""
    group = "NA"
    platform = "ILLUMINA"
    purity = "NA"
    n = 5

    m = re.search("pF", sample_name)
    if m:
        alias = "fatherpool"

    m = re.search("pM", sample_name)
    if m:
        alias = "motherpool"


    line = "\t".join([sample_name,
                      alias,
                      group,
                      platform,
                      purity,
                      str(n)])
    lines.append(line)


with open(samples_tsv, "a") as f:
    lines = [line+"\n" for line in lines]
    f.writelines(lines)
print("Saved (appended) to:", samples_tsv)
