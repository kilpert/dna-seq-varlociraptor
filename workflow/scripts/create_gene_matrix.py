from collections import defaultdict
import argparse
import os

def samplename(filename):
    return os.path.basename(filename).split(".")[0]


parser = argparse.ArgumentParser()
parser.add_argument("genes")
parser.add_argument("cnvpytor", nargs="+")
args = parser.parse_args()

# load genes
gene_positions = defaultdict(lambda: "?:?-?")
with open(args.genes) as f:
    for line in f:
        chrom, start, end, symbol, typ = line.strip().split("\t")
        gene_positions[symbol] = f"{chrom}:{start}-{end}"

d = defaultdict(dict)
for filename in args.cnvpytor :
    with open(filename, "r") as f:
        sample = samplename(filename)
        firstline = True
        for line in f:
            if firstline:
                firstline = False
                continue
            split = line.strip().split("\t")
            chrom, start, end, alt = split[0:4]
            called_genes = set(split[5:])
            for gene in called_genes:
                if gene not in gene_positions:
                    gene_split = gene.split("-")
                    gene_split_corrected = []
                    for i in range(len(gene_split)):
                        if gene_split[i][0].isnumeric():
                            gene_split_corrected[-1] = gene_split_corrected[-1] + "-" + gene_split[i]
                        else:
                            gene_split_corrected.append(gene_split[i])
                    continue
                    for gene in gene_split_corrected:
                        d[gene][sample] = (sample, chrom, start, end, alt)
                else:
                    d[gene][sample] = (sample, chrom, start, end, alt)

samplenames = list(sorted(set([samplename(filename) for filename in args.cnvpytor])))
print("#", "symbol", sep="\t", end="\t")
print(*sorted(samplenames), sep="\t")

def convert_dup_del(x):
    if x == "<DEL>":
        return "-"
    if x == "<DUP>":
        return "+"
    assert(False)

for gene in sorted(d):
    included = [convert_dup_del(d[gene][sample][4]) if sample in d[gene] else "" for sample in samplenames]
    print(gene_positions[gene], end="\t")
    print(gene, end="\t")
    print(*included, sep="\t")


#sort -k4,4 genes.bed | groupBy -g 1,4 -c 4,2,3 -o count,min,max | awk -v OFS='\t' '{print $1, $4, $5, $2, $3}' | cut -f 1-4 > genes.grouped.bed