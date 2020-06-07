import argparse, re
import sys

parser = argparse.ArgumentParser()
parser.add_argument("vcf")
parser.add_argument("--fields", nargs="*", default=[])
parser.add_argument("--fieldnames", nargs="*", default=[])
parser.add_argument("--genotypes", nargs="*", default=[])
args = parser.parse_args()
seperators = [",", "|", ";"]

#output_infos = ["ANN[*][3]", "gnomad_AF"]
filename = args.vcf
output_infos = args.fields
genotype_infos = args.genotypes
fieldnames = args.fieldnames
use_infos = []
use_genotypes = []

assert(len(output_infos) == len(fieldnames))

for x in output_infos:
    name = re.match(r"\w*", x).group() 
    matches = re.findall(r"\[[0-9]+\]|\[\*\]", x)
    indices = []
    for m in matches:
        value = m[1:-1]
        if value.isnumeric():
            value = int(value)
        indices.append(value)
    use_infos.append((name, indices))


for x in genotype_infos:
    name = re.match(r"\w*", x).group()
    rest = x[len(name):]
    if len(rest) > 0:
        assert(rest[0] == "[" and rest[-1] == "]")
        rest = rest[1:-1]
        use_genotypes.append((name, rest.split(",")))
    else:
        use_genotypes.append((name, []))


def parse(indices, seperators, s):
    if len(indices) == 0:
        return s
    current_index = indices[0]
    current_sep = seperators[0]
    split = s.split(current_sep)
    if current_index == "*":
        return current_sep.join(set(parse(indices[1:], seperators[1:], s) for s in split))
    else:
        return split[current_index]


def get_gt_sum(s, gt):
    ret = []
    for g in gt:
        rg = re.escape(g)
        m = re.search("[0-9]*" + rg ,s)
        if m:
            ret.append(int(m.group(0)[:-len(g)]))
    return sum(ret)

# def parse_genotype(obs):
#     #8N-4V-4N+2V+

# print header
print("#chrom", "pos", "ref", "alt", "qual", sep="\t", end="\t")
print(*fieldnames, sep="\t", end="\t")
print()

if filename == "/dev/stdin" or filename == "-":
    f = sys.stdin
else:
    f = open(filename, 'r')

for line in f:
    if line.startswith("#"):
        continue
    split = line.strip().split("\t")
    CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = split[0:9]
    GENOTYPES = split[9:]
    #print(CHROM, POS, REF, ",".join([f"<{x.type}>" for x in variant.ALT]), variant.QUAL, sep="\t", end="")
    print(CHROM, POS, REF, ALT, QUAL, sep="\t", end="")

    info = {}
    for x in INFO.split(";"):
        x_split = x.split("=", 1)
        x0 = x_split[0]
        if len(x_split) == 2:
            x1 = x_split[1]
        else:
            x1 = True
        info[x0] = x1


    for variable, indices in use_infos:
        value = info.get(variable, "")
        print(end="\t")
        print(parse(indices, seperators, value), end="")


    d = {f:i for i,f in enumerate(FORMAT.split(":"))}
    for name, inside in use_genotypes:
        name_index = d[name]
        for g in GENOTYPES:
            print(end="\t")
            value = g.split(":")[name_index]
            if len(inside) > 0:
                value = get_gt_sum(value, inside)
            print(value, end="")
        #print(*[s.data.RC if s.data.RC is not None else "." for s in variant.samples], sep="\t", end="")
    print()
