import argparse, re
import vcf

parser = argparse.ArgumentParser()
parser.add_argument("vcf")
parser.add_argument("--fields", nargs="*", default=[])
parser.add_argument("--fieldnames", nargs="*", default=[])
parser.add_argument("--genotypes", action="store_true", default=False)
args = parser.parse_args()
seperators = [",", "|", ";"]

#output_infos = ["ANN[*][3]", "gnomad_AF"]
filename = args.vcf
output_infos = args.fields
fieldnames = args.fieldnames
use_infos = []

assert(len(output_infos) == len(fieldnames))

for x in output_infos:
    name = re.match(r"\w*", x).group()
    matches = re.findall(r"\[[0-9]+\]|\[\*\]", x)
    indices = []
    x = len(matches)
    for m in matches:
        value = m[1:-1]
        if value.isnumeric():
            value = int(value)
        indices.append(value)
    use_infos.append((name, indices))


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

# print header
print("#chrom", "pos", "ref", "alt", "qual", sep="\t", end="\t")
print(*fieldnames, sep="\t", end="\t")
print()

for variant in vcf.Reader(open(filename, 'r')): # or VCF('some.bcf')
    print(variant.CHROM, variant.POS, variant.REF, ",".join([f"<{x.type}>" for x in variant.ALT]), variant.QUAL, sep="\t", end="")
    info = variant.INFO
    for variable, indices in use_infos:
        value = info.get(variable, "")[0]
        print(end="\t")
        print(parse(indices, seperators, value), end="")

    if args.genotypes:
        print(end="\t")
        print(*[s.data.RC if s.data.RC is not None else "." for s in variant.samples], sep="\t", end="")
    print()
