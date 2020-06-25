from pysam import VariantFile
from collections import defaultdict
import numpy as np
import argparse
import tempfile, os

parser = argparse.ArgumentParser()
parser.add_argument("vcf")
parser.add_argument("cand_mother_field")
parser.add_argument("cand_father_field")
args = parser.parse_args()

transcripts_mother = defaultdict(list)
transcripts_father = defaultdict(list)

p_transcript_compound = defaultdict(lambda: 0)


def P(Q):
    return np.power(10,-Q/10)

def at_least_one(ps):
    return 1 - np.prod(1 - P(ps))

with tempfile.TemporaryDirectory() as tmpdirname:
    with VariantFile(args.vcf, "r") as f, VariantFile(os.path.join(tmpdirname, "tmp.bcf"), "w", header=f.header) as tmp_out:
        for record in f:
            tmp_out.write(record)
            p_mother = record.info[args.cand_mother_field]
            p_father = record.info[args.cand_father_field]
            transcripts = set(ann.split("|")[6] for ann in record.info["ANN"])
            for transcript in transcripts:
                transcripts_mother[transcript].extend(p_mother)
                transcripts_father[transcript].extend(p_father)

        for transcript in set(transcripts_mother) & set(transcripts_father):
            p_mother = np.array(transcripts_mother[transcript], dtype=np.float128)
            p_father = np.array(transcripts_father[transcript], dtype=np.float128)
            p = at_least_one(p_mother) * at_least_one(p_father)
            if p > 0:
                p_transcript_compound[transcript] = -10 * np.log10(p)
            else:
                p_transcript_compound[transcript] = np.nan

    with VariantFile(os.path.join(tmpdirname, "tmp.bcf"), "r") as f, VariantFile("-", "w", header=f.header) as out:
        f.header.info.add("PROB_COMPOUND_HETEROZYGOUS","A","Float","Posterior probability for compound heterozygous variant (PHRED)")
        out.header.info.add("PROB_COMPOUND_HETEROZYGOUS","A","Float","Posterior probability for compound heterozygous variant (PHRED)")

        for record in f:
            probs = []
            new_ann = []
            for ann in record.info["ANN"]:
                transcript = ann.split("|")[6]
                p = p_transcript_compound[transcript]
                if p > 0:
                    new_ann.append(f"{ann}|{p}")
                    probs.append(p)
            record.info["ANN"] = new_ann
            if not new_ann:
                continue

            record.info["PROB_COMPOUND_HETEROZYGOUS"] = float(max(probs))
            out.write(record)