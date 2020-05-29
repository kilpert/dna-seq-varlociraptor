rule vcf_to_tsv:
    input:
        vcf="strling.vcf"
    output:
        tsv="strling.tsv"
    conda:
        "../envs/utils.yaml"
    shell:
        """python workflow/scripts/vcf2tsv.py {input.vcf} --fields ANN[*][1] --fields ANN[*][1] ANN[*][3] --fieldnames region gene --genotypes > {output.tsv}"""

rule create_vcf:
    input:
        expand("results/strling_call/{sample}-genotype.bed", sample=samples["sample_name"].values),
    output:
        vcf="strling.vcf",
        tmp_bed=temp("tmp/tmp.strling.bed"),
        tmp_vcf=temp("tmp/tmp.strling.vcf"),
    params:
        samples="\t".join(samples["sample_name"].values)
    shell:
        """
        python workflow/scripts/merge_bed.py {input} > {output.tmp_bed}
        snpEff -Xmx20g ann hg38 {output.tmp_bed} > {output.tmp_vcf}
        echo "##fileformat=VCFv4.2" > {output.vcf}
        cat {output.tmp_vcf} | grep ^## >> {output.vcf}
        echo -e '##FORMAT=<ID=RC,Number=A,Type=Float,Description="Estimated number of repeats">' >> {output.vcf}
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{params.samples}" >> {output.vcf}
        cat {output.tmp_vcf} | grep -v ^# >> {output.vcf}
        """


rule genotype_to_bed_and_filter:
    input:
        genotype="results/strling_call/{sample}-genotype.txt"
    output:
        bed="results/strling_call/{sample}-genotype.bed"
    shell:
        """
        cat {input.genotype} |
        grep -v "^#" |
        awk 'BEGIN {{FS="\\t"; OFS=FS}} {{ print $1, $2, $3, $4, $6}}' |
        sort > {output.bed}"""

#if (length($4) > 2 && $6 != "nan" &&  $6 >= 30)

rule strling_call:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bin="results/strling/extract/{sample}.bin",
        merged="results/strling/merge/merged-bounds.txt",
        ref="resources/genome.fasta",
    output:
        bounds="results/strling_call/{sample}-bounds.txt",
        genotype="results/strling_call/{sample}-genotype.txt",
        unplaced="results/strling_call/{sample}-unplaced.txt",
    params:
        prefix="results/strling_call/{sample}"
    log:
        "logs/strling/{sample}.call.log"
    conda:
        "../envs/strling.yaml"
    shell:
        "workflow/scripts/strling call -f {input.ref} {input.bam} {input.bin} -o {params.prefix} -b {input.merged} 2> {log}"


rule strling_merge:
    input:
        bam=expand("results/strling/extract/{sample}.bin", sample=samples["sample_name"].values),
        ref="resources/genome.fasta",
    output:
        merged="results/strling/merge/merged-bounds.txt"
    params:
        prefix="results/strling/merge/merged"
    log:
        "logs/strling/merge.log"
    conda:
        "../envs/strling.yaml"
    shell:
        "workflow/scripts/strling merge -f {input.ref} -o {params.prefix} {input.bam} 2> {log}"


rule strling_extract:
    input:
        bam="results/recal/{sample}.sorted.bam",
        ref="resources/genome.fasta",
    output:
        bin="results/strling/extract/{sample}.bin",
    log:
        "logs/strling/{sample}.extract.log"
    conda:
        "../envs/strling.yaml"
    shell:
        "workflow/scripts/strling extract -f {input.ref} {input.bam} {output.bin} 2> {log}"