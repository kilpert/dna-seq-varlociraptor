rule vcf2tsv:
    input:
        "{prefix}.bcf"
    output:
        "{prefix}.tsv"
    conda:
        "../envs/utils.yaml"
    shell:
        "bcftools view {input} | python workflow/scripts/vcf2tsv.py /dev/stdin --fields ANN[*][1] --fields ANN[*][1] ANN[*][3] --fieldnames region gene dp_ af_ obs_ --genotype DP AF OBS[V+,V-] > {output}"