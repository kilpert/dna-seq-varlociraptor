#rule cnvpytor_index:
#    input:
#        "resources/genome.fasta"
#    output:
#        "results/cnvpytor/index/genome.pytor"
#    wrapper:
#        "file:/vol/huge/christo/snakemake-wrappers/bio/cnvpytor/index"

rule cnvpytor:
    threads:
        20
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bam.bai",
        ref="resources/genome.fasta"
    output:
        pytor="results/cnvpytor/{sample}.pytor",
        calls10000="results/cnvpytor/{sample}.10000.tsv",
        calls100000="results/cnvpytor/{sample}.100000.tsv",
    log:
        log="logs/cnvpytor/{sample}.log"
    shell:
        """
        (
        cnvpytor -root {output.pytor} -rd {input.bam} -j {threads} -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 1 16 17 18 19 20 21 22 X
        cnvpytor -root {output.pytor} -his 10000 100000 -j {threads}
        cnvpytor -root {output.pytor} -partition 10000 100000 -j {threads}
        cnvpytor -root {output.pytor} -call 10000 -j {threads} > {output.calls10000}
        cnvpytor -root {output.pytor} -call 100000 -j {threads} > {output.calls100000}
        ) 2> {log}
        """

rule cnvpytor_to_vcf:
    threads:
        4
    input:
        calls="results/cnvpytor/{sample}.{length}.tsv"
    output:
        calls="results/cnvpytor/{sample}.{length}.vcf"
    log:
        log="logs/cnvpytor/tovcf.{sample}.{length}.log"
    shell:
        "(perl workflow/scripts/cnvnator2VCF.pl {input.calls} -reference=ressources/genome.fasta | snpEff -Xmx20G hg38 - > {output.calls}) 2> {log}"

rule extract_genes:
    threads:
        1
    input:
        calls="results/cnvpytor/{sample}.{length}.vcf"
    output:
        calls="results/cnvpytor/{sample}.{length}.genes.txt"
    log:
        log="logs/cnvpytor/extract_genes.{sample}.{length}.log"
    shell:
        """
        SnpSift filter "ANN[ANY].EFFECT != 'intergenic_region' & exists ANN[ANY].GENE" {input.calls} | 
        SnpSift extractFields - CHROM POS END ALT "ANN[*].GENE" |
        tr "\&" "\\t" > {output}
        """

rule build_matrix:
    input:
        genes="resources/genes.grouped.bed",
        calls=expand("results/cnvpytor/{sample}.{length}.genes.txt", sample=samples["sample_name"], length=[10000, 100000])
    output:
        "results/cnvpytor_matrix.tsv"
    shell:
        "python workflow/scripts/create_gene_matrix.py {input.genes} {input.calls} > {output}"
