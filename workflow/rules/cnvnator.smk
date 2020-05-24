rule cnvnator:
    input:
        bam="results/recal/{sample}.sorted.bam",
        bai="results/recal/{sample}.sorted.bam.bai",
        ref="resources/genome.fasta"
    output:
        root="results/cnvnator/{sample}.root"
    conda:
        "../envs/cnvnator.yaml"
    shell:
        """
        workflow/scripts/cnvnator -root {output.root} -tree {input.bam} -chrom $(seq 1 22) X Y
        workflow/scripts/cnvnator -root {output.root} -tree {input.bam} -his 1000 -fasta {input.ref}
        workflow/scripts/cnvnator -root {output.root} -stat 1000 
        workflow/scripts/cnvnator -root {output.root} -partition 1000
        workflow/scripts/cnvnator -root {output.root} -call 1000
        """
