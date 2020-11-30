## person id of the uploaded hg19 candidates
person_id = os.path.basename(config["varvis_hg19_csv"]).split("_")[0]


rule varvis_hg19_csv_to_candidates_vcf:
    input:
        config["varvis_hg19_csv"]
    output:
        "results/candidates/hg19/{person_id}.hg19_candidates.vcf"
    conda:
        "../envs/plink.yaml"
    log:
        "results/log/{person_id}.hg19_candidates.vcf.log"
    shell:
        "bash workflow/scripts/csv2vcf.sh {input} {output} "
        "> {log} "


## rule hg19_vcf_to_bed:
##     input:
##         "results/candidates/hg19/{person_id}.hg19_candidates.vcf"
##     output:
##         "results/candidates/hg19/{person_id}.hg19_candidates.bed"
##     shell:
##         "cat {input} "
##         "| grep -v '^#' "
##         """| awk 'BEGIN{{OFS="\t"}}{{l=length($5);print $1, $2-1, $2+l-1, $3}}' """
##         "> {output} "


rule download_hg19ToHg38_over_chain:
    output:
        "resources/hg19ToHg38.over.chain.gz" # UCSC style chromosome names
    cache:
        True
    shell:
        "curl -L http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz "
        "| pigz -dc "
        "| sed 's/chr//g' " # convert chromosome names to Ensembl style
        "| pigz -9 "
        "> {output}"


rule picard_liftovervcf:
    input:
        vcf = "results/candidates/hg19/{person_id}.hg19_candidates.vcf",
        chain = "resources/hg19ToHg38.over.chain.gz",
        fasta = "resources/genome.fasta"
    output:
        vcf="results/candidates/hg38/{person_id}.hg38_candidates.vcf",
        reject="results/candidates/hg38/{person_id}.hg38_candidates.rejected.vcf"
    log:
        "results/log/{person_id}.picard_liftovervcf.log"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard -Xmx6g LiftoverVcf "
        "I={input.vcf} "
        "CHAIN={input.chain} "
        "R={input.fasta} "
        "O={output.vcf} "
        "REJECT={output.reject} "
        ">{log} 2>&1 "


# rule annotatevarvis:
    # input:
        # "varvis.tsv"
        #   "varlociraptor.vcf"
    # output:
        # "results/parent_annotated/{candidates}.{poolfather}.{poolmothers}.csv"
    # shell:
        # ""
