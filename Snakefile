import csv

configfile: "config.yaml"

with open(config["samples"]) as f:
    SAMPLES = [row["sample"] for row in csv.DictReader(f, delimiter="\t")]

ENV = "envs/rnaseq.yaml"

HISAT2_INDEX = multiext(config["hisat2_index"],
                        ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                        ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")


wildcard_constraints:
    read = "1|2"


rule all:
    input:
        expand("results/fastqc_raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]),
        expand("results/fastqc_trimmed/{sample}_trimmed_{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]),
        expand("results/aligned/{sample}.bam.bai", sample=SAMPLES),

# ----------------------------------------------------------
rule fastqc_raw:
    input:
        config["raw_reads"]
    output:
        html = "results/fastqc_raw/{sample}_{read}_fastqc.html",
        zip  = "results/fastqc_raw/{sample}_{read}_fastqc.zip"
    log: "logs/fastqc_raw/{sample}_{read}.log"
    threads: 2
    resources: mem_mb = config["mem_small"]
    conda: ENV
    shell:
        "fastqc --threads {threads} --outdir results/fastqc_raw {input} > {log} 2>&1"


rule trim:
    input:
        # fill in {read}, keep {sample} as a wildcard
        r1 = expand(config["raw_reads"], read=1, allow_missing=True),
        r2 = expand(config["raw_reads"], read=2, allow_missing=True)
    output:
        r1   = "results/trimmed/{sample}_trimmed_1.fastq.gz",
        r2   = "results/trimmed/{sample}_trimmed_2.fastq.gz",
        html = "results/trimmed/{sample}.fastp.html",
        json = "results/trimmed/{sample}.fastp.json"
    log: "logs/fastp/{sample}.log"
    params: extra = config["fastp_params"]
    threads: config["threads"]
    resources: mem_mb = config["mem_small"]
    conda: ENV
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--detect_adapter_for_pe --thread {threads} {params.extra} "
        "--html {output.html} --json {output.json} > {log} 2>&1"


rule fastqc_trimmed:
    input:
        "results/trimmed/{sample}_trimmed_{read}.fastq.gz"
    output:
        html = "results/fastqc_trimmed/{sample}_trimmed_{read}_fastqc.html",
        zip  = "results/fastqc_trimmed/{sample}_trimmed_{read}_fastqc.zip"
    log: "logs/fastqc_trimmed/{sample}_{read}.log"
    threads: 2
    resources: mem_mb = config["mem_small"]
    conda: ENV
    shell:
        "fastqc --threads {threads} --outdir results/fastqc_trimmed {input} > {log} 2>&1"

# ------------------------------------------------ 
rule hisat2_index:
    # Runs only if the index files are missing (or older than the genome).
    input:
        config["genome_fa"]
    output:
        HISAT2_INDEX
    log: "logs/hisat2_index.log"
    params: prefix = config["hisat2_index"]
    threads: config["threads"]
    resources:
        mem_mb = config["mem_large"],
        runtime = 600
    conda: ENV
    shell:
        "hisat2-build -p {threads} {input} {params.prefix} > {log} 2>&1"


rule hisat2_align:
    input:
        r1  = "results/trimmed/{sample}_trimmed_1.fastq.gz",
        r2  = "results/trimmed/{sample}_trimmed_2.fastq.gz",
        idx = HISAT2_INDEX
    output:
        bam = "results/aligned/{sample}.bam",
        bai = "results/aligned/{sample}.bam.bai"
    log: "logs/hisat2/{sample}.log"          # HISAT2 alignment summary (read by MultiQC)
    params: prefix = config["hisat2_index"]
    threads: config["threads"]
    resources: mem_mb = config["mem_large"]
    conda: ENV
    shell:
        # align -> sort by position -> index; SAM is never written to disk
        "hisat2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} 2> {log} "
        "| samtools sort -@ {threads} -o {output.bam} - "
        "&& samtools index {output.bam}"