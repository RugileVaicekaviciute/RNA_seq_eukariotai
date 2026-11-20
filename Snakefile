import csv

configfile: "config.yaml"

with open(config["samples"]) as f:
    SAMPLES = [row["sample"] for row in csv.DictReader(f, delimiter="\t")]

ENV = "envs/rnaseq.yaml"

wildcard_constraints:
    read = "1|2"


# ---------------------------------------------------------------- final target
rule all:
    input:
        expand("results/fastqc_raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]),
        expand("results/fastqc_trimmed/{sample}_trimmed_{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]),


# ---------------------------------------------------------- 2.1 preprocessing
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
