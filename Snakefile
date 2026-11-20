import csv

configfile: "config.yaml"

# Sample names = "sample" column of the sample table (nothing is hard-coded here)
with open(config["samples"]) as f:
    SAMPLES = [row["sample"] for row in csv.DictReader(f, delimiter="\t")]

HISAT2_INDEX = multiext(config["hisat2_index"],
                        ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                        ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2")

wildcard_constraints:
    read = "1|2"


# ---------------------------------------------------------------- final target
rule all:
    input:
        expand("results/fastqc_raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]),
        expand("results/fastqc_trimmed/{sample}_trimmed_{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]),
        expand("results/aligned/{sample}.bam.bai", sample=SAMPLES),
        "results/featurecounts/counts.txt",
        expand("results/kallisto/{sample}/abundance.tsv", sample=SAMPLES),


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
    conda: "envs/qc.yaml"
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
    conda: "envs/qc.yaml"
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
    conda: "envs/qc.yaml"
    shell:
        "fastqc --threads {threads} --outdir results/fastqc_trimmed {input} > {log} 2>&1"


# ------------------------------------------------ 2.2 alignment-based quantification
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
    conda: "envs/align.yaml"
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
    conda: "envs/align.yaml"
    shell:
        # align -> sort by position -> index; SAM is never written to disk
        "hisat2 -p {threads} -x {params.prefix} -1 {input.r1} -2 {input.r2} 2> {log} "
        "| samtools sort -@ {threads} -o {output.bam} - "
        "&& samtools index {output.bam}"


rule featurecounts:
    # One run over all BAM files -> one gene x sample count table
    input:
        bams = expand("results/aligned/{sample}.bam", sample=SAMPLES),
        gtf  = config["annotation_gtf"]
    output:
        counts  = "results/featurecounts/counts.txt",
        summary = "results/featurecounts/counts.txt.summary"   # read by MultiQC
    log: "logs/featurecounts.log"
    params: strand = config["featurecounts_strand"]
    threads: config["threads"]
    resources: mem_mb = config["mem_large"]
    conda: "envs/quant.yaml"
    shell:
        # -p --countReadPairs: count fragments (read pairs), not single reads
        "featureCounts -p --countReadPairs -s {params.strand} -T {threads} "
        "-a {input.gtf} -o {output.counts} {input.bams} > {log} 2>&1"


# ------------------------------------------------ 2.3 alignment-free quantification
rule kallisto_index:
    # Runs only if the index is missing (kept in ref/ for reuse)
    input:
        config["transcripts_fa"]
    output:
        config["kallisto_index"]
    log: "logs/kallisto_index.log"
    threads: config["threads"]
    resources:
        mem_mb = config["mem_large"],
        runtime = 600
    conda: "envs/quant.yaml"
    shell:
        "kallisto index -t {threads} -i {output} {input} > {log} 2>&1"


rule kallisto_quant:
    input:
        idx = config["kallisto_index"],
        r1  = "results/trimmed/{sample}_trimmed_1.fastq.gz",
        r2  = "results/trimmed/{sample}_trimmed_2.fastq.gz"
    output:
        "results/kallisto/{sample}/abundance.tsv",
        "results/kallisto/{sample}/abundance.h5",
        "results/kallisto/{sample}/run_info.json"
    log: "logs/kallisto/{sample}.log"          # pseudoalignment summary (read by MultiQC)
    threads: config["threads"]
    resources: mem_mb = config["mem_large"]
    conda: "envs/quant.yaml"
    shell:
        "kallisto quant -t {threads} -i {input.idx} -o results/kallisto/{wildcards.sample} "
        "{input.r1} {input.r2} > {log} 2>&1"
