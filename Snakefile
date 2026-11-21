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
# Every result of the workflow (also the input of the MultiQC report)
RESULTS = [
    expand("results/fastqc_raw/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]),
    expand("results/fastqc_trimmed/{sample}_trimmed_{read}_fastqc.zip", sample=SAMPLES, read=[1, 2]),
    expand("results/aligned/{sample}.bam.bai", sample=SAMPLES),
    "results/featurecounts/counts.txt",
    expand("results/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
    expand("results/rseqc/{sample}.{qc}",
           sample=SAMPLES,
           qc=["bam_stat.txt", "infer_experiment.txt", "read_distribution.txt",
               "geneBodyCoverage.txt", "inner_distance_freq.txt"]),
    "results/deeptools/correlation_heatmap.png",
    "results/deeptools/pca.png",
]

rule all:
    input:
        RESULTS,
        "results/multiqc/multiqc_report.html"


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


# ------------------------------------------------ 2.3 alignment quality assessment
rule qc_bed:
    # GTF -> BED12 with one canonical transcript per protein-coding gene.
    # A small, clean gene model keeps RSeQC fast (esp. gene body coverage)
    # and avoids counting the same region many times through isoforms.
    input:
        config["annotation_gtf"]
    output:
        config["qc_bed"]
    log: "logs/qc_bed.log"
    resources: mem_mb = config["mem_small"]
    conda: "envs/rseqc.yaml"
    shell:
        "(grep 'tag \"Ensembl_canonical\"' {input} | grep 'gene_type \"protein_coding\"' > {output}.tmp.gtf "
        "&& gtfToGenePred {output}.tmp.gtf {output}.tmp.gp "
        "&& genePredToBed {output}.tmp.gp {output} "
        "&& rm {output}.tmp.gtf {output}.tmp.gp) > {log} 2>&1"


rule rseqc_bam_stat:                     # a) alignment statistics
    input: "results/aligned/{sample}.bam"
    output: "results/rseqc/{sample}.bam_stat.txt"
    log: "logs/rseqc/{sample}.bam_stat.log"
    resources: mem_mb = config["mem_small"]
    conda: "envs/rseqc.yaml"
    shell: "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer_experiment:             # b) library strandedness
    input:
        bam = "results/aligned/{sample}.bam",
        bed = config["qc_bed"]
    output: "results/rseqc/{sample}.infer_experiment.txt"
    log: "logs/rseqc/{sample}.infer_experiment.log"
    resources: mem_mb = config["mem_small"]
    conda: "envs/rseqc.yaml"
    shell: "infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log}"


rule rseqc_read_distribution:            # c) reads in CDS / UTR / introns / intergenic
    input:
        bam = "results/aligned/{sample}.bam",
        bed = config["qc_bed"]
    output: "results/rseqc/{sample}.read_distribution.txt"
    log: "logs/rseqc/{sample}.read_distribution.log"
    resources: mem_mb = config["mem_small"]
    conda: "envs/rseqc.yaml"
    shell: "read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log}"


rule rseqc_gene_body_coverage:           # d) 5'->3' coverage (RNA degradation check)
    input:
        bam = "results/aligned/{sample}.bam",
        bai = "results/aligned/{sample}.bam.bai",
        bed = config["qc_bed"]
    output: "results/rseqc/{sample}.geneBodyCoverage.txt"
    log: "logs/rseqc/{sample}.gene_body_coverage.log"
    resources:
        mem_mb = config["mem_small"],
        runtime = 1440                   # slowest step of the workflow
    conda: "envs/rseqc.yaml"
    shell: "geneBody_coverage.py -i {input.bam} -r {input.bed} -o results/rseqc/{wildcards.sample} > {log} 2>&1"


rule rseqc_inner_distance:               # e) insert size between read pairs
    input:
        bam = "results/aligned/{sample}.bam",
        bed = config["qc_bed"]
    output: "results/rseqc/{sample}.inner_distance_freq.txt"
    log: "logs/rseqc/{sample}.inner_distance.log"
    resources: mem_mb = config["mem_small"]
    conda: "envs/rseqc.yaml"
    shell: "inner_distance.py -i {input.bam} -r {input.bed} -o results/rseqc/{wildcards.sample} > {log} 2>&1"


rule deeptools_summary:
    # read counts of every sample over the same genes -> input for f) and g)
    input:
        bams = expand("results/aligned/{sample}.bam", sample=SAMPLES),
        bais = expand("results/aligned/{sample}.bam.bai", sample=SAMPLES),
        bed  = config["qc_bed"]
    output: "results/deeptools/multibam_genes.npz"
    log: "logs/deeptools/multibam_summary.log"
    threads: config["threads"]
    resources: mem_mb = config["mem_small"]
    conda: "envs/deeptools.yaml"
    shell:
        "multiBamSummary BED-file --BED {input.bed} --bamfiles {input.bams} --smartLabels "
        "-p {threads} -o {output} > {log} 2>&1"


rule deeptools_correlation:              # f) correlation between samples
    input: "results/deeptools/multibam_genes.npz"
    output:
        png = "results/deeptools/correlation_heatmap.png",
        tab = "results/deeptools/correlation_matrix.tab"    # read by MultiQC
    log: "logs/deeptools/correlation.log"
    resources: mem_mb = config["mem_small"]
    conda: "envs/deeptools.yaml"
    shell:
        "plotCorrelation -in {input} --corMethod spearman --skipZeros --whatToPlot heatmap "
        "--plotNumbers -o {output.png} --outFileCorMatrix {output.tab} > {log} 2>&1"


rule deeptools_pca:                      # g) principal component analysis
    input: "results/deeptools/multibam_genes.npz"
    output:
        png = "results/deeptools/pca.png",
        tab = "results/deeptools/pca.tab"                   # read by MultiQC
    log: "logs/deeptools/pca.log"
    resources: mem_mb = config["mem_small"]
    conda: "envs/deeptools.yaml"
    shell:
        "plotPCA -in {input} -o {output.png} --outFileNameData {output.tab} > {log} 2>&1"


# ------------------------------------------------ 2.4 reporting
rule multiqc:
    # Collects FastQC, fastp, HISAT2, featureCounts, kallisto, RSeQC and
    # deepTools results and logs into one HTML report
    input:
        RESULTS
    output:
        "results/multiqc/multiqc_report.html"
    log: "logs/multiqc.log"
    resources: mem_mb = config["mem_small"]
    conda: "envs/multiqc.yaml"
    shell:
        "multiqc --force --outdir results/multiqc results logs > {log} 2>&1"
