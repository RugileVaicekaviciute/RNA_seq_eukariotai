#!/bin/bash
# Download paired-end FASTQ files for every sample in samples.tsv from ENA into data/raw/.
# Run once on the login node, before Snakemake:
#     bash download_reads.sh
# ENA mirrors SRA and gives ready gzipped FASTQ files + MD5 checksums.
set -euo pipefail

mkdir -p data/raw

tail -n +2 samples.tsv | cut -f1 | while read -r run; do
    # ENA answers with a 2-line table: header + "run_accession  fastq_ftp  fastq_md5"
    # (paired-end: two URLs and two MD5s, separated by ';')
    info=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv" | tail -n 1)

    paste <(cut -f2 <<<"$info" | tr ';' '\n') <(cut -f3 <<<"$info" | tr ';' '\n') |
    while read -r url md5; do
        f=data/raw/$(basename "$url")
        if [ -s "$f" ] && echo "$md5  $f" | md5sum -c --quiet - 2>/dev/null; then
            echo "$f already downloaded, skipping"; continue
        fi
        wget -c -O "$f" "ftp://$url"
        echo "$md5  $f" | md5sum -c -
    done
done

ls -lh data/raw
