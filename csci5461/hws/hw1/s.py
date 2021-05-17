# -*- coding: utf-8 -*-
SAMPLES = ["NIST7035_TAAGGCGA_L001",
           "NIST7035_TAAGGCGA_L002",
           "NIST7086_CGTACTAG_L001",
           "NIST7086_CGTACTAG_L002"
            ]

rule all:
    input:
        "merge/final_merge.bam"



rule bwa:
    input:
        read1 = "{sample}_R1_001.fastq.gz",
        read2 = "{sample}_R2_001.fastq.gz",
        ref = "human_g1k_v37.fasta.gz"
    output:
        "mapped_reads/{sample}.sam"
    shell:
        "bwa mem {input.read1} {input.read2} > {output}"

rule fix:
    input:
        "mapped_reads/{sample}.sam"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "samtools fixmate -O bam {input} {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.sorted.bam")
    shell:
        "samtools sort {input} -o {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.sorted.bam"
    output:
        protected("sorted_reads/{sample}.sorted.bam.bai")
    shell:
        "samtools index {input}"

rule merge:
    input:
        bam = expand("sorted_reads/{sample}.sorted.bam", sample=SAMPLES)
    output:
        "merge/final_merge.bam"
    shell:
        "samtools merge {output} {input.bam}"
