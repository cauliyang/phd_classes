# -*- coding: utf-8 -*-
SAMPLES = [
            "SRR341549",
            "SRR1770413",
            ]
RG = {
        "SRR341549": "O104_H4",
        "SRR1770413":"K12",
        }

rule all:
    input:
        "e_colis.vcf"


rule bwa_map:
    input:
        read1 = "fq/{sample}_1.fastq",
        read2 = "fq/{sample}_2.fastq",
        ref = "E.coli_K12_MG1655.fa",
    params:
        rg = lambda wildcards: RG[wildcards.sample]

    output:
        "mapped_reads/{sample}.raw.sam"
    shell:
        r"bwa mem -R '@RG\tID:{params.rg}\tSM:{params.rg}' {input.ref} {input.read1} {input.read2} -o {output}"

rule sam2bam:
    input:
        "mapped_reads/{sample}.raw.sam"
    output:
        temp("mapped_reads/{sample}.raw.bam")
    shell:
        "samtools view -b {input} > {output}"


rule bam_sort:
    input:
        "mapped_reads/{sample}.raw.bam"
    output:
        "sorted_reads/{sample}.raw.sorted.bam"
    shell:
        "sambamba sort {input} -o {output}"

rule makeup:
    input:
        "sorted_reads/{sample}.raw.sorted.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "sambamba markdup {input} {output}"

rule index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule call:
    input:
        ref = "E.coli_K12_MG1655.fa",
        bam = expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai = expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "e_colis.vcf"
    shell:
      "freebayes -f {input.ref} --ploidy 1  {input.bam} > {output}"


