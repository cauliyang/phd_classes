SAMPLES = ["daughter_SRR089532",
           "father_SRR089545",
           "mother_SRR089553",
            ]

rule all:
    input:
        expand("calls/{sample}.raw.vcf", sample=SAMPLES)



rule bowtie:
    input:
        read1 = "fq/{sample}_1.fastq.gz",
        read2 = "fq/{sample}_2.fastq.gz"
    output:
        temp("mapped_reads/{sample}.sam")
    shell:
        "bowtie2 -x chr17 -1 {input.read1} -2 {input.read2} -S {output}"

rule sam2bam:
    input:
        "mapped_reads/{sample}.sam"
    output:
        temp("mapped_reads/{sample}.bam")
    shell:
        "samtools view -bS {input} > {output}"


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

rule bcftools_call:
    input:
        fa="chr17.fa",
        bam = "sorted_reads/{sample}.sorted.bam"
    output:
        temp("calls/{sample}.mp")
    shell:
        "samtools mpileup -uf {input.fa} {input.bam} -o {output}"

rule bcftools_final:
    input:
        "calls/{sample}.mp"
    output:
        "calls/{sample}.raw.vcf"
    shell:
        "bcftools view -Ov {input} > {output}"
