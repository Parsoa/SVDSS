configfile: "config.yaml"


from os.path import join as pjoin

REF = config["fa"]
BAM = config["bam"]

ODIR = config["out"]
THREADS = config["threads"]

SVDSS_BIN = config["bin"]


rule run:
    input:
        pjoin(ODIR, "calls.vcf"),


rule index:
    input:
        fa=REF,
    output:
        fmd=REF + ".fmd",
    threads: 1
    benchmark:
        pjoin(ODIR, "benchmark", "pp-index.txt")
    shell:
        """
        {SVDSS_BIN} index --reference {input.fa} --index {output.fmd}
        """


rule smooth:
    input:
        fa=REF,
        bam=BAM,
    output:
        bam=pjoin(ODIR, "smoothed.bam"),
    params:
        wd=pjoin(ODIR),
    threads: THREADS
    shell:
        """
        {SVDSS_BIN} smooth --reference {input.fa} --bam {input.bam} --threads {threads} > {output.bam}
        samtools index {output.bam}
        """


rule search:
    input:
        fmd=REF + ".fmd",
        bam=pjoin(ODIR, "smoothed.bam"),
    output:
        sfs=pjoin(ODIR, "specifics.txt"),
    threads: THREADS
    shell:
        """
        {SVDSS_BIN} search --index {input.fmd} --bam {input.bam} --threads {threads} > {output.sfs}
        """


rule call:
    input:
        fa=REF,
        bam=pjoin(ODIR, "smoothed.bam"),
        sfs=pjoin(ODIR, "specifics.txt"),
    output:
        vcf=pjoin(ODIR, "calls.vcf"),
    threads: THREADS
    shell:
        """
        {SVDSS_BIN} call --reference {input.fa} --bam {input.bam} --sfs {input.sfs} --threads {threads} | bcftools sort > {output.vcf}
        """
