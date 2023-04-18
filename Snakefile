configfile: "config.yaml"

from os.path import join as pjoin

REF = config["fa"]
BAM = config["bam"]

ODIR = config["out"] # os.getcwd()
THREADS = config["threads"]

SVDSS_BIN = config["bin"]

rule run:
    input:
        pjoin(ODIR, "svdss.vcf")

rule pp_index:
    input:
        fa = REF
    output:
        fmd = REF + ".fmd"
    threads: 1
    benchmark: pjoin(ODIR, "benchmark", "pp-index.txt")
    shell:
        """
        {SVDSS_BIN} index --reference {input.fa} --index {output.fmd}
        """

rule pp_reconstruct:
    input:
        fa = REF,
        bam = BAM
    output:
        bam = pjoin(ODIR, "smoothed.bam")
    params:
        wd = pjoin(ODIR)
    threads: THREADS
    shell:
        """
        {SVDSS_BIN} smooth --reference {input.fa} --bam {input.bam} --threads {threads} > {output.bam}
        """

rule pp_search:
    input:
        fmd = REF + ".fmd",
        bam = pjoin(ODIR, "smoothed.bam")
    output:
        sfs = pjoin(ODIR, "solution_batch_0.assembled.sfs")
    params:
        wd = pjoin(ODIR)
    threads: THREADS
    shell:
        """
        {SVDSS_BIN} search --index {input.fmd} --bam {input.bam} --threads {threads} --workdir {params.wd} --assemble
        """

rule pp_call:
    input:
        fa = REF,
        bam = pjoin(ODIR, "smoothed.bam"),
        sfs = pjoin(ODIR, "solution_batch_0.assembled.sfs")
    output:
        vcf = pjoin(ODIR, "svdss.vcf")
    params:
        wd = pjoin(ODIR)
    threads: THREADS
    shell:
        """
        {SVDSS_BIN} call --reference {input.fa} --bam {input.bam} --threads {threads} --workdir {params.wd}
        bcftools sort {params.wd}/svs_poa.vcf > {output.vcf}
        """
