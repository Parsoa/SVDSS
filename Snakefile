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
        {SVDSS_BIN} index --fastq {input.fa} --index {output.fmd}
        """

rule pp_reconstruct:
    input:
        fa = REF,
        bam = BAM
    output:
        bam = pjoin(ODIR, "reconstructed.selective.bam")
    params:
        wd = pjoin(ODIR)
    threads: THREADS
    shell:
        """
        {SVDSS_BIN} reconstruct --reference {input.fa} --bam {input.bam} --workdir {params.wd} --threads {threads}
        """

rule pp_sortreconstruct:
    input:
        bam = pjoin(ODIR, "reconstructed.selective.bam")
    output:
        bam = pjoin(ODIR, "reconstructed.selective.sorted.bam")
    params:
        sam = pjoin(ODIR, "reconstructed.sam"),
        tmpd = pjoin(ODIR)
    threads: 1
    shell:
        """
        samtools sort -T {output.bam}.sort-tmp {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule pp_search:
    input:
        fmd = REF + ".fmd",
        bam = pjoin(ODIR, "reconstructed.selective.sorted.bam")
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
        bam = pjoin(ODIR, "reconstructed.selective.sorted.bam"),
        sfs = pjoin(ODIR, "solution_batch_0.assembled.sfs")
    output:
        vcf = pjoin(ODIR, "svdss.vcf")
    params:
        wd = pjoin(ODIR)
    threads: THREADS
    shell:
        """
        n=$(ls {params.wd}/solution_batch_*.assembled.sfs | wc -l)
        {SVDSS_BIN} call --reference {input.fa} --bam {input.bam} --threads {threads} --workdir {params.wd} --batches ${{n}}
        bcftools sort {params.wd}/svs_poa.vcf > {output.vcf}
        """