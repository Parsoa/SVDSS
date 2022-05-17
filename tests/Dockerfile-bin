FROM fedora:35

# RUN apt update -qqy
# RUN apt install -qqy wget bcftools samtools
# RUN apt clean ; apt autoclean -->
RUN dnf install -y wget samtools bcftools

# RUN wget -q https://github.com/Parsoa/SVDSS/releases/download/v1.0.2/SVDSS_linux_x86-64
COPY SVDSS_linux_x86-64 /
RUN chmod +x SVDSS_linux_x86-64

COPY run-svdss.sh /
COPY input /
RUN bash run-svdss.sh ./SVDSS_linux_x86-64 22.fa 22.bam output