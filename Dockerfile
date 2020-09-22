FROM ubuntu:latest

RUN apt update && apt install -y python3-pip git

#RUN git clone https://github.com/Redmar-van-den-Berg/MultiQC.git /multiqc \
#        && cd /multiqc \
#        && git checkout chimera \
#        && pip3 install .

RUN mkdir /pgx 
COPY setup.py multiqc_pgx /pgx/

RUN mkdir /pgx/multiqc_pgx/
COPY multiqc_pgx /pgx/multiqc_pgx/

RUN cd /pgx && pip3 install -e .
