FROM python:3.8

ADD sigProfilerMatrixGenerator_hg19.py .
ADD alexandrov_sig.py .
ADD vcf_alexandrov_docker.sh .
COPY requirements.txt .

RUN pip install -r requirements.txt 
RUN python3 sigProfilerMatrixGenerator_hg19.py
RUN mkdir /vcf_files

ENTRYPOINT ["/bin/bash", "/vcf_alexandrov_docker.sh"]



