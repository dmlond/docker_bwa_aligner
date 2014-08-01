FROM centos:centos6
MAINTAINER Darin London <darin.london@duke.edu>

RUN ["/usr/sbin/useradd", "bwa_user"]
RUN ["/usr/bin/yum", "install", "-y", "wget"]
ADD attribution.txt /home/bwa_user/attribution.txt
RUN ["chown", "bwa_user", "/home/bwa_user/attribution.txt"]
WORKDIR /home/bwa_user
USER bwa_user
RUN ["mkdir", "-p", "data"]
RUN ["wget", "-O", "/home/bwa_user/data/ERR022523_1.fastq.gz", "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022523/ERR022523_1.fastq.gz"] 
RUN ["wget", "-O", "/home/bwa_user/data/ERR022523_2.fastq.gz", "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022523/ERR022523_2.fastq.gz"] 
VOLUME ["/home/bwa_user/data"]
CMD ["cat", "/home/bwa_user/attribution.txt"]
