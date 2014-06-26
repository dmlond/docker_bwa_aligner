# This will install bwa 0.5.9-r16 and samtools 0.1.18 (r982:295)
FROM blalor/centos
MAINTAINER Darin London <darin.london@duke.edu>

RUN ["/usr/bin/yum", "install", "-y", "--nogpgcheck", "wget"]
RUN ["mkdir", "-p", "/data"]
RUN ["wget", "-O", "/data/ERR022523_1.fastq.gz", "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022523/ERR022523_1.fastq.gz"] 
RUN ["wget", "-O", "/data/ERR022523_2.fastq.gz", "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022523/ERR022523_2.fastq.gz"] 
VOLUME ["/data"]
ADD attribution.txt /etc/attribution.txt
CMD ["cat", "/etc/attribution.txt"]
