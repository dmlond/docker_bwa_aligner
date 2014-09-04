# This will install bwa 0.5.9-r16 and samtools 0.1.18 (r982:295)
FROM dmlond/bwa_samtools_base
MAINTAINER Darin London <darin.london@duke.edu>

ADD bwa_aligner.pl /usr/local/bin/bwa_aligner.pl
RUN ["chmod", "777", "/usr/local/bin/bwa_aligner.pl"]
USER bwa_user
ENTRYPOINT ["/usr/local/bin/bwa_aligner.pl"]

