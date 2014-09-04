FROM dmlond/bwa_samtools_base
MAINTAINER Darin London <darin.london@duke.edu>

ADD split_raw.pl /usr/local/bin/split_raw.pl
RUN ["chmod", "777", "/usr/local/bin/split_raw.pl"]
USER bwa_user
ENTRYPOINT ["/usr/local/bin/split_raw.pl"]
