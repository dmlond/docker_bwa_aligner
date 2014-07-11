FROM dmlond/bwa_samtools_base
MAINTAINER Darin London <darin.london@duke.edu>

ADD split_raw.pl /usr/bin/split_raw.pl
RUN ["chmod", "777", "/usr/bin/split_raw.pl"]
ENTRYPOINT ["/usr/bin/split_raw.pl"]
