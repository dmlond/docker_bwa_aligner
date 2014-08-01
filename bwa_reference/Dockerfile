FROM dmlond/bwa_samtools_base
MAINTAINER Darin London <darin.london@duke.edu>

RUN ["/usr/bin/yum", "install", "-y", "--nogpgcheck", "wget"]
ADD add_reference.pl /usr/local/bin/add_reference.pl
RUN ["chmod", "777", "/usr/local/bin/add_reference.pl"]
USER bwa_user
ENTRYPOINT ["/usr/local/bin/add_reference.pl"]
