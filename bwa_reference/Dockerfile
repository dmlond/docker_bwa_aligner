# This will install bwa 0.5.9-r16 and samtools 0.1.18 (r982:295)
FROM dmlond/bwa_samtools_base
MAINTAINER Darin London <darin.london@duke.edu>

RUN ["/usr/bin/yum", "install", "-y", "--nogpgcheck", "wget"]
ADD add_reference.pl /usr/local/bin/add_reference.pl
RUN ["chmod", "777", "/usr/local/bin/add_reference.pl"]
ENTRYPOINT ["/usr/local/bin/add_reference.pl"]
