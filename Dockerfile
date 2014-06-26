# This will install bwa 0.5.9-r16 and samtools 0.1.18 (r982:295)
FROM blalor/centos
MAINTAINER Darin London <darin.london@duke.edu>

RUN ["/usr/bin/yum", "clean", "all"]
RUN ["/usr/bin/yum", "distro-sync", "-q", "-y", "--nogpgcheck"]
RUN ["/usr/bin/yum", "update", "-q", "-y","--nogpgcheck"]
RUN ["/usr/bin/yum", "install", "-y", "bwa", "samtools"]
ADD versions.bash /usr/local/bin/versions.bash
CMD ["/usr/local/bin/versions.bash"]
