# This will install bwa 0.5.9-r16 and samtools 0.1.18 (r982:295)
FROM centos:centos6
MAINTAINER Darin London <darin.london@duke.edu>

RUN ["/usr/sbin/useradd", "bwa_user"]
RUN ["/usr/bin/yum", "clean", "all"]
RUN ["/usr/bin/yum", "distro-sync", "-q", "-y", "--nogpgcheck"]
RUN ["/usr/bin/yum", "update", "-q", "-y","--nogpgcheck"]
ADD http://download.fedoraproject.org/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm /root/
RUN ["rpm", "-ivh", "/root/epel-release-6-8.noarch.rpm"] 
RUN ["/usr/bin/yum", "install", "-y", "wget", "bwa", "samtools"]
ADD versions.bash /usr/local/bin/versions.bash
RUN ["chmod", "777", "/usr/local/bin/versions.bash"]
CMD ["/usr/local/bin/versions.bash"]
