FROM centos:centos6
MAINTAINER Darin London <darin.london@duke.edu>

RUN ["/usr/sbin/useradd", "bwa_user"]
USER bwa_user
RUN ["mkdir", "-p", "/home/bwa_user/bwa_indexed"]
VOLUME ["/home/bwa_user/bwa_indexed"]
CMD ["/bin/true"]
