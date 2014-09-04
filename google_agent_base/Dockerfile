#this must start from dmlond/google_agent_candidate
# which itself must be an ancestor in some way
# from centos:centos6 or centos:centos7
# with epel-release-6.8.noarch.rpm
#
# it is designed to act as an intermediary
# between a non google_agent application wanting
# to become a google_agent wrapped application
# so should start from dmlond/google_agent_candidate
# and produce dmlond/google_agent_candidate which
# can then serve as the basis for an agent container
FROM dmlond/google_agent_candidate
MAINTAINER Darin London <darin.london@duke.edu>
USER root
RUN ["/usr/bin/yum", "install", "-y", "--nogpgcheck", "libyaml", "libyaml-devel", "tar", "make", "gcc", "readline", "readline-devel", "openssl","openssl-devel","libxml2-devel","libxslt","libxslt-devel"]
ADD http://cache.ruby-lang.org/pub/ruby/2.1/ruby-2.1.2.tar.gz /root/
WORKDIR /root
RUN ["tar", "-zxf", "ruby-2.1.2.tar.gz"]
WORKDIR /root/ruby-2.1.2
RUN ["./configure", "--enable-shared", "--disable-install-doc"]
RUN ["make"]
RUN ["make", "install"]
WORKDIR /root/ruby-2.1.2/ext/readline
RUN ["/usr/local/bin/ruby", "extconf.rb"]
RUN ["make"]
RUN ["make", "install"]
WORKDIR /root/ruby-2.1.2/ext/zlib
RUN ["/usr/local/bin/ruby", "extconf.rb"]
RUN ["make"]
RUN ["make", "install"]
WORKDIR /root/ruby-2.1.2/ext/openssl
ENV top_srcdir /root/installs/ruby-2.1.2
RUN ["make"]
RUN ["make", "install"]
WORKDIR /root
RUN ["rm", "-rf", "/root/ruby-2.1.2*"]
RUN ["gem", "install", "nokogiri", "--", "--use-system-libraries"]
RUN ["gem", "install", "spreadsheet_agent"]
