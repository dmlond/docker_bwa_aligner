FROM dmlond/google_split_agent_candidate

ADD split_agent.rb /usr/local/bin/split_agent.rb
RUN ["chmod", "777", "/usr/local/bin/split_agent.rb"]
RUN ["mkdir", "-p", "/home/bwa_user/agent_conf"]
ADD agent.conf.yml /home/bwa_user/agent_conf/agent.conf.yml
USER bwa_user
ENTRYPOINT ["/usr/local/bin/split_agent.rb"]
#ENTRYPOINT ["/bin/bash"]
#CMD []



