FROM dmlond/bwa_aligner_agent_candidate

RUN ["mkdir", "-p", "/home/bwa_user/agent_conf"]
ADD bwa_aligner_agent.rb /usr/local/bin/bwa_aligner_agent.rb
ADD bwa_aligner_runner.rb /usr/local/bin/bwa_aligner_runner.rb
RUN ["chmod", "777", "/usr/local/bin/bwa_aligner_agent.rb", "/usr/local/bin/bwa_aligner_runner.rb"]
ADD agent.conf.yml /home/bwa_user/agent_conf/agent.conf.yml
WORKDIR /home/bwa_user
USER bwa_user
ENTRYPOINT ["/usr/local/bin/bwa_aligner_runner.rb"]
#ENTRYPOINT ["/bin/bash"]
CMD []
