split_agent
==================

This is a dmlond/google_agent_base wrapper agent for
dmlond/split_raw. It takes the same arguments as 
dmlond/split_raw, but uses the Ruby SpreadsheetAgent::Agent
to populate a worksheet configured with agent.conf.yml with
the entries that are produced by the dmlond/split_raw application.
The subset raw files can then be processed by a variety of other
dmlond/google_agent_base wrapper agents, such as dmlond/bwa_aligner_agent.

Configuring dmlond/split_agent
-

Before you can build dmlond/split_agent (see next for
instructions), you must first create a file called
agent.conf.yml in the split_agent build directory.
For reference, you can use the agent.conf.yml.example
file included. Copy agent.conf.yml.example to agent.conf.yml.
This is a YAML file with the account information
required to connect to the Google Drive system,
and a specific Spreadsheet to interact with.

You must change the following fields:

guser: your google account userid
gpass: your password.

It is highly recommended
that you activate [Google 2-Step Verificaiton](https://www.google.com/landing/2step/)
in google, and then create an [App Password](https://support.google.com/accounts/answer/185833?hl=en)
for your agent to use in this field.  Then, if your App Password gets stolen or compromised,
you can delete it and create a new one.

spreadsheet_name: the name of the spreadsheet in google drive that your agent will use.

This must have a worksheet in it called 'alignment'.   It should have the following column names in the first row:

subset, build, reference, ready, bwa_aligner, complete

send_to: this cannot be null, but it is not used.  The agent runs in foreground mode, and prints errors to
STDOUT and STDERR which can be accessed by the docker logs command on a running agent container.

building dmlond/split_agent
-

To build dmlond/split_agent, you have to download dmlond/split_raw from Dockerhub, or build it
from scratch.  If you download dmlond/split_raw, you can 'tag' it with the tag 'dmlond/google_agent_candidate'.
If you build it from scratch, you can either build it with '-t dmlond/google_agent_candidate'
or build it with '-t dmlond/split_raw', and then tag it with 'dmlond/google_agent_candidate'.
Then build dmlond/google_agent_base to produce an image tagged 'dmlond/google_split_agent_candidate'
Then you can build dmlond/split_agent.

Here are the steps starting from a download of split_raw from the dockerhub.  When you
run dmlond/split_raw, docker will download it from Dockerhub, and then it will run and print
the usage instructions for split_raw.  It will then be available to tag appropriately.

```bash
$ sudo docker run dmlond/split_raw
$ sudo docker tag dmlond/split_raw dmlond/google_agent_candidate
$ sudo docker build -t dmlond/google_split_agent_candidate google_agent_base
$ sudo docker build -t dmlond/split_agent split_agent
```

You should consult the README for docker_bwa_aligner for instructions on how
to build split_raw from scratch.

running dmlond/split_agent
-

You run dmlond/split_agent as you would run dmlond/split_raw.
If your agent.conf.yml is set up correctlly, it will run the
split_raw script, and, with each subset file printed by it to STDOUT,
create a row in the 'aligned' worksheet 'subset' field.  All other
fields are left blank for the row.

