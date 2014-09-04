bwa_aligner_agent
==================

This is composed of both a dmlond/google_agent_base wrapper agent for
dmlond/bwa_aligner_raw and a SpreadsheetAgent::Runner
script which monitors the configured spreadsheet worksheet
for new subsets that are ready to align, and runs the bwa_aligner_agent
on them.

Configuring dmlond/bwa_aligner_agent
-

Before you can build dmlond/bwa_aligner_agent (see next for
instructions), you must first create a file called
agent.conf.yml in the bwa_aligner_agent build directory.
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

Note, if you are wanting to dmlond/bwa_aligner_agent in conjunction with
dmlond/split_agent, you should make sure that the agent.conf.yml
file is the same for both before building them.

building dmlond/bwa_aligner_agent
-

To build dmlond/bwa_aligner_agent, you have to download dmlond/bwa_aligner_raw from Dockerhub, or build it
from scratch.  If you download dmlond/bwa_aligner_raw, you can 'tag' it with the tag 'dmlond/google_agent_candidate'.
If you build it from scratch, you can either build it with '-t dmlond/google_agent_candidate'
or build it with '-t dmlond/bwa_aligner_raw', and then tag it with 'dmlond/google_agent_candidate'.
Then build dmlond/google_agent_base to produce an image tagged 'dmlond/google_bwa_aligner_agent_candidate'
Then you can build dmlond/bwa_aligner_agent.

Here are the steps starting from a download of bwa_aligner_raw from the dockerhub.  When you
run dmlond/bwa_aligner_raw, docker will download it from Dockerhub, and then it will run and print
the usage instructions for bwa_aligner_raw.  It will then be available to tag appropriately.

```bash
$ sudo docker run dmlond/bwa_aligner_raw
$ sudo docker tag dmlond/bwa_aligner_raw dmlond/google_agent_candidate
$ sudo docker build -t dmlond/google_bwa_aligner_agent_candidate google_agent_base
$ sudo docker build -t dmlond/bwa_aligner_agent bwa_aligner_agent
```

You should consult the README for docker_bwa_aligner for instructions on how
to build bwa_aligner_raw from scratch.

running dmlond/bwa_aligner_agent
-

You should run dmlond/bwa_aligner_agent in 'deamon' mode with the -d switch,
and provide the volumes that are required by dmlond/bwa_aligner.  The
runner will run continuously, waiting for new subsets to be made ready in the
configured spreadsheet.  To align a subset against a specific build and reference,
fill in the 'build' and 'reference' fields for the subset(s) based on the arguments
you would pass to dmlond/bwa_aligner.  Enter '1' in the 'ready' field, and the runner will
at some point pick up the entry and attempt to run the bwa_aligner_agent script on it.

You can run more than one dmlond/bwa_aligner_agent, sharing the same volumes, to run multiple
alignments in parallel.

While an agent is running, you will see 'r:containerid' in the 'bwa_aligner' field.  At some
point this will change to '1' if it completes, or 'f:contiainerid' if it fails. Sometimes,
the agent completes or fails, but cannot update the spreadsheet.  You should monitor the
logs of each running agent container, and determine whether a run of an agent completed or
failed, and manually update the spreadsheet appropriately.  You can supply the 'containerid'
in the field to docker logs or docker inspect to find out information about the specific
instance of the agent.


