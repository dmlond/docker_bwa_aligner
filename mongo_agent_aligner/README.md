MongoAgent Implementation
=========================

This is an implementation of the [docker_bwa_alignment](https://github.com/dmlond/docker_bwa_aligner).
It uses [crane](https://github.com/michaelsauter/crane) to make running the pipeline
easier.

It consists of four agents to perform different parts of the pipeline:
  - [mongo_agent_alignment](https://github.com/dmlond/mongo_agent_alignment)
  - [mongo_agent_split_raw](https://github.com/dmlond/mongo_agent_split_raw)
  - [mongo_agent_align_subset](https://github.com/dmlond/mongo_agent_align_subset)
  - [mongo_agent_merge_bam](https://github.com/dmlond/mongo_agent_merge_bam)

It also consists of [mongo_agent_merge_monitor](https://github.com/dmlond/mongo_agent_merge_monitor)
to monitor the align_subset jobs for each alignment, and create new merge_bam
jobs when they all complete without errors.

There are two seed tasks that need to be run which download the raw and reference
files and insert/update the mongodb alignment task for the demo alignment:
  - [mongo_agent_seed_alignment](https://github.com/dmlond/mongo_agent_seed_alignment)
  - [mongo_agent_seed_reference](https://github.com/dmlond/mongo_agent_seed_reference)

Finally, there is the [mongo_agent_watch_alignment](https://github.com/dmlond/mongo_agent_watch_alignment)
application which can be run to watch the progress of the alignment by continuously querying the mongodb
to find information about tasks for the four agents.

---Run the analysis
Make a data directory in the same directory as this README.md file
```
$ mkdir data
```

Then run:
```bash
$ crane lift
$ crane logs -f watcher
```

one can run just the seeds:
```bash
$ crane run seeds
```

or agents:
```
$ crane run agents
```

or run an individual part of the pipeline based on their keys in crane.yml
under the containers: heading.
```
$ crane run alignment
```

to clean up after your run:
```
$ crane stop
$ crane rm
```

you can monitor the downloads with (note, if you add -f it
will continue to follow the log until you hit crtl-c, without it
it just prints the last few lines and exits):
```bash
$ crane logs seedraw
$ crane logs seedref
```
See the crane documentation for more details
