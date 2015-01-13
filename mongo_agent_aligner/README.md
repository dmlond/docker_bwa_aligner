MongoAgent Implementation
=========================

This is an implementation of the [docker_bwa_alignment](https://github.com/dmlond/docker_bwa_aligner).
It uses [crane](https://github.com/michaelsauter/crane) to make running the pipeline
easier.

It consists of four agents to perform different parts of the pipeline:
  - mongo_agent_alignment [github](https://github.com/dmlond/mongo_agent_alignment) [docker hub](https://registry.hub.docker.com/u/dmlond/mongo_agent_alignment)
  - mongo_agent_split_raw [github](https://github.com/dmlond/mongo_agent_split_raw) [docker hub](https://registry.hub.docker.com/u/dmlond/mongo_agent_split_raw)
  - mongo_agent_align_subset [github](https://github.com/dmlond/mongo_agent_align_subset) [docker hub](https://registry.hub.docker.com/u/dmlond/mongo_agent_align_subset)
  - mongo_agent_merge_bam [github](https://github.com/dmlond/mongo_agent_merge_bam) [docker hub](https://registry.hub.docker.com/u/dmlond/mongo_agent_merge_bam)

It also consists of mongo_agent_merge_monitor [github](https://github.com/dmlond/mongo_agent_merge_monitor) [docker hub](https://registry.hub.docker.com/u/dmlond/mongo_agent_merge_monitor)
to monitor the align_subset jobs for each alignment, and create new merge_bam
jobs when they all complete without errors.

There are two seed tasks that need to be run which download the raw and reference
files and insert/update the mongodb alignment task for the demo alignment:
  - mongo_agent_seed_alignment [github](https://github.com/dmlond/mongo_agent_seed_alignment) [docker hub](https://registry.hub.docker.com/u/dmlond/mongo_agent_seed_alignment)
  - mongo_agent_seed_reference [github](https://github.com/dmlond/mongo_agent_seed_reference) [docker hub](https://registry.hub.docker.com/u/dmlond/mongo_agent_seed_reference)

Finally, there is the mongo_agent_watch_alignment [github](https://github.com/dmlond/mongo_agent_watch_alignment) [docker hub](https://registry.hub.docker.com/u/dmlond/mongo_agent_watch_alignment)
application which can be run to watch the progress of the alignment by continuously querying the mongodb
to find information about tasks for the four agents.

---Run the analysis with crane
Make a data directory in the same directory as this README.md file
```
mkdir data
chmod 777 data
```

Then run:
```bash
crane lift
crane logs -f watcher
```

one can run just the seeds:
```bash
crane run seeds
```

or agents:
```
crane run agents
```

or run an individual part of the pipeline based on their keys in crane.yml
under the containers: heading.
```
crane run alignment
```

to clean up after your run:
```
crane stop
crane rm
```

you can monitor the downloads with (note, if you use -f  with logs it
will continue to follow the log until you hit crtl-c, without it
it just prints the last few lines and exits):
```bash
crane logs seedraw
crane logs seedref
```
See the crane documentation for more details


###Raw docker run command version with 3 parallel alignsubset agents

```bash
MONGOID=$(docker run --name mongodb -d mongo:latest mongod --smallfiles)
REFVOLID=$(docker run --name refvol dmlond/bwa_reference_volume)
docker run -d --name seedref -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol dmlond/mongo_agent_seed_reference
SEEDRAWID=$(docker run -d --name seedraw -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb -v $(pwd)/data:/home/bwa_user/data dmlond/mongo_agent_seed_alignment)
docker run -d --name alignment -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol --volumes-from seedraw dmlond/mongo_agent_alignment
docker run -d --name splitraw -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol --volumes-from seedraw dmlond/mongo_agent_split_raw
docker run -d --name alignsubset1 -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol --volumes-from seedraw dmlond/mongo_agent_align_subset
docker run -d --name alignsubset2 -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol --volumes-from seedraw dmlond/mongo_agent_align_subset
docker run -d --name alignsubset3 -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol --volumes-from seedraw dmlond/mongo_agent_align_subset
docker run -d --name mergebam -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol --volumes-from seedraw dmlond/mongo_agent_merge_bam
docker run -d --name mergemonitor -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol --volumes-from seedraw dmlond/mongo_agent_merge_monitor
WATCHERID=$(docker run -d --name watcher -e "MONGO_HOST=mongodb" -e "MONGO_DB=agents" -e "QUEUE=test" --link mongodb:mongodb --volumes-from refvol --volumes-from seedraw dmlond/mongo_agent_watch_alignment)
docker logs -f $WATCHERID
```

to clean up, make sure you remove the volumes along with refvol and seedraw.
```bash
docker rm --volumes=true $REFVOLID
docker rm --volumes=true $SEEDRAWID
```
