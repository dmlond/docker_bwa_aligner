MongoAgent Implementation
=========================

This is an implementation of the [docker_bwa_alignment](https://github.com/dmlond/docker_bwa_aligner).
It uses [fig.sh](http://www.fig.sh/) to create the agents and monitor linked to
the bwa_plasmodium_data and bwa_reference_volume volumes, and a mongodb. It
adds the specific Plasmodium reference dataset that the plasmodium_data is
aligned against.  It also encodes the required Environment variables needed
to tie them all together using the same mongodb container into MongoAgents.

---Run the analysis

```bash
$ fig up -d
# wait about 3-5 minutes, use fig logs plasref and fig logs datavol to ensure the
# references and data are downloaded
$ docker ps | grep dmlond/mongo_agent_alignment
# copy the Docker Container ID
$ docker exec -ti $CONTAINERID /bin/bash
container> irb
irb> a = MongoAgent::Agent.new({name: 'alignment_agent', queue: ENV["QUEUE"]})
irb> a.db[a.queue].insert({
  build: 'pf3D7_v2.1.5',
  reference: 'Pf3D7_v2.1.5.fasta.gz',
  raw_file: 'ERR022523_1.fastq.gz',
  agent_name: 'alignment_agent',
  ready:true
})
# you can monitor the process using these commonds
irb> a.db(a.queue).find.count
#  repeat a few times to see how the overall tasks in the queue build up
# monitor the alignment task until it is complete without errors (hopefully)
irb> a.get_tasks({agent_name: 'alignment_agent'}).count
irb> a.get_tasks({agent_name: 'alignment_agent', complete: true}).count
irb> a.get_tasks({agent_name: 'alignment_agent', complete: true, error_encountered: false}).count
# monitor the split_agent task until it is complete without errors
irb> a.get_tasks({agent_name: 'split_agent'}).count
irb> a.get_tasks({agent_name: 'split_agent', complete: true}).count
irb> a.get_tasks({agent_name: 'split_agent', complete: true, error_encountered: false}).count
### do the same with align_subset_agent, then merge_bam_agent
# see the information produced by the merge_agent (you can look at the others as well)
irb> a.get_tasks({agent_name: 'merge_bam_agent', complete: true}).first
irb> exit
container> exit
$
```
