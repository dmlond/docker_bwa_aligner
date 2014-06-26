docker_bwa_aligner
==================

This is a proof of concept docker container based system to run a simple sequence alignment workflow.
The basic parts of this system are:

1. a data volume which exports a directory /data to other containers.
2. a reference volume which exports a dictory /bwa_indexed to other containers
3. a docker containerized application which mounts the reference volume and
downloads a publicly available fasta reference genome file into a subdirectory of /bwa_indexed
named for a specific build, and then indexes the fasta file with bwa and samtools.  To save space
the application can gzip the reference fasta before indexing it. Multiple builds of different genomes
can be stored in the same reference container.
4. a docker containerized application which mounts the data volume and reference volume, and
runs an alignment on single-end or paired end reads in /data against a specific build and bwa indexed
fasta file in /bwa_indexed to produce a sorted bam file named for the read (or the first paired read) in
/data.

Running the workflow
-

This is a simple workflow that aligns publicly available raw sequence reads from the Genome Epidemiology Network (1)
against the publicly available reference Plasmodium falciparum genome (2).  You can run the following, with or without the Build
step discussed below.  If you do not build the containers manually, the docker system will automatically pull the images
from the public Dockerhub repositories mentioned below (mac osx users with boot2docker do not need sudo).

```bash
sudo docker run –name bwa_references dmlond/bwa_reference_volume
sudo docker run -i –volumes-from bwa_references dmlond/bwa_reference -i pf3D7_v2.1.5 ftp://ftp.sanger.ac.uk/pub/pathogens/Plasmodium/falciparum/3D7/3D7.version2.1.5/Pf3D7_v2.1.5.fasta -z
sudo docker run –name plasmodium_data dmlond/bwa_plasmodium_data
sudo docker run -i –volumes-from bwa_references –volumes-from plasmodium_data dmlond/bwa_aligner ERR022523_1.fastq.gz pf3D7_v2.1.5 Pf3D7_v2.1.5.fasta.gz -p ERR022523_2.fastq.gz
mkdir ~/archive
sudo docker run –rm –volumes-from plasmodium_data -v /home/${USER}:/archive dmlond/bwa_samtools_base cp /data/ERR022523_1.fastq.gz.bam /archive/
```

You can run the dmlond/bwa_reference and dmlond/bwa_aligner containers without arguments (or mounted volumes) to get a list of requirements

```bash
sudo docker run dmlond/bwa_reference
sudo docker run dmlond/bwa_aligner
```

Build
-
To build your own version of these containers, you need to install docker.io, clone this repository into docker_bwa_aligner, and
run the following commands (mac osx users with boot2docker do not need sudo):

```bash
sudo docker build -t dmlond/bwa_samtools_base docker_bwa_aligner/
sudo docker build -t  dmlond/bwa_plasmodium_data docker_bwa_aligner/bwa_plasmodium_data
sudo docker build -t dmlond/bwa_reference_volume docker_bwa_aligner/bwa_reference_volume
sudo docker build -t dmlond/bwa_reference docker_bwa_aligner/bwa_reference
sudo docker build -t dmlond/bwa_aligner docker_bwa_aligner/bwa_aligner
```

It is possible to build these containers with a namespace other than dmlond. To do so, you will need
to modify the Dockerfiles in bwa_reference and bwa_aligner and change 'FROM dmlond/bwa_samtools_base'
to the name of your bwa_samtools_base image.  All of the other containers can be built with different names
without affecting the workflow.

DockerHub Repositories
-

All of these repositories are publicly available from dockerhub

* [dmlond/bwa_samtools_base](https://registry.hub.docker.com/u/dmlond/bwa_samtools_base)
* [dmlond/bwa_plasmodium_data](https://registry.hub.docker.com/u/dmlond/bwa_plasmodium_data)
* [dmlond/bwa_reference_volume](https://registry.hub.docker.com/u/dmlond/bwa_reference_volume)
* [dmlond/bwa_reference](https://registry.hub.docker.com/u/dmlond/bwa_reference)
* [dmlond/bwa_aligner](https://registry.hub.docker.com/u/dmlond/bwa_aligner)

Instead of building your own images using this repository, you can simply run the workflow with docker, and it will automatically pull down
these images at run time.


References
-
1. [Genome Epidemiology Network](http://www.malariagen.net/data)
2. [Plasmodium falciparum reference genome] (http://www.nature.com/nature/journal/v419/n6906/abs/nature01097.html)

License
-------
Docker containers that demonstrate a proof of concept bwa alignment workflow
Copyright (c) 2014, Duke University
All rights reserved. Darin London

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the {organization} nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

