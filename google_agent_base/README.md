google_agent_base
==================

Intermediate image specification designed to facilitate the creation of a wrapper application around a normal 
docker application which uses the ruby spreadsheet_agent gem to turn the application into a spreadsheet_agent
in a subsumption architecture(1).  It installs ruby 2.1.2 and spreadsheet_agent.  See the documentation
for spreadsheet_agent(1,2,3) for more information.

To turn a normal application into a spreadsheet_agent:

1. build the normal application starting from centos:centos6 or centos:centos7
with epel-release-6.8.noarch.rpm, using the title 'dmlond/google_agent_candidate'

```bash
$ sudo docker build -t dmlond/google_agent_candidate split_raw/
```

Alternative, if the image of the application that you want to wrap with a google_agent
is already built, you can simply tag it with this title.

```bash
$ sudo docker tag dmlond/split_raw dmlond/google_split_agent_candidate
```

2. create a build context and Dockerfile for the wrapper application which ADDs a
wrapper spreadsheet_agent ruby script and a config file with google connection information.
It can start FROM any title you wish, but it is recommended that the expected image
from which it starts denote that the image is a google_agent_candidate targetted to 
this application.  For example, the split_agent application starts FROM
'google_split_agent_candidate'


3. build google_agent_base with the title targetted to become the agent application
you created the build context for in 2 above.

```bash
$ sudo docker build -t dmlond/google_split_agent_candidate google_agent_base
```

3. build your agent application

```bash
$ sudo docker build -t dmlond/split_agent split_agent

See split_agent and bwa_aligner_agent for examples.

References
-
1. [spreadsheet_agent gem](https://rubygems.org/gems/spreadsheet_agent)
2. [spreadsheet_agent github](https://github.com/dmlond/spreadsheet_agent)
3. [spreadsheet_agent documentation](http://rubydoc.info/gems/spreadsheet_agent)

License
-------
Wrapper context to create a spreadsheet_agent wrapper script from some other application
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

