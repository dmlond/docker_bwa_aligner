#!/usr/local/bin/ruby

require 'spreadsheet_agent'

usage = "bwa_aligner_agent.rb subset_name"
subset_name = ARGV[0]
unless subset_name
  abort usage
end

$stderr.puts "ALIGNING SUBSET #{ subset_name }"

agent = SpreadsheetAgent::Agent.new(
  agent_name: 'bwa_aligner',
  page_name: 'alignment',
  debug: true,
  keys: {"subset" => subset_name},
  max_selves: 1,
  config_file: '/home/bwa_user/agent_conf/agent.conf.yml'
)

agent.process! do |entry|
  no_errors = true
  begin
    command = [
      "/usr/local/bin/bwa_aligner.pl",
      '-s', entry['subset'],
      '-b', entry['build'],
      '-R', entry['reference']
    ]
    $stderr.puts "RUNNING #{ command.join(' ') }"
    unless system *command
      raise "Could not run #{ command }"
    end
  rescue Exception => error
    $stderr.puts "THERE WAS AN ERROR #{ error.message }"
    no_errors = false
  end

  no_errors
end
exit
