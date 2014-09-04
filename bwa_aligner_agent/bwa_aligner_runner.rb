#!/usr/local/bin/ruby

require 'spreadsheet_agent'

@runner = SpreadsheetAgent::Runner.new(
  agent_bin: "/usr/local/bin",
  config_file: '/home/bwa_user/agent_conf/agent.conf.yml'
)

@runner.skip_entry do |entry|
  $stderr.puts "CHECKING #{ entry.inspect }"
  if entry['ready'] == "1"
    if entry['bwa_aligner'].nil? || entry['bwa_aligner'].empty?
      false
    else
      true
    end
  else
    true
  end
end

@max_exceptions = 5
@concurrent_exceptions = 0
while @concurrent_exceptions < @max_exceptions

  begin
    @runner.process! do |entry, page|
      goal_agent = [@runner.agent_bin, "bwa_aligner_agent.rb"].join('/')
      if File.executable? goal_agent
        cmd = [goal_agent]

        @runner.query_fields.each do |query_field|
          if entry[query_field]
            cmd.push entry[query_field]
          end
        end
        command = cmd.join(' ')
        command += '&'
        system command
      else
        raise "#{ goal_agent } not executable!"
      end

      @concurrent_exceptions = 0
      sleep 31
    end
  rescue Exception => e
    $stderr.puts "EXCEPTION #{ e.message }"
    @concurrent_exceptions = @concurrent_exceptions + 1
    $stderr.puts "#{ @concurrent_exceptions } have been experienced"
  end

end
