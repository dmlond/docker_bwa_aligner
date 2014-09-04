#!/usr/local/bin/ruby

@subset_col = 1
require 'spreadsheet_agent'
agent = SpreadsheetAgent::Agent.new(
  agent_name: 'split_agent',
  page_name: 'alignment',
  keys: {subset: 'foo'},
  config_file: '/home/bwa_user/agent_conf/agent.conf.yml'
)
@ws = agent.worksheet
@current_row = @ws.num_rows + 1

command = ["/usr/local/bin/split_raw.pl"] + ARGV
@exit_status = 0
IO.popen(command) do |subset_files|
  subset_files.each do |current_subset_filename|
    @ws[@current_row,@subset_col] = "#{ current_subset_filename.chomp }"
    @ws.save
    @current_row = @current_row + 1
  end
  subset_files.close
  @exit_status = $?.exitstatus
end
exit @exit_status
