#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my $usage = q|USAGE: sudo docker run dmlond/split_raw -f fastq_file [-s size] [-e extension] [-z boolean][-t bp]
 -s: number of entries in each subset fastq file.  default 100000
 -e: the characters to attach to the fastq_file name to create the subset fastq files. default bfq
 -z: whether the input fastq files are compressed (gzip and zip supported). default false. subset files
     will be gzipped whether or not the input is gzipped.  It is highly recommended that the input be compressed.
 -t: if specified, sequence and corresponding quality are truncated to only the first bp basepairs

 This takes a fastq_file in /data, and splits it into gzipped subsets of size entries.  Each subset
 file is stored in /data, and named with the fastq_file name with _subset_number.extension.gz
 appended.
 This container requires a volume to be mounted with a /data directory exported.  This can
 be done with a host directory using the -v docker run flag, or with a data container using
 the --volumes-from container_name docker run flag.
|;

my $data = '/data';
die $usage unless (-d $data);

my $opt = {};
Getopt::Std::getopts('zf:s:e:t:', $opt);

my $fastq_file = $opt->{f} or die $usage;
my $count = $opt->{s};
$count =  100000 unless ($count);
my $ext = $opt->{e};
$ext = 'bfq' unless($ext);
my $compressed = $opt->{z} ? 1 : 0;
my $truncate = $opt->{t};

my $input_file = join('/', $data, $fastq_file);
unless (-f $input_file) {
    die "${fastq_file} not found in ${data}\n";
}
my $prefix = $input_file;
my $fc = 1;
my $c = 0;
my $one_line_seen;
my $command = ($compressed) ? 'zcat': 'cat';

open( my $in_h, '-|', $command, $input_file ) or die "Could not open ${input_file} $!\n";
my $current_file_name = "${prefix}_${fc}.${ext}.gz";
open( my $out_h, '|-', "gzip -c > ${current_file_name}" ) or die $!;
while (my $line = <$in_h>) {
    $one_line_seen = 1;
    if ($c == $count) {
        close $out_h;
        print File::Basename::basename($current_file_name)."\n";
        $fc++;
	$current_file_name = "${prefix}_${fc}.${ext}.gz";
	open($out_h, '|-', "gzip -c > ${current_file_name}") or die $!;
        $c = 0;
    }
    print $out_h $line;
    $line = <$in_h>;
    print $out_h ($truncate) ? substr($line, 0, $truncate)."\n" : $line;
    $line = <$in_h>;
    print $out_h $line;
    $line = <$in_h>;
    print $out_h ($truncate) ? substr($line, 0, $truncate)."\n" : $line;
    $c++;
}
print File::Basename::basename($current_file_name)."\n";
close $in_h;
close $out_h;
unless ($one_line_seen) {
    unlink $current_file_name;
}
exit;
