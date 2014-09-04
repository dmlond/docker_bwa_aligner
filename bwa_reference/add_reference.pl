#!/usr/bin/perl
use strict;
use File::Basename;

my $usage = $0.' [-i or -b] build fasta_url [-z]'."\n-i uses is bwa index, -b uses bwtsw index\nif -z is specified, the fasta is gzipped after being downloaded\n";
unless (-d '/home/bwa_user/bwa_indexed') {
  print STDERR "cannot find /home/bwa_user/bwa_indexed directory, perhaps you need to include the REFERENCE volume?\n";
  exit(1);
}
my %types = (
  '-i' => 'is',
  '-b' => 'bwtsw'
);

my $type = shift or die $usage;
die $usage unless ($types{$type});
my $build = shift or die $usage;
my $fasta_url = shift or die $usage;
my $gzip = shift;

my $root = join('/', '/home/bwa_user', 'bwa_indexed', $build);
my $file_path = join('/', $root, File::Basename::basename($fasta_url));

`mkdir -p ${root}`;
if ($?) {
  print STDERR "Could not create directory ${root} $!\n";
  exit(1);
}
`wget -O ${file_path} ${fasta_url}`;  
if ($?) {
  print STDERR "Could not get ${fasta_url} to ${file_path} $!\n";
  exit(1);
}
if ($gzip) {
  `gzip ${file_path}`;
  if ($?) {
    print STDERR "Could not gzip ${file_path} $!\n";
    exit(1);
  }
  $file_path .= '.gz'
}

`bwa index -a ${ types{$type} } ${file_path}`;
if ($?) {
  print STDERR "Could not index ${file_path} to ${types{$type}}\n";
  exit(1);
}
`samtools faidx ${file_path}`;
if ($?) {
  print STDERR "Could not samtools faidx ${file_path} $!\n";
  exit(1);
}
exit;

