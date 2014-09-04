#!/usr/bin/perl
use strict;
use Getopt::Std;
use File::Basename;

my $usage = q|docker run {data_volume} {reference_volume} dmlond/bwa_aligner {options}

Aligns a fastq file, or set of paired end fastq files, to a reference genome using bwa
and samtools to produce a sorted bam file.

data_volume:  The container requires that the directory /home/bwa_user/data exist, and
  that all file_names passed in with -s and -p be found in this directory.  This can
  be mounted either with the -v or --volumes-from docker run directives.  If mounting
  a data volume container, the container must specify /data as a VOLUME.  If mounting
  with -v, use -v /path/to/hostdir:/data.
  ** The aligner produces an sai file, final bam file, and many temporary bam files
  that are created by samtools when it sorts the final bam file in this directory, 
  so make sure there is enough storage space to handle this based on the size of the
  input files.

reference_volume: The container requires that the directory /home/bwa_user/bwa_indexed
 exist, with a subdirectory named for the build specified with -b.  The reference
 file_name specified with -R must be found in /bwa_indexed/${build}/ along with all
 of the bwa index and samtools faidx index files.  This can be mounted either with
 the -v or --volumes-from docker run directives.  If mounting a data volume
 container, the container must specify /home/bwa_user/bwa_indexed as a VOLUME, and
 it must be owned by bwa_user user and bwa_user group.  If mounting with
 -v, use -v /path/to/hostdir:/bwa_indexed.

options:
  required:
    -s file_name: name of sanger formatted fastq file
      in /data/ (can be gzip compressed or uncompressed)
    -b build: the name of the directory in /bwa_indexed/ where the reference
      file specified with -R is located
    -R file_name: name of the reference fastq file located in
      /bwa_indexed/${build}/. This file must have been indexed with
      both bwa, and samtools faidx, so that all of the accompanying
      index files are located with this file in the same directory.
      (can be gzip compressed or uncompressed)
  optional:
    -p file_name: if reads are paired_end, the second sanger formatted fastq file in the
        pair located in /data/ (can be gzip compressed or uncompressed)
    -o file_name: if specified the bam file produced will be named with this name in
        /data/.  Defaults to the file_name passed in with the -s flag, with .bam appended.
|;

my $opt = {};
Getopt::Std::getopts('s:b:R:p:o:', $opt);

my $raw_seq = $opt->{s} or die $usage;
my $build = $opt->{b} or die $usage;
my $reference = $opt->{R} or die $usage;
my $pair_seq = $opt->{p};
my $final_bam_file = $opt->{o};

my $raw_seq_path = "/home/bwa_user/data/${raw_seq}";
unless (-f $raw_seq_path) {
  print STDERR "Cannot find ${raw_seq_path}\nPerhaps you need to include a DATA volume?\n";
  exit(1);
}

my $pair_seq_path;
if ($pair_seq) {
  $pair_seq_path = "/home/bwa_user/data/${pair_seq}";
  unless (-f $pair_seq_path) {
    print STDERR "Paired end ${pair_seq_path} not found\n";
    exit(1);
  }
}

my $reference_path = "/home/bwa_user/bwa_indexed/${build}/${reference}";
unless (-f $reference_path) {
  print STDERR "Cannot find ${reference_path}\nPerhaps you need to include a REFERENCE volume?\n";
  exit(1);
}

unless (-f "${reference_path}.bwt") {
  print STDERR "Cannot find ${reference_path} bwa index files.\n";
  exit(1);
}

unless (-f "${reference_path}.fai") {
  print STDERR "Cannot find ${reference_path}.fai samtools faidx indexed file.\n";
  exit(1);
}

if ($final_bam_file) {
  $final_bam_file = join('/', '/home/bwa_user/data', $final_bam_file);
} else {
  $final_bam_file = "${raw_seq_path}.bam";
}

if ($pair_seq) {
  paired();
} else {
  single();
}
sort_bam();
exit;

sub single {
  my $output_sai = "${raw_seq_path}.$$.sai";

  print STDERR "ALIGNING\n";
  `bwa aln ${reference_path} ${raw_seq_path} > ${output_sai}`;
  if ($?) {
    print STDERR "PROBLEM ALIGNING $!\n";
    exit(1);
  }
  print STDERR "GENERATING UNSORTED BAM\n";
  `bwa samse ${reference_path} ${output_sai} ${raw_seq_path} | samtools view -bt ${reference_path}.fai - > ${final_bam_file}.unsorted`;
  if ($?) {
   print STDERR "Could not generate bam file $!\n";
   exit(1);
  }
  unlink($output_sai);
}

sub paired {
  my $output_pair1_sai = "${raw_seq_path}.$$.sai";
  my $output_pair2_sai = "${pair_seq_path}.$$.sai";

  print STDERR "ALIGNING PAIR1\n";
  `bwa aln ${reference_path} ${raw_seq_path} > ${output_pair1_sai}`;
  if ($?) {
    print STDERR "PROBLEM ALIGNING PAIR1 $!\n";
    exit(1);
  }

  print STDERR "STDERR ALIGNING PAIR2\n";
  `bwa aln ${reference_path} ${pair_seq_path} > ${output_pair2_sai}`;
  if ($?) {
    print STDERR "PROBLEM ALIGNING PAIR2 $!\n";
    exit(1);
  }

  print STDERR "GENERATING UNSORTED BAM\n";
  `bwa sampe ${reference_path} ${output_pair1_sai} ${output_pair2_sai} ${raw_seq_path} ${pair_seq_path} | samtools view -bt ${reference_path}.fai - > ${final_bam_file}.unsorted`;
  if ($?) {
   print STDERR "Could not generate unsorted bam file $!\n";
   exit(1);
  }
  unlink($output_pair1_sai);
  unlink($output_pair2_sai);
}

sub sort_bam {
  print STDERR "SORTING\n";
  `samtools sort ${final_bam_file}.unsorted ${final_bam_file}`;
  if ($?) {
    print STDERR "Could not generate sorted bam file\n";
    exit(1);
 }
 unlink("${final_bam_file}.unsorted");

 #STUPID SAMTOOLS sort appends .bam!!!!
 `mv ${final_bam_file}.bam ${final_bam_file}`;
 system('samtools', 'flagstat', ${final_bam_file});
}
