#!/usr/bin/perl
use strict;

my $usage = $0." sanger_raw_seq build reference [-p sanger_pair_seq]\naligns both sanger_raw_seq and sanger_pair_seq as paired end if -p is provided\nThis container requires a data volume mapped to the container /data directory, using the --volumes-from or -v docker run flags. This directory must contain the sanger_raw_seq and sanger_pair_seq files.  It will produce a bam file in /data named for the sanger_raw_seq with .bam appended.  It will also produce intermediate sai and unsorted bam files in this directory, but will remove the unsorted bam file before exiting if everything runs smoothly.  Make sure /data maps to adquate storage.\nThe container also requires a reference volume mapped to the container /bwa_indexed directory, either using --volumes-from or -v.  Within this directory, you should organize your bwa_indexed reference fasta files into subdirectories named for the 'build' argument provided.  The aligner expects to align its input to /bwa_indexed/build/reference\n";
my $raw_seq = shift or die $usage;
my $build = shift or die $usage;
my $reference = shift or die $usage;
my $type = shift;
my $pair_seq;
if ($type) {
  $pair_seq = shift or die $usage;
}
my $raw_seq_path = "/data/${raw_seq}";
unless (-f $raw_seq_path) {
  print STDERR "Cannot find ${raw_seq_path}\nPerhaps you need to include a DATA volume?\n";
  exit(1);
}

my $pair_seq_path;
if ($pair_seq) {
  $pair_seq_path = "/data/${pair_seq}";
  unless (-f $pair_seq_path) {
    print STDERR "Paired end ${pair_seq_path} not found\n";
    exit(1);
  }
}

my $reference_path = "/bwa_indexed/${build}/${reference}";
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

if ($type) {
  paired();
} else {
  single();
}
exit;

sub single {
  my $output_sai = "${raw_seq_path}.sai";
  my $output_bam = "${raw_seq_path}.bam";
  print STDERR "ALIGNING\n";
  `bwa aln ${reference_path} ${raw_seq_path} > ${output_sai}`;
  if ($?) {
    print STDERR "PROBLEM ALIGNING $!\n";
    exit(1);
  }
  print STDERR "GENERATING UNSORTED BAM\n";
  `bwa samse ${reference_path} ${output_sai} ${raw_seq_path} | samtools view -bt ${reference_path}.fai - > ${output_bam}.unsorted`;
  if ($?) {
   print STDERR "Could not generate bam file $!\n";
   exit(1);
  }
  sort_bam($output_bam);
}

sub paired {
  my $output_pair1_sai = "${raw_seq_path}.sai";
  my $output_pair2_sai = "${pair_seq_path}.sai";
  my $output_bam = "${raw_seq_path}.bam";

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
  `bwa sampe ${reference_path} ${output_pair1_sai} ${output_pair2_sai} ${raw_seq_path} ${pair_seq_path} | samtools view -bt ${reference_path}.fai - > ${output_bam}.unsorted`;
  if ($?) {
   print STDERR "Could not generate unsorted bam file $!\n";
   exit(1);
  }
  sort_bam($output_bam);
}

sub sort_bam {
  my $output_bam = shift;
  print STDERR "SORTING\n";
  `samtools sort ${output_bam}.unsorted ${output_bam}`;
  if ($?) {
    print STDERR "Could not generate sorted bam file\n";
    exit(1);
 }
 #STUPID SAMTOOLS!!!!
 `mv ${output_bam}.bam ${output_bam}`;
 `rm ${output_bam}.unsorted`;
 system('samtools', 'flagstat', ${output_bam});
}

