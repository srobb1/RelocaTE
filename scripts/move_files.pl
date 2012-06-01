#!/usr/bin/perl -w
use strict;
use File::Spec;
use File::Basename;

my $dir       = shift;
my $final_dir = shift;
my $dir_path  = File::Spec->rel2abs($dir);

$final_dir = File::Spec->rel2abs($final_dir);
my $lowest_dir = basename($dir_path);

my $CP1 =
  "$dir_path/*$lowest_dir*.fq  $final_dir/fq_split_by_chromosome/$lowest_dir/.";
my $CP2 =
"$dir_path/*$lowest_dir*.sam  $final_dir/sam_split_by_chromosome/$lowest_dir/.";
my $CP3 =
"$dir_path/*$lowest_dir*.sorted.bam*  $final_dir/bam_split_by_chromosome/$lowest_dir/.";

`mv $CP1`;
`mv $CP2`;
`mv $CP3`;

