#!/usr/bin/perl -w
use strict;
use File::Spec;
use File::Basename;

## converts sam files to bam files using samtools

my $dir         = shift;
my $genome_path = shift;
my $sample      = shift;
my $dir_path    = File::Spec->rel2abs($dir);
my $lowest_dir  = basename($dir_path);

if ( -e "$dir_path/$sample.$lowest_dir.sam" ) {
  `samtools view -b -S -T $genome_path $dir_path/$sample.$lowest_dir.sam >  $dir_path/$sample.$lowest_dir.bam`;
  `samtools sort  $dir_path/$sample.$lowest_dir.bam  $dir_path/$sample.$lowest_dir.sorted`;
  `samtools index  $dir_path/$sample.$lowest_dir.sorted.bam`;
}
else {
  warn
"$dir_path/$sample.$lowest_dir.sam does not exist, therefore cannot convert to bam file";
}

if ( -e "$dir_path/$sample.$lowest_dir.unPaired.sam" ) {
  `samtools view -b -S -T $genome_path $dir_path/$sample.$lowest_dir.unPaired.sam >  $dir_path/$sample.$lowest_dir.unPaired.bam`;
  `samtools sort  $dir_path/$sample.$lowest_dir.unPaired.bam  $dir_path/$sample.$lowest_dir.unPaired.sorted`;
  `samtools index  $dir_path/$sample.$lowest_dir.unPaired.sorted.bam`;
}
