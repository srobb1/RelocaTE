#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;

## takes 2 mate files and makes sure that the files are matched. if they are not matched, this
## script cleans them up and prints unmatched reads to STDOUT

##$file_1 and _2 should be mate fq files
my $file_1;
my $file_2;

GetOptions(
  '1|file_1:s' => \$file_1,
  '2|file_2:s' => \$file_2,
  'h|help'     => \&getHelp,
);

sub getHelp {
  print "
This script takes 2 paired fastq files and matches mates and prints unPaired reads to STDOUT. 
Please name the fastq files with .fq extention.
usage:
./clean_pairs.pl [-1 fastq file 1] [-2 fastq file 2][-h] 

options:
-1 STR          fq file 1 [required, no default]
-2 STR          fq file 2 [required, no default]
-h              this message
";
  exit 1;
}

if ( !defined $file_1 and !defined $file_2 ) {
  print "\n\nMust provide a two files containing paired reads\n\n";
  &getHelp;
  exit 1;
}

my $file_1_path = File::Spec->rel2abs($file_1);
my $file_2_path = File::Spec->rel2abs($file_2);

my $sample;
##sample_1.fq
##sampe_pair1.fq
if ( $file_1 =~ /(\S+?)_?(\S+)\.(fq|fastq)/ ) {
  $sample = $1;
}

my %pairs;

open INFASTQ_1, "$file_1" or die "problem opening $file_1 $!";
my $fq_1 = $file_1;
my $fq_2 = $file_2;
$fq_1 =~ s/(\.fq|\.fastq)$//;
$fq_2 =~ s/(\.fq|\.fastq)$//;
open OUTFASTQ_1, ">$fq_1.matched.fq";
open OUTFASTQ_2, ">$fq_2.matched.fq";

while ( my $header = <INFASTQ_1> ) {
  chomp $header;
  my $header_to_store;
  if ( $header =~ /(\S+)(\.r|\.f)\s*/ ) {
    $header_to_store = $1;
  }
  elsif ( $header =~ /(\S+)(\/1|\/2)\s*/ ) {
    $header_to_store = $1;
  }
  elsif ( $header =~ /(\S+)/ ) {
    $header_to_store = $1;
  }
  else {
    $header_to_store = $header;
  }
  my $seq = <INFASTQ_1>;
  warn "$file_1 at record $header has a bad format\n" if !defined $seq;
  chomp $seq;
  
  my $qual_header = <INFASTQ_1>;
  warn "$file_1 at record $header has a bad format\n" if !defined $qual_header;
  chomp $qual_header;

  my $qual = <INFASTQ_1>;
  warn "$file_1 at record $header has a bad format\n" if !defined $qual;
  chomp $qual;
  $seq = "$header\n$seq\n$qual_header\n$qual";
  $pairs{$header_to_store} = $seq;

}
if ( $file_2 !~ /(fq|fastq)$/ ) {
  die "$file_2 is not a fastq file\n";
}
open INFASTQ_2, "$file_2" or die $!;
my $no_match_file_to_print;

while ( my $header = <INFASTQ_2> ) {
  chomp $header;
  my $header_to_store;
  if ( $header =~ /(\S+)(\.r|\.f)\s*/ ) {
    $header_to_store = $1;
  }
  elsif ( $header =~ /(\S+)(\/1|\/2)\s*/ ) {
    $header_to_store = $1;
  }
  elsif ( $header =~ /(\S+)/ ) {
    $header_to_store = $1;
  }
  else {
    $header_to_store = $header;
  }
  my $seq_2 = <INFASTQ_2>;
  warn "$file_2 at record $header has a bad format\n" if !defined $seq_2;
  chomp $seq_2;

  my $qual_header = <INFASTQ_2>;
  warn "$file_2 at record $header has a bad format\n" if !defined $qual_header;
  chomp $qual_header;

  my $qual = <INFASTQ_2>;
  warn "$file_2 at record $header has a bad format\n" if !defined $qual;
  chomp $qual;
  $seq_2 = "$header\n$seq_2\n$qual_header\n$qual";

  if ( exists $pairs{$header_to_store} ) {
    my $seq_1 = $pairs{$header_to_store};
    print OUTFASTQ_1 "$seq_1\n";
    print OUTFASTQ_2 "$seq_2\n";

#remove the found keys, so that only the unparied from the first file are in hash at end of loop
    delete $pairs{$header_to_store};
  }
  else {    #print any seq in file2 that was not in file one
    $no_match_file_to_print .= "$seq_2\n";
  }
}

#will save any remaining unpaired seqs to a variable
foreach my $header ( sort keys %pairs ) {
  my $seq = $pairs{$header};
  $no_match_file_to_print .= "$seq\n";
}

##print nomatches to STDOUT
print $no_match_file_to_print if defined $no_match_file_to_print;
