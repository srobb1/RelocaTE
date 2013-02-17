#!/usr/bin/perl -w
use strict;
##converts fq to fa

if (!defined @ARGV){
  die "Do not run directly, to be called by relocaTE.pl\n";
}

my ( $infq, $outfa ) = ( $ARGV[0], $ARGV[1] );

if ( $infq !~ /\.(fq|fastq)$/ ) {
  die
"Was expecting a fq file to convert to a fa file, but I have this instead $infq\n";
}

open INFQ,  $infq     or die "Can't open $infq $!";
open OUTFA, ">$outfa" or die $!;

while ( my $header = <INFQ> ) {
  my $seq         = <INFQ>;
  my $qual_header = <INFQ>;
  my $qual        = <INFQ>;
  if ( substr( $header, 0, 1 ) ne '@' ) {
    die "ERROR: expected \'\@\' but saw $header";
  }
  print OUTFA ">", substr( $header, 1 );
  print OUTFA $seq;
}
close INFQ;
close OUTFA;
