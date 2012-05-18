#!/usr/bin/perl -w
use strict;
##blat parser and fq trimmer
## April 21, 2012: changed the behavior of more than one blat hit to same query seq.
## if on the same strand keep best TE hit
## if on diff strands, this is a likely tandem insert, throw out reads, make a list
## of potential tandem insert reads

my $blat_file          = shift;
my $fq_file_1          = shift;
my $len_cutoff         = shift;
my $mismatch_allowance = shift;

open INFQ, $fq_file_1 or die $!;
my @blat_path = split '/', $blat_file;
my $filename  = pop @blat_path;
my $blat_path = join '/', @blat_path;
pop @blat_path;
my $te_path = join '/', @blat_path;
my $out_fq_path = "$te_path/te_containing_fq";
my $out_fa_path = "$te_path/te_only_read_portions_fa";


open (INBLAT, $blat_file) or die "Please provide a blat output file\n";
open (OUTTANDEM,  ">$out_fq_path/$filename.potential_tandemInserts_containing_reads.list.txt") or die "Can't open $filename.potential_tandemInserts_containing_reads.list.txt or writing\n";

<INBLAT>;    #get rid of blat header lines
<INBLAT>;
<INBLAT>;
<INBLAT>;
<INBLAT>;

my %coord;
while ( my $line = <INBLAT> ) {
  chomp $line;
  my @blat     = split /\t/, $line;
  my $match    = $blat[0];
  my $mismatch = $blat[1];
  my $strand   = $blat[8];
  my $qName    = $blat[9];
  my $qLen     = $blat[10];
  my $qStart   = $blat[11];
  my $qEnd   = $blat[12] - 1; #get all values into 1st base = 0 postion notation
  my $tLen   = $blat[14];
  my $tStart = $blat[15];
  my $tEnd   = $blat[16] - 1; #get all values into 1st base = 0 postion notation
  my $block_qStarts = $blat[19];
  my ($block_qStart) = split ',', $block_qStarts;
  my $addRecord = 0;
  if (exists $coord{$qName}){
    if ($strand eq $coord{$qName}{strand}){ 
      ##if there is are tandem insertions, theses reads
      ##will call many false insertions events.
      ##if on the same strand, the matches are not likely to overlap, 
      ##the two regions would be separated by mismatches or gaps if they were not
      ##tandem insetions
      print OUTTANDEM "$qName\n";
      $qStart = $qStart <  $coord{$qName}{start} ? $qStart :  $coord{$qName}{start}; 
      $qEnd = $qEnd >  $coord{$qName}{end} ? $qEnd :  $coord{$qName}{end}; 
      $tStart = $tStart <  $coord{$qName}{tStart} ? $tStart :  $coord{$qName}{tStart}; 
      $tEnd = $tEnd >  $coord{$qName}{tEnd} ? $tEnd :  $coord{$qName}{tEnd};
      $match = ($coord{$qName}{match} + $match) > $qLen ? $qLen : $coord{$qName}{match} + $match ;
      $addRecord = 1;
    }else {
      ##keep the best match to TE
      if ($coord{$qName}{match} >= $match ){
        $addRecord = 0;
      }else{
        $addRecord = 1;
      }
    }
  }else {
    $addRecord = 1;
  }
  if ($addRecord) {
    $coord{$qName}{match}    = $match;
    $coord{$qName}{len}      = $qLen;
    $coord{$qName}{start}    = $qStart;
    $coord{$qName}{end}      = $qEnd;
    $coord{$qName}{tLen}     = $tLen;
    $coord{$qName}{mismatch} = $mismatch;
    $coord{$qName}{strand}   = $strand;
    $coord{$qName}{tStart}   = $tStart;
    $coord{$qName}{tEnd}     = $tEnd;
  }
}


my $TE = "unspecified";
my $FA = "unspecified";

#$fa.te_$TE.blatout
if ( $filename =~ /(\S+)\.te_(\S+)\.blatout/ ) {
  $FA = $1;
  $TE = $2;
}
##checking to see if these files already exist
my $outfq       = 0;
my $outte5      = 0;
my $outte3      = 0;
if ( -e "$out_fq_path/$FA.te_$TE.ContainingReads.fq" and -s "$out_fq_path/$FA.te_$TE.ContainingReads.fq") {
  $outfq = 1;
}
if ( -e "$out_fa_path/$FA.te_$TE.five_prime.fa" and -s "$out_fa_path/$FA.te_$TE.five_prime.fa") {
  $outte5 = 1;
}
if ( -e "$out_fa_path/$FA.te_$TE.three_prime.fa" and -s "$out_fa_path/$FA.te_$TE.three_prime.fa") {
  $outte3 = 1;
}
open OUTFQ,  ">$out_fq_path/$FA.te_$TE.ContainingReads.fq" if !$outfq;
open OUTTE5, ">$out_fa_path/$FA.te_$TE.five_prime.fa"      if !$outte5;
open OUTTE3, ">$out_fa_path/$FA.te_$TE.three_prime.fa"     if !$outte3;
##go thru each fq record in the fq files. if the name of the seq is in the blat file
##trim the seq
while ( my $line = <INFQ> ) {
  chomp $line;
  my $header = $line;
  $header = substr( $header, 1 );
  $header =~ s/^(\S+)( .*)/$1/;
  chomp( my $seq        = <INFQ> );
  chomp( my $qualHeader = <INFQ> );
  chomp( my $qual       = <INFQ> );

  if ( exists $coord{$header} ) {
    my $start    = $coord{$header}{start};
    my $len      = $coord{$header}{len};
    my $end      = $coord{$header}{end};
    my $tStart   = $coord{$header}{tStart};
    my $tEnd     = $coord{$header}{tEnd};
    my $tLen     = $coord{$header}{tLen};
    my $mismatch = $coord{$header}{mismatch};
    my $match    = $coord{$header}{match};
    my $strand   = $coord{$header}{strand};

    #want to cut and keep anything not matching to database TE
    my ( $trimmed_seq, $trimmed_qual );

    #query read overlaps 5' end of database TE & trimmed seq > cutoff
    if (
      $tStart == 0
      and ( ( $len - ( $match + $mismatch ) ) > $len_cutoff )
      and ( $mismatch / $len ) <= $mismatch_allowance
      )
    {
      my ( $tS, $tE, $qS, $qE ) =
        ( $tStart + 1, $tEnd + 1, $start + 1, $end + 1 );
      ## $te_subseq = portion of the seq that matches to TE
      ## $trimmed_seq = portion of the seq that does not match to TE
      my $te_subseq = substr( $seq, $start, ( $end - $start + 1 ) );
      if ( $strand eq '-' ) {
        ( $te_subseq = reverse $te_subseq ) =~ tr/AaGgTtCcNn/TtCcAaGgNn/;
        ## start at the end of the match and go to end of string
        $trimmed_seq  = substr $seq,  $end + 1;
        $trimmed_qual = substr $qual, $end + 1;
      }
      else {    ## strand is positive
        $trimmed_seq  = substr $seq,  0, $start;
        $trimmed_qual = substr $qual, 0, $start;
      }
      next if length $trimmed_seq < $len_cutoff;
      print OUTTE5
">$header $qS..$qE matches $TE:$tS..$tE mismatches:$mismatch\n$te_subseq\n"
        if !$outte5;
    }
    #query read overlaps 3' end of database TE & trimmed seq > cutoff
    elsif ( ( $tEnd == $tLen - 1 )
      and ( $len - ( $match + $mismatch ) > $len_cutoff )
      and ( $mismatch / $len ) <= $mismatch_allowance
      ) 
    {
      my ( $tS, $tE, $qS, $qE ) =
        ( $tStart + 1, $tEnd + 1, $start + 1, $end + 1 );
      my $te_subseq = substr( $seq, $start, ( $end - $start + 1 ) );
      if ( $strand eq '-' ) {
        $trimmed_seq  = substr $seq,  0, $start;
        $trimmed_qual = substr $qual, 0, $start;
        ( $te_subseq = reverse $te_subseq ) =~ tr/AaGgTtCcNn/TtCcAaGgNn/;
      }
      else {    ## strand is +
        $trimmed_seq  = substr $seq,  $end + 1;
        $trimmed_qual = substr $qual, $end + 1;
      }
      next if length $trimmed_seq < $len_cutoff;
      print OUTTE3
">$header $qS..$qE matches $TE:$tS..$tE mismatches:$mismatch\n$te_subseq\n"
        if !$outte3;
    }
    if ( defined $trimmed_seq ) {
      print "\@$header\n";
      print "$trimmed_seq\n";
      print "$qualHeader\n";
      print "$trimmed_qual\n";
    }
    ##TE containing reads
    ##any read that was in the blat file is written here
    ##in te_containing.fq file
    if ( !$outfq ) {
      print OUTFQ "\@$header\n";
      print OUTFQ "$seq\n";
      print OUTFQ "$qualHeader\n";
      print OUTFQ "$qual\n";
    }
  }
}
