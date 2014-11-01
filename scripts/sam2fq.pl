#!/usr/bin/perl -w
use strict;
use File::Spec;

## converts sam files to fq files

## example usage for multiple directories:
## if directories are Chr1, Chr2, Chr3 --- Chr12
## for i in `seq 1 12` ; do sam2fq.pl Chr$i ; done

## updated on 03/21/2012: realized that if any alignments were unmapped, or unpaired in the sam file
## that was not named unpaired, the sequences were not printed anywhere, they were lost. made a change
## to allow the printing of any unpaired or unmapped seqs to the unpaired.fq file. if the unpaired
## file already exists, it will be appended to. Also changed the output files form _1.fq and _2.fq
## to _p1.fq and _p2.fq

my $dir         = shift;
my $dir_path    = File::Spec->rel2abs($dir);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

##if flag contains 16 it is on the minus strand and is reported as the revcomp in the sam file
my $flag_minus_strand = 16;
#my %flag_minus_strand = (
#  16,  1, 17,  1, 18,  1, 19,  1, 20,  1, 21,  1, 22,  1, 23,  1, 24,  1,
#  25,  1, 26,  1, 27,  1, 28,  1, 29,  1, 30,  1, 31,  1, 48,  1, 49,  1,
#  50,  1, 51,  1, 52,  1, 53,  1, 54,  1, 55,  1, 56,  1, 57,  1, 58,  1,
#  59,  1, 60,  1, 61,  1, 62,  1, 63,  1, 80,  1, 81,  1, 82,  1, 83,  1,
#  84,  1, 85,  1, 86,  1, 87,  1, 88,  1, 89,  1, 90,  1, 91,  1, 92,  1,
#  93,  1, 94,  1, 95,  1, 112, 1, 113, 1, 114, 1, 115, 1, 116, 1, 117, 1,
#  118, 1, 119, 1, 120, 1, 121, 1, 122, 1, 123, 1, 124, 1, 125, 1, 126, 1,
#  127, 1, 144, 1, 145, 1, 146, 1, 147, 1, 148, 1, 149, 1, 150, 1, 151, 1,
#  152, 1, 153, 1, 154, 1, 155, 1, 156, 1, 157, 1, 158, 1, 159, 1, 176, 1,
#  177, 1, 178, 1, 179, 1, 180, 1, 181, 1, 182, 1, 183, 1, 184, 1, 185, 1,
#  186, 1, 187, 1, 188, 1, 189, 1, 190, 1, 191, 1, 208, 1, 209, 1, 210, 1,
#  211, 1, 212, 1, 213, 1, 214, 1, 215, 1, 216, 1, 217, 1, 218, 1, 219, 1,
#  220, 1, 221, 1, 222, 1, 223, 1, 240, 1, 241, 1, 242, 1, 243, 1, 244, 1,
#  245, 1, 246, 1, 247, 1, 248, 1, 249, 1, 250, 1, 251, 1, 252, 1, 253, 1,
#  254, 1, 255, 1
#);

my %dupCheck;
opendir( DIR, $dir_path ) || die "$!";
foreach my $file ( readdir(DIR) ) {
  next unless $file =~ /^(\S+)\.sam$/;
  my $filebase = $1;
  open INSAM, "$dir_path/$file" or die "Can't open $file to convert 2 fq $!";
  my $unpaired = "$dir_path/$filebase" . "_unPaired.fq";
  if ( $filebase =~ /unPaired/ ) {
    $unpaired = "$dir_path/$filebase.fq";
  }
  if ( !-s $unpaired ) {
    open OUTUNPAIRED, ">$unpaired";
  }
  else {
    open OUTUNPAIRED, ">>$unpaired";
  }
  if ( $filebase !~ /unPaired/ ) {
    open OUTFQ_1, ">$dir_path/$filebase" . "_p1.fq";
    open OUTFQ_2, ">$dir_path/$filebase" . "_p2.fq";
  }
  while ( my $line = <INSAM> ) {
    next if $line =~ /^@/;
    chomp $line;
    my @samLine = split /\t/, $line;
    my (
      $name, $flag,  $refName, $three, $four, $five,
      $six,  $seven, $eight,   $seq,   $qual, @rest
    ) = @samLine;
    if ( $flag & $flag_minus_strand ) {
    #if ( exists $flag_minus_strand{$flag} ) {
      my $revseq = reverse uc($seq);
      $revseq =~ tr/ATGC/TACG/;
      $seq = $revseq;
      my $revqual = reverse $qual;
      $qual = $revqual;
    }

    # since we already split the sam files per target
    # do not print read if it already has been printed
    # checking read name plus seq make a unique record
    if ( exists $dupCheck{$name}{$seq} ) {
      next;
    }
    else {
      $dupCheck{$name}{$seq} = 1;
    }
    my $toPrint = "\@$name\n$seq\n\+\n$qual\n";
    if ( $file =~ /unPaired/ ) {
      print OUTUNPAIRED $toPrint;
    }
    elsif ( $flag & 64 ) {    #first mate
    #elsif ( $flag >= 64 and $flag < 128 ) {    #first mate
      print OUTFQ_1 $toPrint;
    }
    elsif ( $flag & 128  ) {    #second mate
    #elsif ( $flag >= 128 and $flag < 192 ) {    #second mate
      print OUTFQ_2 $toPrint;
    }
    elsif ( $name =~ /\.f$/ ) {                 #454 first mate
      print OUTFQ_1 $toPrint;
    }
    elsif ( $name =~ /\.r$/ ) {                 #454 second mate
      print OUTFQ_2 $toPrint;
    }
    elsif ( $name =~ /\.fn$/ ) {                #454 unpaired
      print OUTUNPAIRED $toPrint;
    }
    else {
      print OUTUNPAIRED $toPrint;
    }
  }
  if ( -z $unpaired ) {
    unlink $unpaired;
  }
}

