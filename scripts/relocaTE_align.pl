#!/usr/bin/perl -w
use strict;
## aligns TE trimmed reads to the reference

if ( !defined @ARGV ) {
  die "Do not run directly, to be called by relocaTE.pl\n";
}

my $scripts     = shift;    #full path to scripts directory
my $path        = shift;    #current/top/TE
my $genome_file = shift;
my $regex_file  = shift;
my $TE          = shift;
my $exper       = shift;
my $bowtie2     = shift;
my $bowtie_sam  = 1;
##get the regelar expression patterns for mates and for the TE
##when passed on the command line as an argument, even in single
##quotes I lose special regex characters
open INREGEX, "$regex_file" or die "$!";
my $mate_file_1;
my $mate_file_2;
my $mate_file_unpaired;
my $TSD;

while ( my $line = <INREGEX> ) {
  chomp $line;
  ( $mate_file_1, $mate_file_2, $mate_file_unpaired, $TSD ) = split /\t/, $line;
  $mate_file_1        =~ s/\.(fq|fastq)$//;
  $mate_file_2        =~ s/\.(fq|fastq)$//;
  $mate_file_unpaired =~ s/\.(fq|fastq)$//;
}

my %flanking_fq;
##mate finder
my @files_1;
my @files_2;
my @files_unpaired;
my @flanking_files = <$path/flanking_seq/*flankingReads.fq>;
foreach my $file (@flanking_files) {
  next if -z $file;    ##next file if size is 0
  if ( $file =~ /$mate_file_unpaired/ ) {
    push( @files_unpaired, $file );
  }
  elsif ( $file =~ /$mate_file_1/ ) {
    push( @files_1, $file );
  }
  elsif ( $file =~ /$mate_file_2/ ) {
    push( @files_2, $file );
  }
}
if ( @files_1 and @files_2 ) {
  for ( my $i = 0 ; $i < @files_1 ; $i++ ) {
    my $file_1 = $files_1[$i];
    $file_1 =~ s/$mate_file_1//;
    for ( my $j = 0 ; $j < @files_2 ; $j++ ) {
      my $file_2 = $files_2[$j];
      $file_2 =~ s/$mate_file_2//;
      if ( $file_1 eq $file_2 ) {
        $flanking_fq{$file_1}{1} = $files_1[$i];
        $flanking_fq{$file_1}{2} = $files_2[$j];
        if (@files_unpaired) {
          for ( my $k = 0 ; $k < @files_unpaired ; $k++ ) {
            my $file_unpaired = $files_unpaired[$k];
            if ( $file_1 eq $file_unpaired ) {
              $flanking_fq{$file_1}{unpaired} = $files_unpaired[$k];
              last;
            }
          }
        }

        #if $file_1 eq $file_2 & are finished with unpaired go back to $i loop
        last;
      }
    }
  }
}
else {    ##if only unmatched files are provided
  my @files_singles = <$path/flanking_seq/*flankingReads.fq>;
  foreach my $file ( sort @files_singles ) {
    next if -z $file;    ##next file if size is 0
    $flanking_fq{$file}{unpaired} = $file;
  }
}
##get directory path and file name in separate variables
my @genome_dir = split '/', $genome_file;
my $genome_fa  = pop @genome_dir;
my $genome_dir = join '/', @genome_dir;
$genome_fa =~ /(.+)\.(fa|fasta)$/;
my $target     = $1;
my $target_dir = "$path/$target";
##make new directories for file output
`mkdir -p $path/bowtie_aln`;
my $te_dir_path = $path;
my @bowtie_out_files;

##align each newly created flanking fq files to the genome
##align all files indvidually
##then align again as mates
##will remove any redundant alignments
foreach my $key ( sort keys %flanking_fq ) {
  foreach my $type ( sort keys %{ $flanking_fq{$key} } ) {
    my $flanking_fq = $flanking_fq{$key}{$type};

    #remove and save filename part of path
    my @fq_path = split '/', $flanking_fq;
    my $fq_name = pop @fq_path;
    $fq_name =~ s/\.fq$//;

    if ( !$bowtie2 and $bowtie_sam ) {
      ##bowtie1 with sam output
`bowtie --sam --sam-nohead --sam-nosq -a -m 1 -v 3 -q $genome_file.bowtie_build_index $flanking_fq  1> $path/bowtie_aln/$target.$fq_name.bowtie.single.out 2>> $path/$target.stderr`;
    }
    elsif ($bowtie2) {
      ##bowtie2 -- need to get comparable -a -m1 -v3 arguments
`bowtie2 --sam-nohead --sam-nosq -x $genome_file.bowtie2_build_index -U $flanking_fq  1> $path/bowtie_aln/$target.$fq_name.bowtie.single.out 2>> $path/$target.stderr`;
    }
    else {
      ## bowtie1 with bowtie output
`bowtie --best -q $genome_file.bowtie_build_index $flanking_fq  1> $path/bowtie_aln/$target.$fq_name.bowtie.single.out 2>> $path/$target.stderr`;
    }

    push @bowtie_out_files,
      "$path/bowtie_aln/$target.$fq_name.bowtie.single.out";
  }    #end of foreach my $type ( sort keys %{ $flanking_fq{$key} } )
  if ( exists $flanking_fq{$key}{1} and exists $flanking_fq{$key}{2} ) {
    my $flanking_fq_1 = $flanking_fq{$key}{1};
    my $flanking_fq_2 = $flanking_fq{$key}{2};
    my @fq_path       = split '/', $flanking_fq_1;
    my $fq_name       = pop @fq_path;
    $fq_name =~ s/\.fq$//;
    if ( -s $flanking_fq_1 and -s $flanking_fq_2 ) {

      #clean reads if both flanking.fq are non-zero file size
      #clean means make sure the mate1&2 files are in the same order
      #and that any unmated reads are in the unpaired file
`$scripts/clean_pairs_memory.pl -1 $flanking_fq_1 -2 $flanking_fq_2 1> $te_dir_path/flanking_seq/$fq_name.unPaired.fq 2>> $path/$target.stderr`;
    }    #end of if ( -s $flanking_fq_1 and -s $flanking_fq_2 )
    if (  -s "$flanking_fq_1.matched"
      and -s "$flanking_fq_2.matched" )
    {
      if ( !$bowtie2 and $bowtie_sam ) {
        ##bowtie1 and sam output
`bowtie --sam --sam-nohead --sam-nosq -a -m 1 -v 3 -q $genome_file.bowtie_build_index -1 $flanking_fq_1.matched -2 $flanking_fq_2.matched 1> $path/bowtie_aln/$target.$fq_name.bowtie.mates.out 2>> $path/$target.stderr`;
      }
      elsif ($bowtie2) {
        ##bowtie2
`bowtie2  --sam-nohead  --sam-nosq -x $genome_file.bowtie2_build_index -1 $flanking_fq_1.matched -2 $flanking_fq_2.matched 1> $path/bowtie_aln/$target.$fq_name.bowtie.mates.out 2>> $path/$target.stderr`;
      }
      else {
        ##bowtie1 with bowtie output
`bowtie --best  -q $genome_file.bowtie_build_index -1 $flanking_fq_1.matched -2 $flanking_fq_2.matched 1> $path/bowtie_aln/$target.$fq_name.bowtie.mates.out 2>> $path/$target.stderr`;
      }
      push @bowtie_out_files,
        "$path/bowtie_aln/$target.$fq_name.bowtie.mates.out";
      if ( !$bowtie2 and $bowtie_sam ) {
        ##bowtie1 and sam output
`bowtie --sam --sam-nohead --sam-nosq -a -m 1 -v 3 -q $genome_file.bowtie_build_index $te_dir_path/flanking_seq/$fq_name.unPaired.fq 1> $path/bowtie_aln/$target.$fq_name.bowtie.unPaired.out 2>> $path/$target.stderr`;
      }
      elsif ($bowtie2) {
        ##bowtie2
`bowtie2 --sam-nohead --sam-nosq -x $genome_file.bowtie2_build_index -U $te_dir_path/flanking_seq/$fq_name.unPaired.fq 1> $path/bowtie_aln/$target.$fq_name.bowtie.unPaired.out 2>> $path/$target.stderr`;
      }
      else {
        ##bowtie1 with bowtie output
`bowtie --best -q $genome_file.bowtie_build_index $te_dir_path/flanking_seq/$fq_name.unPaired.fq 1> $path/bowtie_aln/$target.$fq_name.bowtie.unPaired.out 2>> $path/$target.stderr`;
      }
      push @bowtie_out_files,
        "$path/bowtie_aln/$target.$fq_name.bowtie.unPaired.out";
    }  # end of if(-s "$flanking_fq_1.matched" and -s "$flanking_fq_2.matched" )
  }  #end of if ( exists $flanking_fq{$key}{1} and exists $flanking_fq{$key}{2})
}    #end of foreach $key

my $files2merge;
my $filecount = scalar @bowtie_out_files;
if ( $filecount > 50 ) {
  my @big_files_2_merge;
  for ( my $i = 0 ; $i < $filecount ; $i = $i + 50 ) {
    $files2merge = join " ", ( splice( @bowtie_out_files, 0, 50 ) );
    if ( scalar @bowtie_out_files > 1 ) {
      `cat $files2merge > $path/bowtie_aln/$target.$TE.merged.bowtie.$i.temp`;
      push @big_files_2_merge,
        "$path/bowtie_aln/$target.$TE.merged.bowtie.$i.temp";
    }
    else {
      ##if there is only 1 file after processing the 50, then just push onto @big_files_2_merge
      push @big_files_2_merge, $files2merge;
    }
  }
  $files2merge = join ' ', @big_files_2_merge if @big_files_2_merge;
}
else {
  $files2merge = join ' ', @bowtie_out_files if @bowtie_out_files;
}

my $merged_bowtie = "$path/bowtie_aln/$target.$TE.bowtie.out";
if ($files2merge) {
  `cat $files2merge > $merged_bowtie`;
}
else {
  `touch $merged_bowtie`;
}
