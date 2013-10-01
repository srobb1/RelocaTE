#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;
my $all_records_dir = shift;
my $genome_file     = shift;
my $te_fasta        = shift;
my $working_dir     = shift;
my $regex_file      = shift;
my $insert_pos_file = shift;
my $blat_dir = shift;
my %seqs;

##get the regelar expression patterns for mates and for the TE
##when passed on the command line as an argument, even in single
##quotes I lose special regex characters
open INREGEX, "$regex_file" or die "Can't read regex file$!";
my $mate_1_pattern;
my $mate_2_pattern;
while ( my $line = <INREGEX> ) {
  chomp $line;
  ( $mate_1_pattern, $mate_2_pattern ) = split /\t/, $line;
}

## store ranges in TE fasta in %ranges
my %ranges;
open TE_FA, "$te_fasta" or die "Can't open $te_fasta\n";
my $TE ;
while ( my $line = <TE_FA> ) {
  chomp $line;
  if ( $line =~ /^>(\S+)/ ) {
    my $te = $1;
    $TE = $te;
    while ( $line =~ /range=(\S+):(\d+)..(\d+)/ig ) {
      my $range_name = $1;
      my $s          = $2;
      my $e          = $3;

      $ranges{$te}{$range_name} = range_create( $s, $e );
    }
  }
}
my %seq_storage;
if (!defined $blat_dir){
  $blat_dir = "$working_dir/blat_output";
}
my @blat_files = <$blat_dir/*blatout>;
if (scalar @blat_files == 0 ){
  die "Can't find blatout files in $blat_dir\n";
} 
foreach my $blat_file (@blat_files) {
  next if $blat_file =~ /unpaired/i;
  my @file_path = split '/', $blat_file;
  my $file_name = pop @file_path;
  my $FA        = $file_name;
  $FA =~ s/te_(.+).blatout/fa/;
  my $te      = $1;
  my $te_mate = $FA;
  $te_mate =~ s/\.fa//;
  my $prefix = $te_mate.'.fq';
  #my $prefix = $te_mate;
  if ( $prefix =~ /$mate_1_pattern/ ) {
      $prefix =~ s/$mate_1_pattern//;
  }elsif ($prefix =~ /$mate_2_pattern/) {
    $prefix =~ s/$mate_2_pattern//;
  }
  $prefix =~ s/\.fq//;
  my @dbs = `ls $all_records_dir/$prefix*fa`;
  print "$prefix: @dbs\n";
  ## @db should only be p1 and p2 == 2 files
  chomp @dbs;
  ##blat parser
  open INBLAT, $blat_file, or die "Please provide a blat output file\n";

  <INBLAT>;    #get rid of blat header lines
  <INBLAT>;
  <INBLAT>;
  <INBLAT>;
  <INBLAT>;
  while ( my $line = <INBLAT> ) {
    my @line        = split /\t/, $line;
    my $matches     = $line[0];
    my $mismatches  = $line[1];
    my $qBaseInsert = $line[5];
    my $tBaseInsert = $line[7];
    my $strand      = $line[8];
    my $qName       = $line[9];
    my $qLen        = $line[10];
    my $tLen        = $line[14];
    my $tStart      = $line[15] + 1;
    my $tEnd        = $line[16];
    my $id          = $qName;
    my $aln_bp      = $matches + $qBaseInsert + $mismatches;
    
    ### want to find hits to the TE
    ## throw out if gap is too big 
    if ($qBaseInsert > 5 or $tBaseInsert > 5){
      next;
    }
    ## if 50% or more of the read matches to the TE, keep it
    next unless ($matches + $mismatches) >= $qLen*.50 ;
    my $add = 0;
    if ( exists $seqs{$te}{$id}{$te_mate}{blat_hit}{matches} ) {
      my $stored_matches = $seqs{$te}{$id}{$te_mate}{blat_hit}{matches};
      my $stored_mm      = $seqs{$te}{$id}{$te_mate}{blat_hit}{mismatches};
      my $stored_qBI     = $seqs{$te}{$id}{$te_mate}{blat_hit}{qBaseInsert};
      my $stored_aln_bp  = $stored_matches + $stored_qBI + $stored_mm;
      ## if blat hit for this read already exists, is it better?
      if ( ($aln_bp) > ($stored_aln_bp) ) {
        $add = 1;
        delete $seqs{$te}{$id}{$te_mate}{blat_hit};
      }
      else {
        ##don't do anything, don't add new one, don't delete old one
      }
    }
    else {
      ##doesnt exist so add it
      $add = 1;
    }



    if ($add){
      $seqs{$te}{$id}{$te_mate}{blat_hit}{qlen}        = $qLen;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{qBaseInsert} = $qBaseInsert;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{tBaseInsert} = $tBaseInsert;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{strand}      = $strand;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{matches}     = $matches;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{mismatches}  = $mismatches;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{tStart}      = $tStart;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{tEnd}        = $tEnd;
      $seqs{$te}{$id}{$te_mate}{$te}                   = 1;

      foreach my $db (@dbs){
        next if $db =~ /unpaired/i;
        $seq_storage{$te}{$db}{$id} = 1;
      }
    }
  }
}

## retrieve seq and pair of any read that matches to the TE
foreach my $te (keys %seq_storage){
  foreach my $mate_fa ( keys %{$seq_storage{$te}} ) {
    my @mate_fa_path = split '/' , $mate_fa;
    my $file_name = pop @mate_fa_path;
    my $mate = $file_name;
    $mate =~ s/\.fa//;
    my $fastacmd_str = join ',' , sort keys %{$seq_storage{$te}{$mate_fa}};
    my @to_get = sort keys %{$seq_storage{$te}{$mate_fa}};
    my $seq_recs = '';
    if ( @to_get > 500 ) {
        for ( my $i = 0 ; $i < @to_get ; $i = $i + 500 ) {
          $fastacmd_str = join ",", ( splice( @to_get, 0, 500 ) );
          $seq_recs .= `fastacmd -d $mate_fa -s $fastacmd_str`;
      }
     }
     else {
       $seq_recs .= `fastacmd -d $mate_fa -s $fastacmd_str`; 
     }

    if ( defined $seq_recs ) {
      my @seq_recs = split />/, $seq_recs;
      ## get rid of first empty record
      shift @seq_recs; 
      ## each record = id\nseq\nseq\n
      foreach my $rec (@seq_recs) {
        my ( $header, @seq ) = split /\n/, $rec;
        my ($id) = split /\s+/, $header;
        $id =~ s/lcl\|//;
        my $seq = join '', @seq;
        $seqs{$te}{$id}{$mate}{seq} = $seq;
        ##print OUTFA ">$id\n$seq\n";
      }
    }
  }
}
##get reads to align to genome and print to a file for bowtie
foreach my $te ( keys %seqs ) {
  my $bowtie_2_aln = "$working_dir/$te.construcTEr.bowtie2aln.fa";
  open BOWTIEFA, ">$bowtie_2_aln" or die "Can't open $bowtie_2_aln, $!\n";
  foreach my $id ( keys %{ $seqs{$te} } ) {
    ##next if there is only 1 mate
    next if keys %{ $seqs{$te}{$id} } == 1;
    foreach my $mate ( keys %{ $seqs{$te}{$id} } ) {
      ##2 mates
      ##for the one that originally did not align to TE
      if ( exists $seqs{$te}{$id}{$mate}
        and !exists $seqs{$te}{$id}{$mate}{$te} )
      {
        my $seq = $seqs{$te}{$id}{$mate}{seq};
        print BOWTIEFA ">$mate,$id\n$seq\n";
      }
    }
  }
##aln to genome
##create bowtie index
  my $bowtie_out = "$working_dir/$te.construcTEr.bowtie.out";
  if ( !-e "$genome_file.bowtie_build_index.1.ebwt" ) {
    `bowtie-build -f $genome_file $genome_file.bowtie_build_index`;
  }
`bowtie -a -m 1 -v 3 -f $genome_file.bowtie_build_index $bowtie_2_aln  1> $bowtie_out 2> $working_dir/bowtie.stderr`;
#`bowtie --best -a -v 2 -f $genome_file.bowtie_build_index $bowtie_2_aln  1> $bowtie_out 2> $working_dir/bowtie.stderr`;
##parse bowtie out and record the genomic locations of alignments
  my $file_path = File::Spec->rel2abs($bowtie_out);
  open( my $BOWTIE_fh, "<", $file_path ) or die "Can't open $file_path $!\n";
  while ( my $line = <$BOWTIE_fh> ) {
    chomp $line;
    my @line = split /\t/, $line;
    my ( $name, $strand, $target, $start, $seq, $qual, $M, $mismatch ) =
      split /\t/, $line;
    my ( $this_mate, $id ) = split ',', $name;
    ##remove /1 or /2 from the read name
    ##7:12:11277:9907:Y/1     +       Chr1    22134042        TTTTTTATAAATGGATAA      DGGGGGDGGGGFGDGGGG      4       7:A>T,16:C>A
    $id =~ s/(^.+?)\/[1|2](\t.+$)/$1$2/;
    if ( exists $seqs{$te}{$id} ) {
      ##if an aln to genome exsts, delete record that says it aligns to TE
      if ( exists $seqs{$te}{$id}{$this_mate}{$te} ) {
        delete $seqs{$te}{$id}{$this_mate}{$te};
      }
      my $end = $start + ( length $seq ) - 1;
      my $id_str = "$start..$end($mismatch)";
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{$id_str}{start} = $start;
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{$id_str}{strand} = $strand;
      my @mismatches = split( ',', $mismatch );
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{$id_str}{mismatch} = @mismatches;
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{$id_str}{seq} = $seq;
    }
  }
}

my %construcTEr;
foreach my $te ( keys %seqs ) {
  foreach my $id ( keys %{ $seqs{$te} } ) {
    my $te_hit     = 0;
    my $genome_aln = 0;
    foreach my $mate ( keys %{ $seqs{$te}{$id} } ) {
      if ( exists $seqs{$te}{$id}{$mate}{blat_hit} ) {
        ##my $matches         = $seqs{$te}{$id}{$mate}{blat_hit}{matches};
        ##my $q_len           = $seqs{$te}{$id}{$mate}{blat_hit}{qlen};
        ##my $percent_aligned = $matches / $q_len;
        ##if ( $percent_aligned >= .95 ) {
        $te_hit = 1;
        ##}
      }
      if ( exists $seqs{$te}{$id}{$mate}{aln} ) {
        foreach my $target ( sort keys %{ $seqs{$te}{$id}{$mate}{aln} } ) {
          foreach
            my $id_str ( sort keys %{ $seqs{$te}{$id}{$mate}{aln}{$target} } )
          {
            my $mismatches =
              $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{mismatch};
            my $seq_len =
              length $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{seq};
            my $percent_mm = $mismatches / $seq_len;
            if ( $percent_mm < 0.1 ) {
              $genome_aln = 1;
            }
            else {
              delete $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str};
            }
          }
        }
      }
    }
    if ( $te_hit and $genome_aln ) {
      my %storage;
      foreach my $mate ( sort keys %{ $seqs{$te}{$id} } ) {
        #my $fq_obj = $seqs{$te}{$id}{$mate}{fq_rec};
        #my $seq    = get_seq($fq_obj);
        my $seq    =  $seqs{$te}{$id}{$mate}{seq};
        if ( exists $seqs{$te}{$id}{$mate}{blat_hit} ) {
          my $tStart     = $seqs{$te}{$id}{$mate}{blat_hit}{tStart};
          my $qlen       = $seqs{$te}{$id}{$mate}{blat_hit}{qlen};
          my $matches    = $seqs{$te}{$id}{$mate}{blat_hit}{matches};
          my $mismatches = $seqs{$te}{$id}{$mate}{blat_hit}{mismatches};
          my $tEnd       = $seqs{$te}{$id}{$mate}{blat_hit}{tEnd};
          my $strand     = $seqs{$te}{$id}{$mate}{blat_hit}{strand};
          $storage{TE}{$te}{start}      = $tStart;
          $storage{TE}{$te}{end}        = $tEnd;
          $storage{TE}{$te}{mismatches} = $mismatches;
          $storage{TE}{$te}{seq}        = $seq;
          $storage{TE}{$te}{strand}     = $strand;
        }
        elsif ( exists $seqs{$te}{$id}{$mate}{aln} ) {
          foreach my $target ( sort keys %{ $seqs{$te}{$id}{$mate}{aln} } ) {
            foreach
              my $id_str ( sort keys %{ $seqs{$te}{$id}{$mate}{aln}{$target} } )
            {
              my $start = $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{start};
              my $strand =
                $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{strand};
              my $mismatches =
                $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{mismatch};
              my $seq_len =
                length $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{seq};
              my $end = $start + $seq_len - 1;
              $storage{genome}{$target}{start}      = $start;
              $storage{genome}{$target}{end}        = $end;
              $storage{genome}{$target}{mismatches} = $mismatches;
              $storage{genome}{$target}{strand}     = $strand;
              $storage{genome}{$target}{seq}        = $seq;
            }
          }
        }
      }
      foreach my $te ( sort keys %{ $storage{TE} } ) {
        my $te_start      = $storage{TE}{$te}{start};
        my $te_read_seq   = $storage{TE}{$te}{seq};
        my $te_end        = $storage{TE}{$te}{end};
        my $te_strand     = $storage{TE}{$te}{strand};
        my $te_mismatches = $storage{TE}{$te}{mismatches};
        foreach my $target ( sort keys %{ $storage{genome} } ) {
          my $start           = $storage{genome}{$target}{start};
          my $genome_read_seq = $storage{genome}{$target}{seq};
          my $end             = $storage{genome}{$target}{end};
          my $strand          = $storage{genome}{$target}{strand};
          my $mismatches      = $storage{genome}{$target}{mismatches};
          $construcTEr{$target}{"$start..$end"}{$target}{start} = $start;
          $construcTEr{$target}{"$start..$end"}{$target}{end}   = $end;
          $construcTEr{$target}{"$start..$end"}{$target}{mismatch} =
            $mismatches;
          $construcTEr{$target}{"$start..$end"}{$target}{strand} = $strand;
          $construcTEr{$target}{"$start..$end"}{$target}{id}     = $id;
          $construcTEr{$target}{"$start..$end"}{$target}{seq} =
            $genome_read_seq;
          $construcTEr{$target}{"$start..$end"}{$te}{start}    = $te_start;
          $construcTEr{$target}{"$start..$end"}{$te}{end}      = $te_end;
          $construcTEr{$target}{"$start..$end"}{$te}{mismatch} = $te_mismatches;
          $construcTEr{$target}{"$start..$end"}{$te}{strand}   = $te_strand;
          $construcTEr{$target}{"$start..$end"}{$te}{id}       = $id;
          $construcTEr{$target}{"$start..$end"}{$te}{seq}      = $te_read_seq;
        }
      }
    } ## end: if ( $te_hit and $genome_aln )
    else {
      delete $seqs{$te}{$id};
    }
  }
}
my %inserts;
if (-e $insert_pos_file){
open INSERTS, "$insert_pos_file";
while (my $line = <INSERTS>){
  next if $line =~ /^TE.+Exper.+chromosome.+insertion_site.+left_flanking_read_count/;
  my ($TE, $TSD, $sample_desc, $target, $pos) = split /\t/, $line;
  my ($last_TSD_bp) = $pos =~ /\d+\.\.(\d+)/; 
  $inserts{$target}{"$target:$pos"}{pos}=$last_TSD_bp;
  }
}
#print Dumper \%inserts;
open OUTFA , ">$working_dir/$TE.construcTEr.reads.fa";
foreach my $target ( sort keys %construcTEr ) {
  foreach my $range (
    sort { ( split /\.\./, $a )[0] <=> ( split /\.\./, $b )[0] }
    keys %{ $construcTEr{$target} }
    )
  {
    foreach my $name ( keys %{ $construcTEr{$target}{$range} } ) {
      my $start    = $construcTEr{$target}{$range}{$name}{start};
      my $end      = $construcTEr{$target}{$range}{$name}{end};
      my $strand   = $construcTEr{$target}{$range}{$name}{strand};
      my $seq      = $construcTEr{$target}{$range}{$name}{seq};
      my $mismatch = $construcTEr{$target}{$range}{$name}{mismatch};
      my $id       = $construcTEr{$target}{$range}{$name}{id};

      #if ($strand eq '-'){
      #  ( $seq = reverse $seq ) =~ tr/AaGgTtCcNn/TtCcAaGgNn/;
      #  $id = $id.".revcom";
      #}
      ## if record is a TE hit - get range info

      my $range_str = '';
      ## only TEs are in %ranges
      if ( exists $ranges{$name} ) {
        my $read_range = range_create( $start, $end );
        foreach my $range_name ( sort keys %{ $ranges{$name} } ) {
          my $te_range = $ranges{$name}{$range_name};
          my $overlap = range_overlap( $te_range, $read_range );
          if ( $overlap >= 5 ) {
            $range_str = $range_str . ",overlap_with_$range_name=$overlap";
          }
        }
        $range_str =~ s/^,//;
      }
      if ($insert_pos_file){
      my ($r_s, $r_e) = split /\.\./ , $range;
      foreach my $insert_str(keys %{$inserts{$target}}){
        my $pos = $inserts{$target}{$insert_str}{pos};
        if ( (((abs ($r_s -$pos)) < 1000) or ((abs ($r_e -$pos)) < 1000)) 
            and  $mismatch == 0){
          push @{$inserts{$target}{$insert_str}{rec}} , ">$id $name:$start..$end ($strand) mismatches=$mismatch $range_str\n$seq\n";   
        }
      }
      }
      #print OUTFA
#">$id $name:$start..$end ($strand) mismatches=$mismatch $range_str\n$seq\n";
    }
  }
}

## trying to make a output table
if ($insert_pos_file) {
  my $out_dir = "$working_dir/insert_fa";
  `mkdir -p $out_dir`;
  open OUTTABLE , ">$working_dir/$TE.inserts.range.coverage.table.txt";
  my $rangeHeader = "TE\tpos";
  my @rangeHeader = sort keys %{ $ranges{$TE} };
  foreach my $name ( sort keys %{ $ranges{$TE} } ) {
    my $name_range = $ranges{$TE}{$name};
    my $s          = range_get_start($name_range);
    my $e          = range_get_end($name_range);
    $rangeHeader .= "\t$name:$s..$e";
  }
  print OUTTABLE $rangeHeader, "\n";
  foreach my $target ( sort keys %inserts ) {
    foreach my $insert_str ( keys %{ $inserts{$target} } ) {

      my @values;
      my @counts;
      open OUT, ">$out_dir/$insert_str.fa";
      foreach my $seq_rec ( @{ $inserts{$target}{$insert_str}{rec} } ) {
        print OUT "$seq_rec";
        for ( my $i = 0 ; $i < @rangeHeader ; $i++ ) {
          if (!defined $values[$i]){
            $values[$i] = 0;
            $counts[$i] = 0;
          }
          my $n = $rangeHeader[$i];
          if ( $seq_rec =~ /$n=(\d+)/ ) {
            my $value = $1;
            $counts[$i]++;
            if ( $value >= $values[$i] ) {
              $values[$i] = $value;
            }
          }
        }
      }
      my @forTable;
      for (my $i = 0 ; $i <@values ; $i++){
        push @forTable , "$values[$i]($counts[$i])";
      }
      print OUTTABLE "$TE\t$insert_str\t", join ("\t" ,@forTable) , "\n";
    }
    close OUT;
  }
}
#####SUBROUTINES########
sub dir_split {
  my $path = shift;
  my @path = split '/', $path;
  return @path;
}

sub filename_split {
  my $file = shift;
  my @file = split /\./, $file;
  return @file;
}

sub range_create {
## range needs to be s<=e
## range is in 1 base notation
## takes 2 numbers and returns an anonymous array
  my @range = ( $_[0], $_[1] );
  if ( $range[0] !~ /^\d+$/ or $range[1] !~ /^\d+$/ ) {
    die "Ranges provided in the TE fasta must be numbers but \n";
  }
  elsif ( $range[0] == 0 or $range[1] == 0 ) {
    die "Ranges are in 1 base notation not 0, numbers must be > 0\n";
  }
  @range = sort { $a <=> $b } @range;
  return \@range;
}

sub range_check {
  my $range = shift;
  if ( ref $range !~ /ARRAY/ ) {
    die "range functions need to be given an array ref\n";
  }
  if ( scalar @$range != 2 ) {
    die "range functions need to have 2 values\n";
  }

  return $range;

}
sub range_get_start {
 my $r = shift;
 return $$r[0];
}
sub range_get_end {
 my $r = shift;
 return $$r[1];
}
sub range_overlap {
  ## returns the size of the overlap
  my ( $r1, $r2 ) = @_;
  $r1 = range_check($r1);
  $r2 = range_check($r2);

  my ( $s,  $e )  = @$r1;
  my ( $s2, $e2 ) = @$r2;
  ## r     |-----------|
  ## r2 |-->
  ## r2    |
  if ( $s2 <= $s and $e2 >= $s ) {
    if ( $e2 <= $e ) {
      return ( $e2 - $s + 1 );
    }
    elsif ( $e2 > $e ) {
      return ( $e - $s + 1 );
    }
    else {
      die "range: error1\n";
    }
  }
  ## r  |---------------|
  ## r2 |--->
  ## r2    |--->
  elsif ( $s2 >= $s and $s2 <= $e and $e2 >= $s ) {
    if ( $e2 <= $e ) {
      return ( $e2 - $s2 + 1 );
    }
    elsif ( $e2 > $e ) {
      return ( $e - $s2 + 1 );
    }
    else {
      die "range: error2\n";
    }

  }
  else {
    return 0;
  }
}
