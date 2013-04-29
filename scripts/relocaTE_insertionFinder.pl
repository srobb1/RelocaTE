#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Bio::DB::Fasta;

if ( !defined @ARGV ) {
  die "Do not run directly, to be called by relocaTE.pl\n";
}

my $required_reads = 1;    ## rightReads + leftReads needs to be > to this value
my $required_left_reads  = 1;       ## needs to be >= to this value
my $required_right_reads = 1;       ## needs to be >= to this value
my $bowtie               = shift;
my $usr_target           = shift;
my $genome_path          = shift;
my $TE                   = shift;
my $regex_file           = shift;
my $exper                = shift;
my $flank_len   = shift;    ##length of seq flanking insertions to be returned
my $existing_TE = shift;
my $mm_allow    = shift;
my $bowtie2     = shift;
my $bowtie_sam  = 1;        ## change to shift or remove in V2
my %existingTE;
my %existingTE_found;

##if flag contains 16 it is on the minus strand and is reported as the revcomp in the sam file
my %flag_minus_strand = (
  16,  1, 17,  1, 18,  1, 19,  1, 20,  1, 21,  1, 22,  1, 23,  1, 24,  1,
  25,  1, 26,  1, 27,  1, 28,  1, 29,  1, 30,  1, 31,  1, 48,  1, 49,  1,
  50,  1, 51,  1, 52,  1, 53,  1, 54,  1, 55,  1, 56,  1, 57,  1, 58,  1,
  59,  1, 60,  1, 61,  1, 62,  1, 63,  1, 80,  1, 81,  1, 82,  1, 83,  1,
  84,  1, 85,  1, 86,  1, 87,  1, 88,  1, 89,  1, 90,  1, 91,  1, 92,  1,
  93,  1, 94,  1, 95,  1, 112, 1, 113, 1, 114, 1, 115, 1, 116, 1, 117, 1,
  118, 1, 119, 1, 120, 1, 121, 1, 122, 1, 123, 1, 124, 1, 125, 1, 126, 1,
  127, 1, 144, 1, 145, 1, 146, 1, 147, 1, 148, 1, 149, 1, 150, 1, 151, 1,
  152, 1, 153, 1, 154, 1, 155, 1, 156, 1, 157, 1, 158, 1, 159, 1, 176, 1,
  177, 1, 178, 1, 179, 1, 180, 1, 181, 1, 182, 1, 183, 1, 184, 1, 185, 1,
  186, 1, 187, 1, 188, 1, 189, 1, 190, 1, 191, 1, 208, 1, 209, 1, 210, 1,
  211, 1, 212, 1, 213, 1, 214, 1, 215, 1, 216, 1, 217, 1, 218, 1, 219, 1,
  220, 1, 221, 1, 222, 1, 223, 1, 240, 1, 241, 1, 242, 1, 243, 1, 244, 1,
  245, 1, 246, 1, 247, 1, 248, 1, 249, 1, 250, 1, 251, 1, 252, 1, 253, 1,
  254, 1, 255, 1
);

if ( $existing_TE ne 'NONE' ) {
  my $blat = 0;
  open INTE, "$existing_TE" or die $!;
  while ( my $line = <INTE> ) {
    chomp $line;
    if ( $line =~ /psLayout version 3/ ) {
      $blat = 1;
      <INTE>;    ## throw out blat header
      <INTE>;    ## throw out blat header
      <INTE>;    ## throw out blat header
      <INTE>;    ## throw out blat header
      chomp( $line = <INTE> );
    }
    next if $line !~ /\b$usr_target[:\s]/;
    next if $line !~ /\b$TE\b/;
    if ($blat) {
      my @blat     = split /\t/, $line;
      my $match    = $blat[0];
      my $mismatch = $blat[1];
      my $strand   = $blat[8];
      my $qName    = $blat[9];
      my $qLen     = $blat[10];
      my $qStart =
        $blat[11] + 1;    #get all values into 1st base = 1 postion notation
      my $qEnd  = $blat[12];
      my $tName = $blat[13];
      my $tLen  = $blat[14];
      my $tStart =
        $blat[15] + 1;    #get all values into 1st base = 1 postion notation
      my $tEnd = $blat[16];
      my $MM = $mismatch / ( $match + $mismatch );

      if ( $tEnd < $tStart ) {
        ( $tStart, $tEnd ) = ( $tEnd, $tStart );
      }
      if (
        $qStart <= 5
        and ( $qEnd >= ( $qLen - 5 ) )
        ## mismatches should be less than 10% the alignment length
        and $MM <= .1
        ## ref hit should be about the same size as the query
        and ( ( $tEnd - $tStart + 1 ) >= ( $qLen - ( $qLen * .1 ) ) )
        and ( ( $tEnd - $tStart + 1 ) <= ( $qLen + ( $qLen * .1 ) ) )
        ## alignment should be about the same size as the query
        and ( ( $match + $mismatch ) >= ( $qLen - ( $qLen * .1 ) ) )
        and ( ( $match + $mismatch ) <= ( $qLen + ( $qLen * .1 ) ) )
        )
      {
        #print "$qName\t$tName:$tStart..$tEnd $MM --\n";
        $existingTE{$qName}{start}{$tStart} = "$qName\t$tName:$tStart..$tEnd";
        $existingTE{$qName}{end}{$tEnd}     = "$qName\t$tName:$tStart..$tEnd";
        $existingTE_found{"$qName\t$tName:$tStart..$tEnd"}{start} = 0;
        $existingTE_found{"$qName\t$tName:$tStart..$tEnd"}{end}   = 0;
      }
      ##will get rid of this. will only use auto run blat and above.
    }
    else {
      my ( $te, $chr, $start, $end ) = $line =~ /(\S+)\t(\S+):(\d+)\.\.(\d+)/;
      ( $start, $end ) = sort { $a <=> $b } ( $start, $end );
      $existingTE{$te}{start}{$start} = $line;
      $existingTE{$te}{end}{$end}     = $line;
      $existingTE_found{$line}{start} = 0;
      $existingTE_found{$line}{end}   = 0;
    }
  }
}

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
  my @line = split /\t/, $line;
  $TSD = $line[3];
}
my $TSD_pattern = $TSD =~ /[\[.*+?]/ ? 1 : 0;    #does $TSD contain a pattern?

##get chromosome sequence for substr of flanking seq
my $db_obj     = Bio::DB::Fasta->new($genome_path);
my $genome_seq = $db_obj->seq($usr_target);

#remove redundant lines.
my %bowtie;
open BOWTIE, "$bowtie"
  or die "there seems to not be a bowtie file that i can't open $bowtie $!";
while ( my $line = <BOWTIE> ) {
  chomp $line;
  my @line = split /\t/, $line;
  next if $line[2] ne $usr_target;
  my $start = $line[3];
  ###bowtie_output is 0 sam is 1

  #remove /1 or /2 from the read name
  #7:12:11277:9907:Y/1
  $line =~ s/(^.+?)\/[1|2](\t.+$)/$1$2/;
  if ( !exists $bowtie{$line} ) {
    $bowtie{$line} = $start;
  }
}

#make new sorted sam array by sorting on the value of the sort hash
my @sorted_bowtie = sort { $bowtie{$a} <=> $bowtie{$b} } keys %bowtie;

my $last_start = 0;
my $last_end   = 0;
my %teInsertions;
my %teReadClusters;
my $count   = 0;
my @bin     = (0);
my $TSD_len = length $TSD;
foreach my $line (@sorted_bowtie) {
  chomp $line;
  my (
    $name, $flag, $target, $start, $MAPQ, $cigar,
    $MRNM, $MPOS, $TLEN,   $seq,   $qual, @tags
  );
  my ( $strand, $M, $mismatch );
  my ( $len, $end );
  if ( !$bowtie2 and !$bowtie_sam ) {
    ## bowtie1 with bowtie output
    ( $name, $strand, $target, $start, $seq, $qual, $M, $mismatch ) =
      split /\t/, $line;
    next if $M > 0;
    $start = $start + 1;
    my @mismatches = split ',', $mismatch;
    my $mm_count = scalar @mismatches;
    ## only 3 mismatch allowed total
    next if $mm_count > 3;
    $len = length $seq;
    $end = $len + $start - 1;
  }
  elsif ($bowtie2) {
    (
      $name, $flag, $target, $start, $MAPQ, $cigar,
      $MRNM, $MPOS, $TLEN,   $seq,   $qual, @tags
    ) = split /\t/, $line;
    next if $target ne $usr_target;
    $len = length $seq;
    $end = $len + $start - 1;

    $strand = '+';
    if ( exists $flag_minus_strand{$flag} ) {
      $strand = '-';
    }
    ##bowtie2: there is no -v option
    my $tooManyMM = 0;
    next if $MAPQ < 40;    ## higher score means more chance of unique map
    my ( $first_map_score, $second_map_score );
    foreach my $tag (@tags) {
      next unless $tag =~ /XM/;
      if ( $tag =~ /XM:i:(\d+)/ ) {
        $tooManyMM = 1 if $1 > 3;
      }
      elsif ( $tag =~ /AS:i:(\-?\d+)/ ) {
        $first_map_score = $1;
      }
      elsif ( $tag =~ /XS:i:(\-?\d+)/ ) {
        $second_map_score = $1;
      }
    }
    if (  defined $second_map_score
      and $second_map_score != 0
      and $second_map_score == $first_map_score )
    {
      next;    ## if the send alignment is as good as the first, skip it
    }
    next if $tooManyMM;
  }
  else {
    ## bowtie1 with sam output
    ## bowtie1: don't need to filter results since we use bowtie -a -m 1 -v 3, already uniq mapping reads
    (
      $name, $flag, $target, $start, $MAPQ, $cigar,
      $MRNM, $MPOS, $TLEN,   $seq,   $qual, @tags
    ) = split /\t/, $line;
    $len = length $seq;
    $end = $len + $start - 1;
    $strand = '+';
    if ( exists $flag_minus_strand{$flag} ) {
      $strand = '-';
    }
  }
  next if $target ne $usr_target;

  ## test to see if we are still within one insertion event or a different one
  ## is this seq aligned to same region, is it in range
  ## if two sets of overlapping reads are separated by 5bp, these two sets
  ## are considered to be one set. ==> $range_allowance
  my $range_allowance = 0;
  #my $range_allowance = 5;
  my $padded_start    = $bin[0] - $range_allowance;
  my $padded_end      = $bin[-1] + $range_allowance;
  if ( ( $start >= $padded_start and $start <= $padded_end )
    or ( $end >= $padded_start and $end <= $padded_end ) )
  {
    push @bin, $start, $end;
    @bin = sort @bin;
    if ($TSD !~ /UNK/i){
      TSD_check( $count, $seq, $start, $name, $TSD, $strand );
    }else{
      calculate_cluster_depth( $count, $seq, $start, $name, $strand );
    }
  }
  else {
    ## if start and end do not fall within last start and end
    ## we now have a different insertion event
    $count++;
    if ($TSD !~ /UNK/i){
      TSD_check( $count, $seq, $start, $name, $TSD, $strand );
    }else{
      calculate_cluster_depth( $count, $seq, $start, $name, $strand );
    }

    #reset last_start, last_end, @bin
    @bin        = ( $start, $end );
    $last_start = $start;
    $last_end   = $end;
  }
}

if ($TSD =~ /UNK/i){
  ## count depth to find TSD in 
  ## if there are 5 reads (2 right, 3 left) they
  ## should only be a depth of 5 at the TSD
  foreach my $cluster ( sort {$a <=> $b} keys %teReadClusters  ){
    my $read_total = $teReadClusters{$cluster}{read_count};
    my $TSD_len;
    foreach my $chrom_pos ( sort {$a <=> $b} keys %{$teReadClusters{$cluster}{depth}} ){
      my $depth = $teReadClusters{$cluster}{depth}{$chrom_pos};
      if ($depth == $read_total){
        $TSD_len++;
      }
    }
    ## if we have a TSD, then we can proceed to the TSD_check
    if ( $TSD_len ){
      $TSD = '.' x $TSD_len; ## create TSD regex, ex: '....'
      foreach my $name ( keys %{$teReadClusters{$cluster}{read_info}} ){
        my $seq = $teReadClusters{$cluster}{read_info}{$name}{seq};
        my $start = $teReadClusters{$cluster}{read_info}{$name}{seq_start};
        my $strand = $teReadClusters{$cluster}{read_info}{$name}{strand};
        TSD_check( $cluster, $seq, $start, $name, $TSD, $strand );
      }
      ## clean up, we don't need this info anymore
      delete $teReadClusters{$cluster} ;
    }
  }
  
}
##outdir/te/sam/sam_file
my $event = 0;
my @path = split '/', $bowtie;
pop @path;    #throw out filename
pop @path;    #throwout sam dir
my $te_dir = join '/', @path;
my $results_dir = "$te_dir/results";
`mkdir -p $results_dir`;
open OUTFASTA, ">$results_dir/$usr_target.$TE.confident_nonref_genomeflank.fa"
  or die $!;
open OUTALL, ">$results_dir/$usr_target.$TE.all_nonref_insert.txt" or die $!;
open OUTGFF, ">$results_dir/$usr_target.$TE.all_insert.gff"
  or die $!;
open OUTTABLE, ">$results_dir/$usr_target.$TE.confident_nonref_insert.txt"
  or die $!;
open OUTLIST,
  ">$results_dir/$usr_target.$TE.confident_nonref_insert_reads_list.txt"
  or die $!;
print OUTGFF "##gff-version	3\n";
##output in tab delimited table
my $tableHeader =
"TE\tTSD\tExper\tchromosome\tinsertion_site\tstrand\tleft_flanking_read_count\tright_flanking_read_count\tleft_flanking_seq\tright_flanking_seq\n";
print OUTTABLE $tableHeader;

my $note;
foreach my $insertionEvent ( sort { $a <=> $b } keys %teInsertions ) {
  foreach my $foundTSD ( sort keys %{ $teInsertions{$insertionEvent} } ) {
    foreach my $start (
      sort { $a <=> $b }
      keys %{ $teInsertions{$insertionEvent}{$foundTSD} }
      )
    {
      my $TE_orient_foward =
        exists $teInsertions{$insertionEvent}{$foundTSD}{$start}{TE_orient}{'+'}
        ? $teInsertions{$insertionEvent}{$foundTSD}{$start}{TE_orient}{'+'}
        : 0;
      my $TE_orient_reverse =
        exists $teInsertions{$insertionEvent}{$foundTSD}{$start}{TE_orient}{'-'}
        ? $teInsertions{$insertionEvent}{$foundTSD}{$start}{TE_orient}{'-'}
        : 0;
      my $TE_orient = $TE_orient_foward > $TE_orient_reverse ? "+" : "-";
      my $start_count =
        $teInsertions{$insertionEvent}{$foundTSD}{$start}{count};
      my $left_count = $teInsertions{$insertionEvent}{$foundTSD}{$start}{left};
      my $right_count =
        $teInsertions{$insertionEvent}{$foundTSD}{$start}{right};
      my @reads = @{ $teInsertions{$insertionEvent}{$foundTSD}{$start}{reads} };

      if (  ( defined $left_count and defined $right_count )
        and ( $left_count >= $required_left_reads )
        and ( $right_count >= $required_right_reads )
        and ( ( $right_count + $left_count ) > $required_reads ) )
      {
        $event++;
        my $coor                  = $start + ( $TSD_len - 1 );
        my $zero_base_coor        = $coor - 1;
        my $sub_string_start      = $zero_base_coor - $flank_len + 1;
        my $seq_start             = $coor - $flank_len + 1;
        my $seq_end               = $coor + $flank_len;
        my $left_flanking_ref_seq = substr $genome_seq, $sub_string_start,
          $flank_len;
        my $right_flanking_ref_seq = substr $genome_seq,
          $zero_base_coor + 1, $flank_len;
        $note = "Non-reference, not found in reference";
        my $coor_start = $coor - length($foundTSD) + 1;
        my $tableLine =
"$TE\t$foundTSD\t$exper\t$usr_target:$coor_start..$coor\t$TE_orient\t$left_count\t$right_count\t$left_flanking_ref_seq\t$right_flanking_ref_seq\n";
        print OUTTABLE $tableLine;
        print OUTGFF
"$usr_target\t$exper\ttransposable_element_insertion_site\t$coor_start\t$coor\t.\t$TE_orient\t.\tID=$TE.te_insertion_site.$usr_target.$coor;Note=$note;left_flanking_read_count=$left_count;right_flanking_read_count=$right_count;left_flanking_seq=$left_flanking_ref_seq;right_flanking_seq=$right_flanking_ref_seq;TSD=$foundTSD\n";
        print OUTFASTA
">$exper.$usr_target:$coor_start..$coor TSD=$foundTSD $usr_target:$seq_start..$seq_end\n$left_flanking_ref_seq$right_flanking_ref_seq\n";
        if ($right_count>0 and $left_count>0){
          print OUTALL
"$TE\t$foundTSD\t$exper\t$usr_target\t$coor_start..$coor\t$TE_orient\tT:$start_count\tR:$right_count\tL:$left_count\n";
        }else {
          print OUTALL
"$TE\tinsufficient_data\t$exper\t$usr_target\t$coor_start..$coor\t$TE_orient\tT:$start_count\tR:$right_count\tL:$left_count\n";
        }
        print OUTLIST "$usr_target:$coor_start..$coor\t", join( ",", @reads ),
          "\n";
      }
      else {
        my $coor = $start + ( $TSD_len - 1 );
        my $coor_start = $coor - length($foundTSD) + 1;
        $left_count  = defined $left_count  ? $left_count  : 0;
        $right_count = defined $right_count ? $right_count : 0;
        if ($right_count>0 and $left_count>0){
          print OUTALL
"$TE\t$foundTSD\t$exper\t$usr_target\t$coor_start..$coor\t$TE_orient\tT:$start_count\tR:$right_count\tL:$left_count\n";
        }else {
          print OUTALL
"$TE\tinsufficient_data\t$exper\t$usr_target\t$coor_start..$coor\t$TE_orient\tT:$start_count\tR:$right_count\tL:$left_count\n";
        }
      }
    }
  }
}

my $allexisting = "$results_dir/$exper.$TE.all_reference.txt";
open ALLEXISTING, ">>$allexisting" or die $!;
print ALLEXISTING
  "strain\tTE\texistingTE_coor\treads_align_2_start\treads_align_2_end\n"
  if -s "$allexisting" < 10;
foreach my $found ( keys %existingTE_found ) {
  my $end_count   = $existingTE_found{$found}{end};
  my $start_count = $existingTE_found{$found}{start};
  my ( $ref, $s, $e ) = $found =~ /(\S+)\:(\d+)\.\.(\d+)/;
  if ( $start_count > 0 or $end_count > 0 ) {
    $note = "Shared, in ref and reads";
  }
  else {
    $note = "Reference-only, only present in reference";
  }
  print OUTGFF
"$ref\t$exper\ttransposable_element_insertion_site\t$s\t$e\t.\t.\t.\tID=$TE.te_insertion_site.$ref.$s..$e;TE_Name=$TE;Note=$note;left_flanking_read_count=$start_count;right_flanking_read_count=$end_count\n";
  print ALLEXISTING "$exper\t$found\t$start_count\t$end_count\n";
}
sub calculate_cluster_depth {
  my ( $event, $seq, $seq_start, $read_name, $strand ) = @_;
  #my ($te_5prime, $te_3prime) = (0,0);
  ## make a note of which end of the TE has been identified
  #if ( $read_name =~ /:5$/){
  #  $teReadClusters{$event}{five_prime}++;
  #}elsif ( $read_name =~ /:3$/ ){
  #  $teReadClusters{$event}{three_prime}++;
  #}
  $teReadClusters{$event}{read_count}++;
  $teReadClusters{$event}{read_info}{$read_name}{seq}=$seq;
  $teReadClusters{$event}{read_info}{$read_name}{seq_start}=$seq_start;
  $teReadClusters{$event}{read_info}{$read_name}{strand}=$strand;
  for (my $i = $seq_start ; $i < $seq_start + (length $seq) ; $i++){
     $teReadClusters{$event}{depth}{$i}++;
  }
}
sub TSD_check {
  ##$seq is entire trimmd read, not just the TSD portion of the read
  ##$start is the first postition of the entire read match to ref
  my ( $event, $seq, $seq_start, $read_name, $tsd, $strand ) = @_;
  $seq = uc $seq;
  my $rev_seq = reverse $seq;
  $rev_seq =~ tr /ATGCN/TACGN/;
  my $result = 0;
  my $TSD;
  my $TE_orient = 0;
  my $pos;
  my $start;    ## first base of TSD
  ##start means that the TE was removed from the start of the read
  ##5 means the trimmed end mapps to the 5prime end of the TE
  ##3 means the trimmed end mapps to the 3prime end of the TE

  if ( $strand eq '+' ) {
    if ( $read_name =~ /start:[35]$/
      and ( $seq =~ /^($tsd)/ or $rev_seq =~ /($tsd)$/ ) )
    {
      $result    = 1;
      $TSD       = $1;
      $pos       = "right";
      $TE_orient = ( substr( $read_name, -1 ) == 5 ) ? '-' : '+';
      $start     = $seq_start;
    }
    ##end means that the TE was removed from the end of the read
    elsif ( $read_name =~ /end:[35]$/
      and ( $seq =~ /($tsd)$/ or $rev_seq =~ /^($tsd)/ ) )
    {
      $result    = 1;
      $TSD       = $1;
      $pos       = "left";
      $TE_orient = ( substr( $read_name, -1 ) == 5 ) ? '+' : '-';
      $start     = $seq_start + ( ( length $seq ) - ( length $TSD ) );
    }
  }
  elsif ( $strand eq '-' ) {
    if ( $read_name =~ /start:[53]$/
      and ( $seq =~ /($tsd)$/ or $rev_seq =~ /^($tsd)/ ) )
    {
      $result    = 1;
      $TSD       = $1;
      $pos       = "left";
      $TE_orient = ( substr( $read_name, -1 ) == 5 ) ? '+' : '-';
      $start     = $seq_start + ( ( length $seq ) - ( length $TSD ) );
    }
    ##end means that the TE was removed from the end of the read
    elsif ( $read_name =~ /end:[53]$/
      and ( $seq =~ /^($tsd)/ or $rev_seq =~ /($tsd)$/ ) )
    {
      $result    = 1;
      $TSD       = $1;
      $pos       = "right";
      $TE_orient = ( substr( $read_name, -1 ) == 5 ) ? '-' : '+';
      $start     = $seq_start;
    }

  }
  ## existingTEs are TEs present in Ref and reads
  if ( $result and $TE_orient ) {
    my ( $tir1_end, $tir2_end );
    if ($pos eq 'left'){
      $tir1_end = $seq_start + length $seq;
    }elsif ($pos eq 'right'){
      $tir2_end = $seq_start - 1;
    }
    if ( defined $tir1_end and exists $existingTE{$TE}{start}{$tir1_end} ) {
      my $te_id = $existingTE{$TE}{start}{$tir1_end};
      $existingTE_found{$te_id}{start}++;
    }
    elsif ( defined $tir2_end and exists $existingTE{$TE}{end}{$tir2_end} ) {
      my $te_id = $existingTE{$TE}{end}{$tir2_end};
      $existingTE_found{$te_id}{end}++;
    }
    else {
      ## insert is not a ref insertion, it is a non-ref
      $teInsertions{$event}{$TSD}{$start}{count}++;
      $teInsertions{$event}{$TSD}{$start}{$pos}++;
      $teInsertions{$event}{$TSD}{$start}{TE_orient}{$TE_orient}++;
      $read_name =~ s/:start|:end//;
      push @{ $teInsertions{$event}{$TSD}{$start}{reads} }, $read_name;
    }
  }
}
