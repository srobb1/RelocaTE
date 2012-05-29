#!/usr/bin/perl -w

## 04212012: added the ability to use a tab-delimited list of existing TE locations 
## to identify if the reference TE-insertion is present in the reads 
## ----
## 02172012: changed the number of flanking reads that are needed for a 
## insert to be called from a total of 2 to (left>1 and right>1) therefore (total>3).  
## Also changed the max allowed mis-matches to 3

use strict;

my $required_reads = 3; ## rightReads + leftReads needs to be > or = to this value
my $required_left_reads = 1; ## needs to be > to this value
my $required_right_reads = 1; ## needs to be > to this value
my $bowtie  = shift;
my $usr_target  = shift;
my $genome_path = shift;
my $TE 		= shift;
my $regex_file	= shift;
my $exper	= shift;
my $flank_len   = shift; ##length of seq flanking insertions to be returned
my $existing_TE        = shift;
my %existingTE;
my %existingTE_found;
if ($existing_TE ne 'NONE'){
  open INTE, "$existing_TE" or die $!;
  while (my $line = <INTE>){
    chomp $line;
    next if $line !~ /$usr_target/;
    next if $line !~ /$TE/;
    my ($te , $chr, $start, $end) = $line =~ /(\S+)\t(\S+):(\d+)\.\.(\d+)/;
    ($start , $end) = sort {$a <=> $b} ($start, $end);
    $existingTE{$te}{start}{$start} = $line;
    $existingTE{$te}{end}{$end} = $line;
  }
}

my $genome_seq;
my $id;

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
my $TSD_pattern = $TSD =~/[\[.*+?]/ ? 1 : 0; #does $TSD contain a pattern?



##get chromosome sequence for substr of flanking seq
open GENOME, "$genome_path" or die "Cannot open genome fasta: $genome_path $!";
my $i = 0;
while ( my $line = <GENOME> ) {
    last if $i > 1;
    chomp $line;
    if ( $line =~ />(\S+)/ ) {
        $i++;
        $id = $1;
        if ( $id ne $usr_target ) {
            warn "-t $usr_target is not found in genome fasta.";
            warn
              "$id is found in genome fasta and will be used to parse results.";
            $usr_target = $id;
        }
    }
    else {
        $genome_seq .= $line;
    }
}
#remove redundant lines.
open BOWTIE, "$bowtie" or die "there seems to not be a bowtie file that i can open $!";
my %bowtie;
while (my $line = <BOWTIE>){
        chomp $line;
        my @line = split /\t/ , $line;
        next if $line[2] ne $usr_target;
        my $start = $line[3];
        #remove /1 or /2 from the read name
        #7:12:11277:9907:Y/1     +       Chr1    22134042        TTTTTTATAAATGGATAA      DGGGGGDGGGGFGDGGGG      4       7:A>T,16:C>A
        $line =~ s/(^.+?)\/[1|2](\t.+$)/$1$2/;
	if (! exists $bowtie{$line}){
		$bowtie{$line}=$start;
	}
}
#make new sorted sam array by sorting on the value of the sort hash
my @sorted_bowtie = sort {  $bowtie{$a} <=> $bowtie{$b} } keys %bowtie;

my $last_start = 0;
my $last_end   = 0;
my %teInsertions;
my $count = 0;
my @bin   = (0);
my $TSD_len = length $TSD;
foreach my $line (@sorted_bowtie) {
    chomp $line;
    my (
        $name,  $strand, $target, $start, $seq,
        $qual,  $M, $mismatch
    ) = split /\t/, $line;
    $start = $start + 1 ; #0 offset
    my $len = length $seq;
    ## format offset:reference-base>read-base
    my @mismatches = split ',' , $mismatch;
    my $mm_count = scalar @mismatches ;
    ## only 3 mismatch allowed total
    next if $mm_count > 3 ;
    my $end = $len + $start - 1;
    next if $target ne $usr_target;

    ## test to see if we are still within one insertion event or a different one
    ## is this seq aligned to same region, is it in range
    ## if two sets of overlapping reads are separated by 5bp, these two sets
    ## are considered to be one set. ==> $range_allowance
    my $range_allowance = 5;
    my $padded_start = $bin[0]  - $range_allowance;
    my $padded_end   = $bin[-1] + $range_allowance;
    if (   ( $start >= $padded_start and $start <= $padded_end )
        or ( $end >= $padded_start and $end <= $padded_end ) )
    {
        push @bin, $start, $end;
        @bin = sort @bin;
        my $first_TSD = uc( substr $seq, 0, $TSD_len );
        my $last_TSD = uc( substr $seq, (-1 * $TSD_len) );
        TSD_check( $count, $first_TSD, $start,   'left',  $name , $TSD);
        TSD_check( $count, $last_TSD,  $end - ($TSD_len - 1), 'right', $name , $TSD);
    }
    else {
        ## if start and end do not fall within last start and end
        ## we now have a different insertion event
        $count++;
        my $first_TSD = uc( substr $seq, 0, $TSD_len );
        my $last_TSD = uc( substr $seq, (-1 * $TSD_len) );
        TSD_check( $count, $first_TSD, $start,   'left',  $name , $TSD);
        TSD_check( $count, $last_TSD,  $end - ($TSD_len - 1), 'right', $name ,$TSD);

        #reset last_start, last_end, @bin
        @bin        = ( $start, $end );
        $last_start = $start;
        $last_end   = $end;
    }
}

my $event = 0;
my @path = split '/' , $bowtie;
my $bowtie_file_name = pop @path; #throw out filename
my $path = join '/' , @path;
my $ref_dir = pop @path;
my $te_dir = pop @path;
my $top_dir_path = join '/' , @path;
my $results_dir = "$top_dir_path/results";
`mkdir -p $results_dir`;
open OUTFASTA, ">$results_dir/$usr_target.$TE.te_insertion_sites.fa"         or die $!;
open OUTALL,   ">$results_dir/$usr_target.$TE.te_insertion.all.txt"         or die $!;
open OUTGFF,   ">$results_dir/$usr_target.$TE.te_insertion_sites.gff"        or die $!;
open OUTTABLE, ">$results_dir/$usr_target.$TE.te_insertion_sites.table.txt"  or die $!;
open OUTLIST,  ">$results_dir/$usr_target.$TE.te_insertion_sites.reads.list" or die $!;
print OUTGFF "##gff-version	3\n";
##output in tab delimited table
my $tableHeader = "TE\tTSD\tExper\tchromosome\tinsertion_site\tleft_flanking_read_count\tright_flanking_read_count\tleft_flanking_seq\tright_flanking_seq\n";
print OUTTABLE $tableHeader;

foreach my $insertionEvent ( sort { $a <=> $b } keys %teInsertions ) {
  foreach my $foundTSD (sort keys %{$teInsertions{$insertionEvent}}){
    foreach my $start (
        sort { $a <=> $b }
        keys %{ $teInsertions{$insertionEvent}{$foundTSD} }
      )
    {
        my $start_count = $teInsertions{$insertionEvent}{$foundTSD}{$start}{count};
        my $left_count  = $teInsertions{$insertionEvent}{$foundTSD}{$start}{left};
        my $right_count = $teInsertions{$insertionEvent}{$foundTSD}{$start}{right};
        my @reads       = @{ $teInsertions{$insertionEvent}{$foundTSD}{$start}{reads} };

        if ( ( defined $left_count and defined $right_count) 
              and ($left_count  >  $required_left_reads)
              and ($right_count >  $required_right_reads )
              and (($right_count + $left_count) >= $required_reads)
            ) {
            $event++;
            my $coor                  = $start + ($TSD_len - 1);
            my $zero_base_coor        = $coor - 1;
            my $sub_string_start      = $zero_base_coor - $flank_len + 1;
            my $seq_start             = $coor - $flank_len + 1;
            my $seq_end               = $coor + $flank_len;
            my $left_flanking_ref_seq = substr $genome_seq, $sub_string_start,
              $flank_len;
            my $right_flanking_ref_seq = substr $genome_seq,
              $zero_base_coor + 1, $flank_len;
my $tableLine = "$TE\t$foundTSD\t$exper\t$usr_target\t$coor\t$left_count\t$right_count\t$left_flanking_ref_seq\t$right_flanking_ref_seq\n";
            print OUTTABLE $tableLine;
            print OUTGFF
"$usr_target\t$exper\ttransposable_element_insertion_site\t$coor\t$coor\t.\t.\t.\tID=$TE.te_insertion_site.$usr_target.$coor;left_flanking_read_count=$left_count;right_flanking_read_count=$right_count;left_flanking_seq=$left_flanking_ref_seq;right_flanking_seq=$right_flanking_ref_seq;TSD=$foundTSD\n";
            print OUTFASTA
">$exper.$usr_target.$coor TSD=$foundTSD $usr_target:$seq_start..$seq_end\n$left_flanking_ref_seq$right_flanking_ref_seq\n";
            print OUTALL
"$TE\t$foundTSD\t$exper\t$usr_target\t$coor\tC:$start_count\tR:$right_count\tL:$left_count\n";
            print OUTLIST "$usr_target:$coor\t", join( ",", @reads ), "\n";
        }
        else {
            my $coor = $start + ($TSD_len - 1);
            $left_count  = defined $left_count  ? $left_count  : 0;
            $right_count = defined $right_count ? $right_count : 0;
            print OUTALL
"$TE\t$foundTSD\t$exper\t$usr_target\t$coor\tC:$start_count\tR:$right_count\tL:$left_count\n";
        }
    }
  }
}
print OUTALL
"
total confident insertions identified by a right AND left mping flanking sequence (C>=$required_reads,R>$required_right_reads,L>$required_left_reads)= $event
Note:C=total read count, R=right hand read count, L=left hand read count\n" if $event > 0;


if (scalar ( keys %existingTE_found ) > 0){
  open FOUND , ">>$results_dir/$exper.existing.$TE.found.txt" or die $!;
  print FOUND "strain\texistingTE_coor\treads_align_2_start\treads_align_2_end\n" if -s "$results_dir/$exper.existing.$TE.found.txt" < 10;;
  foreach my $found (keys %existingTE_found){
    my $end_count =  exists $existingTE_found{$found}{end} ? $existingTE_found{$found}{end} : 0;
    my $start_count = exists $existingTE_found{$found}{start} ? $existingTE_found{$found}{start} : 0; 
     
    print FOUND "$exper\t$found\t$start_count\t$end_count\n";
  }
}

sub TSD_check {
    my ( $event, $seq, $start, $pos, $read, $tsd ) = @_;
    my $rev_seq =  reverse $seq;
    $rev_seq =~ tr /AaTtGgCc/TtAaCcGg/;
    my $result = 0;
    my $TSD;
    if ($TSD_pattern and ($seq =~ /($tsd)/ or $rev_seq =~ /($tsd)/ ) ) {
	$result = 1;
        $TSD = $1;
    }elsif (!$TSD_pattern and ($seq eq $tsd or $rev_seq eq $tsd) ) {
    	$result = 1;
        $TSD = $tsd;
    }
    if ( $result ) {
        my ($tir1_end, $tir2_end) = ( ($start + (length $seq)) , ($start - 1));
        if (exists $existingTE{$TE}{start}{$tir1_end}){
          my $te_id = $existingTE{$TE}{start}{$tir1_end};
          $existingTE_found{$te_id}{start}++;
        }
        elsif (exists $existingTE{$TE}{end}{$tir2_end}){
          my $te_id = $existingTE{$TE}{end}{$tir2_end};
          $existingTE_found{$te_id}{end}++;
        }else{
          $teInsertions{$event}{$TSD}{$start}{count}++;
          $teInsertions{$event}{$TSD}{$start}{left}++  if $pos eq 'left';
          $teInsertions{$event}{$TSD}{$start}{right}++ if $pos eq 'right';
          push @{ $teInsertions{$event}{$TSD}{$start}{reads} }, $read;
        }
   }
}
