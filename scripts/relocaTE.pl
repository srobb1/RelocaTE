#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use Cwd;
use strict;

##change $scripts to location of relocaTE scripts
my $scripts = '~/bin/relocaTE_editing';

if ( !defined @ARGV ) {
  &getHelp();
}
my $genomeFasta = 'NONE';
my $te_fasta;
my $target             = 'NONE';
my $len_cutoff         = 10;
my $mismatch_allowance = 0;
my $fq_dir;
my $exper = 'not.given';
my $mate_file_1        = '_p1';
my $mate_file_2        = '_p2';
my $mate_file_unpaired = '.unPaired';
my $workingdir;
my $outdir     = 'outdir_teSearch';
my $parallel   = 1;
my $qsub_array = 1;
my ( $blat_minScore, $blat_tileSize ) = ( 10, 7 );
my $flanking_seq_len = 100;
my $existing_TE      = 'NONE';
GetOptions(
  'p|parallel:i'         => \$parallel,
  'a|qsub_array:i'       => \$qsub_array,
  'e|exper:s'            => \$exper,
  'w|workingdir:s'       => \$workingdir,
  'o|outdir:s'           => \$outdir,
  'd|fq_dir:s'           => \$fq_dir,
  'g|genomeFasta:s'      => \$genomeFasta,
  't|te_fasta:s'         => \$te_fasta,
  'l|len_cutoff:i'       => \$len_cutoff,
  'm|mismatch:f'         => \$mismatch_allowance,
  '1|mate_1_id:s'        => \$mate_file_1,
  '2|mate_2_id:s'        => \$mate_file_2,
  'u|unpaired_id:s'      => \$mate_file_unpaired,
  'bm|blat_minScore:i'   => \$blat_minScore,
  'bt|blat_tileSize:i'   => \$blat_tileSize,
  'f|flanking_seq_len:i' => \$flanking_seq_len,
  'x|existing_TE:s'      => \$existing_TE,
  'h|help'               => \&getHelp,
);
my $current_dir;
$qsub_array = 0 if $parallel == 0;

if ( defined $workingdir and -d $workingdir ) {
  $current_dir = File::Spec->rel2abs($workingdir);
  $current_dir =~ s/\/$//;
}
else {
  $current_dir = cwd();
}
my $mapping = 1;

if ( !defined $genomeFasta ) {
  print "\n\nPlease provide reference genome by using -g Genome fasta path\n";
  &getHelp();
}
elsif ( $genomeFasta eq 'NONE' ) {
  print
"You did not provide a genome fasta, if you proceed only reads containing the TE will be found, no mapping of insertions will be performed\n";
  print "Proceed without mapping?\n";
  my $answer;
  while ( $answer = <STDIN> ) {

    #Exit if it was just spaces (or just an enter)
    last if $answer =~ /^\s*|\n$/;
  }
  if ( $answer =~ /n/i ) {
    &getHelp();
  }
  else {
    $mapping = 0;
  }
}
elsif ( !-e $genomeFasta ) {
  print "$genomeFasta does not exist. Check file name.\n";
  &getHelp();
}
my $genome_path;
if ( -e $genomeFasta ) {
  $genome_path = File::Spec->rel2abs($genomeFasta);
}
if ( !defined $te_fasta ) {
  print
"\n\nPlease provide fasta file containing transposable elements by using -t TE fasta path\n";
  &getHelp();
}
elsif ( !-e $te_fasta ) {
  print "$te_fasta does not exist. Check file name.\n";
  &getHelp();
}
else {
  my $first_line = `head -n1 $te_fasta`;
  if ( $first_line !~ /^>\S+\s+TSD=\S+/ ) {
    die
"The TE_fasta:$te_fasta does not have the proper format:\n>TE_NAME TSD=TSD\nSEQUENCE\n";
  }
}
my @fq_files;
my %fq_files;
if ( !defined $fq_dir ) {
  print "\n\nPlease provide a directory of paired fastq files\n";
  &getHelp();
}
elsif ( $fq_dir eq 'SKIP' ) {
  ##skip all other steps for processing the raw fq files
}
elsif ( !-d $fq_dir ) {
  print
"\n\nCheck the spelling or location of $fq_dir, Please provide a directory of paired fastq files\n";
  &getHelp();
}
else {
  my $fq_path = File::Spec->rel2abs($fq_dir);
  @fq_files = <$fq_path/*fq>;
  my @fastq_files = <$fq_path/*fastq>;
  push @fq_files, @fastq_files;
  if ( scalar @fq_files == 0 ) {
    print "Must provide at least 1 short read file\n";
    &getHelp();
  }
}
my $existing_TE_path;
if ( $existing_TE ne 'NONE' ) {
  if ( !-e $existing_TE ) {
    print "The existing_TE file:$existing_TE, you provided can be not found\n";
    &getHelp();
  }
  else {
    $existing_TE_path = File::Spec->rel2abs($existing_TE);
    my $line = `head $existing_TE`;
    if ( $line !~ /\S+\t\S+:\d+\.\.\d+/ ) {
      print "The existing_TE file is not in the appropriate format:
   
mping   Chr12:839604..840033
mping   Chr11:23200534..23200105

TE_name<tab>ref_seqname:first_Base_Of_TIR1..Last_base_of_TIR2

   ";
      &getHelp();
    }
  }
}

sub getHelp {
  print ' 
usage:
./relocaTE.pl [-t TE_fasta_file][-g chromosome_genome_fasta][-d dir_of_fq][-e short_sample_name][-h] 

options:

**required:
-g STR          genome fasta file path. If not provided will only align reads to TE and remove TE seq from short reads. [no default]
-t STR          fasta containing nucleotide sequences of transposable elements with TSD=xxx in the desc. [no default]
-d STR          directory of paired and unpaired fastq files (paired _1.fq & _2.fq) (.fq or .fastq is acceptable)  [no default]

**recommended:
-e STR          Short sample name, will be used in the output files to create IDs for the insert (ex. A123) [not.given]
-o STR          name for directory to contain output directories and files, will be created for the run (ex. 04222012_A123) [outdir_teSearch]

**optional:
-p INT          run each genome sequence separetly, parallel. The alternative (0) would be to run one after the other (int, 0=false or 1=true) [1] 
-a INT          if running in parallel, run jobs as a qsub PBS array when many are creaated (see: man qsub option -t).(int, 0=false or 1=true) [1] 
-w STR          base working directory, needs to exist, will not create, full path [cwd] 
-l INT          len cutoff for the TE trimmed reads to be aligned [10] 
-m FRACTION     mismatch allowance for alignment to TE (ex 0.1) [0] 
-1 STR		regular expression to identify mate 1 paired files ex: file_p1.fq or file_1.noNumbers.fq [\'_p1\']
-2 STR          regular expression to identify mate 2 paired files ex: file_p2.fq or file_2.noNumbers.fq [\'_p2\']
-u STR          regular expression to identify unpaired files ex: file.unPaired.fq[\'.unPaired\'] 
-bm INT		blat minScore value, in comparison of reads to TE sequence [10]
-bt INT		blat tileSize value, in comparison of reads to TE sequence  [7]
-f INT		length of the sequence flanking the found insertion to be returned. This sequence is taken from the reference genome [100]
-x STR		tab-delimited file containing the coordinates of existing TE.
-h              this message

SAMPLE Existing TE (the two columns are tab-delimited)
mping   Chr12:839604..840033
mping   Chr12:1045463..1045892

SAMPLE TE FASTA
>mping TSD=TTA
GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATG
ATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTT
TCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGT
CCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAATGATCCCAGCAA
CTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGT
TTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATT
GTGACTGGCC

Must contain "TSD=", can be a Perl regular express.  
  Example: any 4 characters: TSD=....
  Example: A or T followed by GCC: TSD=(A|T)GCC 
  Example: CGA followed by any character then an A then CT or G: TSD=CGA.A(CT|G) 

';

  exit 1;
}
$outdir =~ s/\/$//;
my $te_path = File::Spec->rel2abs($te_fasta);
my $top_dir = $outdir;
if ($qsub_array) {
  `mkdir -p $current_dir/$top_dir/shellscripts`;
  open QSUBARRAY, ">$current_dir/$top_dir/run_these_jobs.sh"
    or die "Can't open $current_dir/$top_dir/run_these_jobs.sh\n";
}
##split genome file into individual fasta files
my @genome_fastas;
if ($mapping) {
  open( INFASTA, "$genome_path" ) || die "$!\n";
  my $i      = 0;
  my $exists = 0;
  while ( my $line = <INFASTA> ) {
    if ($exists) {
      next unless $line =~ /^>(\S+)/;
    }
    if ( $line =~ /^>(\S+)/ ) {
      my $id = $1;
      if ( $id =~ /\|/ ) {
        $id   =~ s/\|/_/g;
        $line =~ s/\|/_/g;
      }
      if ( $i > 0 ) {
        close(OUTFASTA);
        $i = 0;
      }
      my @genome_dir = split '/', $genome_path;
      pop @genome_dir;
      my $genome_dir = join '/', @genome_dir;
      my $new_file = "$genome_dir/$id.fa";
      push @genome_fastas, $new_file;
      if ( -e $new_file ) {
        $exists = 1;
        next;
      }
      else {
        $exists = 0;
      }
      open( OUTFASTA, ">$new_file" ) or die "$!\n";
      print OUTFASTA $line;
      $i++;
    }
    elsif ( $line !~ /^>/ ) {    ##should be sequence
      print OUTFASTA $line;
    }
    else {
      die "Your genome fasta file is in a unexpected format. 
I was expecting a line of seqeunce but found something else:
$line\n";
    }
  }
  close(INFASTA);
  close(OUTFASTA);

  #create bowtie index
  my $cmd;
  if ( !-e "$genome_path.bowtie_build_index.1.ebwt" ) {
    $cmd = "bowtie-build -f $genome_path $genome_path.bowtie_build_index";
  }
  $genome_path =~ /.+\/(.+)\.fa$/;
  my $ref = $1;
  if ( $parallel and defined $cmd ) {
    my $shell_dir = "$current_dir/$top_dir/shellscripts/step_1";
    `mkdir -p $shell_dir`;
    open OUTSH, ">$shell_dir/$ref.formatGenome.step_1.sh";
    print OUTSH "$cmd\n";
    close OUTSH;
    `chmod +x $shell_dir/$ref.formatGenome.step_1.sh`;
  }
  elsif ( $parallel and !defined $cmd ) {
    my $step2_file =
"$current_dir/$top_dir/shellscripts/step_1_not_needed_genomefasta_already_formatted";
    my $shell_dir = "$current_dir/$top_dir/shellscripts";
    `mkdir -p $shell_dir`;
    `touch $step2_file` if $parallel;
  }
  elsif ( defined $cmd ) {
    ##run it now
    `$cmd`;
  }
  if ($qsub_array) {
    if (
      !-e "$current_dir/$top_dir/shellscripts/step_1_not_needed_genomefasta_already_formatted"
      and $qsub_array )
    {
      print QSUBARRAY
"qsub $current_dir/$top_dir/shellscripts/step_1/$ref.formatGenome.step_1.sh\n";
    }
  }
}    ##end if($mapping)

my @fq;
my @fa;

#convert fq files to fa for blat
open QSUBARRAY2, ">$current_dir/$top_dir/shellscripts/run.step_2.sh"
  if $qsub_array;
my $fq_count = 0;
if ( $fq_dir ne 'SKIP' ) {
  foreach my $fq (@fq_files) {
    my $fq_path = File::Spec->rel2abs($fq);
    push @fq, $fq_path;
    my $fa = $fq;
    if ( $fa =~ s/\.(fq|fastq)$/.fa/ ) {
      push @fa, $fa;
      if ( !-e $fa ) {
        my $cmd = "$scripts/relocaTE_fq2fa.pl $fq_path $fa";
        if ($parallel) {
          my @fq_path   = split '/', $fq_path;
          my $fq_name   = pop @fq_path;
          my $shell_dir = "$current_dir/$top_dir/shellscripts/step_2";
          `mkdir -p $shell_dir`;
          my $outsh = ">$shell_dir/$fq_count." . "fq2fa.sh";
          open OUTSH, ">$outsh";
          print OUTSH "$cmd\n";
        }
        else {
          `$cmd`;
        }
      }
      else {
        my $shell_dir = "$current_dir/$top_dir/shellscripts";
        `mkdir -p $shell_dir`;
        my $step2_file =
"$current_dir/$top_dir/shellscripts/step_2_not_needed_fq_already_converted_2_fa";
        `touch $step2_file` if $parallel;
      }
    }
    else {
      print
"$fq does not seem to be a fastq based on the file extension. It should be fq or fastq\n";
      &getHelp();
    }
    $fq_count++;
  }
  if (
    !-e "$current_dir/$top_dir/shellscripts/step_2_not_needed_fq_already_converted_2_fa"
    and $qsub_array )
  {
    print QSUBARRAY "qsub -t 0-", $fq_count - 1,
      " $current_dir/$top_dir/shellscripts/run.step_2.sh\n";
    print QSUBARRAY2
      "sh $current_dir/$top_dir/shellscripts/step_2/\$PBS_ARRAYID.fq2fa.sh";
  }
  elsif ($qsub_array) {
    unlink "$current_dir/$top_dir/shellscripts/run.step_2.sh";
  }
}    ##end if $fq_dir ne 'SKIP'
close QSUBARRAY2;
##split TE fasta into single record fastas
my @te_fastas;
my %TSD;

##split the TE fasta of many seqs into individual files
open( INFASTA, "$te_fasta" ) || die "$!\n";
my $i = 0;
while ( my $line = <INFASTA> ) {
  if ( $line =~ /^>(\S+)\s+TSD=(\S+)/ ) {
    my $id = $1;
    $TSD{$id} = $2;
    if ( $i > 0 ) {
      close(OUTFASTA);
      $i = 0;
    }
    my $te_file = "$id.fa";
    $te_file =~ s/\|/_/g;
    ##create new dir for files: workingDir/outdir/TE/
    my $te_dir = "$current_dir/$top_dir/$id";
    push @te_fastas, "$te_dir/$te_file";
    `mkdir -p $te_dir`;
    open( OUTFASTA, ">$te_dir/$te_file" ) or die "$!\n";
    print OUTFASTA $line;
    $i++;
  }
  elsif ( $line =~ /^>/ and $line !~ /TSD=/ ) {
    die
"The TE_fasta:$te_fasta does not have the proper format:\n>TE_NAME TSD=TSD\nSEQUENCE\n";
  }
  else {    ##should be sequence
    print OUTFASTA $line;
  }
}
close(INFASTA);
close(OUTFASTA);

#foreach TE fasta blat against target chromosome and parse and find insertion sites
foreach my $te_path (@te_fastas) {
  my @path     = split '/', $te_path;
  my $te_fasta = pop @path;
  my $path     = join '/', @path;
  my $TE       = $te_fasta;
  $TE =~ s/\.fa//;
  `mkdir -p $path/blat_output`;
  `mkdir -p $path/flanking_seq`;
  `mkdir -p $path/te_containing_fq`;
  `mkdir -p $path/te_only_read_portions_fa`;

  #blat fa files against te.fa
  my @flanking_fq;
  my $fq_file_count = scalar @fq;

  open QSUBARRAY3, ">$current_dir/$top_dir/shellscripts/$TE.run.step_3.sh"
    if $qsub_array;
  open QSUBARRAY4, ">$current_dir/$top_dir/shellscripts/$TE.run.step_5.sh"
    if $qsub_array;
  for ( my $i = 0 ; $i < $fq_file_count ; $i++ ) {
    my $fa = $fa[$i];
    my $fq = $fq[$i];

    #remove and save filename part of path
    my @fa_path = split '/', $fa;
    my $fa_name = pop @fa_path;
    $fa_name =~ s/\.fa$//;
    my $shell_dir = "$current_dir/$top_dir/shellscripts/step_3/$TE";
    if ($parallel) {
      `mkdir -p $shell_dir`;
      open OUTSH, ">$shell_dir/$i.$TE.blat.sh";
    }

    #use pre-existing blatout files
    if ( !-e "$path/blat_output/$fa_name.te_$TE.blatout" ) {
      my $cmd =
"blat -minScore=$blat_minScore -tileSize=$blat_tileSize $te_path $fa $path/blat_output/$fa_name.te_$TE.blatout";
      print OUTSH "$cmd\n" if $parallel;
      `$cmd` if !$parallel;
    }

    #use pre-esixting te_containing_fq files
    my $te_Containing_fq =
      "$path/te_containing_fq/$fa_name.te_$TE.ContainingReads.fq";
    if ( -e $te_Containing_fq ) {
      $fq = $te_Containing_fq;
    }
    my $cmd =
"perl $scripts/relocaTE_trim.pl $path/blat_output/$fa_name.te_$TE.blatout $fq $len_cutoff $mismatch_allowance > $path/flanking_seq/$fa_name.te_$TE.flankingReads.fq";
    if ($parallel) {
      print OUTSH "$cmd\n";
      close OUTSH;
      `chmod +x $shell_dir/*blat.sh`;
    }
    else {
      `$cmd` if !$parallel;
    }
  }
  if ($qsub_array) {
    print QSUBARRAY "qsub -t 0-", $fq_file_count - 1,
      " $current_dir/$top_dir/shellscripts/$TE.run.step_3.sh\n";
    print QSUBARRAY3
"sh $current_dir/$top_dir/shellscripts/step_3/$TE/\$PBS_ARRAYID.$TE.blat.sh";
  }
  ##if a genome file was provided, align seqs to genome
  ##if no genome file was provided, will only blat and trim reads of te seq
  if ($mapping) {
    my $param_path = "$current_dir/$top_dir/$TE";
    my $outregex   = "$param_path/regex.txt";
    open OUTREGEX, ">$outregex" or die $!;
    print OUTREGEX "$mate_file_1\t$mate_file_2\t$mate_file_unpaired\t$TSD{$TE}";
    my $cmd =
"$scripts/relocaTE_align.pl $scripts $param_path $genome_path $outregex $TE $exper";
    if ( !$parallel ) {
      `$cmd`;
    }
    else {
      my $shell_dir = "$current_dir/$top_dir/shellscripts/step_4/$TE";
      $genome_path =~ /.+\/(.+)\.fa$/;
      my $ref = $1;
      `mkdir -p $shell_dir`;
      open OUTSH, ">$shell_dir/$ref.$TE.align.sh";
      print OUTSH "$cmd\n";
      close OUTSH;
      `chmod +x $shell_dir/$ref.$TE.align.sh`;
      print QSUBARRAY "qsub $shell_dir/$ref.$TE.align.sh\n";
    }

    my $genome_count = 0;
    foreach my $genome_file (@genome_fastas) {
      $genome_file =~ /.+\/(.+)\.fa$/;
      my $target = $1;
      $genome_path =~ /.+\/(.+)\.fa$/;
      my $ref           = $1;
      my $merged_bowtie = "$path/$ref/bowtie_aln/$ref.$TE.bowtie.out";
      my $cmd =
"$scripts/relocaTE_insertionFinder.pl $merged_bowtie $target $genome_file $TE $outregex $exper $flanking_seq_len $existing_TE_path";
      if ( !$parallel ) {
        `$cmd`;
      }
      else {
        my $shell_dir = "$current_dir/$top_dir/shellscripts/step_5/$TE";
        `mkdir -p $shell_dir`;
        open OUTSH, ">$shell_dir/$genome_count.$TE.findSites.sh";
        print OUTSH "$cmd\n";
        close OUTSH;
        `chmod +x $shell_dir/*findSites.sh`;
      }
      $genome_count++;
    }
    if ($qsub_array) {
      print QSUBARRAY "qsub -t 0-", $genome_count - 1,
        " $current_dir/$top_dir/shellscripts/$TE.run.step_5.sh\n";
      print QSUBARRAY4
"sh $current_dir/$top_dir/shellscripts/step_5/$TE/\$PBS_ARRAYID.$TE.findSites.sh";
    }
  }
  if ($qsub_array) {
    close QSUBARRAY3;
    close QSUBARRAY4;
  }
}

##cat all '.te_insertion_sites.table.txt' results into one file
foreach my $te_path (@te_fastas) {
  my @path     = split '/', $te_path;
  my $te_fasta = pop @path;
  my $path     = join '/', @path;
  my $TE       = $te_fasta;
  $TE =~ s/\.fa//;
    if ($parallel){
      `mkdir -p $current_dir/$top_dir/shellscripts/step_6/$TE`;
      open FINISH , ">$current_dir/$top_dir/shellscripts/step_6/$TE/step_6.$TE.finishing.sh";
      print FINISH 
"echo \"TE\tTSD\tExper\tchromosome\tinsertion_site\tleft_flanking_read_count\tright_flanking_read_count\tleft_flanking_seq\tright_flanking_seq\" > $path/results/temp
for i in \`ls $path/results/*.$TE.te_insertion_sites.table.txt\` ; do grep -v flanking_read_count \$i >> $path/results/temp ; done
mv $path/results/temp $path/results/all.$TE.te_insertion_sites.table.txt\n";
      `chmod +x $current_dir/$top_dir/shellscripts/step_6/$TE/step_6.$TE.finishing.sh`;
    }
    if ($qsub_array){
      print QSUBARRAY "qsub $current_dir/$top_dir/shellscripts/step_6/$TE/step_6.$TE.finishing.sh\n"; 
    }
    if (!$parallel and !$qsub_array){
      ##do it now
      `echo \"TE\tTSD\tExper\tchromosome\tinsertion_site\tleft_flanking_read_count\tright_flanking_read_count\tleft_flanking_seq\tright_flanking_seq\" > $path/results/temp`;
        my @files = `ls $path/results/*.$TE.te_insertion_sites.table.txt`;
        foreach my $file (@files){
          chomp $file;
         `grep -v flanking_read_count $file  >> $path/results/temp`;
        } 
       `mv $path/results/temp $path/results/all.$TE.te_insertion_sites.table.txt`;
     }
    close FINISH;
}

if ($qsub_array) {
  close QSUBARRAY;
}
