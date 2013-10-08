#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use Cwd;
use FindBin qw($RealBin);
use File::Path qw(make_path remove_tree);
use strict;

my $scripts = $RealBin;

if ( !defined @ARGV ) {
  &getHelp();
}
my $genome_fasta = 'NONE';
my $te_fasta;
my $target             = 'NONE';
my $len_cutoff         = 10;
my $mismatch_allowance = 0;
my $fq_dir;
my $exper              = 'not.given';
my $mate_file_1        = '_p1';
my $mate_file_2        = '_p2';
my $mate_file_unpaired = '.unPaired';
my $workingdir;
my $outdir     = 'outdir_teSearch';
my $parallel   = 0;
my $qsub_array = 0;
my $qsub_q;
my ( $blat_minScore, $blat_tileSize ) = ( 10, 7 );
my $flanking_seq_len = 100;
my $existing_TE      = 'NONE';
my $bowtie2          = 0;
my $nonLTR           = 0;
my $relax_align      = 0;
my $relax_reference      = 0;
GetOptions(
  'p|parallel:i'         => \$parallel,
  'a|qsub_array:i'       => \$qsub_array,
  'q|qsub_q:s'           => \$qsub_q,
  'e|exper:s'            => \$exper,
  'w|workingdir:s'       => \$workingdir,
  'o|outdir:s'           => \$outdir,
  'd|fq_dir:s'           => \$fq_dir,
  'g|genome_fasta:s'     => \$genome_fasta,
  't|te_fasta:s'         => \$te_fasta,
  'l|len_cutoff:i'       => \$len_cutoff,
  'm|mismatch:f'         => \$mismatch_allowance,
  '1|mate_1_id:s'        => \$mate_file_1,
  '2|mate_2_id:s'        => \$mate_file_2,
  'u|unpaired_id:s'      => \$mate_file_unpaired,
  'bm|blat_minScore:i'   => \$blat_minScore,
  'bt|blat_tileSize:i'   => \$blat_tileSize,
  'f|flanking_seq_len:i' => \$flanking_seq_len,
  'r|reference_ins:s'   => \$existing_TE,
  'b2|bowtie2:i'        => \$bowtie2,
  'ra|relax_align:i'    => \$relax_align,
  'rr|relax_reference:i'    => \$relax_reference,
  'h|help' => \&getHelp,
);
my $current_dir;

$parallel = 1 if $qsub_array == 1;
if ( defined $qsub_q ) {
  $qsub_q = "-q $qsub_q";
}
else {
  $qsub_q = '';
}

if ( defined $workingdir and -d $workingdir ) {
  $current_dir = File::Spec->rel2abs($workingdir);
  $current_dir =~ s/\/$//;
}
else {
  $current_dir = cwd();
}
my $mapping = 1;

if ( !defined $genome_fasta ) {
  print "\n\nPlease provide reference genome by using -g Genome fasta path\n";
  die "\nuse -h option to get help\n";
}
elsif ( $genome_fasta eq 'NONE' ) {
  print
"A reference genome fasta was NOT provided. Proceeding without a reference will result in only the reads containing the TE being identified, no mapping of insertions will be performed\n";
  print "Proceed without mapping? (y|n) \n";
  my $answer = <STDIN>;
  if ( $answer =~ /n/i ) {
    &getHelp();
  }
  elsif ( $answer =~ /y/i ) {
    $mapping = 0;
  }
  print "Great, proceeding without aligning to a reference genome.\n";
}
elsif ( !-e $genome_fasta ) {
  print "$genome_fasta does not exist. Check file name.\n";
  die "\nuse -h option to get help\n";
}
my $genome_path;
if ( -e $genome_fasta ) {
  $genome_path = File::Spec->rel2abs($genome_fasta);
}
if ( !defined $te_fasta ) {
  print
"\nPlease provide fasta file containing transposable elements by using -t TE fasta path

SAMPLE TE FASTA:
>mping TSD=TTA
GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATG
ATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTT
TCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGT
CCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAATGATCCCAGCAA
CTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGT
TTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATT
GTGACTGGCC

FASTA header must contain \"TSD=\", can be a Perl regular expression.  
  Example: these exact characters TTA: TSD=TTA 
  Example: any 4 characters: TSD=....
  Example: A or T followed by GCC: TSD=(A|T)GCC 
  Example: CGA followed by any character then an A then CT or G: TSD=CGA.A(CT|G) 
\n";
  die "\nuse -h option to get help\n";
}
elsif ( !-e $te_fasta ) {
  print "$te_fasta does not exist. Check file name.\n";
  die "\nuse -h option to get help\n";
}
else {
  open INFILE, $te_fasta or die "Can't open $te_fasta\n";
  my $first_line = <INFILE>;
  close INFILE;
  if ( $first_line !~ /^>\S+\s+TSD=\S+/ ) {
    die
"The TE_fasta:$te_fasta does not have the proper format:\n>TE_NAME TSD=TSD\nSEQUENCE\n";
  }
}
my @fq_files;
my %fq_files;
if ( !defined $fq_dir ) {
  print "\n\nPlease provide a directory of paired fastq files\n";
  die "\nuse -h option to get help\n";
}
elsif ( $fq_dir eq 'SKIP' ) {
  ##skip all other steps for processing the raw fq files
}
elsif ( !-d $fq_dir ) {
  print
"\n\nCheck the spelling or location of $fq_dir, Please provide a directory of paired fastq files\n";
  die "\nuse -h option to get help\n";
}
else {
  my $fq_path = File::Spec->rel2abs($fq_dir);
  @fq_files = <$fq_path/*.fq>;
  my @fastq_files = <$fq_path/*.fastq>;
  push @fq_files, @fastq_files;
  if ( scalar @fq_files == 0 ) {
    print "Must provide at least 1 short read file\n";
    die "\nuse -h option to get help\n";
  }
}
my $existing_TE_path = 'NONE';
my $existing_blat    = 0;
if ( $existing_TE ne 'NONE' ) {
  if ( $existing_TE eq '1' ) {
    ##run blat
    $existing_blat = 1;
  }
  elsif ( !-e $existing_TE ) {
    print "$existing_TE does not exist\n";
    print "Please use -r 1 or provide a file that exists\n";
    die "\nuse -h option to get help\n";
  }
  else {
    $existing_TE_path = File::Spec->rel2abs($existing_TE);
    open INFILE, $existing_TE or die "Can't open $existing_TE\n";
    my $first_line = <INFILE>;
    close INFILE;
    if ( $first_line !~ /\S+\t\S+:\d+\.\.\d+/ ) {
      print "The existing_TE file is not in the appropriate format:
   
SAMPLE Reference TEs (the two columns are tab-delimited):
mping   Chr12:839604..840033
mping   Chr11:23200534..23200105

TE_name<tab>ref_seqname:first_Base_Of_TIR1..Last_base_of_TIR2

or (recommended) use \'-r 1\' for RelocaTE to find your TE in the reference
   ";
      die "\nuse -h option to get help\n";
    }
  }
}

sub getHelp {
  print ' 
RelocaTE live:
usage:
./relocaTE.pl [-t TE_fasta_file][-g chromosome_genome_fasta][-d dir_of_fq][-e short_sample_name][-h] 

options:

**required:
-t |--te_fasta		file		fasta containing nucleotide sequences of transposable elements with 
					TSD=xxx in the desc. [no default]
-d |--fq_dir		dir		directory of paired and unpaired fastq files (paired _p1.fq & _p2.fq)
					(.fq or .fastq is acceptable)  [no default]

**recommended: 
-g |--genome_fasta	file		genome (reference) fasta file path. If not provided will only align 
					reads to TE and remove TE seq from short reads. [no default]
-e |--exper 		STR		Short sample name, will be used in the output files to create IDs for
					the insert (ex. A123) [not.given]
-o |--outdir 		STR		name for directory to contain output directories and files, will be
					created for the run (ex. 04222012_A123) [outdir_teSearch]

**optional:
-p |--parallel 		INT		Break down the single big job of relocaTE into as many smaller jobs as
					possible. The alternative (0) would be to run one after the other
					(int, 0=false or 1=true) [0] 
-q |--qsub_q 		STR		same as qsub -q option, not required [no default]
-a |--qsus_array	INT		if \'-a 1\' , create qsub PBS array jobs to run the many shell scripts
					created in the \'-a 1\' option. (see: man qsub option -t).(
					0=false or 1=true) [0] 
-w |--workingdir	dir		base working directory, needs to exist, will not be creates, full path
					required [cwd] 
-l |--len		INT		len cutoff for the TE trimmed reads to be aligned [10] 
-m |--mismatch		FRACTION	mismatch allowance for alignment to TE (ex 0.1) [0] 
-1 |--mate_1_id		STR		string to uniquely identify mate 1 paired files ex: file_p1.fq [_p1]
-2 |--mate_2_id		STR		pattern to uniquely identify mate 2 paired files ex: file_p2.fq [_p2]
-u |--unpaired_id	STR		pattern to uniquely identify unpaired files ex: file.unPaired.fq [.unPaired] 
-bm|--blat_minScore	INT		blat minScore value, used by blat in the comparison of reads to TE sequence [10]
-bt|--blat_tileSize	INT		blat tileSize value, used by blat in the comparison of reads to TE sequence  [7]
-f |--flanking_seq	INT		length of the sequence flanking the found insertion to be returned. This
					sequence is taken from the reference genome [100]
-r |--reference_ins	STR		To identify reference and shared insertions (reference and reads)
					choose option-1 or option-2. 
					option-1) (recommended) use \'-r 1\' to have RelocaTE find the location of your TE in the 
					reference.
					option-2) input the file name of a tab-delimited file containing the coordinates
					of TE insertions pre-existing in the reference sequence. [no default]
-b2 |--bowtie2	        INT             to use bowtie2 use \'-b2 1\' else for bowtie use \'-b2 0\' [0]
-h  |--help				this message


See documentation for more information. http://srobb1.github.com/RelocaTE/

';
  exit 1;
}
if ( $outdir eq '' or $outdir =~ /^\s+$/ or !defined $outdir ) {
  die "your -o option has an incorrect value, it needs to be something
\nuse -h option to get help\n";
}
else {
  $outdir =~ s/\/$//;
}
my $te_path = File::Spec->rel2abs($te_fasta);
my @outdir = split /\//, $outdir;
$outdir = pop @outdir;
my $top_dir = $outdir;
my @depend;
my $shellscripts = "$current_dir/$top_dir/shellscripts";
if ($qsub_array) {
  mkdir "$current_dir/$top_dir";
  mkdir "$shellscripts";
  open QSUBARRAY, ">$current_dir/$top_dir/run_these_jobs.sh"
    or die "Can't open $current_dir/$top_dir/run_these_jobs.sh\n";
}
elsif ($parallel) {
  mkdir "$current_dir/$top_dir";
  mkdir "$shellscripts";
  open PARALLEL, ">$current_dir/$top_dir/run_these_jobs.sh"
    or die "Can't open $current_dir/$top_dir/run_these_jobs.sh\n";
}
else {
  mkdir "$current_dir/$top_dir";
}
my $qsub_formatGenome_cmd = 0;
## get names of each ref sequecne
my @genome_seqs;
if ( $mapping > 0 ) {
  open( INFASTA, "$genome_path" ) || die "$!\n";
  while ( my $line = <INFASTA> ) {
    next unless $line =~ /^>(\S+)/;
    if ( $line =~ /^>(\S+)/ ) {
      my $id = $1;
      if ( $id =~ /\|/ ) {
        $id =~ s/\|/_/g;
      }
      push @genome_seqs, $id;
    }
    else {
      die "Your genome FASTA file is in a unexpected format. 
>SEQNAME
SEQUENCE
>SEQNAME2
SEQUENCE2\n";
    }
  }
  close(INFASTA);

  #create bowtie index
  my $cmd;
  if ( !$bowtie2 and !-e "$genome_path.bowtie_build_index.1.ebwt" ) {
    $cmd =
"bowtie-build -f $genome_path $genome_path.bowtie_build_index 12> $current_dir/$top_dir/bowtie-build.out";
    $qsub_formatGenome_cmd = 1;
  }
  elsif ( $bowtie2 and !-e "$genome_path.bowtie2_build_index.1.bt2" ) {
    $cmd =
"bowtie2-build -f $genome_path $genome_path.bowtie2_build_index 12> $current_dir/$top_dir/bowtie-build2.out";
    $qsub_formatGenome_cmd = 1;
  }
  my $ref = 'ref';
  if ( $genome_path =~ /(?:.+\/)?(.+)\.(fa|fasta)$/ ) {
    $ref = $1;
  }
  if ( $parallel and defined $cmd ) {
    my $shell_dir = "$shellscripts/step_1";
    mkdir $shell_dir;
    open OUTSH, ">$shell_dir/step_1.$ref.formatGenome.sh";
    print OUTSH "$cmd\n";
    close OUTSH;
    chmod 0755, "$shell_dir/step_1.$ref.formatGenome.sh";
    print PARALLEL "sh $shell_dir/step_1.$ref.formatGenome.sh\n"
      if !$qsub_array;
  }
  elsif ( $parallel and !defined $cmd ) {
    my $step1_file =
      "$shellscripts/step_1_not_needed_genome_fasta_already_formatted";
    my $shell_dir = "$shellscripts";
    mkdir $shell_dir;
    if ($parallel) {
      open STEP1, ">$step1_file" or die "Can't Open $step1_file\n";
      print STEP1 '';
      close STEP1;
    }
  }
  elsif ( defined $cmd ) {
    ##run it now
    print "Formatting the reference genome: $genome_path\n";
    system($cmd);
  }
  if ($qsub_array) {
    if ( !-e "$shellscripts/step_1_not_needed_genome_fasta_already_formatted" )
    {
      my $job = "$shellscripts/step_1/step_1.$ref.formatGenome.sh";
      print QSUBARRAY
        "STEP1=\`qsub -e $shellscripts -o $shellscripts $qsub_q $job\`
echo \$STEP1\n";

    }
  }
}    ##end if($mapping)
my $top_blat_output_dir = "$current_dir/$top_dir/blat_output";
mkdir $top_blat_output_dir;
my $nonLTR_blat_params = '';
if ( $nonLTR ){
  #$nonLTR_blat_params =  '-noTrimA -stepSize=5';
  $nonLTR_blat_params =  '-noTrimA';
}
##run existing TE blat against ref if the file does not exsit
my $existingTE_blatout = "$top_blat_output_dir/existingTE.blatout";
my $qsub_existingTE_cmd = 0;
my $existing_blat_cmd =
"blat $genome_path $nonLTR_blat_params $te_path $existingTE_blatout 1> $top_blat_output_dir/existingTE.blat.stdout";
#"blat -noTrimA $genome_path $te_path $current_dir/$top_dir/existingTE.blatout 1> $current_dir/$top_dir/existingTE.blat.stdout";
if ($existing_blat) {
  ##if running blat set existing_TE_path to blatout
  $existing_TE_path = $existingTE_blatout; 
  #$existing_TE_path = "$current_dir/$top_dir/existingTE.blatout";
  if ( $parallel
    and !-e $existingTE_blatout )
  {
    my $shell_dir = "$shellscripts";
    if ( !-d $shell_dir ) {
      mkdir $shell_dir;
    }
    open OUTSH, ">$shell_dir/step_0.existingTE_blat.sh"
      or die "Can't open $shell_dir/step_0.existingTE_blat.sh for writing $!\n";
    print PARALLEL "sh $shell_dir/step_0.existingTE_blat.sh\n" if !$qsub_array;
    print OUTSH "$existing_blat_cmd\n";
    if ($qsub_array) {
      $qsub_existingTE_cmd = 1;
      print QSUBARRAY
"EXISTINGTE=`qsub -e $shellscripts -o $shellscripts $qsub_q $shellscripts/step_0.existingTE_blat.sh`
echo \$EXISTINGTE\n";
    }
    close OUTSH;
  }
  elsif ( !-e $existingTE_blatout ) {
    ## do it now
    print "finding TEs ($te_path) in the reference genome ($genome_path)\n";
    system($existing_blat_cmd);
  }
}

my @fq;
my @fa;

#convert fq files to fa for blat
open QSUBARRAY2, ">$shellscripts/step_2.fq2fa.sh"
  if $qsub_array;
my $fq_count = 0;
my @convert2fa;
if ( $fq_dir ne 'SKIP' ) {
  foreach my $fq (@fq_files) {
    my $fq_path = File::Spec->rel2abs($fq);
    push @fq, $fq_path;
    my $fa = $fq;
    if ( $fa =~ s/\.(fq|fastq)$/.fa/ ) {
      push @fa, $fa;
      if ( !-e $fa ) {
        push @convert2fa , $fa;
        my $cmd = "perl $scripts/relocaTE_fq2fa.pl $fq_path $fa";
        if ($parallel) {
          my @fq_path   = split '/', $fq_path;
          my $fq_name   = pop @fq_path;
          my $shell_dir = "$shellscripts/step_2";
          my $file_count = @convert2fa -1;
          mkdir $shell_dir;
          my $outsh = "$shell_dir/$file_count." . "fq2fa.sh";
          open OUTSH, ">$outsh";
          print PARALLEL "sh $outsh\n" if !$qsub_array;
          print OUTSH "$cmd\n";
        }
        else {
          ##run it now
          print "Converting $fq_path to fasta for blat\n";
          system($cmd);
        }
      }
    }
    else {
      print
"$fq does not seem to be a fastq based on the file extension. It should be fq or fastq\n";
      &getHelp();
    }
    $fq_count++;
  }

  if ( !@convert2fa ) {
    my $shell_dir = "$shellscripts";

    mkdir $shell_dir;
    my $step2_file =
      "$shellscripts/step_2_not_needed_fq_already_converted_2_fa";

    if ($parallel) {
      open STEP2, ">$step2_file" or die "Can't Open $step2_file\n";
      print STEP2 '';
      close STEP2;
    }
  }


  #if ( !-e "$shellscripts/step_2_not_needed_fq_already_converted_2_fa"
  if ( @convert2fa and $qsub_array )
  {
    my $end = @convert2fa - 1;
    my $job = "$shellscripts/step_2.fq2fa.sh";
    if ( !@depend ) {
      print QSUBARRAY
        "STEP2=\`qsub -e $shellscripts -o $shellscripts $qsub_q -t 0-$end $job\`
echo \$STEP2\n";
      @depend = ( "STEP2", "afterokarray" );
    }
    else {
      my ( $last_job, $afterok ) = @depend;
      @depend = ( "STEP2", "afterokarray" );
      my $jobName = $depend[0];
      print QSUBARRAY
"$jobName=`qsub -e $shellscripts -o $shellscripts $qsub_q -t 0-$end -W depend=$afterok:\$$last_job $job`
echo \$$jobName\n";
    }
    print QSUBARRAY2 "sh $shellscripts/step_2/\$PBS_ARRAYID.fq2fa.sh";
  }
  elsif ($qsub_array) {
    unlink "$shellscripts/step_2.fq2fa.sh";
  }
}    ##end if $fq_dir ne 'SKIP'
close QSUBARRAY2;

##split the TE fasta of many seqs into individual files
my @te_fastas;
my %TSD;
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

    mkdir $te_dir;
    my $results_dir = "$te_dir/results";
    remove_tree ($results_dir) if -e $results_dir;
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
my $depend = 1 if @depend;
foreach my $te_path (@te_fastas) {
  if ($depend) {
    @depend = ( "STEP2", "afterokarray" );
  }
  else {
    @depend = ();
  }
  my @path     = split '/', $te_path;
  my $te_fasta = pop @path;
  my $path     = join '/', @path;
  my $TE       = $te_fasta;
  $TE =~ s/\.fa//;


  mkdir "$path/blat_output";
  mkdir "$path/flanking_seq";
  mkdir "$path/te_containing_fq";
  mkdir "$path/te_only_read_portions_fa";

  #blat fa files against te.fa
  my @flanking_fq;
  my $fq_file_count = scalar @fq;

  open QSUBARRAY3, ">$shellscripts/step_3.$TE.blat.sh"
    if $qsub_array;
  open QSUBARRAY4, ">$shellscripts/step_5.$TE.finder.sh"
    if $qsub_array;
  for ( my $i = 0 ; $i < $fq_file_count ; $i++ ) {
    my $fa = $fa[$i];
    my $fq = $fq[$i];

    #remove and save filename part of path
    my @fa_path = split '/', $fa;
    my $fa_name = pop @fa_path;
    $fa_name =~ s/\.fa$//;
    my $shell_dir = "$shellscripts/step_3/$TE";
    if ($parallel) {
      make_path( $shell_dir, { mode => 0755 } );
      open OUTSH, ">$shell_dir/$i.$TE.blat.sh"
        or die "Can't open $shell_dir/$i.$TE.blat.sh $!\n";
      print PARALLEL "sh $shell_dir/$i.$TE.blat.sh\n" if !$qsub_array;
    }

    #use pre-existing blatout files
    if ( !-e "$path/blat_output/$fa_name.te_$TE.blatout" ) {
      my $cmd =
"blat $nonLTR_blat_params -minScore=$blat_minScore -tileSize=$blat_tileSize $te_path $fa $path/blat_output/$fa_name.te_$TE.blatout 1>> $path/blat_output/blat.out";
#"blat -noTrimA -minScore=$blat_minScore -tileSize=$blat_tileSize $te_path $fa $path/blat_output/$fa_name.te_$TE.blatout 1>> $path/blat_output/blat.out";
      print OUTSH "$cmd\n" if $parallel;
      print "Finding reads in $fa_name that contain sequence of $TE\n"
        if !$parallel;
      system($cmd) if !$parallel;
    }

    #use pre-existing te_containing_fq files
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
      chmod 0755, "$shell_dir/*blat.sh";
    }
    else {
      ##run it now
      print "Trimming $fq reads of $TE sequence\n" if !$parallel;
      system($cmd) if !$parallel;
    }
  }
  if ($qsub_array) {
    my $end  = $fq_file_count - 1;
    my $job  = "$shellscripts/step_3.$TE.blat.sh";
    my $desc = $TE;
    $desc =~ s/\W/_/;
    if ( !@depend ) {
      print QSUBARRAY
"STEP_3_$desc=\`qsub -e $shellscripts -o $shellscripts $qsub_q -t 0-$end $job\`
echo \$STEP_3_$desc\n";
      @depend = ( "STEP_3_$desc", "afterokarray" );
    }
    else {
      my ( $last_job, $afterok ) = @depend;
      @depend = ( "STEP_3_$desc", "afterokarray" );
      my $jobName = $depend[0];
      print QSUBARRAY
"$jobName=`qsub -e $shellscripts -o $shellscripts $qsub_q -t 0-$end -W depend=$afterok:\$$last_job $job`
echo \$$jobName\n";

    }
    print QSUBARRAY3 "sh $shellscripts/step_3/$TE/\$PBS_ARRAYID.$TE.blat.sh";
  }
  ##if a genome file was provided, align seqs to genome
  ##if no genome file was provided, will only blat and trim reads of te seq
  if ($mapping) {
    my $param_path = "$current_dir/$top_dir/$TE";
    my $outregex   = "$param_path/regex.txt";
    open OUTREGEX, ">$outregex" or die $!;
    print OUTREGEX "$mate_file_1\t$mate_file_2\t$mate_file_unpaired\t$TSD{$TE}";
    my $cmd =
"$scripts/relocaTE_align.pl $scripts $param_path $genome_path $outregex $TE $exper $bowtie2 $relax_align";
#"$scripts/relocaTE_align.pl $scripts $param_path $genome_path $outregex $TE $exper $bowtie2";
    if ( !$parallel ) {
      ## run now
      print "Aligning $TE trimmed reads to the reference ($genome_path)\n";
      system($cmd);
    }
    else {
      my $shell_dir = "$shellscripts/step_4/$TE";
      $genome_path =~ /.+\/(.+)\.(fa|fasta)$/;
      my $ref = $1;

      #`mkdir -p $shell_dir`;
      make_path( $shell_dir, { mode => 0755 } );

      #mkdir $shell_dir;
      open OUTSH, ">$shell_dir/step_4.$ref.$TE.align.sh";
      print OUTSH "$cmd\n";
      print PARALLEL "sh $shell_dir/step_4.$ref.$TE.align.sh\n" if !$qsub_array;
      close OUTSH;

      #`chmod +x $shell_dir/step_4.$ref.$TE.align.sh`;
      chmod 0755, "$shell_dir/step_4.$ref.$TE.align.sh";
      if ($qsub_array) {
        my $existing_depend = '';
        if ( $qsub_formatGenome_cmd ) {
          $existing_depend = "-W depend=afterok:\$STEP1" if !@depend;
          $existing_depend = ",afterok:\$STEP1" if @depend;
          #$existing_depend = ",depend=afterok:\$STEP1" if @depend;
        }
        my $job  = "$shell_dir/step_4.$ref.$TE.align.sh";
        my $desc = $TE;
        $desc =~ s/\W/_/;
        if ( !@depend ) {
          print QSUBARRAY
"STEP_4_$desc=\`qsub -e $shellscripts -o $shellscripts $qsub_q $existing_depend $job\`
echo \$STEP_4_$desc\n";
          @depend = ( "STEP_4_$desc", "afterok" );
        }
        else {
          my ( $last_job, $afterok ) = @depend;
          @depend = ( "STEP_4_$desc", "afterok" );
          my $jobName = $depend[0];
          print QSUBARRAY
"$jobName=`qsub -e $shellscripts -o $shellscripts $qsub_q -W depend=$afterok:\$$last_job","$existing_depend $job`
echo \$$jobName\n";
        }
      }
    }

    my $genome_count = 0;
    foreach my $seq_id (@genome_seqs) {
      $genome_path =~ /.+\/(.+)\.(fa|fasta)$/;
      my $ref           = $1;
      my $merged_bowtie = "$path/bowtie_aln/$ref.$TE.bowtie.out";
      my $cmd =
"$scripts/relocaTE_insertionFinder.pl $merged_bowtie $seq_id $genome_path $TE $outregex $exper $flanking_seq_len $existing_TE_path $mismatch_allowance $bowtie2 $relax_reference $relax_align";
      if ( !$parallel ) {
        ##run it now
        print "Finding $TE insertions in $seq_id\n";
        system($cmd);
      }
      else {
        my $shell_dir = "$shellscripts/step_5/$TE";
        make_path( $shell_dir, { mode => 0755 } );
        open OUTSH, ">$shell_dir/$genome_count.$TE.findSites.sh";
        print OUTSH "$cmd\n";
        close OUTSH;
        print PARALLEL "sh $shell_dir/$genome_count.$TE.findSites.sh\n"
          if !$qsub_array;
        chmod 0755, "$shell_dir/$genome_count.$TE.findSites.sh";
      }
      $genome_count++;
    }
    if ($qsub_array) {
      my $end             = $genome_count - 1;
      my $job             = "$shellscripts/step_5.$TE.finder.sh";
      my $existing_depend = '';
      if ($qsub_existingTE_cmd) {
        $existing_depend = "-W depend=afterok:\$EXISTINGTE" if !@depend;
        $existing_depend = ":\$EXISTINGTE" if @depend;
      }
      my $desc = $TE;
      $desc =~ s/\W/_/;
      if ( !@depend ) {
        print QSUBARRAY
"STEP_5_$desc=\`qsub -e $shellscripts -o $shellscripts $qsub_q -t 0-$end $existing_depend $job\`
echo \$STEP_5_$desc\n";
        @depend = ( "STEP_5_$desc", "afterokarray" );
      }
      else {
        my ( $last_job, $afterok ) = @depend;
        @depend = ( "STEP_5_$desc", "afterokarray" );
        my $jobName = $depend[0];
        print QSUBARRAY
"$jobName=`qsub -e $shellscripts -o $shellscripts $qsub_q -t 0-$end -W depend=$afterok:\$$last_job",
          "$existing_depend $job`
echo \$$jobName\n";
      }
      print QSUBARRAY4
        "sh $shellscripts/step_5/$TE/\$PBS_ARRAYID.$TE.findSites.sh";
    }
  }
  if ($qsub_array) {
    close QSUBARRAY3;
    close QSUBARRAY4;
  }
}
## Finished, clean up, cat files
##cat all '.te_insertion_sites.table.txt' results into one file
foreach my $te_path (@te_fastas) {
  my @path     = split '/', $te_path;
  my $te_fasta = pop @path;
  my $TE       = $te_fasta;
  my $path     = join '/', @path;
  $TE =~ s/\.fa//;

  pop @path;
  my $pre_path     = join '/', @path;

  if ($parallel) {

    my $shell_dir = "$shellscripts/step_6/$TE";
    make_path( $shell_dir, { mode => 0755 } );
    open FINISH, ">$shellscripts/step_6/$TE/step_6.$TE.finishing.sh";
    print PARALLEL "sh $shellscripts/step_6/$TE/step_6.$TE.finishing.sh\n"
      if !$qsub_array;
## toDo: Need to put all of this in a perl script then have the
## finishing shell script execute that perl script
    print FINISH "
`mkdir -p $path/results/all_files`

#combine confident insertions to one file
echo \"TE\tTSD\tExper\tinsertion_site\tstrand\tleft_flanking_read_count\tright_flanking_read_count\tleft_flanking_seq\tright_flanking_seq\" > $path/results/temp
for i in \`ls $path/results/*.$TE.confident_nonref_insert.txt\` ; do grep -v flanking_read_count \$i >> $path/results/temp ; done
mv $path/results/temp $path/results/$exper.$TE.confident_nonref.txt
mv $path/results/*.$TE.confident_nonref_insert.txt $path/results/all_files

#combine all insertions to one file
echo \"TE\tTSD\tExper\tchromosome\tinsertion_site\tstrand\tcombined_read_count\tright_flanking_read_count\tleft_flanking_read_count\" > $path/results/temp2
for i in \`ls $path/results/*.$TE.all_nonref_insert.txt\` ; do grep -v total \$i | grep -v Note >> $path/results/temp2 ; done
mv $path/results/temp2 $path/results/$exper.$TE.all_nonref.txt
mv $path/results/*.$TE.all_nonref_insert.txt $path/results/all_files

#combine confident insertions ref seqs to one file
for i in \`ls $path/results/*.$TE.confident_nonref_genomeflank.fa\` ; do cat \$i  >> $path/results/temp3 ; done
mv $path/results/temp3 $path/results/$exper.$TE.confident_nonref_genomeflanks.fa
mv $path/results/*.$TE.confident_nonref_genomeflank.fa $path/results/all_files

#combine confident insertions gff to one file
echo \"##gff-version 3\" > $path/results/temp4
for i in \`ls $path/results/*.$TE.all_insert.gff\` ; do grep -v gff \$i  >> $path/results/temp4 ; done        
mv $path/results/temp4 $path/results/$exper.$TE.all_inserts.gff
mv $path/results/*.$TE.all_insert.gff $path/results/all_files

#combine confident insertions reads to one file
for i in \`ls $path/results/*.$TE.confident_nonref_insert_reads_list.txt\` ; do cat \$i  >> $path/results/temp5 ; done
mv $path/results/temp5 $path/results/$exper.$TE.confident_nonref_reads_list.txt
mv $path/results/*.$TE.confident_nonref_insert_reads_list.txt $path/results/all_files

# move other outfiles somewhere else
if [ -e $pre_path/bowtie-build.out ] ; then 
  mv $pre_path/bowtie-build.out $path/bowtie_aln/.
fi

";
    `chmod +x $shellscripts/step_6/$TE/step_6.$TE.finishing.sh`;
  }
  if ($qsub_array) {
    my $job  = "$shellscripts/step_6/$TE/step_6.$TE.finishing.sh";
    my $desc = $TE;
    $desc =~ s/\W/_/;
    if ( !@depend ) {
      my $jobName = "STEP_6_$desc";
      print QSUBARRAY
        "$jobName=\`qsub -e $shellscripts -o $shellscripts $qsub_q $job\`
echo \$$jobName\n";
      @depend = ( "STEP_6_$desc", "afterok" );
    }
    else {
      my ( $last_job, $afterok ) = ( "STEP_5_$desc", "afterokarray" );
      @depend = ( "STEP_6_$desc", "afterok" );
      my $jobName = $depend[0];
      print QSUBARRAY
"$jobName=`qsub -e $shellscripts -o $shellscripts $qsub_q -W depend=$afterok:\$$last_job $job`
echo \$$jobName\n";
    }
  }
  if ( !$parallel and !$qsub_array ) {
    ##do it now
    ##combine and delete individual chr files for confident sites
    print "Finishing and cleaning up\n";
`echo \"TE\tTSD\tExper\tinsertion_site\tstrand\tleft_flanking_read_count\tright_flanking_read_count\tleft_flanking_seq\tright_flanking_seq\" > $path/results/temp`;
    my @files = `ls $path/results/*.$TE.confident_nonref_insert.txt`;
    foreach my $file (@files) {
      chomp $file;
      `grep -v flanking_read_count $file  >> $path/results/temp`;
      unlink $file;
    }
    `mv $path/results/temp $path/results/$exper.$TE.confident_nonref.txt`;

    ##combine and delete individual chr files for all sites
`echo \"TE\tTSD\tExper\tchromosome\tinsertion_site\tstrand\tcombined_read_count\tright_flanking_read_count\tleft_flanking_read_count\" > $path/results/temp2`;
    @files = `ls $path/results/*.$TE.all_nonref_insert.txt`;
    foreach my $file (@files) {
      chomp $file;
      `grep -v total $file | grep -v Note  >> $path/results/temp2`;
      unlink $file;
    }
    `mv $path/results/temp2 $path/results/$exper.$TE.all_nonref.txt`;

    ##combine and delete individual chr fasta files
    @files = `ls $path/results/*.$TE.confident_nonref_genomeflank.fa`;
    foreach my $file (@files) {
      chomp $file;
      `cat $file >> $path/results/temp3`;
      unlink $file;
    }
`mv $path/results/temp3 $path/results/$exper.$TE.confident_nonref_genomeflanks.fa`;

    ##combine and delete individual chr gff files
    `echo \"##gff-version 3\" > $path/results/temp4`;
    @files = `ls $path/results/*.$TE.all_insert.gff`;
    foreach my $file (@files) {
      chomp $file;
      `grep -v gff $file >> $path/results/temp4`;
      unlink $file;
    }
    `mv $path/results/temp4 $path/results/$exper.$TE.all_inserts.gff`;

    ##combine and delete individual chr reads list
    @files = `ls $path/results/*.$TE.confident_nonref_insert_reads_list.txt`;
    foreach my $file (@files) {
      chomp $file;
      `cat $file >> $path/results/temp5`;
      unlink $file;
    }
    # move other outfiles somewhere else
if (-e "$pre_path/bowtie-build.out" ){ 
   `mv $pre_path/bowtie-build.out $path/bowtie_aln/.`;
}
  if ( -d  "$pre_path/shellscripts"){
     `rm -rf $pre_path/shellscripts`;
  }

`mv $path/results/temp5 $path/results/$exper.$TE.confident_nonref_reads_list.txt`;

    # move other outfiles somewhere else
if (-e "$pre_path/bowtie-build.out" ){ 
   `mv $pre_path/bowtie-build.out $path/bowtie_aln/.`;
}
  if ( -d  "$pre_path/shellscripts"){
     `rm -rf $pre_path/shellscripts`;
   }


    print "$TE results are found in $path/results\n";
  }
  close FINISH;
}

if ($qsub_array) {
  close QSUBARRAY;
## this would happen before IO was finished on the file
  #  system ("qsub $qsub_q $current_dir/$top_dir/run_these_jobs.sh");
  print "$current_dir/$top_dir/run_these_jobs.sh was created
-- Run this script, \'sh $current_dir/$top_dir/run_these_jobs.sh\' 
-- This script will submit all jobs to the queue in the appropriate order.
-- Be sure to check the error files in $current_dir/$top_dir/shellscripts. They should all be file size 0.
\n";
}
elsif ($parallel) {

  #system (sort "$current_dir/$top_dir/run_these_jobs.sh");
  print
    "Run each command line statement in $current_dir/$top_dir/run_these_jobs.sh
--Run these in order (step_1,step_2,step_3, so on) for each TE.
--For example, all the step_3 scripts for a specific TE should be successfully completed (finished without errors) 
before running a step_4 script of the same TE.
--All scripts of the same step can be run in parallel (at the same time).
\n";
}
