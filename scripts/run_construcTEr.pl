#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use Cwd;
use strict;
use Data::Dumper;
use FindBin qw($RealBin);

##change $scripts to location of relocaTE scripts
my $scripts = $RealBin;


if ( !defined @ARGV ) {
  &getHelp();
}
my $genomeFasta;
my $te_fasta;
my $fq_dir;
my $te_fq_dir;
my $workingdir;
my $outdir   = 'construcTEr';
my $parallel = 1;
my $arrayJob = 1;
my $mate_file_1        = '_1\D*?fq';
my $mate_file_2        = '_2\D*?fq';
my $blat_dir;
#my $unpaired        = 'unPaired\D*?fq';
my $insert_file = 0;
GetOptions(
  'p|parallel:i'    => \$parallel,
  'a|arrayJob:i'    => \$arrayJob,
  'i|insert_file:s' => \$insert_file,    
  'b|blat_dir:s' => \$blat_dir,    
  'w|workingdir:s'  => \$workingdir,
  'o|outdir:s'      => \$outdir,
  'd|fq_dir:s'      => \$fq_dir,
  'g|genomeFasta:s' => \$genomeFasta,
  't|te_fasta:s'    => \$te_fasta,
  '1|mate_1_id:s'   => \$mate_file_1,
  '2|mate_2_id:s'   => \$mate_file_2,
#  'u|unpaired_id:s'   => \$unpaired,
  'h|help'          => \&getHelp,
);

=cut
print "
-p $parallel
-i $insert_file
-w $workingdir
-o $outdir
-d $fq_dir
-t $te_fasta
-g $genomeFasta
-1 $mate_file_1
-2 $mate_file_2

";
=cut

my $current_dir;

if ( defined $workingdir and -d $workingdir ) {
  $current_dir = File::Spec->rel2abs($workingdir);
  $current_dir =~ s/\/$//;
}
else {
  $current_dir = cwd();
}
if ( !defined $genomeFasta ) {
  print "\n\nPlease provide reference genome by using -g Genome fasta path\n";
  &getHelp();
} elsif ( !-e $genomeFasta ) {
  print "$genomeFasta does not exist. Check file name.\n";
  &getHelp();
}
if ( !defined $te_fasta ) {
  print
"\n\nPlease provide fasta file containing transposable elements by using -t TE fasta path\n";
  &getHelp();
}elsif ( !-e $te_fasta ) {
  print "$te_fasta does not exist. Check file name.\n";
  &getHelp();
}
else {
  my $first_line = `head -n1 $te_fasta`;
  if ( $first_line !~ /^>\S/ ) {
    die
"The TE_fasta:$te_fasta does not have the proper fasta format\n";
  }
}
if ( $insert_file and !-e $insert_file) {
  print "$insert_file does not exist. Check file name.\n";
  &getHelp();
}
if ( defined $blat_dir and -d $blat_dir ) {
  $blat_dir = File::Spec->rel2abs($blat_dir);
  $blat_dir =~ s/\/$//;
}else {
  print "$blat_dir does not exist, need to run blats\n";
}
if ( !defined $fq_dir ) {
  print "\n\nPlease provide a directory of complete set of paired fastq files\n";
  &getHelp();
}
elsif ( !-d $fq_dir ) {
  print
"\n\nCheck the spelling or location of $fq_dir, Please provide a directory of paired fastq files\n";
  &getHelp();
}

sub getHelp {
  print ' 
usage:
./construcTEr.pl [-t TE_fasta_file][-g chromosome_genome_fasta][-d dir_of_fq][-h] 

options:

**required:
-g STR          genome fasta file path. If not provided will only align reads to TE and remove TE seq from short reads. [no default]
-t STR          fasta containing nucleotide sequences of transposable elements with TSD=xxx in the desc. [no default]
-d STR          directory of paired and unpaired fastq files (.fq or .fastq is acceptable)  [no default]
-1 STR          regular expression to identify mate 1 paired files [\'_1\D*?fq\']
-2 STR          regular expression to identify mate 2 paired files [\'_2\D*?fq\']

**recommended:
## maybe, i removed this: -f STR		directory of TE subdirectory of relocaTE output ex: /home/usr/relocaTE_output/myTE/ [no default]

**optional:
-i STR          name of the file that contains insert positions,relocaTE:all.$te.te_insertion_sites.table.txt (TE<tab>.<tab>ref<tab>pos)
-o STR          name for directory to contain output directories and files, will be created for the run (ex. 04222012_A123) [construcTEr]
-p INT          run each genome sequence separetly, parallel. The alternative (0) would be to run one after the other (int, 0=false or 1=true) [1] 
-a INT          create an array job script for the shellscripts created with -p 1? (int, 0=false or 1=true) [1] 
-w STR          base working directory, needs to exist, will not create, full path [cwd] 
-h              this message

SAMPLE TE FASTA: with any number (0 or more) custom defined ranges
>mping TSD=TTA range=TIR1:1..15 range=TIR2:416..430 range=other:yourStart..yourEnd
GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATG
ATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTT
TCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGT
CCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAATGATCCCAGCAA
CTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGT
TTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATT
GTGACTGGCC

';

  exit 1;
}
$outdir =~ s/\/$//;
my $te_path = File::Spec->rel2abs($te_fasta);
my $top_dir = $outdir;

## get fq files
$fq_dir = File::Spec->rel2abs($fq_dir);
my @fq_files = <$fq_dir/*.fq>;
push @fq_files , <$fq_dir/*.fastq>;
## create bowtie index
my $genome_file = File::Spec->rel2abs($genomeFasta);
if ( !-e "$genome_file.bowtie_build_index.1.ebwt" ) {
        `bowtie-build -f $genome_file $genome_file.bowtie_build_index`;
}

## convert fq files to fa for blat and seq retrieval
my @fa;
my $shell_script_count = 0;
my $shell_dir; 
if ($parallel){
  $shell_script_count++ ;
  $shell_dir = "$current_dir/$top_dir/shellscripts/step_$shell_script_count";
  `mkdir -p $shell_dir`;
}
if ($arrayJob and $parallel){
  open ARRAY , ">$current_dir/$top_dir/run.these.jobs.sh" or die "Can't open $current_dir/$top_dir/run.these.jobs.sh for writing\n";
}
my $file_count = 0;
my $format = 0;
my $fq2fa = 0;
foreach my $fq (@fq_files) {
  next if $fq =~ /unparied/i;
  my $cmd;
  my $fq_path = File::Spec->rel2abs($fq);
  my $fa = $fq;
  my @fq_path   = split '/', $fq_path;
  my $fq_name   = pop @fq_path;
  if ( $fa =~ s/\.(fq|fastq)$/.fa/ ) {
    push @fa, $fa;
    if ( !-e $fa ) {
        $fq2fa = 1;
        $cmd = "$scripts/relocaTE_fq2fa.pl $fq_path $fa";
        if ($parallel){
          my $outsh = "$shell_dir/$file_count" . "fq2fa.sh";
          open OUTSH, ">$outsh";
          print OUTSH "$cmd\n";
          close OUTSH;
         }else {
            `$cmd`;
         }
     }
     if (!-e "$fa.nin"){
        $format = 1;
        $cmd = "formatdb -i $fa -p F -o T";
       if ($parallel){
         my $outsh = "$shell_dir/$file_count" . ".fa4formatdb.sh";
         open OUTSH, ">$outsh" or die "Can't open $outsh for writting $!\n";
         print OUTSH "$cmd\n";
         close OUTSH;
       }else {
        `$cmd`;
       }
     }
  }
  else {
    print
"$fq does not seem to be a fastq based on the file extension. It should be fq or fastq\n";
    &getHelp();
  }
  $file_count++;
}


if ($parallel and $arrayJob and ($format or $fq2fa)){
  #$shell_script_count++;
  my $sh = "$shell_dir/run_construcTEr.step_" . $shell_script_count . ".sh";
  open SH , ">$sh";
  print SH "sh $shell_dir/\$PBS_ARRAYID.fa4formatdb.sh\n";
  print ARRAY "qsub -t 0-" , $file_count-1 , " $shell_dir/run_construcTEr.step_", $shell_script_count ,".sh\n";
  close SH;
}

##split TE fasta into single record fastas
my @te_fastas;

##split the TE fasta of many seqs into individual files
open( INFASTA, "$te_fasta" ) || die "$!\n";
my $i = 0;
while ( my $line = <INFASTA> ) {
  if ( $line =~ /^>(\S+)/ ) {
    my $id = $1;
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
  elsif ( $line =~ /^>/  ) {
    die
"The TE_fasta:$te_fasta does not have the proper format:\n>TE_NAME\nnSEQUENCE\n";
  }
  else {    ##should be sequence
    print OUTFASTA $line;
  }
}
close(INFASTA);
close(OUTFASTA);


#foreach TE fasta blat against target chromosome and parse and find insertion sites
if (!defined $blat_dir){
if ($parallel){
  $shell_script_count++;
}
foreach my $te_path (@te_fastas) {
  my @path     = split '/', $te_path;
  my $te_fasta = pop @path;
  my $path     = join '/', @path;
  my $TE       = $te_fasta;
  $TE =~ s/\.fa//;
  `mkdir -p $path/blat_output`;
   if ($parallel) {
    $shell_dir = "$current_dir/$top_dir/shellscripts/step_$shell_script_count/$TE";
    `mkdir -p $shell_dir`;
   }

  #blat fa files against te.fa
  my $file_count = 0;
  foreach my $fa ( @fa ) {
    #remove and save filename part of path
    my @fa_path = split '/', $fa;
    my $fa_name = pop @fa_path;
    $fa_name =~ s/\.fa$//;
    if ($parallel) {
      #$shell_dir = "$current_dir/$top_dir/shellscripts/step_$shell_script_count/$TE";
     #`mkdir -p $shell_dir`;
      open OUTSH, ">$shell_dir/$file_count.blat.sh";
    }
    if ( !-e "$path/blat_output/$fa_name.te_$TE.blatout" ) {
      my $cmd =
"blat -minScore=10 -tileSize=7 $te_path $fa $path/blat_output/$fa_name.te_$TE.blatout";
      print OUTSH "$cmd\n" if $parallel;
      `$cmd` if !$parallel;
    }
    $file_count++;
  }
  if ($parallel and $arrayJob){
    #$shell_script_count++;
    my $sh = "$shell_dir/run_construcTEr.$TE.step_" . $shell_script_count .".sh";
    open SH , ">$sh";
    print SH "sh $shell_dir/\$PBS_ARRAYID.blat.sh\n";
    print ARRAY "qsub -t 0-" , $file_count-1 , " $shell_dir/run_construcTEr.$TE.step_", $shell_script_count ,".sh\n";
    close SH;
  }
}
}

if ($parallel){
  $shell_script_count++;
}
foreach my $te_path (@te_fastas) {
  my @te_path = split '/' , $te_path;
  pop @te_path;
  my $TE = pop @te_path; 
  my $outregex   = "$current_dir/$top_dir/$TE/$TE.regex.txt";
  open OUTREGEX, ">$outregex" or die $!;
  print OUTREGEX  "$mate_file_1\t$mate_file_2";
  my $cmd =  "$scripts/construcTEr.pl $fq_dir $genome_file $te_path $current_dir/$top_dir/$TE $outregex $insert_file $blat_dir";
  if ($parallel) {
    $shell_dir = "$current_dir/$top_dir/shellscripts/step_$shell_script_count/$TE";
    `mkdir -p $shell_dir`;
    open OUTSH, ">$shell_dir/$TE.construcTEr.sh" or die "Can't open $shell_dir/$TE.construcTEr.sh";
    print OUTSH $cmd;
  }else{
    `$cmd`;
  }
  if ($parallel and $arrayJob){
    print ARRAY "qsub $shell_dir/$TE.construcTEr.sh";
  }
}
