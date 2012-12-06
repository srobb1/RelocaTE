#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;
use File::Basename;
use FindBin qw($RealBin);
## trims and filters reads, aligns to genome, splits by target.
## produces fq, sam and bam files for each individual target
my $scripts_dir = $RealBin;
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

my $dir = '.';
my $genomeFasta;
my ( $minLength, $minQuality, $minPercent ) = ( 50, 20, 80 );
my $Q             = 33;
my $insertLength  = 500;
my $mate_1_id     = "_p1";
my $mate_2_id     = "_p2";
my $unpaired_id     = ".unPaired";
my $split         = 0;
my $filter_trim   = 1;
my $tempDir       = $current_dir;
my $bin_per_chrom = 1;
my $prefix = 'none';
my $bwa_quality = 10;

GetOptions(
  'd|dir:s'           => \$dir,
  '1|mate_1_id:s'     => \$mate_1_id,
  '2|mate_2_id:s'     => \$mate_2_id,
  'u|unpaired_id:s'     => \$unpaired_id,
  'g|genomeFasta:s'   => \$genomeFasta,
  'l|minLength:i'     => \$minLength,
  'q|minQuality:i'    => \$minQuality,
  'p|minPercent:i'    => \$minPercent,
  's|quality:i'       => \$Q,
  'i|insertLength:i'  => \$insertLength,
  'x|split:i'         => \$split,
  'f|filter_trim:i'   => \$filter_trim,
  't|tempDir:s'       => \$tempDir,
  'b|bin_per_chrom:i' => \$bin_per_chrom,
  'r|prefix:s'        => \$prefix,
  'bq|bwa_quality:i'  => \$bwa_quality,
  'h|help'            => \&getHelp,
);

$prefix  = !defined $prefix ? ''       : $prefix . '.';
if ( !defined $genomeFasta ) {
  print "\n\nPlease provide reference genome by using -g Genome fasta path\n";
  &getHelp();
}

sub getHelp {
  print "
usage:
./raw_paired_reads_2_split_by_target.pl [-d fq_file_directory] [-1 mate_pair_file_1_id][-2 mate_pair_file_2_id][-l minLength] [-q minQuality] [-p minPercent] [-s quality_offset] [-i insert_size] [-h] 

options:
-d STR		directory of original raw fq files (must be .fq or .fastq) [.]
-r STR		prefix for merged bam and fq files (ex. A123) [none]
-g STR		genome fasta file path [no default]
-l INT		min length for fastq_quality_trimmer [50]
-q INT		min quality score for fastq_quality_trimmer [20]
-p INT		min percent for fastq_quality_filter [80]
-s INT		quality score offset type Sanger(33) or Illumina(64) [33]
-i INT		insert library length [500]
-1 STR		file containing mate 1 id (ex reads_1.fq). If not paired end use 'NONE'. [_p1]
-2 STR		file containing mate 2 id (ex reads_2.fq). If not paried end use 'NONE'. [_p2]
-u STR		file containing unparied (ex reads.unPaired.fq). If no id use 'NONE' (ex reads.fq) . [.unPaired]
-b INT	        split the sam and bam files and organize by chromosome yes=1 no=0 [1]	
-x INT	        split fq file into smaller files (1,000,000/file) yes=1 no=0 [0]	
-f INT	        run fastq_quality_filter and fastq_quality_trimmer yes=1 no=0 [1]	
-t STR		location to create temp directories [/dev/shm]
-h 		this message
";

  exit 1;
}

my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
  localtime(time);
$year += 1900;
my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );

my $date = "$year-$abbr[$mon]-$mday";    ## gives 2012-Dec-17

unless ( -e $genomeFasta ) {
  print "$genomeFasta does not exist. Check file name.\n";
  &getHelp();
}

my $log_file    = "$current_dir/$date-log.txt";
`touch $log_file`;

open JOBS_SH, ">$current_dir/$date.process_raw_reads.sh";
print JOBS_SH "sh $current_dir/p\$PBS_ARRAYID.process_raw_reads.sh";
open ARRAY_SH, ">$current_dir/$date.arrayJob.sh";

my $genome_path = File::Spec->rel2abs($genomeFasta);
my $bwa_c_switch = '';
if ($Q eq 64){
  $bwa_c_switch = '-c';
}
my $bwa_index = "bwa index $bwa_c_switch -a bwtsw $genome_path";
unless ( -e "$genome_path.rsa" ) {
  open GINDEX, ">$current_dir/genome_indexing.sh";
  print GINDEX "#!/bin/bash\n\n";
  print GINDEX "echo \"$bwa_index\"\n";
  print GINDEX "$bwa_index\n";
  print GINDEX "chmod ago+rwx $genome_path*\n";
  close GINDEX;
  print "\n\n !! Make sure that you run \"sh genome_indexing.sh\" before you process the fastq files !!\n\n";
}
my $dir_path = File::Spec->rel2abs($dir);
my ($mate_1_path , $mate_2_path);
my (@filelist_1, @filelist_2);
my $paired = 1;
if ($mate_1_id =~ /NONE/i){
  ##no paired end data, single end only
  $mate_1_path = "$dir_path/*.f*q";
  @filelist_1  = < $mate_1_path >;
  @filelist_2  = ();
  $paired = 0;
  $mate_1_id = '';
  $mate_2_id= '';
}else{
  ##check to make sure that files can be found with the mate patterns provided
  $mate_1_path = "$dir_path/*$mate_1_id.f*q";
  $mate_2_path = "$dir_path/*$mate_2_id.f*q";
  @filelist_1  = < $mate_1_path >;
  @filelist_2  = < $mate_2_path >;
}
if ($unpaired_id =~ /NONE/i){
  $unpaired_id = '';
}
if ( !@filelist_1 ) {
  print
"Cannot find any files in $dir_path that are similar to your mate_file_1 pattern $mate_1_id\n";
  &getHelp();
}
if ( $paired  and !@filelist_2 ) {
  print
"Cannot find any files in $dir_path that are similar to your mate_file_2 pattern $mate_2_id\n";
  &getHelp();
}
my $bwa_I_switch = '';
if ($Q eq 64){
  $bwa_I_switch = '-I';  
}

my %files;
##split into smaller files
if ($split) {
  mkdir "$current_dir/split_by_number_fq"
    unless -d "$current_dir/split_by_number_fq";
  opendir( DIR, $dir_path ) || die "$!";
  foreach my $file ( readdir(DIR) ) {
    my ( $volume, $directories, $filename ) = File::Spec->splitpath($file);
    next unless ( $filename =~ /\.(fastq|fq)$/ );
    print "\n\nI am spliting your fastq files right now into many smaller fastq files.
We only have to do this one time.
This might take a little bit:\n\n";
`$scripts_dir/fastq_split.pl -s 1000000 -o split_by_number_fq/ $dir_path/$file`;
  }
}
if ($split) {
  $dir_path = "$current_dir/split_by_number_fq";
}
opendir( DIR, $dir_path ) || die "Can't Open $dir_path $!";
my $fq_ext;
foreach my $file ( readdir(DIR) ) {
  my ( $volume, $directories, $filename ) = File::Spec->splitpath($file);
  my ( $filename_base, $sampleName, $pairID, $suffix );
  if ( $paired ){
    next unless ( $filename =~ /((\S+)($mate_1_id|$mate_2_id))\.(fastq|fq)$/ );
    ( $filename_base, $sampleName, $pairID, $suffix ) = ( $1, $2, $3, $4 );
  }
  else {
    next unless $filename =~ /(fastq|fq)$/;
    if ($filename =~ /(\S+)\.(fastq|fq)$/){ 
      ( $filename_base , $suffix ) = ( $1, $2);
      $sampleName = $filename_base;
    }
  }
  push @${ $files{$sampleName} }, $filename_base;
  $fq_ext = $suffix;
}

my $ext = ".fq";

my $count = 0;
foreach my $sample ( sort keys %files ) {
  open OUTFILE, ">$current_dir/p$count.process_raw_reads.sh";
  print OUTFILE "#!/bin/bash\n\n";
  print OUTFILE "compute_node=`hostname`\n";
  print OUTFILE "sample=$sample\n";
  print OUTFILE "echo \"$sample \$compute_node \" >> $log_file\n";
  print OUTFILE "umask 002\n";
  my (
    @trim_filter, @clean,               @aln,
    @sam,         @split_sam_by_target, @sam2fq,
    @sam2bam,     @mergeBam,            @catfq
  );

  print OUTFILE "
if [ -d $tempDir ]; then
 tmp_dir=`mktemp --tmpdir=$tempDir -d`
else
 tmp_dir=`mktemp --tmpdir=/tmp -d`
fi\n";
  print OUTFILE "cd \$tmp_dir\n";

  #for unpaired, unlabeled in name files reads.fq
  print OUTFILE "for i in \`ls $dir_path/$sample$unpaired_id.f*q\` ; do ln -s \$i \$tmp_dir/$sample.unPaired.fq ; done\n" if $paired ;

  #foreach single file write the trim and filter and the aln commands
  foreach my $file ( sort @${ $files{$sample} } ) {
    if ($filter_trim) {
      push @trim_filter,
"fastq_quality_trimmer -Q$Q -l $minLength -t $minQuality -i $dir_path/$file.$fq_ext |fastq_quality_filter -Q$Q -q $minQuality -p $minPercent -v -o \$tmp_dir/$file$ext";
      push @aln,
"bwa aln -t 8 $bwa_I_switch $genome_path \$tmp_dir/$file$ext > \$tmp_dir/$file.sai";
    }
    else {
      print OUTFILE "ln -s $dir_path/$file.$fq_ext \$tmp_dir/$file$ext\n";
      push @aln,
"bwa aln -t 8 -q $bwa_quality $bwa_I_switch $genome_path \$tmp_dir/$file$ext > \$tmp_dir/$file.sai" if !$paired;
      push @aln,
"bwa aln -t 8 -q $bwa_quality $bwa_I_switch $genome_path \$tmp_dir/$file$ext > \$tmp_dir/$file.sai" if $paired;
    }
  }

  #foreach potential paired sample write the following commands
  my ( $pair1, $pair2 );
  my $pairs = @${ $files{$sample} };
  if ( $pairs == 2 ) {
    ( $pair1, $pair2 ) = sort @${ $files{$sample} };
    push @clean,
"$scripts_dir/clean_pairs.pl -1 \$tmp_dir/$pair1$ext -2 \$tmp_dir/$pair2$ext >> \$tmp_dir/$sample.unPaired.fq";
    push @clean, "mv \$tmp_dir/$pair1.matched$ext \$tmp_dir/$pair1$ext";
    push @clean, "mv \$tmp_dir/$pair2.matched$ext \$tmp_dir/$pair2$ext";
    ##after cleaning the 2 paired files a file of unPaired reads is generated
    ##run bwa aln on this file
    ##and bwa samse
    push @clean,
"if [ -s \$tmp_dir/$sample.unPaired.fq ] ; then bwa aln -t 8 $bwa_I_switch $genome_path \$tmp_dir/$sample.unPaired.fq > \$tmp_dir/$sample.unPaired.sai ; fi";
    push @clean,
"if [ -e \$tmp_dir/$sample.unPaired.sai ] ; then bwa samse  $genome_path \$tmp_dir/$sample.unPaired.sai \$tmp_dir/$sample.unPaired.fq   > \$tmp_dir/$sample.unPaired.sam ; fi" if $paired;
    push @split_sam_by_target,
"if [ -e \$tmp_dir/$sample.unPaired.sam ] ; then $scripts_dir/splitSam_byTarget.pl -s \$tmp_dir/$sample.unPaired.sam ; fi"
      if $bin_per_chrom;
    push @sam,
"bwa sampe -a $insertLength $genome_path \$tmp_dir/$pair1.sai \$tmp_dir/$pair2.sai \$tmp_dir/$pair1$ext \$tmp_dir/$pair2$ext  > \$tmp_dir/$sample.sam";
    if ($bin_per_chrom) {
      push @split_sam_by_target,
        "$scripts_dir/splitSam_byTarget.pl -s \$tmp_dir/$sample.sam";
      push @sam2fq,
"for i in `ls \$tmp_dir/split_by_target` ; do $scripts_dir/sam2fq.pl \$tmp_dir/split_by_target/\$i ; done";
      push @sam2bam,
"for i in `ls \$tmp_dir/split_by_target` ; do $scripts_dir/sam2bam.pl \$tmp_dir/split_by_target/\$i $genome_path $sample; done";
    }
  }
  elsif ( $pairs == 1 ) {
    ( $pair1, $pair2 ) = sort @${ $files{$sample} };
    push @sam,
"bwa samse $genome_path \$tmp_dir/$pair1.sai \$tmp_dir/$pair1$ext   > \$tmp_dir/$sample.sam";
    push @split_sam_by_target,
      "$scripts_dir/splitSam_byTarget.pl -s \$tmp_dir/$sample.sam" if $bin_per_chrom;
  }
  else {
    warn
      "error: $sample has $pairs. This sample should have 2 pairs or just 1.\n";
  }

  foreach my $trim_filter (@trim_filter) {
    print OUTFILE "echo \"$sample start trim\"\n";
    print OUTFILE "$trim_filter\n\n";
    print OUTFILE "echo \"$sample end trim\"\n\n";
  }
  foreach my $clean (@clean) {
    print OUTFILE "echo \"$sample start clean\"\n";
    print OUTFILE "$clean\n\n";
    print OUTFILE "echo \"$sample end clean\"\n\n";
  }
  foreach my $aln (@aln) {
    print OUTFILE "echo \"$sample start bwa aln\"\n";
    print OUTFILE "$aln\n\n";
    print OUTFILE "echo \"$sample end bwa aln\"\n\n";
  }
  foreach my $sam (@sam) {
    print OUTFILE "echo \"$sample start bwa sam\"\n";
    print OUTFILE "$sam\n\n";
    print OUTFILE "echo \"$sample end bwa sam\"\n\n";
  }
  foreach my $line (@split_sam_by_target) {
    print OUTFILE "echo \"$sample start split by target\"\n";
    print OUTFILE "$line\n\n";
    print OUTFILE "echo \"$sample end split by target\"\n\n";
  }
  foreach my $line (@sam2fq) {
    print OUTFILE "echo \"$sample start sam2fq\"\n";
    print OUTFILE "$line\n\n";
    print OUTFILE "echo \"$sample end sam2fq\"\n\n";
  }
  foreach my $line (@sam2bam) {
    print OUTFILE "echo \"$sample start sam2bam\"\n";
    print OUTFILE "$line\n\n";
    print OUTFILE "echo \"$sample end sam2bam\"\n\n";
  }
  print OUTFILE
"if [ ! -d \"$current_dir/bam_for_all_reads\" ] ; then mkdir -m 0775 $current_dir/bam_for_all_reads ; fi\n";

print OUTFILE "if [ -e \$tmp_dir/$sample.unPaired.sam ] ; then 
samtools view -h -b -S -T $genome_path \$tmp_dir/$sample.unPaired.sam >  \$tmp_dir/$sample.unPaired.bam
samtools sort  \$tmp_dir/$sample.unPaired.bam  \$tmp_dir/$sample.unPaired.sorted
samtools index  \$tmp_dir/$sample.unPaired.sorted.bam
cp \$tmp_dir/$sample.unPaired.sorted.bam* $current_dir/bam_for_all_reads/. ; fi \n";

  print OUTFILE
"if [ ! -f \"\$tmp_dir/$sample.sam\" ] ; then echo \" $sample.sam does not exist\" ;fi\n";
  print OUTFILE "echo \"$sample start samtools\"\n";
  print OUTFILE "echo \"$sample samtools:view\"\n\n";
  print OUTFILE
"samtools view -h -b -S -T $genome_path \$tmp_dir/$sample.sam >  \$tmp_dir/$sample.bam\n";
  print OUTFILE "echo \"$sample start samtools:sort\"\n";
  print OUTFILE
    "samtools sort  \$tmp_dir/$sample.bam  \$tmp_dir/$sample.sorted\n";
  print OUTFILE "echo \"$sample samtools:index\"\n";
  print OUTFILE "samtools index  \$tmp_dir/$sample.sorted.bam\n";
  print OUTFILE "echo \"$sample end samtools\"\n\n";
  print OUTFILE
"echo \"$sample start cp \$tmp_dir/$sample.sorted.bam $current_dir/bam_for_all_reads \"\n";
  print OUTFILE
    "cp \$tmp_dir/$sample.sorted.bam* $current_dir/bam_for_all_reads/.\n";
  print OUTFILE
"echo \"$sample end cp \$tmp_dir/$sample.sorted.bam $current_dir/bam_for_all_reads \"\n\n";
  print OUTFILE
"if [ -f \$tmp_dir/$sample.sorted.bam ] && [ ! -f $current_dir/bam_for_all_reads/$sample.sorted.bam ] ; then echo \" $sample.sorted.bam did not copy to $current_dir/bam_for_all_reads\" ; fi\n";

  if ($filter_trim) {
    print OUTFILE "mkdir -m 0775 -p $current_dir/fq_split_by_number_filtered\n";
    print OUTFILE
"cp \$tmp_dir/*.fq \$tmp_dir/*unPaired.fq $current_dir/fq_split_by_number_filtered\n";
  }
  if ($bin_per_chrom) {
    print OUTFILE "##make directories
mkdir -m 0775 -p $current_dir/fq_split_by_chromosome
mkdir -m 0775 -p $current_dir/bam_split_by_chromosome
for i in `ls \$tmp_dir/split_by_target` ; do touch \$tmp_dir/list.txt ; echo \$i >> \$tmp_dir/list.txt ; done
for i in `cat \$tmp_dir/list.txt` ; do mkdir -m 0775 -p $current_dir/bam_split_by_chromosome/\$i ; done
for i in `cat \$tmp_dir/list.txt` ; do mkdir -m 0775 -p $current_dir/fq_split_by_chromosome/\$i ; done

##move bam to user dir
for i in `cat \$tmp_dir/list.txt` ; do cp  \$tmp_dir/split_by_target/\$i/*\$i*sorted.bam $current_dir/bam_split_by_chromosome/\$i/. ; done

##move fq files
for i in `cat \$tmp_dir/list.txt` ; do cp  \$tmp_dir/split_by_target/\$i/*fq $current_dir/fq_split_by_chromosome/\$i/. ; done
#for i in `cat \$tmp_dir/list.txt` ; do  if [ -e \$tmp_dir/split_by_target/\$i/*\$i.unPaired.fq ] ; then cp  \$tmp_dir/split_by_target/\$i/*\$i.unPaired.fq $current_dir/fq_split_by_chromosome/\$i. ; fi ; done

##move list.txt
if [ ! -e $current_dir/$date-targets.list.txt ] ; then cp \$tmp_dir/list.txt $current_dir/$date-targets.list.txt ; fi 

";
  }
  print OUTFILE "cd $current_dir\n";
  print OUTFILE "rm -rf \$tmp_dir\n";
  print OUTFILE "echo \"$sample end_of_commands\" >> $log_file\n";

  $count++;
}
$count--;
print ARRAY_SH "qsub -t 0-$count $current_dir/$date.process_raw_reads.sh";

open DOLAST_SH , ">$current_dir/$date-do_last.sh";
$prefix =~ s/\.$//;
print DOLAST_SH "
for i in `cat $current_dir/$date-targets.list.txt` ; do $scripts_dir/cat_fq_shell.pl $current_dir/fq_split_by_chromosome/\$i $prefix $tempDir ; done
for i in `cat $current_dir/$date-targets.list.txt` ; do $scripts_dir/merge_bam_shell.pl $current_dir/bam_split_by_chromosome/\$i $prefix $tempDir ; done


";

print "

Ok, here we go ...

  -- To process the raw fastq files either:
	Run each of the newly created p#.process.raw.reads.sh scripts
        \"sh p0.process_raw_reads.sh\"

        --OR-- if you have PBS

        Run the date.arrayJob.sh script.
        \"sh $date.arrayJob.sh\"  

        You might have to add the -q queueName to this script.

  -- Check for errors. 
       Errors often arise if you do not have enough disk space available
       Each process outputs error messages to the screen or to qsub error files if using PBS qsub

  -- Now you can combine all the small files back into 1 file per chromosome or ref seq
   	
       Run \"sh $date-do-last.sh\"
       This will create more shell scritps (like the ones listed below) that you can run one or more at a time.
        - Chr1.mergeBam.sh
        - Chr2.cat_fq.sh
        - and many more
  
  -- If there are no errors in the combination of bam files or of fq files, you can 
         delete all the directories inside the bam_split_by_target and fq_split_by_target directories

  *** If this is your first time, make sure to copy this info somewhere so that 
  *** you can refer back to it


";
