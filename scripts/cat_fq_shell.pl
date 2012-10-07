#!/usr/bin/perl -w
use strict;
use File::Spec;
use FindBin qw($RealBin);

##change $scripts to location of relocaTE scripts
my $scripts = $RealBin;

#combines all the fq files in one directory. Create shell scripts that can be run individually or in parallel

## for multiple directories try this:
## if directories are Chr1, Chr2, Chr3 ... Chr12
## for i in `seq 1 12` ; do cat_fq_shell.pl Chr$i prefix tempdir; done

## updated 03/31/2012: changed the output files from _1.fq to _p1.fq and _2.fq to _p2.fq

if ( !defined @ARGV ) {
  die
"for a single directory: 
            cat_fq_shell.pl dir_of_fq_files prefix [tempDir:/scratch] [clean 1|0 default:1]
for multiple directories:
            for i in `ls` ; do cat_fq_shell.pl \$i prefix [tempDir:/scratch] [clean 1|0 default:1]; done
or if you have numbered ref seqs, ex Chr1, Chr2 ... Chr12:
            for i in `seq 1 12` ; do cat_fq_shell.pl Chr\$i prefix [tempDir:/scratch] [clean 1|0 default:1]; done

Note: \'clean\' means that the resulting fq files will be matched, and any unparied seqs will be put into an unParied.fq file

Then run each shell script individually or in parallel\n\n";

}

my $dir     = shift;
my $prefix  = shift;
my $tempDir = shift;
my $clean   = shift;
$prefix = !defined $prefix ? '' : $prefix . '.';
$clean  = !defined $clean  ? 1  : 0;
my $dir_path = File::Spec->rel2abs($dir);

$tempDir = defined $tempDir ? $tempDir : '/scratch';
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

my @dirs       = split '/', $dir_path;
my $lowest_dir = pop @dirs;
my $one_up     = join '/', @dirs;

open SH, ">$current_dir/$lowest_dir.cat_fq.sh";
print SH "#!/bin/bash\n\n";

print SH "echo mktemp\n";
print SH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
print SH "echo \"tmp_dir=\$tmp_dir\"\n";
print SH "cd \$tmp_dir\n";

my $mate_1   = "$lowest_dir" . "_p1.fq";
my $mate_2   = "$lowest_dir" . "_p2.fq";
my $unpaired = "$lowest_dir" . ".unPaired.fq";

print SH "cat $dir_path/*_p1.fq > \$tmp_dir/$mate_1\n";
print SH "cat $dir_path/*_p2.fq > \$tmp_dir/$mate_2\n";
print SH
"if [ -s $dir_path/$unpaired ] ; then cat $dir_path/*unPaired.fq >  \$tmp_dir/$unpaired.tmp ; fi\n";

if ($clean) {
  print SH
"$scripts/clean_pairs.pl -1 \$tmp_dir/$mate_1 -2 \$tmp_dir/$mate_2 > \$tmp_dir/$unpaired.tmp2\n";
  print SH
"if [ -e \$tmp_dir/$unpaired.tmp ] ; then cat \$tmp_dir/$unpaired.tmp \$tmp_dir/$unpaired.tmp2 > \$tmp_dir/$unpaired ; else mv \$tmp_dir/$unpaired.tmp2 \$tmp_dir/$unpaired ; fi\n";
  $mate_1 = "$lowest_dir" . "_p1.matched.fq";
  $mate_2 = "$lowest_dir" . "_p2.matched.fq";
}
print SH "cp \$tmp_dir/$mate_1 $one_up/$prefix$mate_1\n";
print SH "cp \$tmp_dir/$mate_2 $one_up/$prefix$mate_2\n";
print SH
"if [ -s \$tmp_dir/$unpaired  ] ; then  cp \$tmp_dir/$unpaired $one_up/$prefix$unpaired ; fi\n";

print SH "cd $current_dir\n";
print SH "rm -rf \$tmp_dir\n";
