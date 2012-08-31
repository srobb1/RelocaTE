#!/usr/bin/perl -w
##modified 11-3-2011
use strict;
use File::Spec;

#provide a directory of bam files to merge

if ( !defined @ARGV ) {
 die 
"for a single directory: 
            merge_bam_and_create_merged_sam_shell.pl dir_of_bam_files prefix [tempDir:/scratch]
for multiple directories:
            for i in `ls` ; do merge_bam_and_create_merged_sam_shell.pl \$i prefix [tempDir:/scratch]; done
or if you have numbered ref seqs, ex Chr1, Chr2 ... Chr12:
            for i in `seq 1 12` ; do merge_bam_and_create_merged_sam_shell.pl Chr\$i prefix [tempDir:/scratch]; done

Then run each shell script individually or in parallel\n\n";
}
my $dir     = shift;
my $prefix  = shift;
my $tempDir = shift;
$prefix  = !defined $prefix ? ''       : $prefix . '.';
$tempDir = defined $tempDir ? $tempDir : '/scratch';
my $dir_path    = File::Spec->rel2abs($dir);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

my @dirs       = split '/', $dir_path;
my $lowest_dir = pop @dirs;
my $one_up     = join '/', @dirs;
pop @dirs;
my $two_up = join '/', @dirs;

open SH, ">$current_dir/$lowest_dir.mergeBam.sh";
print SH "#!/bin/bash\n\n";

print SH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
print SH "cd \$tmp_dir\n";

#merge bam files
print SH "samtools merge \$tmp_dir/$lowest_dir.merged.bam $dir_path/*bam\n";

#sort and index merged bam
print SH
"samtools sort  \$tmp_dir/$lowest_dir.merged.bam  \$tmp_dir/$lowest_dir.merged.sorted\n";
print SH "samtools index  \$tmp_dir/$lowest_dir.merged.sorted.bam\n";

=cut
#convert merged and sorted bam to a sam
print SH
"samtools view -h \$tmp_dir/$lowest_dir.merged.sorted.bam -o \$tmp_dir/$lowest_dir.merged.sorted.sam\n";
print SH
"if [ ! -d \"$two_up/sam_split_by_chromosome\" ]; then mkdir $two_up/sam_split_by_chromosome ; fi\n";
print SH
"for i in `ls *.sam*` ; do cp \$tmp_dir/\$i $two_up/sam_split_by_chromosome/$prefix\$i ; done\n";
=cut

print SH
  "for i in `ls *merged.sorted.bam*` ; do cp \$tmp_dir/\$i $one_up/$prefix\$i ; done\n";
print SH "cd $current_dir\n";
print SH "rm -rf \$tmp_dir\n";
