RelocaTE: A tool to identify the locations of transposable element insertion events that are present in DNA short read data but absent in the reference genome sequence.

RelocaTE is a collection of scripts in which short reads (paired or unpaired), a fasta containing the sequences of transposable elements and a reference genome sequence are the input and the output is a series of files containing the locations (relative to the reference genome) of TE insertions in the short reads. These insertions are insertions that are present only in the short reads and not present in the reference genome. If a tab-delimited file containing the coordinates of TEs in the reference is provided a list of the number of reads that support the presence of existing TE insertions is produced.

CharacTErizer: A companion tool compares the numbers of reads that flank the TE sequence and contain genomic sequence to the number of reads that span a predicted insertion site with no gaps. These spanners contain no TE sequence. The ratio of spanners to flankers is used to classify the insertion as homozygous, heterozygous, new (somatic) or other.

<a href="#cmd">RelocaTE Command Line Options</a>
- <a href="#t">-t Str:  TE FASTA File</a>
- <a href="#d">-d Str: directory of fq files</a>
- <a href="#g">-g Str: reference genome fasta file</a>
- <a href="#e">-e Str: Sample identifier</a>
- <a href="#o">-o Str: output directory name</a>
- <a href="#1">-1 Str: unique mate/pair 1 string</a>
- <a href="#2">-2 Str: unique mate/pair 2 string</a>
- <a href="#u">-u Str: unique unPaired string</a>
- <a href="#p">-p 0|1: split into many jobs</a>
- <a href="#a">-a 0|1: create PBS array job script</a>
- <a href="#w">-w Str: working directory name</a>
- <a href="#l">-l n:  min length cutoff for alignment to reference</a>
- <a href="#bm">-bm n: blat minScore for alignment to TE</a>
- <a href="#m">-m n<=0:  mismatches allowed in blat alignment to TE</a>
- <a href="#bt">-bt n: blat tileSize for alignmetn to TE</a>
- <a href="#f">-f n: length of the insertion site flanking seq to be returned </a>
- <a href="#x">-x Str: file name of file with locations of existing TE insertions in reference</a>


###<a name="cmd">RelocaTE Command Line Options</a>:
####<a name="t">-t TE Fasta File</t>

Required. No default value.

The file name of the fasta file containing the nucleotide sequence of one or many transposable elements.  The sequence should include the complete terminal inverted repeats (TIRs) [or LTR] but not include the target site in the sequence proper. The target site should be provided in the description portion of the fasta file in the following format, TSD=xyz. The TSD will be searched for in the both the forward and reverse strand. [[after testing the use of no TSD, write weather TSD= can be used]]

SAMPLE TE FASTA:
<pre>
>mping TSD=TTA
GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAA
ATGATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGC
GTCGTTTCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGAT
CTCTTGCGTCCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAA
TGATCCCAGCAACTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTC
CACTGTGGGGATTGTTTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGT
GCAATGACACTAGCCATTGTGACTGGCC
</pre>
The fasta must contain TSD=
<br>This can be a Perl regular expression.
<br>
<pre>
Example: these exact characters TTA: TSD=TTA
Example: any 4 characters: TSD=....
Example: A or T followed by GCC: TSD=(A|T)GCC
Example: CGA followed by any character then an A then CT or G: TSD=CGA.A(CT|G)
</pre>
####<a name="d">-d directory of fq files</a>

Required. No default value.

The name of the directory of paired and unpaired fastq files (paired _p1.fq & _p2.fq). Both the .fq and .fastq extensions are accepted as extensions of fastq files. If something different is used RelocaTE will not recognize those files as being fastq files. 


####<a name="g">-g reference genome fasta file</a>

Optional, Recommended. No default value.

The file name of the fasta file containing the genome sequence of the reference. If it is provided, the locations of the insertion events will be reported relative to the reference genome coordinates. 

If the genome sequence is not provided a series of files will be generated. One set will contain the intact reads that align to the TE. The second and third set of files will be made up of trimmed reads.  The second set will be only the trimmed portion of the reads in the first set that align to the TE. The third set will contain the trimmed portion of the reads that do not align to the TE, therefore the portion of the reads that should align to the genome sequence not containing a TE insertion.


####<a name="e">-e Sample identifier</a>

Optional, Recommended. The default value is not.given

A short string for sample name. This string will be used in the output files to create IDs for the insert (ex. A123)

####-o output directory name

Optional, Recommended. The default value is outdir_teSearch

A short string for the output directory name. This string will be used to create a directory to contain the output files and directories in the current working directory. The complete path is not required, only the desired name for the directory. 

####<a name="1>-1 unique mate/pair 1 string</a>

Optional, Recommended. The default value is _p1

A string to identify mate 1 paired files. Should contain the unique text and any text between the unique text and the fq extension. This string will be used in a regular expression to identify the files as a mate 1 file, so the string should not be found in the mate 2 file or the unpaired files

Example:
<pre>
If the files are named as such: file_1.fq
The string would be:            _1

File:                           file_1.noNumbers.fq
String:                         _1.noNumbers

File:                           file_1_1.fq (and mate = file_1_2.fq)
String:                         _1

Issue with last scenario:       _1 will recognize both mates.
Suggestion:                     rename files to file_1_p1.fq and file_1_p2.fq. 
                                Now the string _p1 can be used to uniquely 
                                identify all _p1 files and no _p2 files.
</pre>
####<a name="2">-2 unique mate/pair 2 string</>

Optional, Recommended. The default value is _p2

See -1 for a more in depth explanation.

Example:
<pre>
File:                            file_p2.fq
String:                          _p2
</pre>
####<a name="u">-u unique unPaired string</a>

Optional, Recommended. The default value is .unPaired.

See -1 for a more in depth explanation.

Example:
<pre>
File:                           file.unParied.fq
String:                         .unParied
</pre>

####<a name="p">-p n</a>

Optional. Default value is 1.

n is 0 or 1.

0: means only one large job will be ran.<br>
1: many shell scripts will be generated for the user to run<br>

Break down the single big job of relocaTE into as many smaller jobs as possible. If selected this option will cause the creation of shell scripts which can be manually ran or submitted to a queue. This enables the jobs to be run in parallel. The folders of shell scripts should be run as ordered. Step_1 needs to run and be complete before Step_2 jobs can be proper started.  If the genome fasta had already been split and indexed this job will be skipped. 

####<a name="a">-a n</a>

Optional. Default value is 1.

n is 0 or 1.

0: Nothing will be done<br>
1: if -p 1 then array jobs will be created<br>

If -p 1 then -a can be 1, otherwise it will be 0. If 1, in addition to the shell scripts generated from -p 1, qsub PBS array job scripts will be made for easier submittion of the shell scripts to the queue. See `man qsub option -t` for more information.

Submit each array job one at a time, waiting for the previous job to be completed before submitting the next.

See run_these_jobs.sh for the array jobs.

####<a name="w">-w working directory name</a>

Optional. Default value is the current working directory.

If a directory different form the cwd is given it needs to exist, will not create. Provide the full path. 

####<a name="l">-l n</a>

Optional. Default value is 10.

n is a value for the length cutoff. This is the minimum length that a needs to be after the removal of TE sequence. This trimmed read will be aligned to the genome. When selecting a custom value consider these points:
  - value can not be greater than the read length
  - How many bps are needed to limit false alignments to the reference?
  - How many bps are needed to recognize the TE? 
  - The answer to the above two questions should not total more than the read length.
 

####<a name="bm">-bm n</a>

Optional. Default value is 10.

n is used for the blat minScore value for the comparison of reads to the TE sequence.

Excerpt directly from Blat manual:

>
> -minScore=N sets minimum score.  This is the matches minus the 
>               mismatches minus some sort of gap penalty.

####<a name="m">-m n\<=0</a>

Optional. Default value is 0.

Any number less than or equal to 0.

Fraction of the bps that aligned to the TE that are allowed to not be an exact match. For example, if 10 bp align to the TE and the allowance is 0.1, 1 bp can be a mismatch.

 
####<a name="bt">-bt n</a>
 
 Optional. Default value is 7.

n is used for the blat tileSize value for the comparison of reads to the the TE sequence.

Excerpt directly from Blat manual:
> -tileSize=N sets the size of match that triggers an alignment.  
>               Usually between 8 and 12
>               Default is 11 for DNA and 5 for protein.
 
####<a name="f">-f n</a>
 
 Optional. Default value is 100.

n is the length of the sequence flanking the found insertion site to be returned in an output fasta file and in the output gff file. This sequence is taken from the reference genome.

####<a name="x">-x File</a>

Optional. No default value.

This is the file name of a tab-delimited file containing the coordinates of existing TE in the reference.  If this file is provided a new file will be generated containing a list of these existing insertions found in the reads. The number of reads supporting the start and the end of the insertion will be reported. 

The format is two columns, neither column have any white space. The first colum is the TE name. Then a tab separates the first column from the second column. And the second column contains the reference sequence name as described in the reference fasta, the starting bp, a colon, double periods, and finally the ending bp.


SAMPLE Existing TE File:
<pre>
mping   Chr12:839604..840033
mping   Chr12:1045463..1045892
</pre>  
  
	

###Quick Start Guide:

1.&nbsp;&nbsp;Get the sequence of your TE, including the TIRs.  Create a fasta file with your sequence, TE name and the TSD. The TSD= is required. With DNA transposons, by definition, during an insertion event the target site is duplicated. Therefore the target site will be used to identify an insertion event.  The reverse complement of each read containing portions of the ends of the provided TE will also be searched for the TSD to identify insertion events. A specific sequence of nucleotides can be used or a perl regular expression. 
For example: 
- if the desired TSD is TT followed by an A or G followed by a C or GT the regular expression would be:
<br>`TSD=TT[AG](C|GT).` 
- Also a very general pattern can be used: 
<br>`any 4 bp TSD=....`  

Example TE Fasta File: 
<pre>
>mping TSD=TTA
GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAA
AATGATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCA
GCGTCGTTTCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCG
GATCTCTTGCGTCCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTA
CAAATGATCCCAGCAACTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAAC
CCCTCCACTGTGGGGATTGTTTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGT
AGCCGTGCAATGACACTAGCCATTGTGACTGGCC
</pre>

2.&nbsp;&nbsp;Get the short reads together.  
&nbsp;&nbsp;&nbsp; - Dividing your reads files into many files of 1 million reads is a nice way to improve the speed of the analysis.  A script included in the package can be used to do this. It is called fastq_split.pl.

<pre>
perl fastq_split.pl 
 Please Provide: 
 -s int The number of sequences per file [1_000_000]
 -o dir The name of the directory to output the files
 followed by a list of files to be split

 Usage:
 ./fastq_split -o split_fq ~/somedir/someRandom.fq
 ./fastq_split -o split_fq ~/somedir/*fq
</pre>

&nbsp;&nbsp;&nbsp; - &nbsp;Check the name of your files. If they are paired they should have some indication of being pair 1 and pair 2.<br>
&nbsp;&nbsp;&nbsp; - &nbsp;Determine the pattern to indicate the 1st and 2nd pair<br>
	Example:
<pre>
Flowcell25_lane1_pair1.fastq and Flowcell25_lane1_pair2.fastq
_pair1 would match the first paired file
_pair2 would match the second paired file
</pre>

3.&nbsp;&nbsp;Get the reference genome fasta. This will be one file containing all the reference sequences of your organism.<br>
4.&nbsp;&nbsp;If you want to determine if you short reads contain the same TE insertions as the reference, create a file with the existing TE insertions found in the reference.

This file name is tab-delimited containing the coordinates of existing TE in the reference.  If this file is provided a new file will be generated containing a list of these existing insertions found in the reads. The number of reads supporting the start and the end of the insertion will be reported. 

The format is two columns, neither column have any white space. The first colum is the TE name. Then a tab separates the first column from the second column. And the second column contains the reference sequence name as described in the reference fasta, a :, the starting bp,  .. , and finally the ending bp.

SAMPLE:
<pre>
mping   Chr12:839604..840033
mping   Chr12:1045463..1045892
</pre>

5.&nbsp;&nbsp;Do you have a queue system available to you, or do you have multiple processors? If so, you can select the –p 1 option. This will create a series of shell scripts that you can submit to the queue or run on multiple processors.<br> 
&nbsp;&nbsp;&nbsp; - These should be run in order. For example all of step_1 should be run and completed before the shell scripts of step_2 are run, and so on.

6.&nbsp;&nbsp;Do you have a queue system available and is it PBS? If so, you can select the –a 1 option to generate array job scripts. These jobs can be submitted to the queue instead of submitting hundreds of individual shell scripts. The array job will submit the hundreds of scripts for you.  
&nbsp;&nbsp;&nbsp; - Make sure to submit the array jobs in order and wait for the job to complete before submitting the next one.

7.&nbsp;&nbsp;You are now ready to run relocaTE.pl with your data. If you run the program without any command line options, it will print out a list of the options and short descriptions.

	
###CharacTErizer:
<pre>
usage:
./characterizer.pl [-s relocaTE table output file][-b bam file or dir of bam files][-g reference genome fasta file][-h] 

options:
-s file         relocaTE output file: *.te_insertion.all.txt [no default]
-b dir/file     bam file of the orginal reads aligned to reference (before TE was trimmed) or directory of bam files [no default]
-g file         reference genome fasta file [no default]
-h              this help message
</pre>

###Requirements:
- blat
- bowtie
- bioperl
- samtools
 
###What does relocaTE.pl actually do?
  1. it splits the supplied reference genome fasta into individual files, one file for each sequence.
  2. if not already done, creates a bowtie index for the complete reference fasta.
  3. if not already done, it converts the fq files to fa files.
  4. it splits the supplied TE fasta into individual files, one file for each sequence.
  5. if not already done, it runs blat: one job for every read file for every TE file. these jobs can be writted to shell scripts if -p 1, and an array job script will be generated if -a 1
  6. if not already done, it runs relocaTE_trim.pl: one job for every blat out file. shell scripts and array jobs will be created if -p 1 and -a 1
  7. it runs relocaTE_align.pl: one job for the one reference fasta. a shell script created if -p 1 and -a 1
  8. it runs relocaTE_insertionFinder.pl: one job for every TE for every sequence of the reference fasta. shell scripts and array jobs will be created if -p 1 and -a 1.
  9. it will concatenate the results of each reference sequence into one file: one job for every TE. shell scripts and array jobs will be created if -p 1 and -a 1.

For more information see documentation: http://docs........