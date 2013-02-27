This is a collection of sample data and command line arguments to run relocaTE and characTErizer.

The sample scritps will only work on a macOSX or linux operating system.
The sample scritps should be executed using a terminal.

These scripts expect that the following programs are installed and included in your path (See http://srobb1.github.com/RelocaTE/ ffor information on where to find these tools):
	1. Samtools
	2. Bowtie (version 1)
	3. BWA
	4. BioPerl
   	5. Blat
	6. Blast: formatdb and fastacmd

STEPS A-C
A. RelocaTE
B. BAM files
C. CharacTErizer

A. To find TE insertions in the included fastq reads run RelocaTE:

  Start by selecting the correct run_relocaTE shell script.
    - if you have access to a PBS queue, use:
  	run_relocaTE_qsub.sh
    - if you do not have access to a PBS queue but do can run jobs in parallel, use:
	run_relocaTE_shell.sh
    - if you do not have access to a PBS queue and can not or do not want to run jobs in parallel or have no idea what any of this means, use:
	run_relocaTE.sh

  You can run this script by changing into the dirctory that contains this script and typing
 	sh run_relocaTE_qsub.sh
	or
        sh run_relocaTE_shell.sh
        or
        sh run_relocaTE.sh

  These scritps will run relocaTE directly or through a series of shell scripts. 
	If you are running
        - run_relocaTE_qsub.sh:
		1) follow the instructions that are printed to the screen (run run_these_jobs.sh)
      		2) once complete, view the results in 02052012_sample/mping/results
        or
        - run_relocaTE_shell.sh
		1) follow the instructions that are printed to the screen (run run_these_jobs.sh)
      		2) once complete, view the results in 02052012_sample/mping/results
        or
        - run_relocaTE.sh
      		1) once complete, view the results in 02052012_sample/mping/results


B. Find Spanners to help classify the insertions (homozygous, heterozygous, etc) by generating a BAM file of the reads not trimmed of TE to the reference.
  A BAM file is included in the sample data set, but one can be genreted by running the included script:
	create_bam.sh by:
	- changing directory into the directory of this script
 	- and typing "sh create_bam.sh"


C. To classify the insertions (homozygous, heterozygous, etc) run_characTErizer.sh 
	- change directory into the directory of this script
	- type "sh run_characTErizer.sh"
        - once complete, view the resulting files in the directory you ran characTErizer.pl
		1) sample.inserts_characTErized.gff: GFF file of the classified insertions including excisions 
		2) sample.inserts_characTErized.txt: Text file of the classified insertions including excisions
		3) excisions_with_footprint.vcfinfo: additional information on the insertions that have been classified as exicision events  
 
