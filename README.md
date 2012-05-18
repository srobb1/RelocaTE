RelocaTE
========

Find the locations of TEs using the TSD in unassembled short reads by comparing to a closely related reference genome assembly

run relocaTE.pl without any options to see information on how to run.

Requirements:
blat
bowtie
bioperl
samtools

What does relocaTE.pl actually do?
  1. it splits the supplied reference genome fasta into individual files, one file for each sequence.
  2. if not already done, creates a bowtie index for the complete reference fasta.
  3. if not already done, it converts the fq files to fa files.
  4. it splits the supplied TE fasta into individual files, one file for each sequence.
  5. if not already done, it runs blat: one job for every read file for every TE file. these jobs can be writted to shell scripts if -p 1, and an array job script will be generated if -a 1
  6. if not already done, it runs relocaTE_trim.pl: one job for every blat out file. shell scripts and array jobs will be created if -p 1 and -a 1
  7. it runs relocaTE_align.pl: one job for the one reference fasta. a shell script created if -p 1 and -a 1
  8. it runs relocaTE_insertionFinder.pl: one job for every TE for every sequence of the reference fasta. shell scripts and array jobs will be created if -p 1 and -a 1.
  9. it will concatenate the results of each reference sequence into one file: one job for every TE. shell scripts and array jobs will be created if -p 1 and -a 1.


--Sofia
