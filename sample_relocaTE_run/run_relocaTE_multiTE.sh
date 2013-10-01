## run this command to run all steps sequentially
## only use this if you do not have access to cluster or qsub
## This is to show the output when using a FASTA with more than one TE.
## note: the added TE will not find any non-ref insertions 

../scripts/relocaTE.pl -t data/TEs.fa -g data/MSUr7.sample.fa -d data/fq -e sample -o 02052012_sample -1 _p1 -2 _p2 -r 1

