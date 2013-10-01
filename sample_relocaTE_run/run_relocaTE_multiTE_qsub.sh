## run this command to run relocaTE
## use this command only if you have access to a PBS queue
## to check run this command qstat. If it returns info about jobs in the queue you can run this command
## This is to show the output when using a FASTA with more than one TE.
## note: the added TE will not find any non-ref insertions

../scripts/relocaTE.pl -t data/TEs.fa -g data/MSUr7.sample.fa -d data/fq -e sample -o 02052012_sample -1 _p1 -2 _p2 -r 1 -p 1 -a 1  

