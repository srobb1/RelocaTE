## run this command to run relocaTE
## run this only if you do not have access to a PBS queue
## to check, run this command, 'qstat'. if you see info about jobs in the queue add '-a 1' to this command
## or use run_relocaTE_qsub.sh instead.
## if you do not see info about jobs in the queue, use this command:

../scripts/relocaTE.pl -t data/mping.fa -g data/MSUr7.sample.fa -d data/fq -e sample -o 02052012_sample -1 _p1 -2 _p2 -r 1 -p 1 

