#!/bin/tcsh
#
#$ -S /bin/tcsh -cwd
#$ -o error.out -j y
#$ -l h_vmem=2G
#$ -l mem_free=2G
#$ -e error.log

#java -Xmx7g -jar /data/ohler/Dina/Drosophila/Research/code/EMwithMC/target/EMwithMC-1.0-SNAPSHOT-jar-with-dependencies.jar $1 $2 $3 $4 $5

python /data/ohler/Dina/Drosophila/Research/code/kmerCount/src/getKmer_Libsvm.py $1 $2 $3 $4
