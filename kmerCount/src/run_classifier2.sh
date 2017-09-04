#!/bin/tcsh
#
#$ -S /bin/tcsh -cwd
#$ -o RError2.out -j y
#$ -l h_vmem=20G
#$ -l mem_free=20G
#$ -e error.log



python cvClassifierROC2.py $1 $2 $3


