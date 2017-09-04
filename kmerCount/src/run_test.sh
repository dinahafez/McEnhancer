#!/bin/tcsh
#
#$ -S /bin/tcsh -cwd
#$ -o R.out -j y
#$ -l h_vmem=20G
#$ -l mem_free=20G
#$ -e error.log



python pyhton_test.py 
