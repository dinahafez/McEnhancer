#!/bin/bash
#
#$ -S /bin/tcsh -cwd
#$ -o error2.out -j y
#$ -l h_vmem=10G
#$ -l mem_free=10G
#$ -e error2.log

dir=$1;
itr=$2;
pos=$3;
order=$4;

discount=$5;
discountn=$6;
dif=$7;
dir_stage2="ModelStage2"

for i in `seq 1 18`;
do  
	./mc_per_cluster.sh $dir $order $pos $i $itr $discount $discountn $dif;
	
	#./mc_per_cluster_stage2.sh $dir_stage2 $order $pos $i $itr $discount $discountn $dif;
done
	
