#!/bin/bash
#
#$ -S /bin/tcsh -cwd
#$ -o error2.out -j y
#$ -l h_vmem=10G
#$ -l mem_free=10G
#$ -e error2.log

#Pre: all files for initialization for each cluster have been generated, labeled, unlabeled, and background files, and saved in result directory

codeVersion=$1"/";
MC_order=$2;	 #Mc order
positive=$3;
clusterNo=$4;	 #cluster no to process
itr=$5;
kmer="5";

discount=$6;
discountn=$7;
dif=$8;

echo $positive
echo $clusterNo
echo $codeVersion

if [ $clusterNo != $positive ]
then
echo "MC order is "$MC_order;				
#create directories 
MC_order_path=$codeVersion; 
parentDir=$codeVersion 
fullPath=$MC_order_path;
if [ ! -d "$fullPath" ]
then
	mkdir -p $fullPath;
fi

fullPath=$fullPath"MC_"$positive"_output/";

if [ ! -d "$fullPath" ]
then
	mkdir  -p $fullPath;
fi


#run Makov Chain 
dhsLabeledFileName=$parentDir$positive"_initialize_single.fa";   
dhsUnLabeledFileName=$parentDir$positive"_unlabeled_single.fa";
dhsBackgroundFileName=$parentDir$clusterNo"_initialize_single.fa";     
dhsLabeledPositiveFileName=$fullPath$positive"vs"$clusterNo"_pos.fa";
echo $dhsLabeledPositiveFileName;
if [ ! -s $dhsLabeledPositiveFileName ]
then
	echo "perl InterpolatedMC.pl $dhsLabeledFileName $dhsUnLabeledFileName $dhsBackgroundFileName $dhsLabeledPositiveFileName $MC_order $itr;"
	perl InterpolatedMC.pl $dhsLabeledFileName $dhsUnLabeledFileName $dhsBackgroundFileName $dhsLabeledPositiveFileName $MC_order $itr $discount $discountn $dif;
	 
fi
fi

