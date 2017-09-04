#!/bin/bash
#
#$ -S /bin/tcsh -cwd
#$ -o error2.out -j y
#$ -l h_vmem=10G
#$ -l mem_free=10G
#$ -e error2.log

#Pre: all files for initialization for each cluster have been generated, labeled, unlabeled, and background files, and saved in result directory
dataDrive="/data/ohler/Dina/Drosophila/Research/data/ClustersDHS/";
resultDrive="/data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection/ModelStage1/MC_multi_output/";

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
MC_order_path=$codeVersion; #"srilm_"$MC_order"order/";  
parentDir=$resultDrive;
fullPath=$parentDir$MC_order_path;
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
dhsLabeledFileName=$parentDir$positive"_pos_10_single.fa";   #initialize_redfly_uniq_single.fa";
dhsUnLabeledFileName=$parentDir$positive"_unlabeled_single.fa";
#Markov chain background could be other or ubiquitous
	dhsBackgroundFileName=$parentDir$clusterNo"_pos_10_single.fa";     #$parentDir$clusterNo"_neg.fa";
 

dhsLabeledPositiveFileName=$fullPath$positive"vs"$clusterNo"_pos.fa";
echo $dhsLabeledPositiveFileName;
if [ ! -s $dhsLabeledPositiveFileName ]
then
	echo "perl InterpolatedMC.pl $dhsLabeledFileName $dhsUnLabeledFileName $dhsBackgroundFileName $dhsLabeledPositiveFileName $MC_order $itr;"
	perl InterpolatedMC.pl $dhsLabeledFileName $dhsUnLabeledFileName $dhsBackgroundFileName $dhsLabeledPositiveFileName $MC_order $itr $discount $discountn $dif;
	 
fi
fi

