#!/bin/bash
#
#$ -S /bin/tcsh -cwd
#$ -o error2.out -j y
#$ -l h_vmem=10G
#$ -l mem_free=10G
#$ -e error2.log

#Pre: all files for initialization for each cluster have been generated, labeled, unlabeled, and background files, and saved in result directory
resultDrive="../result/ModelStage1/"; # "/data/ohler/Dina/Drosophila/Research/results_idr_0.2/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/ModelSelection/";

clusterFile2=$1	 #cluster no to process

kmer="5";         #no of kmers
background="ubiq";   #could be ubiq or other


DIR=$(dirname "${clusterFile2}")"/"

clusterFile=$(basename "$clusterFile2")

IFS='_' read -a array <<< "${clusterFile}";
clusterNo=${array[0]};
FileName="${clusterFile%%.*}";



if [ $clusterNo == "1" ] || [ $clusterNo == "2" ] || [ $clusterNo == "3" ] || [ $clusterNo == "4" ] || [ $clusterNo == "5" ] || [ $clusterNo == "6" ] || [ $clusterNo == "7" ] || [ $clusterNo == "8" ] || [ $clusterNo == "9" ] || [ $clusterNo == "10" ] ; then 
	echo "non"
else   

parentDir=$resultDrive;
fullPath=$parentDir$DIR;


ubiqFileCheck=$fullPath"all_non_nervous_10_merged_single.fa";


#get kmer for positive, negative and ubiquious files, to prepare for classification using logreg
posFile=$FileName;   
negFile="all_pos_not_"$clusterNo; 
ubiqFile="all_non_nervous_10_merged_single";

echo "Generating kmers ";

if [ $background == "other" ]
then
	if [ ! -s $fullPath$kmer"_mer/"$posFile".tab" ]
	then 
		python kmerCount/src/getKmer_normalize_byGene_Libsvm.py $fullPath $posFile $fullPath $kmer
	
	fi
	if [ ! -s $fullPath$kmer"_mer/"$negFile".tab" ]
	then 
		perl /data/ohler/Dina/Drosophila/Research/code/pipeline/getNegativeFiles_MC.pl $fullPath $clusterNo;
		python kmerCount/src/getKmer_normalize_byGene_Libsvm.py $fullPath $posFile $fullPath $kmer

		
	fi
else
	if [ ! -s $fullPath$kmer"_mer/"$posFile".tab" ]
	then 
		python kmerCount/src/getKmer_normalize_byGene_Libsvm.py $fullPath $posFile $fullPath $kmer
	fi

fi

echo "Generating kmers finished";
#classify between positive vs all_not_cluster, as well as between positive vs ubiquitous

posClassifier=$fullPath$kmer"_mer/"$posFile".tabParallel";
negClassifier=$fullPath$kmer"_mer/"$negFile".tabParallel";
ubiqClassifier=$fullPath$kmer"_mer/"$ubiqFile".tabParallel"; 

echo "classifing";

if [ $background == "other" ]
then 
	while [ ! -e $negClassifier -o ! -e $posClassifier ];
	do	
		echo "classification cannot be started coz the kmers are still not ready--will sleep";
		sleep 10   #sleep for 30 seconds
	done
	#./run_classifier.sh $fullPath $posFile $negFile $kmer
	echo "running classifier";

	ClassifyFile=$fullPath$kmer"_mer/"$posFile"_"$negFile"_Classifier.ROC";
	echo $ClassifyFile;
	if [ ! -s $ClassifyFile ]
	then	
		java -Xmx7g -jar /data/ohler/Dina/Drosophila/Research/code/ClassifyKmerCounts/target/ClassifyKmerCounts-1.2-SNAPSHOT-jar-with-dependencies.jar $fullPath $posFile $negFile $kmer
	fi
else
	echo "running classifier";

ClassifyFile=$fullPath$kmer"_mer/"$posFile"_"$ubiqFile"_Classifier.ROC";
echo $ClassifyFile;
	if [ ! -s $ClassifyFile ]
	then	
		java -Xmx7g -jar /data/ohler/Dina/Drosophila/Research/code/ClassifyKmerCounts/target/ClassifyKmerCounts-1.2-SNAPSHOT-jar-with-dependencies.jar $fullPath $posFile $ubiqFile $kmer
	fi
fi

fi
