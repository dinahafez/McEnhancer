import string;
import sys;
import numpy as np;

inputDirectory = "/data/ohler/Dina/Drosophila/Research/results/BaseLine_assign_DHS_to_closest_gene/ClustersDHS/" ;# sys.argv[1];
seqFile = "111_pos_single" ;#sys.argv[2];
resultDirectory = "/data/ohler/Dina/Drosophila/Research/results/BaseLine_assign_DHS_to_closest_gene/ClustersDHS/";#sys.argv[3];
kmerLength = int(5);#   sys.argv[4]);

def countKmers(currSeq,kmerDict,numKmers):
	# Count kmers here
	kmerCountArray = np.zeros(numKmers+1);
	totalKmers = float(len(currSeq) - kmerLength+1);
	for i in range(0,len(currSeq)-kmerLength+1):
		currKmer = currSeq[i:i+kmerLength];
		kmerCountArray[kmerDict[currKmer]] = kmerCountArray[kmerDict[currKmer]] + float(1)/totalKmers;
		

	return kmerCountArray;

def main():
	# read in kmerDict
	kmerDict = dict();
	#kmerDictFileName = "../data/kmerList/all" + str(kmerLength) + "mer.list";
	kmerDictFileName = "/data/ohler/Dina/Drosophila/Research/code/kmerCount/data/kmerList/all" + str(kmerLength) + "mer.list";
	

	kmerDictFile = open(kmerDictFileName,'r');
	kmerDictFile.readline();
	maxNumber = -1;
	for line in kmerDictFile:
		elems = string.split(line.rstrip());
		kmerDict[elems[0]] = int(elems[1]);
		if (int(elems[1]) > maxNumber):
			maxNumber = int(elems[1]);
	print maxNumber;
	kmerDictFile.close();


	# Open and write file header
	#kmerCountFileName = "../results/" + seqFile + "_" + str(kmerLength) + "mer.tab";
	kmerCountFileName = resultDirectory + str(kmerLength) +"_mer/" + seqFile + ".tab";

	kmerCountFile = open(kmerCountFileName,'w');
	kmerCountFile.write("sid");
	for i in range(0,maxNumber+1):
		kmerCountFile.write("\t"); kmerCountFile.write("kmer"+str(i));
	kmerCountFile.write("\n");

	# Read sequence file and get kmers
	#sequenceFileName = "../data/sequence/" + seqFile + ".fa";
	sequenceFileName = inputDirectory + seqFile + ".fa";

	sequenceFile = open(sequenceFileName,'r');
	seqID = -1;
	for line in sequenceFile:
		if (line[0] == ">"):
			seqID = line[1:].rstrip();
		else:
			currKmerArray = countKmers(line.rstrip(), kmerDict, maxNumber);
			kmerCountFile.write(seqID);
			for kmerCount in currKmerArray:
				kmerCountFile.write("\t"); kmerCountFile.write(str(kmerCount));
			kmerCountFile.write("\n");
	kmerCountFile.close();
	sequenceFile.close();

main();
