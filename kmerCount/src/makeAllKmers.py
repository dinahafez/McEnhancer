import string;
import sys;

kmerLength = int(sys.argv[1]);

def makeString(currPos):
	nucleotideList = ['A','C','G','T'];
	finalReturnList = list();
	if (currPos < kmerLength):
		for nucleotide in nucleotideList:
			#print nucleotide,"\t",currPos;
			currReturnList = makeString(currPos + 1);
			#print currReturnList
			for i in range(0,len(currReturnList)):
				currReturnList[i] = nucleotide + currReturnList[i];
			finalReturnList = finalReturnList + currReturnList;
		return finalReturnList;
	else:
		#print "HERE";
		return nucleotideList;
	
def revcomp(currString):
	revNucleotideDict = dict();
	revNucleotideDict["A"] = "T";
	revNucleotideDict["C"] = "G";
	revNucleotideDict["G"] = "C";
	revNucleotideDict["T"] = "A";
	
	currRevString = "";
	for nucleotide in currString:
		currRevString = currRevString + revNucleotideDict[nucleotide];
	return currRevString[::-1];

def main():

	kmerFileName = "../data/all" + str(kmerLength) + "mer.list";	
	kmerFile = open(kmerFileName,'w');
	kmerFile.write("kmer\tindex\n");

	kmerList = makeString(1);	
	print len(kmerList);
	
	kmerDict = dict();
	currIndex = 0;	

	for kmer in kmerList:
		revKmer = revcomp(kmer);
		if (kmer not in kmerDict and revKmer not in kmerDict):
			kmerDict[kmer] = currIndex;
			kmerDict[revKmer] = currIndex;
			currIndex = currIndex + 1;

	print len(kmerDict);

	for kmer in kmerDict:
		kmerFile.write(kmer); kmerFile.write("\t");
		kmerFile.write(str(kmerDict[kmer])); kmerFile.write("\n");
	
	kmerFile.close();

	indexKmerDict = dict();
	for kmer in kmerDict:
		currIndex = kmerDict[kmer];
		if currIndex not in indexKmerDict:
			indexKmerDict[currIndex] = list();
		indexKmerDict[currIndex].append(kmer);
	print len(indexKmerDict);


main();


