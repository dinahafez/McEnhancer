import string;
import sys;
import numpy as np;
import subprocess;

inputDirectory = sys.argv[1]; #"/data/ohler/Dina/Drosophila/Research/results/BaseLine_assign_DHS_to_closest_gene/ClustersDHS/" ;# sys.argv[1];
seqFile = sys.argv[2]; #"1_pos_single" ;#sys.argv[2];
resultDirectory = sys.argv[3]; #"/data/ohler/Dina/Drosophila/Research/results/BaseLine_assign_DHS_to_closest_gene/ClustersDHS/";#sys.argv[3];
kmerLengthString = sys.argv[4];
kmerLength = int(kmerLengthString);

def countKmers(currSeq,kmerDict,numKmers):
    # Count kmers here
    kmerCountArray = np.zeros(numKmers+1);
    totalKmers = float(len(currSeq) - kmerLength+1);
    for i in range(0,len(currSeq)-kmerLength+1):
        currKmer = currSeq[i:i+kmerLength];
        kmerCountArray[kmerDict[currKmer]] = kmerCountArray[kmerDict[currKmer]] + float(1);#/totalKmers;

    return kmerCountArray;

def main():
    # read in kmerDict
    kmerDict = dict();
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
    #singleIndex = seqFile.index(".fa");
    outFile = seqFile;
    kmerCountFileName = resultDirectory + str(kmerLength) +"_mer/" + outFile + "_count.tab";
    print kmerCountFileName;
    kmerCountFile = open(kmerCountFileName,'w');
    kmerCountFile.write("sid");
    for i in range(0,maxNumber+1):
        kmerCountFile.write("\t"); kmerCountFile.write("kmer"+str(i+maxNumber+1));
    kmerCountFile.write("\n");

    # Read sequence file and get kmers
    
    sequenceFileName = inputDirectory + seqFile + ".fa";

    sequenceFile = open(sequenceFileName,'r');
    seqID = -1;
    for line in sequenceFile:
        if (line[0] == ">"):
            seqID = line[1:].rstrip();
            geneName = seqID.split(';')[0];
        else:
            seq = line.rstrip();
            currKmerArray = countKmers(seq, kmerDict, maxNumber);
            kmerCountFile.write(geneName+"_"+str(len(seq)));
            for kmerCountIndex in range(len(currKmerArray)):
            #for kmerCount in currKmerArray:
                if (currKmerArray[kmerCountIndex] != 0):
                    kmerCountFile.write("\t"); 
                    ind = kmerCountIndex + 1 + maxNumber + 1;
                    kmerCountFile.write(str(ind));
                    kmerCountFile.write(":"); 
                    kmerCountFile.write(str(currKmerArray[kmerCountIndex]));
                    
            kmerCountFile.write("\n");
    kmerCountFile.close();
    sequenceFile.close();
    
    
    args = resultDirectory +" "+ seqFile + " "+kmerLengthString;
    pipe = subprocess.call(["perl", "/data/ohler/Dina/Drosophila/Research/code/kmerCount/src/normalize_byGene.pl", resultDirectory , seqFile,kmerLengthString])
        
    paralerFileName = resultDirectory + str(kmerLength) +"_mer/" + outFile + ".tabParallel";
    paralerFile = open(paralerFileName,'w');
    paralerFile.close();
    
main();
