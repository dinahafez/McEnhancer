import os;
import string;
import numpy as np;
import sys;
from random import shuffle;

args = sys.argv;

cvFold = 4;
iterations = 10;

#tfType = args[1];
#hssRegion = args[2];
	
kmers = '5';#args[1];
cluster = '11'; #args[2];
condition = 'after_MC'; #args[3];

#sequenceList1 = args[2];
#sequenceList2 = args[3];

sequenceList1 = cluster+"_pos_single";
sequenceList2 =  cluster+"_pos_single"; #"all_pos_not_"+cluster+"_single";
#sequenceList2 = "ubiquitous"; 
#condition = "before_MC";

#dataDrive = "/nfs/igsp/ohler_fc/ohlerlab/anirudh/mm9DHSClassifier/";
resultDrive = "/data/ohler/Dina/Drosophila/Research/results/MarkovChain_5order_classification_python/";

sequenceSetDict = dict();
#sequenceSetDict[0] = dataDrive + "input/" + tfType + "/" + hssRegion + "/" + sequenceList1 + ".tab";
#sequenceSetDict[1] = dataDrive + "input/" + tfType + "/" + hssRegion + "/" +sequenceList2 + ".tab";
sequenceSetDict[0] = resultDrive + condition + "/" + kmers + "_mer/" +  sequenceList1 + ".tab";
sequenceSetDict[1] = resultDrive + condition + "/" + kmers + "_mer/" + sequenceList2 + ".tab";

classifierStatsFileName = resultDrive + condition + "/" + kmers + "_mer/"+ sequenceList1 + "_" + sequenceList2 + "_Classifier.ROC";
classifierStatsFile = open(classifierStatsFileName,'w');
classifierStatsFile.write("List1\tList2\tROC\n");

regCoeffFileNamePrefix = resultDrive + condition + "/" + kmers + "_mer/"+ sequenceList1 + "_" + sequenceList2 + "_Coeff";

comparisonDict = dict();
for i in range(0,len(sequenceSetDict)):
	comparisonDict[i] = list();
comparisonDict[0].append(1);

compareNameList = list();
compareNameList.append(sequenceList1);compareNameList.append(sequenceList2);

for i in range(0,len(sequenceSetDict)):
	compareFile1Name = sequenceSetDict[i];
	compareFile1 = open(compareFile1Name,'r');
	file1HeaderLine = compareFile1.readline();  #sid kmer0 kmer1 ... kmern
	compareList1 = list();
	for line in compareFile1:   #id count1 count2 ... countn
		elems = string.split(line.rstrip());
		compareList1.append(elems[1:]);
		
	compareFile1.close();
	rmFileCmdStr = "rm " + sequenceSetDict[i];
	#os.system(rmFileCmdStr);
	numFeatures = len(compareList1[0]);

	cvList1Sizes = np.zeros(cvFold);
	compareList1Size = len(compareList1);
	for k in range(0,cvFold):
		cvList1Sizes[k] = int(compareList1Size/cvFold);
	listRemainder = compareList1Size % cvFold;
	for k in range(0,listRemainder):
		cvList1Sizes[k] = cvList1Sizes[k] + 1;

	currComparisonList = comparisonDict[i];
	for j in currComparisonList:

		currRegCoeffList =  list();
		for regCoeffIndex in range(0,numFeatures+1):
			currRegCoeffList.append(list());

		compareFile2Name = sequenceSetDict[j];
		compareFile2 = open(compareFile2Name,'r');		
		file2HeaderLine = compareFile2.readline();
		compareList2 = list();
		for line in compareFile2:
			elems = string.split(line.rstrip());
			compareList2.append(elems[1:]);
		compareFile2.close();
		rmFileCmdStr = "rm " + sequenceSetDict[j];
		#os.system(rmFileCmdStr);
	
		cvList2Sizes = np.zeros(cvFold);
		compareList2Size = len(compareList2);
		for k in range(0,cvFold):
			cvList2Sizes[k] = int(compareList2Size/cvFold);
		listRemainder = compareList2Size % cvFold;
		for k in range(0,listRemainder):
			cvList2Sizes[k] = cvList2Sizes[k] + 1;
		
		compiledAuROC = 0.0;
		
		for iter in range(0,iterations):
			shuffle(compareList1);
			shuffle(compareList2);
			for k in range(0,cvFold):
				numTrainGenes = 0;
				numLambdaGenes = 0
				numTestGenes = 0;

				currTrainValues = list();
				currTrainClasses = list();
				currLambdaValues = list();
				currLambdaClasses = list();
				currTestValues = list();
				currTestClasses = list();

				arrayPointer = 0;

				totalPosTrainGeneList = list();
				cvPosTrainGeneList = list();
				cvPosLambdaGeneList = list();
				cvPosTestGeneList = list();

				for m in range(0,cvFold):
					if (m!=k):
						for n in range(0,int(cvList1Sizes[m])):
							totalPosTrainGeneList.append(arrayPointer+n);
					else:
						for n in range(0,int(cvList1Sizes[m])):
							cvPosTestGeneList.append(arrayPointer + n);
					arrayPointer = arrayPointer + int(cvList1Sizes[m]);

				shuffle(totalPosTrainGeneList);
				totalTrainNum = len(totalPosTrainGeneList);
				cvTrainNum = (totalTrainNum * 5)/6;
				cvLambdaNum = totalTrainNum - cvTrainNum;

				numTrainGenes = numTrainGenes + cvTrainNum;
				numLambdaGenes = numLambdaGenes + cvLambdaNum;
				numTestGenes = numTestGenes + len(cvPosTestGeneList);
				

				for n in range(0,totalTrainNum):
					if (n<cvTrainNum):
						cvPosTrainGeneList.append(totalPosTrainGeneList[n]);
					else:
						cvPosLambdaGeneList.append(totalPosTrainGeneList[n]);

				arrayPointer = 0;
				totalNegTrainGeneList = list();
				cvNegTrainGeneList = list();
				cvNegLambdaGeneList = list();
				cvNegTestGeneList = list();

				for m in range(0,cvFold):
					if (m!=k):
						for n in range(0,int(cvList2Sizes[m])):
							totalNegTrainGeneList.append(arrayPointer+n);
					else:
						for n in range(0,int(cvList2Sizes[m])):
							cvNegTestGeneList.append(arrayPointer + n);
					arrayPointer = arrayPointer + int(cvList2Sizes[m]);

				shuffle(totalNegTrainGeneList);
				totalTrainNum = len(totalNegTrainGeneList);
				cvTrainNum = (totalTrainNum * 5)/6;
				cvLambdaNum = totalTrainNum - cvTrainNum;

				numTrainGenes = numTrainGenes + cvTrainNum;
				numLambdaGenes = numLambdaGenes + cvLambdaNum;
				numTestGenes = numTestGenes + len(cvNegTestGeneList);
				
				for n in range(0,totalTrainNum):
					if (n<cvTrainNum):
						cvNegTrainGeneList.append(totalNegTrainGeneList[n]);
					else:
						cvNegLambdaGeneList.append(totalNegTrainGeneList[n]);

				for f in range(0,numFeatures):
					for n in range(0,len(cvPosTrainGeneList)):
						currFeatureList = compareList1[cvPosTrainGeneList[n]];
						currTrainValues.append(currFeatureList[f]);				
					for n in range(0,len(cvNegTrainGeneList)):
						currFeatureList = compareList2[cvNegTrainGeneList[n]];
						currTrainValues.append(currFeatureList[f]);				
					for n in range(0,len(cvPosLambdaGeneList)):
						currFeatureList = compareList1[cvPosLambdaGeneList[n]];
						currLambdaValues.append(currFeatureList[f]);				
					for n in range(0,len(cvNegLambdaGeneList)):
						currFeatureList = compareList2[cvNegLambdaGeneList[n]];
						currLambdaValues.append(currFeatureList[f]);				
					for n in range(0,len(cvPosTestGeneList)):
						currFeatureList = compareList1[cvPosTestGeneList[n]];
						currTestValues.append(currFeatureList[f]);				
					for n in range(0,len(cvNegTestGeneList)):
						currFeatureList = compareList2[cvNegTestGeneList[n]];
						currTestValues.append(currFeatureList[f]);				
				for n in range(0,len(cvPosTrainGeneList)):
					currTrainClasses.append("1");				
				for n in range(0,len(cvPosLambdaGeneList)):
					currLambdaClasses.append("1");				
				for n in range(0,len(cvPosTestGeneList)):
					currTestClasses.append("1");				

				for n in range(0,len(cvNegTrainGeneList)):
					currTrainClasses.append("-1");				
				for n in range(0,len(cvNegLambdaGeneList)):
					currLambdaClasses.append("-1");				
				for n in range(0,len(cvNegTestGeneList)):
					currTestClasses.append("-1");				

				arrayPointer = 0;
				
				trainFeatureOutfileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + "currCVFeatures.train";
				trainFeatureOutfile = open(trainFeatureOutfileName,'w');
				trainFeatureOutfile.write("%%MatrixMarket matrix array real general\n");
				trainFeatureOutfile.write(str(numTrainGenes));trainFeatureOutfile.write("\t");
				trainFeatureOutfile.write(str(numFeatures));trainFeatureOutfile.write("\n");
				for entry in currTrainValues:
					trainFeatureOutfile.write(entry);trainFeatureOutfile.write("\n");
				trainFeatureOutfile.close();
				
				trainClassOutfileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + "currCVClass.train";
				trainClassOutfile = open(trainClassOutfileName,'w');
				trainClassOutfile.write("%%MatrixMarket matrix array real general\n");
				trainClassOutfile.write(str(numTrainGenes));trainClassOutfile.write("\t1\n");
				for entry in currTrainClasses:
					trainClassOutfile.write(entry);trainClassOutfile.write("\n");
				trainClassOutfile.close();
				
				lambdaFeatureOutfileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + "currCVFeatures.lambda";
				lambdaFeatureOutfile = open(lambdaFeatureOutfileName,'w');
				lambdaFeatureOutfile.write("%%MatrixMarket matrix array real general\n");
				lambdaFeatureOutfile.write(str(numLambdaGenes));lambdaFeatureOutfile.write("\t");
				lambdaFeatureOutfile.write(str(numFeatures));lambdaFeatureOutfile.write("\n");
				for entry in currLambdaValues:
					lambdaFeatureOutfile.write(entry);lambdaFeatureOutfile.write("\n");
				lambdaFeatureOutfile.close();
				
				lambdaClassOutfileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + "currCVClass.lambda";
				lambdaClassOutfile = open(lambdaClassOutfileName,'w');
				lambdaClassOutfile.write("%%MatrixMarket matrix array real general\n");
				lambdaClassOutfile.write(str(numLambdaGenes));lambdaClassOutfile.write("\t1\n");
				for entry in currLambdaClasses:
					lambdaClassOutfile.write(entry);lambdaClassOutfile.write("\n");
				lambdaClassOutfile.close();
				
				testFeatureOutfileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + "currCVFeatures.test";
				testFeatureOutfile = open(testFeatureOutfileName,'w');
				testFeatureOutfile.write("%%MatrixMarket matrix array real general\n");
				testFeatureOutfile.write(str(numTestGenes));testFeatureOutfile.write("\t");
				testFeatureOutfile.write(str(numFeatures));testFeatureOutfile.write("\n");
				for entry in currTestValues:
					testFeatureOutfile.write(entry);testFeatureOutfile.write("\n");
				testFeatureOutfile.close();
				
				testClassOutfileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + "currCVClass.test";
				testClassOutfile = open(testClassOutfileName,'w');
				testClassOutfile.write("%%MatrixMarket matrix array real general\n");
				testClassOutfile.write(str(numTestGenes));testClassOutfile.write("\t1\n");
				for entry in currTestClasses:
					testClassOutfile.write(entry);testClassOutfile.write("\n");
				testClassOutfile.close();
				
				modelFileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + ".model";
				classifyResultFileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + ".result";
				labelPredFileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + "cvTest.preds";
				performanceFileName = resultDrive + condition + "/" + kmers + "_mer/tempProcessingFiles/" + sequenceList1 + "_" + sequenceList2 + "currCV.ROC";
				lambdaOptimal = -1;
				auROCMax = 0;
				lambdaInitial = 0.001;
				for lambdaMult in range(0,1):
					lambdaThis = lambdaInitial + (float(lambdaMult)*(lambdaInitial));
					#feature file, labels, lambda, model output file
					runTrainCmd = "/data/ohler/Dina/Drosophila/Research/code/l1_logreg/l1_logreg_train   " + trainFeatureOutfileName + " " + trainClassOutfileName + " " + str(lambdaThis)  + " " + modelFileName +" > train.out";
					os.system(runTrainCmd);
					# model fiel, test feautres, results
					runClassifyCmd = "/data/ohler/Dina/Drosophila/Research/code/l1_logreg/l1_logreg_classify -q -p " +modelFileName + " " + lambdaFeatureOutfileName + " " + classifyResultFileName +">classify.out";
					os.system(runClassifyCmd);
				
					classifyResultList = list();
					classifyResultFile = open(classifyResultFileName,'r');
					for dummyLine in range(0,7):
						classifyResultFile.readline();
					for line in classifyResultFile:
						classifyResultList.append(line.strip());
					classifyResultFile.close();

					labelPredFile = open(labelPredFileName,'w');
					labelPredFile.write('labels\tpredictions\n');
					for lambdaGeneCounter in range(0,len(currLambdaClasses)):
						labelPredFile.write(currLambdaClasses[lambdaGeneCounter]); 				
						labelPredFile.write("\t"); 				
						labelPredFile.write(classifyResultList[lambdaGeneCounter]); 				
						labelPredFile.write("\n"); 				
					labelPredFile.close();

					runAuROCCmd = "R CMD BATCH  --silent --args " +labelPredFileName+" "+ performanceFileName + "< getAuROC_test.r > test.out ";
					#runAuROCCmd = "R CMD BATCH --no-save \"--args "+ labelPredFileName + " "+performanceFileName+ "\"  getAuROC.r  test.out "
					print runAuROCCmd;
					os.system(runAuROCCmd);

					performanceFile = open(performanceFileName,'r');
					currAuROC = float((performanceFile.readline()).rstrip());
					performanceFile.close();
					if (currAuROC > auROCMax):
						auROCMax = currAuROC;
						lambdaOptimal = lambdaThis;
				
				
				runTrainCmd = "/data/ohler/Dina/Drosophila/Research/code/l1_logreg/l1_logreg_train -q -s " + trainFeatureOutfileName + " " + trainClassOutfileName + " " + str(lambdaOptimal)  + " " + modelFileName +" >train.out";
				os.system(runTrainCmd);

				modelFile = open(modelFileName,'r');
				for lineCount in range(0,9):
					modelFile.readline();
				for regCoeffIndex in range(0,numFeatures+1):
					currRegCoeffList[regCoeffIndex].append(float((modelFile.readline()).strip()));
				modelFile.close();	
				

				runClassifyCmd = "/data/ohler/Dina/Drosophila/Research/code/l1_logreg/l1_logreg_classify -q -p " + modelFileName + " " + testFeatureOutfileName + " " + classifyResultFileName + " >classify_AllTFs.out";
				os.system(runClassifyCmd);
			
				classifyResultList = list();
				classifyResultFile = open(classifyResultFileName,'r');
				for dummyLine in range(0,7):
					classifyResultFile.readline();
				for line in classifyResultFile:
					classifyResultList.append(line.strip());
				classifyResultFile.close();

				labelPredFile = open(labelPredFileName,'w');
				labelPredFile.write('labels\tpredictions\n');
				for testGeneCounter in range(0,len(currTestClasses)):
					labelPredFile.write(currTestClasses[testGeneCounter]); 				
					labelPredFile.write("\t"); 				
					labelPredFile.write(classifyResultList[testGeneCounter]); 				
					labelPredFile.write("\n"); 				
				labelPredFile.close();

			
				runAuROCCmd = "R --vanilla --args  < getAuROC_test.r > test.out";
			#	runAuROCCmd = "R CMD BATCH --no-save \"--args "+ labelPredFileName + " "+performanceFileName+ "\"  getAuROC.r  test.out "
				print runAuROCCmd;				
				os.system(runAuROCCmd);

				performanceFile = open(performanceFileName,'r');
				currAuROC = float((performanceFile.readline()).rstrip());
				performanceFile.close();
				compiledAuROC = compiledAuROC + currAuROC;


		compiledAuROC = compiledAuROC/(cvFold*iterations);
		compareName1 = compareNameList[i];
		compareName2 = compareNameList[j];
		print compareName1,"\t",compareName2;
		currClassifierOutput =  compareName1 + "\t" + compareName2;	
		currClassifierOutput  = currClassifierOutput + "\t" + str(compiledAuROC);
			
		#regCoeffFileNamePrefix = resultDrive + "before_MC/" + kmers + "_mer/" ;# +"/Coeff";
		currRegCoeffFileName = regCoeffFileNamePrefix;  #   + "_" + compareName1 + "_" + compareName2 + ".coeffs";
		currRegCoeffFile = open(currRegCoeffFileName,'w');
		for coeffList in currRegCoeffList:
			for coeff in coeffList:
				currRegCoeffFile.write(str(coeff));currRegCoeffFile.write("\t");
			currRegCoeffFile.write("\n");
		currRegCoeffFile.close();
		
		classifierStatsFile.write(currClassifierOutput);classifierStatsFile.write("\n");

classifierStatsFile.close();
