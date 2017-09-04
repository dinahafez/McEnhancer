library("ROCR")
args <- commandArgs(trailingOnly = TRUE)
inputFile=args[1];
outputFile=args[2];
inputFile;
outputFile;
currPreds <- read.table(inputFile, header = TRUE);
attach(currPreds);
pred <- prediction(currPreds$predictions,currPreds$labels);
perf <- performance(pred,"auc");
auROC <- perf@y.values[[1]];
write(auROC,file = outputFile,ncolumns = 1);

