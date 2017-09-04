library("ROCR")
inputFile = "/data/ohler/Dina/Drosophila/Research/results/MarkovChain_classification/before_MC/5_mer/tempProcessingFiles/11_pos_all_pos_not_11cvTest.preds";
outputFile = "/data/ohler/Dina/Drosophila/Research/results/MarkovChain_classification/before_MC/5_mer/tempProcessingFiles/11_pos_all_pos_not_11cvTest.ROC";
inputFile;
outputFile;
currPreds <- read.table(inputFile, header = TRUE);
attach(currPreds);
pred <- prediction(currPreds$predictions,currPreds$labels);
perf <- performance(pred,"auc");
auROC <- perf@y.values[[1]];
write(auROC,file = outputFile,ncolumns = 1);

