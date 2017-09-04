library("ROCR")
#args <- commandArgs(trailingOnly = TRUE)
#inputFile=args[1];
#outputFile=args[2];
inputFile="/data/ohler/Dina/Drosophila/Research/results/MarkovChain_5order_classification/after_MC/6_mer/tempProcessingFiles/11_pos_single20.pos_ubiquitous.pred";
#outputFile;
currPreds <- read.table(inputFile, header = TRUE);

pred <- prediction(currPreds$predictions,currPreds$labels);
perf <- performance(pred,"auc");
auROC <- perf@y.values[[1]];
auROC
#write(auROC,file = outputFile,ncolumns = 1);

