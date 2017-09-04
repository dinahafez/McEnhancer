import os;
import string;
import numpy as np;
import sys;
from random import shuffle;

#runClassifyCmd = "R CMD BATCH --no-save getAuROC_test.r test.out ";  #working
runClassifyCmd = "R CMD BATCH --no-save \"--args /data/ohler/Dina/Drosophila/Research/results/MarkovChain_classification/before_MC/5_mer/tempProcessingFiles/12_pos_all_pos_not_12cvTest.preds /data/ohler/Dina/Drosophila/Research/results/MarkovChain_classification/before_MC/5_mer/tempProcessingFiles/12_pos_all_pos_not_12currCV.ROC\"  getAuROC.r  test.out "
#working !!!!
os.system(runClassifyCmd);
