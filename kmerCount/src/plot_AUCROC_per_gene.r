setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/ClustersDHS/MC_5order/MC_ubiq_output/Cluster_11/")
allAUC=NULL;
sumpercluster = 0;
countpercluster=0;

#pdf("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/Graphs/MC_5order_classify_5kmer_ubiquitous")

for(i in 11:39)
{
  sumpercluster = 0;
  countpercluster=0;
  clusterAUC = NULL;
  
  dir = paste("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V2/ClustersDHS/MC_5order/MC_ubiq_output/Cluster_",i, sep  = "");
  dir1 = paste(dir, "/5_mer/", sep="");
  setwd(dir1);
        files = list.files(pattern="*.ROC")
  for (afile in files)
  {
    auc = read.table(afile, header=FALSE)$V3
    clusterAUC = c(clusterAUC,auc)
    sumpercluster = sumpercluster + auc;
    countpercluster = countpercluster+1;
    
  }  
  allAUC = cbind(allAUC,clusterAUC)
 
}
boxplot(allAUC,names= seq (11,39),outline=F, ylim=c(0.5,1))

#title("Initialize redfly MC 5 order 5mer UNMERGED: one against ubiquitous")
#dev.off()