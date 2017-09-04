setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_Eileen_known_enhancers_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
#pdf("/data/ohler/Dina/Drosophila/Research/results/Graphs/Initialization_Stark_Redfly_uniq_genes_only_model_selection.pdf")
xe=NULL;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
  auc = read.table(file3, header=FALSE)$V3
  allAUC = c(allAUC,auc)
  sum = sum+auc;
  total=total+1;
  xe=c(xe,as.numeric(i));
  }
  else
  {
    allAUC = c(allAUC,-1)
  }
}

plot(allAUC, ylim=c(0.5,1), xaxt="n");
axis(1,at=xe-10, labels=xe);
average = sum/total;
abline(h=average ,col='red',lty=3)

title("Initialization Eileen known enhancers: unique genes only, one feature vector")
#dev.off()