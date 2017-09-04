average_1=NULL;
average_2=NULL;
####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V3/ClustersDHS/MC_3order/MC_ubiq_output/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
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
  }
}
average = sum/total;
average_1=c(average_1,average);
all_1=allAUC
####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V3/ClustersDHS/MC_3order/MC_other_output_against_all_pos_not/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
 # file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_1=c(average_1,average);
####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V3/ClustersDHS/MC_3order/MC_other_output_against_25/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_25_initialize_redfly_uniq_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_1=c(average_1,average);

####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V3/ClustersDHS/MC_3order/MC_other_output_against_35/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_35_initialize_redfly_uniq_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_1=c(average_1,average);
####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V3/ClustersDHS/MC_3order/MC_other_output_against_37/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_37_initialize_redfly_uniq_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_1=c(average_1,average);



####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/closestGene_uniq/ClustersDHS/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_single_ubiquitous_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_2=c(average_2,average);
all_2 = allAUC;
####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/closestGene_uniq/ClustersDHS/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_single_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  #file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_2=c(average_2,average);
####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/closestGene_uniq/ClustersDHS/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_single_25_pos_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_2=c(average_2,average);
####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/closestGene_uniq/ClustersDHS/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_single_35_pos_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_2=c(average_2,average);
####################################################
setwd("/data/ohler/Dina/Drosophila/Research/results/closestGene_uniq/ClustersDHS/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_single_37_pos_single_Classifier.ROC", sep="")
  print (file3)
  if (file.exists(file3))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
  }
}
average = sum/total;
average_2=c(average_2,average);


#pdf("/data/ohler/Dina/Drosophila/Research/results/Initialization_redfly_V3/Graphs/compare_MC_3order_baseline")
plot(average_1, ylim=c(0,1), xaxt="n", col="red",pch = 19)
points(average_2, col="blue", pch=19)
axis(1,at=c(1,2,3,4,5),labels=c("ubiq","all_other","cluster_25", "cluster_35", "cluster_37"))
#abline(h=average ,col='red',lty=3)
title("MC vs baseline -- uniq genes")
legend("bottomright", inset=.05, title= "aucROC",c("MC","baseline"), fill=c("red","blue"), horiz=TRUE)
#dev.off()
