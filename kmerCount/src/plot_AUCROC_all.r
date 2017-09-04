average_1=NULL;
average_2=NULL;

Stark_dir="/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/5_mer/";
Stark_dir_stage2="/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/MC_2nd_3order/MC_ubiq_output/5_mer/";
####################################################
####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//5_mer/")
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_Eileen_known_enhancers_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output/5_mer/")
allAUC=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  #if ((i != 27 )&&(i!=34)&&(i!=39))
  {
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  
  
  file_stark=paste(Stark_dir,file3,sep="")
  if (file.exists(file3))
    #if (file.exists(file_stark) && (i != 27))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC = c(allAUC,auc)
    sum = sum+auc;
    total=total+1;
    print (file3)
  }
}
}
average = sum/total;
average_1=c(average_1,average);

####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//MC_2nd_3order/MC_ubiq_output/5_mer/")
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output_0.1/5_mer/")
allAUC2=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
 # if ((i != 27 )&&(i!=34)&&(i!=39))
  {
    file = paste(i,"_pos_all_pos_not_", sep="")
    file2 = paste(file,i,sep="")
    file3 = paste(file2,"_Classifier.ROC",sep="")
    file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
        
    
  #if (file.exists(file3))
  file_stark=paste(Stark_dir_stage2,file3,sep="")
  if (file.exists(file3))
    #if (file.exists(file_stark) && (i != 27))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC2 = c(allAUC2,auc)
    sum = sum+auc;
    total=total+1;
    print (file3)
  }
}
}
average = sum/total;
average_1=c(average_1,average);

####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//5_mer/")
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output_0.1/5_mer_one_vector/")
allAUC3=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  if ((i != 27 )&&(i!=34)&&(i!=39))
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  #print (file3)
  
  file_stark=paste(Stark_dir,file3,sep="")
  if (file.exists(file3))
    #if (file.exists(file_stark) && (i != 27))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC3 = c(allAUC3,auc)
    sum = sum+auc;
    total=total+1;
  }
}
}
average = sum/total;
average_1=c(average_1,average);

####################################################
#setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_redFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output//MC_2nd_3order/MC_ubiq_output/5_mer/")
setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/MC_ubiq_output_0.1/5_mer/")
allAUC4=NULL;
sum = 0;
total=0;
for(i in 11:39)
{
  if ((i != 27 )&&(i!=34)&&(i!=39))
{
  file = paste(i,"_pos_all_pos_not_", sep="")
  file2 = paste(file,i,sep="")
  file3 = paste(file2,"_Classifier.ROC",sep="")
  file3 = paste(i,"_pos_ubiquitous_single_Classifier.ROC", sep="")
  #print (file3)
  
  
  #if (file.exists(file3))
  file_stark=paste(Stark_dir_stage2,file3,sep="")
  if (file.exists(file3))
    #if (file.exists(file_stark) && (i != 27))
  {
    auc = read.table(file3, header=FALSE)$V3
    allAUC4 = c(allAUC4,auc)
    sum = sum+auc;
    total=total+1;
  }
}
}
average = sum/total;
average_1=c(average_1,average);

#########PLOT
#pdf("/data/ohler/Dina/Drosophila/Research/results/Graphs/Initialization_both_cutoff_0.1_multiple_vs_one_DHS_feature_vector.pdf")

par(ann=FALSE)
plot(jitter(rep(0,length(allAUC)),amount=0.25), allAUC,
     xlim=range(-0.5,3.5), ylim=range(0.5,1),
     axes=FALSE,frame.plot=TRUE)
points(jitter(rep(1,length(allAUC2)), amount=0.25), allAUC2, col=2)
points(jitter(rep(2,length(allAUC3)), amount=0.25), allAUC3, col=1)
points(jitter(rep(3,length(allAUC4)), amount=0.25), allAUC4, col=2)

##Add in the y-axis
axis(2, seq(0,1,by=0.1))

##Add in the x-axis labels
mtext("DHS_one_vec", side = 1, at=-0.25)
mtext("DHS_multiple_vec", side = 1, at=1)
mtext("DHS_one_filt", side = 1, at=2)
mtext("DHS_multiple_filt", side = 1, at=3.0)
##Add in the means
segments(-0.25, mean(allAUC), 0.25, mean(allAUC))
segments(0.75, mean(allAUC2), 1.25, mean(allAUC2))
segments(1.75, mean(allAUC3), 2.25, mean(allAUC3))
segments(2.75, mean(allAUC4), 3.25, mean(allAUC4))
title("Initialization both, cutoff 0.1, multiple vs one DHS feature vector")
#dev.off()
