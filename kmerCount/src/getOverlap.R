setwd("/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/");
par(mfrow=c(4,3)) 
for(i in 28:39)
{
  #file = paste("/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_RedFly_uniqGenes/ClustersDHS/MC_3order/Overlap/",i, sep="")
  file2 = paste(i,"_overlap.bed",sep="")
  print (file2)
  if ( file.info(file2)$size != 0 )
  {
      overlap = read.table(file2, h=F)
      figfile = paste("/data/ohler/Dina/Drosophila/Research/results/Initialization_Stark_uniqGenes/ClustersDHS/MC_3order/Overlap/histograms/",i, sep="")
      fig = paste(figfile,"_overlap.pdf", sep="")
      #pdf(fig)
      t = paste("cluster ",i)
      hist(overlap$V5, breaks = 100, main=t, xlab="#common DHSs")
      
      t = paste("cluster ",i)
      title(t)
      #dev.off()
  }
}