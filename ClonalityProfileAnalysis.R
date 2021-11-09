#  Loading libraries

suppressPackageStartupMessages({
  library(MAST);  library(immunarch);  library(plyr); library(data.table);  library(reshape);  library(scater); library(ComplexHeatmap);
  library(edgeR);  library(ROCR);  library(scRNA.seq.funcs);  library(gridExtra);  library(NMF);  library(RColorBrewer);  library(ggplot2);
  library(dplyr);  library(MASS);  library(pscl);  library(ape);  library(factoextra);  library(NbClust);  library(igraph);  library(stringdist);
  library(tcR);  library(purrr);  library(tidyr);  library(dendextend);  library(doParallel);  library(foreach);  library(dtplyr);  library(tidyverse);
  library(ggpubr);  library(corrplot);  library(e1071);  library(bootstrap);  library(DAAG);  library(RColorBrewer);  library(ppcor);  library(factoextra);
  library(mixOmics);  library(glmnet);  library(psych);  library(caret);  library(randomForest);  library(PRROC);  library(viridis); library(grid);
  library(ggplotify);  library(ROCR);   library(jpeg);   library(lme4);   library(glmmTMB);   library(reshape);   library(protr);  library(msa);  library(Biostrings);  library(ggseqlogo)
})

chains<-c("TRD","TRG","TRA","TRB","IGH","IGL","IGK")
for (c in chains) {
  print(c)
  irep.ra<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'_WK0.rds')) # This file including only AIRR-seq data from RA patients at baseline should be generated from the "Chain".rds file that can be downloaded from the Rheumatology Research Group website.
  irep.ctrl<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'_CTRL.rds')) # This file including only AIRR-seq data from RA patients at baseline should be generated from the "Chain".rds file that can be downloaded from the Rheumatology Research Group website.
  
  meta_col<-irep.ra$meta
  samples_col<-meta_col$Sample;  data_col<-irep.ra$data[match(samples_col,names(irep.ra$data))];
  irep.ra<-list(data_col,meta_col);  names(irep.ra)<-c("data","meta")
  table(names(irep.ra$data)==as.character(irep.ra$meta$Sample))
  
  irep.all<-c(irep.ra$data,irep.ctrl$data)
  names(irep.all)<-paste0("Sample-",seq(1,length(names(irep.all))))
  pr<-as.data.frame(pubRep(irep.all, "aa", .verbose=F));  rownames(irep.all)<-irep.all$CDR3.aa; rownames(pr)<-pr$CDR3.aa
  
  count<-pr[,c(-1,-2)];
  freq<-as.data.frame(apply(count, 2, function(x) (x/sum(x, na.rm=T))))
  freq$Clone<-rownames(freq)
  data2plot<-gather(freq, condition, measurement, colnames(freq)[1]:colnames(freq)[ncol(freq)-1], factor_key=T)
  
  samples<-colnames(count)
  df.det.clones<-data.frame("Sample"="test", "Detected.Clones"=0, stringsAsFactors=F)
  counter=0
  for (s in samples) {
    detected.clones<-length(which(count[,s]>0))
    counter=counter+1
    df.det.clones[counter,]<-c(s,detected.clones)
  }
  
  colnames(data2plot)[2]<-"Sample"
  data2plot.final<-merge(df.det.clones, data2plot, by="Sample")
  data2plot.final$Detected.Clones<-paste0("Clones (n) = ",data2plot.final$Detected.Clones)
  data2plot.final$Detected.Clones<-as.factor(data2plot.final$Detected.Clones)
  
  # Considering only individual detected clones - different scale
  p<-ggdensity(data2plot.final, x="measurement", rug=T, color="Sample", fill="Sample") +
    labs(title=paste0("\n",c,"\n"), x=paste0("\nClone Frequency (%)\n"), y="\nDensity\n") +
    facet_wrap(~Sample, scales="free") + geom_text(data=data2plot.final, mapping=aes(x=Inf*0.5, y=Inf*0.5, hjust=1, vjust=1.9, label=Detected.Clones), size=2) +
    theme(legend.position="none", axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(hjust=0.5, size=16, face="bold"))
  jpeg(paste0("Output_Data/ClonalityProfile_",c,".jpeg"), res=300, width=6500, height=6500)
  plot(p)
  dev.off()
}
