####  Loading libraries  ####

suppressPackageStartupMessages({
  library(MAST);  library(immunarch);  library(plyr); library(data.table);  library(reshape);  library(scater); library(ComplexHeatmap);
  library(edgeR);  library(ROCR);  library(scRNA.seq.funcs);  library(gridExtra);  library(NMF);  library(RColorBrewer);  library(ggplot2);
  library(dplyr);  library(MASS);  library(pscl);  library(ape);  library(factoextra);  library(NbClust);  library(igraph);  library(stringdist);
  library(tcR);  library(purrr);  library(tidyr);  library(dendextend);  library(doParallel);  library(foreach);  library(dtplyr);  library(tidyverse);
  library(ggpubr);  library(corrplot);  library(e1071);  library(bootstrap);  library(DAAG);  library(RColorBrewer);  library(ppcor);  library(factoextra);
  library(mixOmics);  library(glmnet);  library(psych);  library(caret);  library(randomForest);  library(PRROC);  library(viridis); library(grid);
  library(ggplotify);  library(ROCR);   library(jpeg);   library(lme4);   library(glmmTMB);   library(reshape);   library(protr);  library(msa);  library(Biostrings);  library(ggseqlogo)
})

# Clones defined by CDR3aa & Data after sample+replicate merging

x1<-readRDS("Data2analyze/Metadata_TRG.rds")
x2<-read.table("Data2analyze/Donor2RNA_Codes.txt", header=T)
ra.wk0<-as.character(x1[x1$Week=="wk0", "Sample"]);     ra.wk12<-as.character(x2$RNA2);     ra.all<-c(ra.wk0,ra.wk12);
ctrl<-as.character(x1[x1$IMID=="CTRL", "Sample"]);      all.wk0<-c(ra.wk0,ctrl);            all.wk12<-c(ra.wk12,ctrl);         all<-c(ra.all,ctrl)

chains<-c("TRA","TRB","TRD","TRG","IGH","IGL","IGK")
vars<-c("Reads","Reads.Sample","Unique","Unique.Sample")

df.airrseq<-data.frame("Chain"="test","AIRRseq.Feature"="test","RA.All"=0,"RA.Baseline"=0,"RA.Week12"=0,"Controls"=0,"Controls.RA.Baseline"=0,"Controls.RA.Week12"=0,"All"=0, stringsAsFactors=F)
count=0
for (c in chains) {
  print(c)
  chain.data<-readRDS(paste0("Data2analyze/Clones_By_CDR3aa/",c,".rds"))
  for (v in vars) {
    chain.data.ra.all<-chain.data$data[(names(chain.data$data) %in% ra.all)]
    chain.data.ra.wk0<-chain.data$data[(names(chain.data$data) %in% ra.wk0)]
    chain.data.ra.wk12<-chain.data$data[(names(chain.data$data) %in% ra.wk12)]
    chain.data.ctrl<-chain.data$data[(names(chain.data$data) %in% ctrl)]
    chain.data.all.wk0<-chain.data$data[(names(chain.data$data) %in% all.wk0)]
    chain.data.all.wk12<-chain.data$data[(names(chain.data$data) %in% all.wk12)]
    chain.data.all<-chain.data$data[(names(chain.data$data) %in% all)]
    if (v=="Reads") {
      out.name<-"Reads (N)"
      ra.all.data<-sum(do.call(rbind, chain.data.ra.all)$Clones)
      ra.wk0.data<-sum(do.call(rbind, chain.data.ra.wk0)$Clones)
      ra.wk12.data<-sum(do.call(rbind, chain.data.ra.wk12)$Clones)
      all.wk0.data<-sum(do.call(rbind, chain.data.all.wk0)$Clones)
      all.wk12.data<-sum(do.call(rbind, chain.data.all.wk12)$Clones)
      ctrl.data<-sum(do.call(rbind, chain.data.ctrl)$Clones)
      all.data<-sum(do.call(rbind, chain.data.all)$Clones)
    }
    if (v=="Reads.Sample") {
      out.name<-"Reads per sample (m)"
      ra.all.data<-round(sum(do.call(rbind, chain.data.ra.all)$Clones)/length(chain.data.ra.all),2)
      ra.wk0.data<-round(sum(do.call(rbind, chain.data.ra.wk0)$Clones)/length(chain.data.ra.wk0),2)
      ra.wk12.data<-round(sum(do.call(rbind, chain.data.ra.wk12)$Clones)/length(chain.data.ra.wk12),2)
      all.wk0.data<-round(sum(do.call(rbind, chain.data.all.wk0)$Clones)/length(chain.data.all.wk0),2)
      all.wk12.data<-round(sum(do.call(rbind, chain.data.all.wk12)$Clones)/length(chain.data.all.wk12),2)
      ctrl.data<-round(sum(do.call(rbind, chain.data.ctrl)$Clones)/length(chain.data.ctrl),2)
      all.data<-round(sum(do.call(rbind, chain.data.all)$Clones)/length(chain.data.all),2)
    }
    if (v=="Unique") {
      out.name<-"Unique clones (N)"
      ra.all.data<-length(unique(do.call(rbind, chain.data.ra.all)$CDR3.aa))
      ra.wk0.data<-length(unique(do.call(rbind, chain.data.ra.wk0)$CDR3.aa))
      ra.wk12.data<-length(unique(do.call(rbind, chain.data.ra.wk12)$CDR3.aa))
      all.wk0.data<-length(unique(do.call(rbind, chain.data.all.wk0)$CDR3.aa))
      all.wk12.data<-length(unique(do.call(rbind, chain.data.all.wk12)$CDR3.aa))
      ctrl.data<-length(unique(do.call(rbind, chain.data.ctrl)$CDR3.aa))
      all.data<-length(unique(do.call(rbind, chain.data.all)$CDR3.aa))
    }
    if (v=="Unique.Sample") {
      out.name<-"Unique clones per sample (m)"
      ra.all.data<-round(length(unique(do.call(rbind, chain.data.ra.all)$CDR3.aa))/length(chain.data.ra.all),2)
      ra.wk0.data<-round(length(unique(do.call(rbind, chain.data.ra.wk0)$CDR3.aa))/length(chain.data.ra.wk0),2)
      ra.wk12.data<-round(length(unique(do.call(rbind, chain.data.ra.wk12)$CDR3.aa))/length(chain.data.ra.wk12),2)
      all.wk0.data<-round(length(unique(do.call(rbind, chain.data.all.wk0)$CDR3.aa))/length(chain.data.all.wk0),2)
      all.wk12.data<-round(length(unique(do.call(rbind, chain.data.all.wk12)$CDR3.aa))/length(chain.data.all.wk12),2)
      ctrl.data<-round(length(unique(do.call(rbind, chain.data.ctrl)$CDR3.aa))/length(chain.data.ctrl),2)
      all.data<-round(length(unique(do.call(rbind, chain.data.all)$CDR3.aa))/length(chain.data.all),2)
    }
    count=count+1
    df.airrseq[count,]<-c(c,out.name,ra.all.data,ra.wk0.data,ra.wk12.data,ctrl.data,all.wk0.data,all.wk12.data,all.data)
  }
}

write.csv(df.airrseq, "Output_Data/Table_AIRRseq.csv", row.names=F)