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


# Principal component analysis

df.guV.sub<-readRDS("Data2analyze/GeneUsageV.rds")
df.guJ.sub<-readRDS("Data2analyze/GeneUsageJ.rds")
df.guVJ.sub<-readRDS("Data2analyze/GeneUsageVJ.rds")

pheno.name="IMID"
week2discard="wk12"
pheno2pred="RA"
phenoOther="CTRL"
grr.meta<-read.csv("Data2analyze/Metadata.csv", sep="\t");  rownames(grr.meta)<-grr.meta$Sample;
grr.meta$RespGM<-grr.meta$Resp
grr.meta$RespGM[(grr.meta$RespGM=="MOD" & !is.na(grr.meta$RespGM))]<-"GOOD"
grr.meta$Resp[grr.meta$Resp=="MOD"]<-NA
col.week<-which(colnames(grr.meta)=="Week")
col.age<-which(colnames(grr.meta)=="Age")
col.gender<-which(colnames(grr.meta)=="Gender")
col.pheno.name<-which(colnames(grr.meta) %in% pheno.name)
cols2keep<-c(col.week, col.age, col.gender, col.pheno.name)
meta00<-grr.meta[,cols2keep];   colnames(meta00)[4]<-"Pheno"
meta0<-meta00[(meta00$Week!=week2discard & !(is.na(meta00$Pheno)) & (meta00$Pheno==pheno2pred | meta00$Pheno==phenoOther)),]
dupl.sample.exc<-c("IXTCB01476","IXTCB01395","IXTCB01434")
meta<-meta0[!(rownames(meta0) %in% dupl.sample.exc),]

filt.gu.005<-function(df.gu, meta) {
  df.guS<-merge(meta, df.gu, by="row.names")
  df.guS.filt<-df.guS[!is.na(df.guS$Row.names),];   rownames(df.guS.filt)<-df.guS.filt$Row.names
  df.guS.filt<-df.guS.filt[,-1]
  
  filt.cutoff<-round(nrow(df.guS.filt)*0.05)
  cols2select<-vector();  c=0;
  for (j in 1:ncol(df.guS.filt)) {
    if (length(which(df.guS.filt[,j]==0))<filt.cutoff) {c=c+1; cols2select[c]<-j}
  }
  
  df.guS.filt2<-df.guS.filt[,cols2select]
  df.guS2analyze<-df.guS.filt2[,grep("TRA|TRB|TRD|TRG|IGH|IGL|IGK",colnames(df.guS.filt2))]
  return(df.guS2analyze)
}

df.guV2analyze.sub.filt005<-filt.gu.005(df.guV.sub, meta)
df.guJ2analyze.sub.filt005<-filt.gu.005(df.guJ.sub, meta)
df.guVJ2analyze.sub.filt005<-filt.gu.005(df.guVJ.sub, meta)

pheno.data<-meta$Pheno

make.pc1pc2.plot<-function(gene.seg, df.gu, pheno.name, pheno.data) {
  gs<-gene.seg
  pd<-pheno.data
  df.gu<-df.gu
  pheno<-pheno.name
  data2print<-as.matrix(df.gu)
  #M<-cor(data2print);  #res1<-cor.mtest(data2print, conf.level=.95)
  res.pca<-prcomp(data2print, scale=T)
  df.pca<-data.frame("IMID"=pd,"PC1"=res.pca$x[,1], "PC2"=res.pca$x[,2])
  var.pc1<-paste0(round(get_eigenvalue(res.pca)[1,2],2),"%")
  var.pc2<-paste0(round(get_eigenvalue(res.pca)[2,2],2),"%")
  g1<-ggplot(df.pca, aes(x=PC1, y=PC2, color=IMID)) + geom_point() + theme_classic() +
    xlim(min(df.pca$PC1)*1.2,max(df.pca$PC1)*1.2) + ylim(min(df.pca$PC2)*1.2,max(df.pca$PC2)*1.2) +
    labs(x=paste0("\nPC1 (",var.pc1,")\n"), y=paste0("\nPC2 (",var.pc2,")\n"), title=paste0("\n",gene.seg," segment\n")) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) + 
    scale_color_manual(values=c("#7DBB7D","#F69191")) + geom_hline(yintercept=0, linetype="dashed", color = "gray") +
    geom_vline(xintercept=0, linetype="dashed", color = "gray") +
    theme(legend.title=element_blank(), legend.position="top", plot.title=element_text(hjust=0.5, face="bold", size=16), legend.key.size=unit(0.5, "cm"), text=element_text(size=12))
  return(g1)
}

gene.seg<-"V";   df.gu.final<-df.guV2analyze.sub.filt005;  pheno.name<-pheno.name;
jpeg(paste0("Output_Data/PCA_",gene.seg,".jpeg"), width=1900, height=1900, res=300)
make.pc1pc2.plot(gene.seg, df.gu.final, pheno.name, pheno.data)
dev.off()

gene.seg<-"J";   df.gu.final<-df.guJ2analyze.sub.filt005;  pheno.name<-pheno.name;
jpeg(paste0("Output_Data/PCA_",gene.seg,".jpeg"), width=1900, height=1900, res=300)
make.pc1pc2.plot(gene.seg, df.gu.final, pheno.name, pheno.data)
dev.off()

gene.seg<-"VJ";   df.gu.final<-df.guVJ2analyze.sub.filt005;  pheno.name<-pheno.name;
jpeg(paste0("Output_Data/PCA_",gene.seg,".jpeg"), width=1900, height=1900, res=300)
make.pc1pc2.plot(gene.seg, df.gu.final, pheno.name, pheno.data)
dev.off()

# RA vs Controls

files<-list(df.guV2analyze.sub.filt005, df.guJ2analyze.sub.filt005, df.guVJ2analyze.sub.filt005)
names(files)<-c("V_SubFam_Filt005", "J_SubFam_Filt005", "VJ_SubFam_Filt005")
p1<-rownames(meta[(meta$Pheno==pheno2pred & meta$Week!=week2discard),])
p2<-rownames(meta[(meta$Pheno==phenoOther & meta$Week!=week2discard),])

reads.sample<-readRDS("Data2analyze/ChainReadsCount.rds")
reads.sample$Reads.All<-reads.sample$TRA+reads.sample$TRB+reads.sample$TRD+reads.sample$TRG+reads.sample$IGH+reads.sample$IGL+reads.sample$IGK

# Distribution after filtering out segments with freq<5%
long=melt(df.guV2analyze.sub.filt005)
ggplot(long, aes (value)) + geom_density() + facet_wrap(~variable, scales="free")

# Statistical test
res.all<-list()
count=0
for (n in names(files)) {
  gs<-strsplit(n,'_')[[1]][1];   fam<-strsplit(n,'_')[[1]][2];   filt<-strsplit(n,'_')[[1]][3]
  df<-files[[n]]
  df.out<-data.frame("Segment"=rep("NA",ncol(df)), "Beta"=rep("NA",ncol(df)), "Pvalue"=rep("NA",ncol(df)),"Pvalue.Wilcox"=rep("NA",ncol(df)), stringsAsFactors=F)
  for (i in 1:ncol(df)) {
    seg<-colnames(df)[i]
    p1.seg<-df[p1,seg]
    p2.seg<-df[p2,seg]
    p1.ids<-p1
    p2.ids<-p2
    p1.pheno<-rep("RA",length(p1.seg))
    p2.pheno<-rep("CTRL",length(p2.seg))
    df.seg<-data.frame("Pheno"=c(p1.pheno,p2.pheno),"values"=c(p1.seg,p2.seg), stringsAsFactors=F)
    rownames(df.seg)<-c(p1.ids,p2.ids)
    df.seg.merged<-merge(df.seg,meta,by="row.names")
    colnames(df.seg.merged)[1]<-"Sample"
    df.seg.final<-merge(df.seg.merged,reads.sample[,c("Sample","Reads.All")], by="Sample")
    df.seg.final$Pheno.x[df.seg.final$Pheno.x=="RA"]<-1
    df.seg.final$Pheno.x[df.seg.final$Pheno.x=="CTRL"]<-0
    df.seg.final$Pheno.x<-as.numeric(df.seg.final$Pheno.x)
    model<-glm(Pheno.x~values+Age+Reads.All, data=df.seg.final, family="binomial")
    pval<-coef(summary(model))["values",'Pr(>|z|)']
    beta<-coef(summary(model))["values",'Estimate']
    pval.wc<-wilcox.test(p1.seg, p2.seg)$p.value
    df.out[i,1]<-seg
    df.out[i,2]<-beta
    df.out[i,3]<-pval
    df.out[i,4]<-pval.wc
  }
  df.out$Pvalue<-as.numeric(df.out$Pvalue)
  df.out.ord<-df.out[order(df.out$Pvalue),]
  df.out.ord$FDR<-p.adjust(df.out.ord$Pvalue, method="fdr")
  write.csv(df.out.ord[,c(1,2,3,5)], paste0("Output_Data/SegmentUsage_RACTRL_",gs,".csv"), row.names=F)
}
  
# Statistical test
res.all<-list()
count=0
phenos<-c("RespGM","ActB","ACPA","RF")
metadata<-read.csv("Data2analyze/Metadata.csv", sep="\t")
metadata$RespGM<-metadata$Resp
metadata$RespGM[metadata$RespGM=="MOD"]<-"GOOD"
metadata<-metadata[metadata$IMID=="RA",]

for (ph in phenos) {
  df.out<-data.frame("Segment"="NA", "Beta"=10, "Pvalue"=10,"Pvalue.Wilcox"=10, stringsAsFactors=F)
  pheno.col<-which(colnames(metadata)==ph)
  colnames(metadata)[pheno.col]<-"Pheno"
  phs<-sort(as.character(unique(metadata$Pheno[!is.na(metadata$Pheno)])))
  p1.name<-phs[1]
  p2.name<-phs[2]
  p1<-as.character(metadata[metadata$Pheno==p1.name & !is.na(metadata$Pheno),"Sample"])
  p2<-as.character(metadata[metadata$Pheno==p2.name & !is.na(metadata$Pheno),"Sample"])
  print(p1)
  print(p2)
  print(p1.name)
  print(p2.name)
  for (n in names(files)) {
    gs<-strsplit(n,'_')[[1]][1];   fam<-strsplit(n,'_')[[1]][2];   filt<-strsplit(n,'_')[[1]][3]
    df<-files[[n]]
    for (i in 1:ncol(df)) {
      seg<-colnames(df)[i]
      p1.seg<-df[p1,seg]
      p2.seg<-df[p2,seg]
      p1.ids<-p1
      p2.ids<-p2
      p1.pheno<-rep(p1.name,length(p1.seg))
      p2.pheno<-rep(p2.name,length(p2.seg))
      df.seg<-data.frame("Pheno"=c(p1.pheno,p2.pheno),"values"=c(p1.seg,p2.seg), stringsAsFactors=F)
      rownames(df.seg)<-c(p1.ids,p2.ids)
      df.seg.merged<-merge(df.seg,meta,by="row.names")
      colnames(df.seg.merged)[1]<-"Sample"
      df.seg.final<-merge(df.seg.merged,reads.sample[,c("Sample","Reads.All")], by="Sample")
      df.seg.final$Pheno.x[df.seg.final$Pheno.x==p1.name]<-1
      df.seg.final$Pheno.x[df.seg.final$Pheno.x==p2.name]<-0
      df.seg.final$Pheno.x<-as.numeric(df.seg.final$Pheno.x)
      model<-glm(Pheno.x~values+Age+Reads.All, data=df.seg.final, family="binomial")
      pval<-coef(summary(model))["values",'Pr(>|z|)']
      beta<-coef(summary(model))["values",'Estimate']
      pval.wc<-wilcox.test(p1.seg, p2.seg)$p.value
      df.out[i,1]<-seg
      df.out[i,2]<-beta
      df.out[i,3]<-pval
      df.out[i,4]<-pval.wc
    }
    df.out$Pvalue<-as.numeric(df.out$Pvalue)
    df.out.ord<-df.out[order(df.out$Pvalue),]
    df.out.ord$FDR<-p.adjust(df.out.ord$Pvalue, method="fdr")
    df.out.ord$Pvalue.Wilcox<-as.numeric(df.out.ord$Pvalue.Wilcox)
    df.out.ord<-df.out.ord[order(df.out.ord$Pvalue.Wilcox),]
    df.out.ord$FDR.Wilcox<-p.adjust(df.out.ord$Pvalue.Wilcox, method="fdr")
    write.csv(df.out.ord[,c(1,2,3,4,5,6)], paste0("Output_Data/SegmentUsage_",ph,"_",gs,".csv"), row.names=F) 
  }
}


# Longitudinal
  
df.guV.sub<-readRDS("Data2analyze/GeneUsageV.rds")
df.guJ.sub<-readRDS("Data2analyze/GeneUsageJ.rds")
df.guVJ.sub<-readRDS("Data2analyze/GeneUsageVJ.rds")

  
  pheno.name="Week"
  week2discard="CTRL"
  pheno2pred="wk0"
  phenoOther="wk12"
  grr.meta<-read.csv("Data2analyze/Metadata.csv", sep="\t");  rownames(grr.meta)<-grr.meta$Sample;
  grr.meta$RespGM<-grr.meta$Resp
  grr.meta$RespGM[(grr.meta$RespGM=="MOD" & !is.na(grr.meta$RespGM))]<-"GOOD"
  grr.meta$Resp[grr.meta$Resp=="MOD"]<-NA
  col.week<-which(colnames(grr.meta)=="Week")
  col.age<-which(colnames(grr.meta)=="Age")
  col.gender<-which(colnames(grr.meta)=="Gender")
  col.pheno.name<-which(colnames(grr.meta) %in% pheno.name)
  cols2keep<-c(col.week, col.age, col.gender, col.pheno.name)
  meta00<-grr.meta[,cols2keep];   colnames(meta00)[4]<-"Pheno"
  meta<-meta00[(meta00$Week!=week2discard & !(is.na(meta00$Pheno)) & (meta00$Pheno==pheno2pred | meta00$Pheno==phenoOther)),]
  
  filt.gu.005<-function(df.gu, meta) {
    df.guS<-merge(meta, df.gu, by="row.names")
    df.guS.filt<-df.guS[!is.na(df.guS$Row.names),];   rownames(df.guS.filt)<-df.guS.filt$Row.names
    df.guS.filt<-df.guS.filt[,-1]
    
    filt.cutoff<-round(nrow(df.guS.filt)*0.05)
    cols2select<-vector();  c=0;
    for (j in 1:ncol(df.guS.filt)) {
      if (length(which(df.guS.filt[,j]==0))<filt.cutoff) {c=c+1; cols2select[c]<-j}
    }
    
    df.guS.filt2<-df.guS.filt[,cols2select]
    df.guS2analyze<-df.guS.filt2[,grep("TRA|TRB|TRD|TRG|IGH|IGL|IGK",colnames(df.guS.filt2))]
    return(df.guS2analyze)
  }
  
  df.guV2analyze.sub.filt005<-filt.gu.005(df.guV.sub, meta)
  df.guJ2analyze.sub.filt005<-filt.gu.005(df.guJ.sub, meta)
  df.guVJ2analyze.sub.filt005<-filt.gu.005(df.guVJ.sub, meta)
  
  files<-list(df.guV2analyze.sub.filt005, df.guJ2analyze.sub.filt005, df.guVJ2analyze.sub.filt005)
  names(files)<-c("V_SubFam_Filt005", "J_SubFam_Filt005", "VJ_SubFam_Filt005")
  p1<-rownames(meta[(meta$Pheno==pheno2pred & meta$Week!=week2discard),])[rownames(meta[(meta$Pheno==pheno2pred & meta$Week!=week2discard),]) %in% as.character(donor.data.filt$RNA1)]
  p2<-rownames(meta[(meta$Pheno==phenoOther & meta$Week!=week2discard),])[rownames(meta[(meta$Pheno==phenoOther & meta$Week!=week2discard),]) %in% as.character(donor.data.filt$RNA2)]
  
  reads.sample.wk12<-readRDS("Data2analyze/ChainReadsCountWeek12.rds")
  reads.sample.wk0<-readRDS("Data2analyze/ChainReadsCount.rds")
  reads.sample<-rbind(reads.sample.wk0,reads.sample.wk12)
  reads.sample$Reads.All<-reads.sample$TRA+reads.sample$TRB+reads.sample$TRD+reads.sample$TRG+reads.sample$IGH+reads.sample$IGL+reads.sample$IGK
  
  # Statistical test
  res.all<-list()
  count=0
  for (n in names(files)) {
    gs<-strsplit(n,'_')[[1]][1];   fam<-strsplit(n,'_')[[1]][2];   filt<-strsplit(n,'_')[[1]][3]
    df<-files[[n]]
    df.out<-data.frame("Segment"=rep("NA",ncol(df)), "Beta"=rep("NA",ncol(df)), "Pvalue"=rep("NA",ncol(df)), stringsAsFactors=F)
    for (i in 1:ncol(df)) {
      seg<-colnames(df)[i]
      p1.seg<-df[p1,seg]
      p2.seg<-df[p2,seg]
      p1.ids<-p1
      p2.ids<-p2
      p1.pheno<-rep("Baseline",length(p1.seg))
      p2.pheno<-rep("Treated",length(p2.seg))   
      df.seg<-data.frame("Pheno"=c(p1.pheno,p2.pheno),"values"=c(p1.seg,p2.seg), stringsAsFactors=F)
      rownames(df.seg)<-c(p1.ids,p2.ids)   
      df.seg.merged<-merge(df.seg,meta,by="row.names")
      colnames(df.seg.merged)[1]<-"Sample"
      donor2merge<-data.frame("Sample"=c(as.character(donor.data.filt$RNA1),as.character(donor.data.filt$RNA2)), "Donor"=c(as.character(donor.data.filt$Donor),as.character(donor.data.filt$Donor)))
      df.seg.merged.v2<-merge(df.seg.merged,donor2merge,by="Sample")
      df.seg.final<-merge(df.seg.merged.v2,reads.sample[,c("Sample","Reads.All")], by="Sample")
      df.seg.final$Pheno.x[df.seg.final$Pheno.x=="Treated"]<-1
      df.seg.final$Pheno.x[df.seg.final$Pheno.x=="Baseline"]<-0
      df.seg.final$Pheno.x<-as.numeric(df.seg.final$Pheno.x)
      model1<-glmer(Pheno.x~values+Age+Reads.All+(1|Donor),  family=binomial, data=df.seg.final)
      pval<-summary(model1)$coefficients["values","Pr(>|z|)"]
      beta<-summary(model1)$coefficients["values","Estimate"]
      df.out[i,1]<-seg
      df.out[i,2]<-beta
      df.out[i,3]<-pval
    }
    df.out$Pvalue<-as.numeric(df.out$Pvalue)
    df.out.ord<-df.out[order(df.out$Pvalue),]
    df.out.ord$FDR<-p.adjust(df.out.ord$Pvalue, method="fdr")
    write.csv(df.out.ord[,c(1,2,3,4)], paste0("Output_Data/SegmentUsage_longitudinal_",gs,".csv"), row.names=F)
  }


