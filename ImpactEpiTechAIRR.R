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

# Correlation SeqDepth - Plate, Integrity, Concentration

tech.var<-read.table("Data2analyze/Technical_Var.txt", header=T)
data<-readRDS("Data2analyze/CloneDet.rds")
data<-data[(data$Sample %in% c(ra.wk0,ctrl)),]

dfm<-merge(data,tech.var,by="Sample")
dfm$Reads<-as.numeric(dfm$Reads)
dfm$Plate<-as.factor(dfm$Plate)

seqdepth.plate<-ggplot(dfm, aes(x = fct_reorder(Plate, Reads, .fun=median, .desc=T), y = Reads)) +  geom_boxplot(aes(fill = fct_reorder(Plate, Reads, .fun = median, .desc =T))) +
  geom_jitter(position=position_jitter(0.2)) + facet_wrap(~Chain, scales="free") +  scale_fill_manual(values=brewer.pal(5,name="PuOr")) +
  labs(title="\nPlate vs. Sequencing Depth\n", x="\nPlate\n", y="\nSequencing depth (N reads)\n") + theme_classic() +
  theme(legend.position="none", plot.title=element_text(hjust=0.5, size=14, face="bold")) +
  stat_compare_means(method="kruskal.test", label="p.format", label.x.npc=0.4, hide.ns=T)

seqdepth.integrity<-ggplot(dfm, aes(x=Reads, y=RNA.Integrity)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nRNA Integrity vs. Sequencing Depth\n", x="\nSequencing Depth (N, reads)\n", y="\n\nRNA Integrity (RIN)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.1)

seqdepth.concentration<-ggplot(dfm, aes(x=Reads, y=RNA.Concentration)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nRNA Concentration vs. Sequencing Depth\n", x="\nSequencing Depth (N, reads)\n", y="\n\nRNA Concentration (ng/ul)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.1)

integrity.plate<-ggplot(dfm[dfm$Chain=="TRA",], aes(x = fct_reorder(Plate, RNA.Integrity, .fun=median, .desc=T), y = RNA.Integrity)) +  geom_boxplot(aes(fill = fct_reorder(Plate, RNA.Integrity, .fun = median, .desc =T))) +
  geom_jitter(position=position_jitter(0.2)) +  scale_fill_manual(values=brewer.pal(5,name="PuOr")) +
  labs(title="\nPlate vs. RNA Integrity\n", x="\nPlate\n", y="\nRNA Integrity (RIN)\n") + theme_classic() +
  theme(legend.position="none", plot.title=element_text(hjust=0.5, size=14, face="bold")) +
  stat_compare_means(method="kruskal.test", label="p.format", label.x.npc=0.4, hide.ns=T)

integrity.concentration<-ggplot(dfm[dfm$Chain=="TRA",], aes(x=RNA.Integrity, y=RNA.Concentration)) + geom_point(color="dodgerblue4", size=0.6) +
  theme_bw()+labs(title="\nRNA Concentration vs. RNA Integrity\n", x="\nRNA Integrity (RIN)\n", y="\n\nRNA Concentration (ng/ul)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.8)

concentration.plate<-ggplot(dfm[dfm$Chain=="TRA",], aes(x = fct_reorder(Plate, RNA.Concentration, .fun=median, .desc=T), y = RNA.Concentration)) +  geom_boxplot(aes(fill = fct_reorder(Plate, RNA.Concentration, .fun = median, .desc =T))) +
  geom_jitter(position=position_jitter(0.2)) +  scale_fill_manual(values=brewer.pal(5,name="PuOr")) +
  labs(title="\nPlate vs. RNA Concentration\n", x="\nPlate\n", y="\nRNA Concentration (ng/ul)\n") + theme_classic() +
  theme(legend.position="none", plot.title=element_text(hjust=0.5, size=14, face="bold")) +
  stat_compare_means(method="kruskal.test", label="p.format", label.x.npc=0.4, hide.ns=T)

pdf("Output_Data/TechVar_Correlation.pdf", width=11, height=12)
plot(seqdepth.concentration)
plot(seqdepth.integrity)
plot(seqdepth.plate)
plot(concentration.plate)
plot(integrity.concentration)
plot(integrity.plate)
dev.off()

# Diversity

df.div.chain<-readRDS("Data2analyze/DivBaselineControls.rds")
df.div.chain.long<-melt(df.div.chain, id.vars=colnames(df.div.chain)[1:14], variable.name="Diversity.Measure")

pdf("Output_Data/Diversity_TechVar.pdf", width=11, height=12)
divs<-c("D50","D20","Gini.Simp.Logit","Inv.Simp","Shannon")
for (d in divs) {
  data2analyze<-df.div.chain.long[df.div.chain.long$Diversity.Measure==d,]
  div.reads<-ggplot(data2analyze, aes(x=value, y=Reads.Chain)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
    theme_bw()+labs(title=paste0("\n",d," vs. Sequencing Depth\n"), x=paste0("\n",d,"\n"), y="\n\nSequencing Depth (N reads)\n") +
    theme(plot.title=element_text(hjust=0.5, face="bold")) + 
    stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.6, label.y.npc=0.7)
  plot(div.reads)
}
vars<-c("RNA.Integrity","RNA.Concentration")
divs<-c("D50","D20","Gini.Simp.Logit","Inv.Simp","Shannon")
for (v in vars) {
  for (d in divs) {
    data2analyze<-df.div.chain.long[(df.div.chain.long$Diversity.Measure==d),]
    data2analyze$value.norm<-data2analyze$value/data2analyze$Reads.Chain
    data2analyze<-merge(data2analyze,tech.var,by="Sample")
    colnames(data2analyze)[which(colnames(data2analyze)==v)]<-"Var"
    data2analyze$Var<-as.numeric(data2analyze$Var)
    div.adj<-ggplot(data2analyze, aes(x=value.norm, y=Var)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
      theme_bw()+labs(title=paste0("\n",d," vs. ",v,"\n"), x=paste0("\n",d," (Adj. by sequencing depth)\n"), y=paste0("\n\n",v,"\n")) +
      theme(plot.title=element_text(hjust=0.5, face="bold")) + 
      stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.6, label.y.npc=0.7)
    plot(div.adj)
  }
}
dev.off()

pdf("Output_Data/Diversity_EpiVar.pdf", width=16, height=9)
divs<-c("D50","D20","Gini.Simp.Logit","Inv.Simp","Shannon")
for (d in divs) {
  data2analyze<-df.div.chain.long[df.div.chain.long$Diversity.Measure==d,]
  data2analyze$value.norm<-data2analyze$value/data2analyze$Reads.Chain
  div.gender<-ggboxplot(data2analyze, x="Chain", y="value.norm", color="Gender", add="jitter", fill="Gender") +
    labs(title=paste0("\n",d," vs. Gender\n"), y=paste0("\n",d," (Adj. by sequencing depth)\n"), x=paste0("\n")) +
    scale_color_manual(values=c("black","black")) + scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
    theme(legend.position="top", plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.key.size=unit(1.3, "cm")) +
    stat_compare_means(aes(group=Gender), method="kruskal.test", label="p.format", label.y.npc=1.00, hide.ns=T)
  div.age<-ggplot(data2analyze, aes(x=value.norm, y=Age)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
    theme_bw()+labs(title=paste0("\n",d," vs. Age\n"), x=paste0("\n",d," (Adj. by sequencing depth)\n"), y=paste0("\n\nAge\n")) +
    theme(plot.title=element_text(hjust=0.5, face="bold")) + 
    stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.5, label.y.npc=0.9)
  grid.arrange(div.gender,div.age, ncol=2)
}
dev.off()


# Unique Clone Number

data<-readRDS("Data2analyze/CloneDet.rds")
data<-data[(data$Sample %in% c(ra.wk0,ctrl)),]

ucn.reads<-ggplot(data, aes(x=Reads, y=Clones)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nUNIQUE CLONES vs. Sequecing depth\n", x="\nSequencing Depth (N, reads)\n", y="\n\nUnique Clones (N)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.1)

data$Clones.Norm<-data$Clone/data$Reads
dfm<-merge(data,tech.var,by="Sample")

ucn.integrity.adj<-ggplot(dfm, aes(x=RNA.Integrity, y=Clones.Norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nUNIQUE CLONES vs. RNA integrity\n", x="\nRNA Integrity (RIN)\n", y="\n\nUnique Clones (Adj. by sequencing depth)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.1)

ucn.concentration.adj<-ggplot(dfm, aes(x=RNA.Concentration, y=Clones.Norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nUNIQUE CLONES vs. RNA concentration\n", x="\nRNA Concentration (ug/ml)\n", y="\n\nUnique Clones (Adj. by sequencing depth)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.1)

pdf("Output_Data/UCN_TechVar.pdf", width=15, height=13)
plot(ucn.reads)
plot(ucn.integrity.adj)
plot(ucn.concentration.adj)
dev.off()

ucn.age.adj<-ggplot(data, aes(x=Age, y=Clones.Norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
             theme_bw()+labs(title="\nUNIQUE CLONES vs. Age\n", x="\nAge (years)\n", y="\n\nUnique Clones (Adj. by sequencing depth)\n") +
             theme(plot.title=element_text(hjust=0.5, face="bold")) + 
             stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.1)

ucn.gender.adj<-ggboxplot(data, x="Chain", y="Clones.Norm", color="Gender", add="jitter", fill="Gender") +
                labs(title="\nUNIQUE CLONES vs. Gender\n", x="", y="\nUnique Clones (Adj. by sequencing depth)\n") + scale_color_manual(values=c("black","black")) +
                scale_fill_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position="right", plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.key.size=unit(1.3, "cm")) +
                stat_compare_means(aes(group=Gender), method="wilcox.test", label="p.format", label.y=max(data$Clones.Norm)*1.1, hide.ns=T)

pdf("Output_Data/UCN_EpiVar.pdf", width=16, height=9)
grid.arrange(ucn.gender.adj,ucn.age.adj, ncol=2)
dev.off()

# Chain Usage

reads.sample.prop<-readRDS("Data2analyze/ChainReadsProp.rds")
reads.sample.prop<-reads.sample.prop[reads.sample.prop$Sample %in% c(ra.wk0,ctrl),]
df.prop.long<-melt(reads.sample.prop);  colnames(df.prop.long)[2]<-"Chain"
metadata<-readRDS("Data2analyze/Metadata_IGH.rds")
df.prop.long.pheno<-merge(metadata, df.prop.long, by="Sample")

reads.sample<-readRDS("Data2analyze/ChainReadsCount.rds")
rs.sum<-reads.sample[,-1]
rs.sum$Reads.All<-rowSums(rs.sum)
rs.sum$Sample<-rownames(rs.sum)
colnames(rs.sum)<-c("Reads.TRA", "Reads.TRB", "Reads.TRG", "Reads.TRD", "Reads.IGH", "Reads.IGL", "Reads.IGK", "Reads.All", "Sample")

df.sums.long1<-merge(df.prop.long.pheno, rs.sum, by="Sample")
df.sums.long<-merge(df.sums.long1,tech.var,by="Sample")

cu.reads<-ggplot(df.sums.long, aes(x=Reads.All, y=value)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nCHAIN USAGE vs. Sequencing Depth\n", x="\nSequencing Depth (N, reads)\n", y="\n\nExpression (% reads)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.9)

df.sums.long$value.norm<-df.sums.long$value/df.sums.long$Reads.All

cu.integrity.adj<-ggplot(df.sums.long, aes(x=RNA.Integrity, y=value.norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nCHAIN USAGE vs. RNA Integrity\n", x="\nRNA Integrity (RIN)\n", y="\n\nExpression (% reads, Adj. by Seq Depth)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

cu.concentration.adj<-ggplot(df.sums.long, aes(x=RNA.Concentration, y=value.norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nCHAIN USAGE vs. RNA Concentration\n", x="\nRNA Concentration (ng/ul)\n", y="\n\nExpression (% reads, Adj. by Seq Depth)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

pdf("Output_Data/CU_TechVar.pdf", width=15, height=13)
plot(cu.reads)
plot(cu.integrity.adj)
plot(cu.concentration.adj)
dev.off()

cu.gender.adj<-ggboxplot(df.sums.long, x="Chain", y="value.norm", color="Gender", add="jitter", fill="Gender") +
  labs(title="\nCHAIN USAGE vs. Gender\n", x="", y="\nExpression (% reads, Adj. by Seq Depth)\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position="right", plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.key.size=unit(1.3, "cm")) +
  stat_compare_means(aes(group=Gender), method="wilcox.test", label="p.format", hide.ns=T)

cu.age.adj<-ggplot(df.sums.long, aes(x=Age, y=value.norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~Chain, scales="free") +
  theme_bw()+labs(title="\nCHAIN USAGE vs. Age\n", x="\nAge (years)\n", y="\n\nExpression (% reads, Adj. by Seq Depth)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

pdf("Output_Data/CU_EpiVar.pdf", width=16, height=9)
grid.arrange(cu.gender.adj,cu.age.adj, ncol=2)
dev.off()

# BCR-related variables: isotype percentage

df.is.perc<-readRDS("Data2analyze/IsotypePerc.rds")
m<-data.frame("IgG"=df.is.perc[df.is.perc$Isotype=="IgG","Perc"],
              "IgE"=df.is.perc[df.is.perc$Isotype=="IgE","Perc"],
              "IgA"=df.is.perc[df.is.perc$Isotype=="IgA","Perc"],
              "IgD"=df.is.perc[df.is.perc$Isotype=="IgD","Perc"],
              "IgM"=df.is.perc[df.is.perc$Isotype=="IgM","Perc"],
              "Reads.All"=df.is.perc[df.is.perc$Isotype=="IgM","Reads.All"],
              "Reads.Chain"=df.is.perc[df.is.perc$Isotype=="IgM","Reads.Chain"],
              "Age"=df.is.perc[df.is.perc$Isotype=="IgM","Age"],
              "Gender"=as.numeric(df.is.perc[df.is.perc$Isotype=="IgM","Gender"]),
              "IMID"=as.numeric(df.is.perc[df.is.perc$Isotype=="IgM","IMID"]),
              "Sample"=as.character(df.is.perc[df.is.perc$Isotype=="IgM","Sample"]))

dfm<-merge(m,tech.var,by="Sample")
dfm.long<-gather(dfm, condition, measurement, IgG:IgM, factor_key=T)

iso.reads<-ggplot(dfm.long, aes(x=Reads.Chain, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nISOTYPE PERCENTAGE vs. Sequencing Depth\n", x="\nSequencing Depth (N, reads)\n", y="\n\nIsotype (%)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.9)

iso.integrity.unadj<-ggplot(dfm.long, aes(x=RNA.Integrity, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nISOTYPE PERCENTAGE vs. RNA Integrity\n", x="\nAge (years)\n", y="\n\nRNA Integrity (RIN)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

iso.concentration.unadj<-ggplot(dfm.long, aes(x=RNA.Concentration, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nISOTYPE PERCENTAGE vs. RNA Concentration\n", x="\nAge (years)\n", y="\n\nRNA Concentration (ng/ul)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

pdf("Output_Data/BCRisoperc_TechVar.pdf", width=13, height=8)
plot(iso.reads)
plot(iso.integrity.unadj)
plot(iso.concentration.unadj)
dev.off()

iso.age.unadj<-ggplot(dfm.long, aes(x=Age, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nISOTYPE PERCENTAGE vs. Age\n", x="\nAge (years)\n", y="\n\nIsotype (%)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

dfm.long$Gender<-as.character(dfm.long$Gender)
dfm.long$Gender[dfm.long$Gender=="1"]<-"FEMALE"
dfm.long$Gender[dfm.long$Gender=="2"]<-"MALE"
iso.gender.unadj<-ggboxplot(dfm.long, x="condition", y="measurement", color="Gender", add="jitter", fill="Gender") +
  labs(title="\nISOTYPE PERCENTAGE vs. Gender\n", x="", y="\nIsotype (% reads)\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position="right", plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.key.size=unit(1.3, "cm")) +
  stat_compare_means(aes(group=Gender), method="wilcox.test", label="p.format", hide.ns=T)

pdf("Output_Data/BCRisoperc_EpiVar.pdf", width=16, height=7)
grid.arrange(iso.gender.unadj,iso.age.unadj, ncol=2)
dev.off()

# BCR-related variables: percentage mutated IgD/M

df.IgDMmut.perc<-readRDS("Data2analyze/IgDMmutPerc.rds")
m<-df.IgDMmut.perc[,c(1,3,4,7,8,9)]
dfm<-merge(m,tech.var,by="Sample")

isomut.reads<-ggplot(dfm, aes(x=Reads.Chain, y=Perc.IgDM.Mut)) + geom_point(color="dodgerblue4", size=0.6) +
  theme_bw()+labs(title="\nPERCENTAGE IgD/M MUTATED vs. Sequencing Depth\n", x="\nSequencing Depth (N, reads)\n", y="\n\nIgD/M mutated (%)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.9)

isomut.integrity.unadj<-ggplot(dfm, aes(x=RNA.Integrity, y=Perc.IgDM.Mut)) + geom_point(color="dodgerblue4", size=0.6) +
  theme_bw()+labs(title="\nPERCENTAGE IgD/M MUTATED vs. RNA Integrity\n", x="\nRNA Integrity (RIN)\n", y="\n\nIgD/M mutated (%)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.9)

isomut.concentration.unadj<-ggplot(dfm, aes(x=RNA.Concentration, y=Perc.IgDM.Mut)) + geom_point(color="dodgerblue4", size=0.6) +
  theme_bw()+labs(title="\nPERCENTAGE IgD/M MUTATED vs. RNA Concentration\n", x="\nRNA Concentration (ng/ul)\n", y="\n\nIgD/M mutated (%)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.9)

pdf("Output_Data/BCRmut_TechVar.pdf", width=7, height=7)
plot(isomut.reads)
plot(isomut.integrity.unadj)
plot(isomut.concentration.unadj)
dev.off()

isomut.gender.unadj<-ggboxplot(dfm, x="Gender", y="Perc.IgDM.Mut", color="Gender", add="jitter", fill="Gender") +
  labs(title="\nPERCENTAGE IgD/M MUTATED vs. Gender\n", x="", y="\nIgD/M mutated (%)\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position="right", plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.key.size=unit(1.3, "cm")) +
  stat_compare_means(aes(group=Gender), method="wilcox.test", label="p.format", hide.ns=T)

isomut.age.unadj<-ggplot(dfm, aes(x=Age, y=Perc.IgDM.Mut)) + geom_point(color="dodgerblue4", size=0.6) +
  theme_bw()+labs(title="\nPERCENTAGE IgD/M MUTATED vs. AGE\n", x="\nAge (years)\n", y="\n\nIgD/M mutated (%)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.9)

pdf("Output_Data/BCRmut_EpiVar.pdf", width=16, height=7)
grid.arrange(isomut.gender.unadj, isomut.age.unadj, ncol=2)
dev.off()

# BCR-related variables: CSR

df.csr<-readRDS("Data2analyze/CSRindex.rds")
m<-data.frame("IgE-IgG"=df.csr[df.csr$CSR.Name=="IgE-IgG","CSR.Index"],
              "IgDM-IgG"=df.csr[df.csr$CSR.Name=="IgDM-IgG","CSR.Index"],
              "IgDM-IgE"=df.csr[df.csr$CSR.Name=="IgDM-IgE","CSR.Index"],
              "IgA-IgDM"=df.csr[df.csr$CSR.Name=="IgA-IgDM","CSR.Index"],
              "IgA-IgG"=df.csr[df.csr$CSR.Name=="IgA-IgG","CSR.Index"],
              "IgA-IgE"=df.csr[df.csr$CSR.Name=="IgA-IgE","CSR.Index"],
              "Reads.All"=df.csr[df.csr$CSR.Name=="IgA-IgE","Reads.All"],
              "Reads.Chain"=df.csr[df.csr$CSR.Name=="IgA-IgE","Reads.Chain"],
              "Age"=df.csr[df.csr$CSR.Name=="IgA-IgE","Age"],
              "Gender"=df.csr[df.csr$CSR.Name=="IgA-IgE","Gender"],
              "IMID"=df.csr[df.csr$CSR.Name=="IgA-IgE","IMID"],
              "Sample"=df.csr[df.csr$CSR.Name=="IgA-IgE","Sample"])

dfm<-merge(m,tech.var,by="Sample")
dfm.long<-gather(dfm, condition, measurement, IgE.IgG:IgA.IgE, factor_key=T)
dfm.long$value.norm<-dfm.long$measurement/dfm.long$Reads.Chain

csr.reads<-ggplot(dfm.long, aes(x=Reads.Chain, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nCLASS SWITCH RECOMBINATION vs. Sequencing Depth\n", x="\nSequencing Depth (N, reads)\n", y="\n\nClass Switch Recombination Index\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.9)

csr.integrity.adj<-ggplot(dfm.long, aes(x=RNA.Integrity, y=value.norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nCLASS SWITCH RECOMBINATION  vs. RNA Integrity\n", x="\nAge (years)\n", y="\n\nClass Switch Recombination Index (Adj. by Seq Depth)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

csr.concentration.adj<-ggplot(dfm.long, aes(x=RNA.Concentration, y=value.norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nCLASS SWITCH RECOMBINATION vs. RNA Concentration\n", x="\nAge (years)\n", y="\n\nClass Switch Recombination Index (Adj. by Seq. Depth)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

pdf("Output_Data/BCRcsr_TechVar.pdf", width=10, height=7)
plot(csr.reads)
plot(csr.integrity.adj)
plot(csr.concentration.adj)
dev.off()

csr.age.adj<-ggplot(dfm.long, aes(x=Age, y=value.norm)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nCLASS SWITCH RECOMBINATION vs. Age\n", x="\nAge (years)\n", y="\n\nClass Switch Recombination Index (Adj. by Seq Depth)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

dfm.long$Gender<-as.character(dfm.long$Gender)
csr.gender.adj<-ggboxplot(dfm.long, x="condition", y="value.norm", color="Gender", add="jitter", fill="Gender") +
  labs(title="\nCLASS SWITCH RECOMBINATION vs. Gender\n", x="", y="\nClass Switch Recombination Index (Adj. by Seq Depth)\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position="right", plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.key.size=unit(1.3, "cm")) +
  stat_compare_means(aes(group=Gender), method="wilcox.test", label="p.format", hide.ns=T)

pdf("Output_Data/BCRcsr_EpiVar.pdf", width=16, height=7)
grid.arrange(csr.gender.adj, csr.age.adj, ncol=2)
dev.off()

# BCR-related variables: Somatic Hypermutations
df.shm<-readRDS("Data2analyze/SHM.rds")
df.shm.igh<-df.shm[df.shm$Chain=="IGH",]
m.igh<-data.frame("IgG34"=df.shm.igh[df.shm.igh$Isotype=="IgG34","SHM.Perc"],
                  "IgG12"=df.shm.igh[df.shm.igh$Isotype=="IgG12","SHM.Perc"],
                  "IgE"=df.shm.igh[df.shm.igh$Isotype=="IgE","SHM.Perc"],
                  "IgA"=df.shm.igh[df.shm.igh$Isotype=="IgA","SHM.Perc"],
                  "IgD"=df.shm.igh[df.shm.igh$Isotype=="IgD","SHM.Perc"],
                  "IgM"=df.shm.igh[df.shm.igh$Isotype=="IgM","SHM.Perc"],
                  "Reads.All"=df.shm.igh[df.shm.igh$Isotype=="IgM","Reads.All"],
                  "Reads.Chain"=df.shm.igh[df.shm.igh$Isotype=="IgM","Reads.Chain"],
                  "Age"=df.shm.igh[df.shm.igh$Isotype=="IgM","Age"],
                  "Gender"=df.shm.igh[df.shm.igh$Isotype=="IgM","Gender"],
                  "IMID"=df.shm.igh[df.shm.igh$Isotype=="IgM","IMID"],
                  "Sample"=df.shm.igh[df.shm.igh$Isotype=="IgM","Sample"])

df.shm.igk<-df.shm[df.shm$Chain=="IGK",]
m.igk<-data.frame("IGK"=df.shm.igk[df.shm.igk$Isotype=="IGK","SHM.Perc"],
                  "Reads.All"=df.shm.igk[df.shm.igk$Isotype=="IGK","Reads.All"],
                  "Reads.Chain"=df.shm.igk[df.shm.igk$Isotype=="IGK","Reads.Chain"],
                  "Age"=df.shm.igk[df.shm.igk$Isotype=="IGK","Age"],
                  "Gender"=df.shm.igk[df.shm.igk$Isotype=="IGK","Gender"],
                  "IMID"=df.shm.igk[df.shm.igk$Isotype=="IGK","IMID"],
                  "Sample"=df.shm.igk[df.shm.igk$Isotype=="IGK","Sample"])

df.shm.igl<-df.shm[df.shm$Chain=="IGL",]
m.igl<-data.frame("IGL"=df.shm.igl[df.shm.igl$Isotype=="IGL","SHM.Perc"],
                  "Reads.All"=df.shm.igl[df.shm.igl$Isotype=="IGL","Reads.All"],
                  "Reads.Chain"=df.shm.igl[df.shm.igl$Isotype=="IGL","Reads.Chain"],
                  "Age"=df.shm.igl[df.shm.igl$Isotype=="IGL","Age"],
                  "Gender"=df.shm.igl[df.shm.igl$Isotype=="IGL","Gender"],
                  "IMID"=df.shm.igl[df.shm.igl$Isotype=="IGL","IMID"],
                  "Sample"=df.shm.igl[df.shm.igl$Isotype=="IGL","Sample"])

m.igl2<-m.igl[,c("Reads.All","Reads.Chain","Age","Gender","IMID","Sample","IGL")]
colnames(m.igl2)[7]<-"measurement"
m.igl2$condition<-"IGL"
m.igl2<-m.igl2[,c(1,2,3,4,5,6,8,7)]

m.igk2<-m.igk[,c("Reads.All","Reads.Chain","Age","Gender","IMID","Sample","IGK")]
colnames(m.igk2)[7]<-"measurement"
m.igk2$condition<-"IGK"
m.igk2<-m.igk2[,c(1,2,3,4,5,6,8,7)]

m.igh2<-m.igh[,c("Reads.All","Reads.Chain","Age","Gender","IMID","Sample","IgG34","IgG12","IgE","IgA","IgD","IgM")]
m.igh2<-gather(m.igh2, condition, measurement, IgG34:IgM, factor_key=T)


m<-rbind(m.igh2,m.igl2,m.igk2)

dfm<-merge(m,tech.var,by="Sample")
dfm.long<-dfm
dfm.long$value.norm<-dfm.long$measurement/dfm.long$Reads.Chain

shm.reads<-ggplot(dfm.long, aes(x=Reads.Chain, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nSOMATIC HYPERMUTATIONS vs. Sequencing Depth\n", x="\nSequencing Depth (N, reads)\n", y="\n\nSomatic Hypermutation Index\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) +
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.4, label.y.npc=0.9)

shm.integrity.unadj<-ggplot(dfm.long, aes(x=RNA.Integrity, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nSOMATIC HYPERMUTATIONS  vs. RNA Integrity\n", x="\nAge (years)\n", y="\n\nSomatic Hypermutation\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

shm.concentration.unadj<-ggplot(dfm.long, aes(x=RNA.Concentration, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nSOMATIC HYPERMUTATIONS vs. RNA Concentration\n", x="\nAge (years)\n", y="\n\nSomatic Hypermutation Index\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

pdf("Output_Data/BCRshm_TechVar.pdf", width=10, height=10)
plot(shm.reads)
plot(shm.integrity.unadj)
plot(shm.concentration.unadj)
dev.off()

dfm.long$Gender<-as.character(dfm.long$Gender)
shm.gender.unadj<-ggboxplot(dfm.long, x="condition", y="measurement", color="Gender", add="jitter", fill="Gender") +
  labs(title="\nSOMATIC HYPERMUTATIONS vs. Gender\n", x="", y="\nSomatic Hypermutation Index\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position="right", plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.key.size=unit(1.3, "cm")) +
  stat_compare_means(aes(group=Gender), method="wilcox.test", label="p.format", hide.ns=T)

shm.age.unadj<-ggplot(dfm.long, aes(x=Age, y=measurement)) + geom_point(color="dodgerblue4", size=0.6) + facet_wrap(~condition, scales="free") +
  theme_bw()+labs(title="\nSOMATIC HYPERMUTATIONS vs. Age\n", x="\nAge (years)\n", y="\n\nSomatic Hypermutation Index\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold")) + 
  stat_cor(aes(label = paste(..rr.label.., sep="~`,`~")),  label.x.npc=0.05, label.y.npc=0.85)

pdf("Output_Data/BCRshm_EpiVar.pdf", width=20, height=10)
grid.arrange(shm.gender.unadj, shm.age.unadj, ncol=2)
dev.off()
