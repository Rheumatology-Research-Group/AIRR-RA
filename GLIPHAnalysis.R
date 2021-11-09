# Loading libraries

suppressPackageStartupMessages({
  library(MAST);  library(immunarch);  library(plyr); library(data.table);  library(reshape);  library(scater); library(ComplexHeatmap);
  library(edgeR);  library(ROCR);  library(scRNA.seq.funcs);  library(gridExtra);  library(NMF);  library(RColorBrewer);  library(pheatmap);
  library(stringdist); library(Peptides);  library(ggpubr);   library(igraph);  library(stringdist);
  library(tcR);  library(purrr);  library(tidyr);  library(dendextend);  library(doParallel);  library(foreach);  library(motifStack)
  library(protr);  library(msa);  library(Biostrings);  library(ggseqlogo);  library(HDMD);  library(viridis);   library(readr)
})

# Creating GLIPH2 Data Format

dir.create("Data2analyze/GLIPH2")

tcr<-readRDS("Data2analyze/Clones_By_CDR3aa/TRB.rds")

# Keeping best V segment
samples<-names(tcr$data)
data<-list()
counter=0
for (s in samples) {
  s.data<-tcr$data[[s]]
  cdr3.b<-s.data$CDR3.aa
  trb.vbest<-unlist(map(strsplit(s.data$V.name,","),1))
  trb.jbest<-unlist(map(strsplit(s.data$J.name,","),1))
  cdr3.a<-rep("NA",length(cdr3.b))
  pheno<-tcr$meta[tcr$meta$Sample==s,"IMID"][[1]]
  ind<-rep(paste0(s,":",pheno),length(cdr3.b))
  count<-s.data$Clones
  df.out<-data.frame("CDR3b"=cdr3.b, "TRBV"=trb.vbest, "TRBJ"=trb.jbest, "CDR3a"=cdr3.a, "Pheno"=ind, "Count"=count)
  counter=counter+1
  print(counter)
  data[[counter]]<-df.out
}
names(data)<-samples
tcrb<-do.call(rbind,data)
saveRDS(tcrb, "Data2analyze/GLIPH2/Data_TCRB.rds")

# Filtering individuals of interest

hla<-readRDS("Data2analyze/GenoHLA.rds")

## RA week0
pheno.name="IMID";  pheno2pred="RA";   week2keep="wk0";   week2discard="wk12";
pheno.info<-read.csv("Data2analyze/Metadata.csv", sep="\t")
pheno.info2<-pheno.info[(pheno.info$Week!=week2discard), c("Sample",pheno.name)];  colnames(pheno.info2)<-c("Sample","Pheno")
#df.pheno<-pheno.info2[(pheno.info2$Pheno==pheno2pred),]
codes<-as.character(pheno.info2[!(is.na(pheno.info2$Pheno)),"Sample"])

samples2keep<-codes[(codes %in% as.character(unlist(map(strsplit(hla,"---"),1))))]

samples<-unlist(map(strsplit(hla,"---"),1))
genos<-unlist(map(strsplit(hla,"---"),2))
df.hla<-data.frame("Samples"=samples,"Genotypes"=genos)
df.hla.filt<-df.hla[(as.character(df.hla$Samples) %in% samples2keep),]
df.hla.filt$Genotypes<-gsub("xxx","\t",df.hla.filt$Genotypes)
colnames(df.hla.filt)[1]<-"Sample"
df.hla.merged<-merge(df.hla.filt,pheno.info2, by="Sample")
hla2analyze<-paste(paste0(df.hla.merged$Sample,":",df.hla.merged$Pheno),df.hla.merged$Genotypes,sep="\t")

tcrb$Pheno2filt<-as.character(unlist(map(strsplit(as.character(tcrb$Pheno),":"),1)))
tcrb.filt<-tcrb[(tcrb$Pheno2filt %in% samples2keep),]
tcrb2analyze<-tcrb.filt[,c(1,2,3,4,5,6)]

write.table(tcrb2analyze, "Data2analyze/GLIPH2/RAvsCTRL_TCRB.txt", row.names=F, col.names=F, quote=F, sep="\t")
write.table(hla2analyze, "Data2analyze/GLIPH2/RAvsCTRL_HLA.txt", quote=F, row.names=F, col.names=F, sep="\t")
  

# Part1: Summarizing results
# Files that are processed by this script can be accessed through the "Download" tab at the Rheumatology Research Group website

# RA

f1<-read.table("Data2analyze/RAweek0_Output_CD4_cluster.csv", sep=",", header=T)
f1$FDR<-p.adjust(f1$Fisher_score, metho="fdr")
df2save.cd4<-f1[,c(2,3,30,4,5,13,14,15)]
colnames(df2save.cd4)<-c("AminoAcid.Pattern","Enrichment.Pval","Enrichment.FDR","Individuals.N","UniqueClones.N","Clones","V.GeneSegment","J.GeneSegment")
write.csv(df2save.cd4,"Output_Data/GLIPH_RAcd4.csv", row.names=F)
f1<-read.table("Data2analyze/RAweek0_Output_CD8_cluster.csv", sep=",", header=T)
f1$FDR<-p.adjust(f1$Fisher_score, metho="fdr")
df2save.cd8<-f1[,c(2,3,30,4,5,13,14,15)]
colnames(df2save.cd8)<-c("AminoAcid.Pattern","Enrichment.Pval","Enrichment.FDR","Individuals.N","UniqueClones.N","Clones","V.GeneSegment","J.GeneSegment")
write.csv(df2save.cd8,"Output_Data/GLIPH_RAcd8.csv", row.names=F)
ra.cd4<-df2save.cd4[df2save.cd4$Enrichment.FDR<0.05 & df2save.cd4$Individuals.N>=3 & df2save.cd4$UniqueClones.N>=3,]
ra.cd8<-df2save.cd8[df2save.cd8$Enrichment.FDR<0.05 & df2save.cd8$Individuals.N>=3 & df2save.cd8$UniqueClones.N>=3,]
ra.enrichCD4<-unique(ra.cd4$AminoAcid.Pattern)
ra.enrichCD8<-unique(ra.cd8$AminoAcid.Pattern)

length(unique(ra.cd4$AminoAcid.Pattern))
length(unique(ra.cd8$AminoAcid.Pattern))

cd4cd8<-unique(as.character(ra.cd4$AminoAcid.Pattern)[as.character(ra.cd4$AminoAcid.Pattern) %in% as.character(ra.cd8$AminoAcid.Pattern)])
cd8cd4<-unique(as.character(ra.cd8$AminoAcid.Pattern)[as.character(ra.cd8$AminoAcid.Pattern) %in% as.character(ra.cd4$AminoAcid.Pattern)])

ra.enrich<-unique(c(ra.cd4$AminoAcid.Pattern,ra.cd8$AminoAcid.Pattern))

# CTRL

f1<-read.table("Data2analyze/CTRL_Output_CD4_cluster.csv", sep=",", header=T)
f1$FDR<-p.adjust(f1$Fisher_score, metho="fdr")
df2save.cd4<-f1[,c(2,3,30,4,5,13,14,15)]
colnames(df2save.cd4)<-c("AminoAcid.Pattern","Enrichment.Pval","Enrichment.FDR","Individuals.N","UniqueClones.N","Clones","V.GeneSegment","J.GeneSegment")
write.csv(df2save.cd4,"Output_Data/GLIPH_CTRLcd4.csv", row.names=F)
f1<-read.table("Data2analyze/CTRL_Output_CD8_cluster.csv", sep=",", header=T)
f1$FDR<-p.adjust(f1$Fisher_score, metho="fdr")
df2save.cd8<-f1[,c(2,3,30,4,5,13,14,15)]
colnames(df2save.cd8)<-c("AminoAcid.Pattern","Enrichment.Pval","Enrichment.FDR","Individuals.N","UniqueClones.N","Clones","V.GeneSegment","J.GeneSegment")
write.csv(df2save.cd8,"Output_Data/GLIPH_CTRLcd8.csv", row.names=F)
ctrl.cd4<-df2save.cd4[df2save.cd4$Enrichment.FDR<0.05 & df2save.cd4$Individuals.N>=3 & df2save.cd4$UniqueClones.N>=3,]
ctrl.cd8<-df2save.cd8[df2save.cd8$Enrichment.FDR<0.05 & df2save.cd8$Individuals.N>=3 & df2save.cd8$UniqueClones.N>=3,]
ctrl.enrichCD4<-unique(ctrl.cd4$AminoAcid.Pattern)
ctrl.enrichCD8<-unique(ctrl.cd8$AminoAcid.Pattern)

# RA-specific

ra.specific.cd4<-as.character(ra.enrichCD4[!(ra.enrichCD4 %in% ctrl.enrichCD4)])
ra.specific.cd8<-as.character(ra.enrichCD8[!(ra.enrichCD8 %in% ctrl.enrichCD8)])
ra.inctrl.cd4<-as.character(ra.enrichCD4[(ra.enrichCD4 %in% ctrl.enrichCD4)])
ra.inctrl.cd8<-as.character(ra.enrichCD8[(ra.enrichCD8 %in% ctrl.enrichCD8)])
ra.cd4.filt<-ra.cd4[(as.character(ra.cd4$AminoAcid.Pattern) %in% ra.specific.cd4),]
ra.cd8.filt<-ra.cd8[(as.character(ra.cd8$AminoAcid.Pattern) %in% ra.specific.cd8),]
ra.cd4.filt$Cell.Type<-rep("CD4.Tcell",nrow(ra.cd4.filt))
ra.cd8.filt$Cell.Type<-rep("CD8.Tcell",nrow(ra.cd8.filt))

# Getting HLA data

ra.cd4.data<-read.table("Data2analyze/RAweek0_Output_CD4_cluster.csv", sep=",", header=T)
ra.cd8.data<-read.table("Data2analyze/RAweek0_Output_CD8_cluster.csv", sep=",", header=T)
ra.cd4.hla<-unique(ra.cd4.data[,c("pattern","hla_score")])
ra.cd8.hla<-unique(ra.cd8.data[,c("pattern","hla_score")])

colnames(ra.cd4.hla)<-c("AminoAcid.Pattern","HLA.Pval")
colnames(ra.cd8.hla)<-c("AminoAcid.Pattern","HLA.Pval")

ra.cd4.filt.hla<-merge(ra.cd4.filt,ra.cd4.hla,by="AminoAcid.Pattern")
ra.cd8.filt.hla<-merge(ra.cd8.filt,ra.cd8.hla,by="AminoAcid.Pattern")
ra.merged<-unique(rbind(ra.cd4.filt.hla,ra.cd8.filt.hla)[,c(1,4,5,9,2,3,10)])
ra.merged<-ra.merged[order(ra.merged$Enrichment.Pval),]
ra.merged.hla.filt<-ra.merged[ra.merged$HLA.Pval<0.05,]
ra.merged.hla.filt$Pval.Comb<-(-log10(ra.merged.hla.filt$Enrichment.Pval)+(-log10(ra.merged.hla.filt$HLA.Pval)))
data.final<-ra.merged.hla.filt[order(ra.merged.hla.filt$Pval.Comb, decreasing=T),]
data.final$HLA.Allele<-"NA"
data.final$HLA.Allele<-c("DRB1*11", "DRB1*15", "DQB1*06", "DQB1*03", "DQA1*03, DRB1-04", "DQB1*02", "DQB1*03", "DQA1*03, DRB1-04", "DQB1*03", "DQA1*03, DRB1-04",
                         "DRB1*11","DRB1*13","DQA1*03, DRB1-04","DQA1*03, DRB1-04","DQA1*05","DRB1*11","DQA1*05","DQA1*03, DRB1-04","DRB1*11","DQA1*03, DRB1-04",
                         "DRB1*11","DQA1*05","DQA1*05","DRB1*13","DQA1*03, DRB1-04","DRB1*07","DRB1*15","DRB1*15","DQA1*05","DRB1*11",
                         "DQA1*05","DRB1*15","DRB1*15","DRB1*15","DQB1*06","DQB1*03","DQB1*03","DRB1*03","DQB1*03","DQA1*03",
                         "DRB1*11","DRB1*11","DQB1*03","DRB1*11","DQA1*05","DRB1*11","DQB1*02","DQB1*03","DRB1*03","DQB1*03",
                         "DQA1*05","DQA1*05","DQB1*03","DRB1*15","DQA1*02","DRB1*11","DQB1*06","DQB1*06","DRB1*13","DQB1*02",
                         "DQB1*02","DQA1*03, DRB1-04, DRB1*11","DRB1*15","DQB1*02","DRB1*11")

data.final<-data.final[order(data.final$HLA.Pval),-8]

write.table(data.final,"Output_Data/GLIPH_Main.csv", quote=F, row.names=F, col.names=T, sep=";")


# Part2: Detailed processing of GLIPH output
# Like in previous step, files that are processed by this script can be accessed through the "Download" tab at the Rheumatology Research Group website

# Is there any cluster that is significantly enriched in RA compared to controls? ####

# CD4
patterns<-system("awk -F',' '{print $2}' RAweek0_Output_CD4_cluster.csv", intern=T)
pvals<-system("awk -F',' '{print $3}' RAweek0_Output_CD4_cluster.csv", intern=T)
df.ra.cd4<-data.frame("Cluster.Pattern"=patterns, "Cluster.Enrichment.Pval"=pvals)
df.ra.cd4<-unique(df.ra.cd4[-1,])
df.ra.cd4<-df.ra.cd4[df.ra.cd4$Cluster.Pattern!="",]
df.ra.cd4$Cluster.Enrichment.Pval<-as.numeric(as.character(df.ra.cd4$Cluster.Enrichment.Pval))
df.ra.cd4$Cluster.Enrichment.FDR<-p.adjust(df.ra.cd4$Cluster.Enrichment.Pval, method="fdr")
nrow(df.ra.cd4[df.ra.cd4$Cluster.Enrichment.FDR<0.05,])
clusters.ra.cd4<-as.character(df.ra.cd4[df.ra.cd4$Cluster.Enrichment.FDR<0.05,"Cluster.Pattern"])
out1<-df.ra.cd4[df.ra.cd4$Cluster.Enrichment.FDR<0.05,]
saveRDS(out1,"Output_Data/Sig_RA_CD4.rds")

# CD8
patterns<-system("awk -F',' '{print $2}' RAweek0_Output_CD8_cluster.csv", intern=T)
pvals<-system("awk -F',' '{print $3}' RAweek0_Output_CD8_cluster.csv", intern=T)
df.ra.cd8<-data.frame("Cluster.Pattern"=patterns, "Cluster.Enrichment.Pval"=pvals)
df.ra.cd8<-unique(df.ra.cd8[-1,])
df.ra.cd8<-df.ra.cd8[df.ra.cd8$Cluster.Pattern!="",]
df.ra.cd8$Cluster.Enrichment.Pval<-as.numeric(as.character(df.ra.cd8$Cluster.Enrichment.Pval))
df.ra.cd8$Cluster.Enrichment.FDR<-p.adjust(df.ra.cd8$Cluster.Enrichment.Pval, method="fdr")
nrow(df.ra.cd8[df.ra.cd8$Cluster.Enrichment.FDR<0.05,])
clusters.ra.cd8<-as.character(df.ra.cd8[df.ra.cd8$Cluster.Enrichment.FDR<0.05,"Cluster.Pattern"])
out2<-df.ra.cd8[df.ra.cd8$Cluster.Enrichment.FDR<0.05,]
saveRDS(out2,"Output_Data/Sig_RA_CD8.rds")

# Common CD4-CD8 patterns
clusters.ra.cd8[clusters.ra.cd8 %in% clusters.ra.cd4]


# Which of these enriched clusters are non-detected in controls?

# CD4
patterns<-system("awk -F',' '{print $2}' CTRL_Output_CD4_cluster.csv", intern=T)
pvals<-system("awk -F',' '{print $3}' CTRL_Output_CD4_cluster.csv", intern=T)
df.ctrl.cd4<-data.frame("Cluster.Pattern"=patterns, "Cluster.Enrichment.Pval"=pvals)
df.ctrl.cd4<-unique(df.ctrl.cd4[-1,])
df.ctrl.cd4<-df.ctrl.cd4[df.ctrl.cd4$Cluster.Pattern!="",]
df.ctrl.cd4$Cluster.Enrichment.Pval<-as.numeric(as.character(df.ctrl.cd4$Cluster.Enrichment.Pval))
df.ctrl.cd4$Cluster.Enrichment.FDR<-p.adjust(df.ctrl.cd4$Cluster.Enrichment.Pval, method="fdr")
nrow(df.ctrl.cd4[df.ctrl.cd4$Cluster.Enrichment.FDR<0.05,])
clusters.ctrl.cd4<-as.character(df.ctrl.cd4[df.ctrl.cd4$Cluster.Enrichment.FDR<0.05,"Cluster.Pattern"])
clusters.ctrl.cd4.nominal<-as.character(df.ctrl.cd4[df.ctrl.cd4$Cluster.Enrichment.Pval<0.05,"Cluster.Pattern"])

ra.specific.cd4.clusters<-clusters.ra.cd4[!(clusters.ra.cd4 %in% clusters.ctrl.cd4.nominal)]
ractrl.shared.cd4.clusters<-clusters.ra.cd4[(clusters.ra.cd4 %in% clusters.ctrl.cd4.nominal)]

# CD8
patterns<-system("awk -F',' '{print $2}' CTRL_Output_CD8_cluster.csv", intern=T)
pvals<-system("awk -F',' '{print $3}' CTRL_Output_CD8_cluster.csv", intern=T)
df.ctrl.cd8<-data.frame("Cluster.Pattern"=patterns, "Cluster.Enrichment.Pval"=pvals)
df.ctrl.cd8<-unique(df.ctrl.cd8[-1,])
df.ctrl.cd8<-df.ctrl.cd8[df.ctrl.cd8$Cluster.Pattern!="",]
df.ctrl.cd8$Cluster.Enrichment.Pval<-as.numeric(as.character(df.ctrl.cd8$Cluster.Enrichment.Pval))
df.ctrl.cd8$Cluster.Enrichment.FDR<-p.adjust(df.ctrl.cd8$Cluster.Enrichment.Pval, method="fdr")
nrow(df.ctrl.cd8[df.ctrl.cd8$Cluster.Enrichment.FDR<0.05,])
clusters.ctrl.cd8<-as.character(df.ctrl.cd8[df.ctrl.cd8$Cluster.Enrichment.FDR<0.05,"Cluster.Pattern"])
clusters.ctrl.cd8.nominal<-as.character(df.ctrl.cd8[df.ctrl.cd8$Cluster.Enrichment.Pval<0.05,"Cluster.Pattern"])

ra.specific.cd8.clusters<-clusters.ra.cd8[!(clusters.ra.cd8 %in% clusters.ctrl.cd8.nominal)]
ractrl.shared.cd8.clusters<-clusters.ra.cd8[(clusters.ra.cd8 %in% clusters.ctrl.cd8.nominal)]


# Is any the RA-specific CD4- and CD8 clusters likely to be presented by a particular HLA?
ra.specific.cd4.clusters<-clusters.ra.cd4[!(clusters.ra.cd4 %in% clusters.ctrl.cd4.nominal)]
ra.specific.cd8.clusters<-clusters.ra.cd8[!(clusters.ra.cd8 %in% clusters.ctrl.cd8.nominal)]

ra.cd4.res<-read.csv("RAweek0_Output_CD4_cluster.csv", sep=",", head=T)
ra.cd8.res<-read.csv("RAweek0_Output_CD8_cluster.csv", sep=",", head=T)

ra.specific.cd4.clusters.hla<-unique(ra.cd4.res[(ra.cd4.res$pattern %in% ra.specific.cd4.clusters),c("pattern","hla_score")])
ra.specific.cd8.clusters.hla<-unique(ra.cd8.res[(ra.cd8.res$pattern %in% ra.specific.cd8.clusters),c("pattern","hla_score")])

ra.specific.cd4.clusters.hla<-ra.specific.cd4.clusters.hla[order(ra.specific.cd4.clusters.hla$hla_score),]
ra.specific.cd8.clusters.hla<-ra.specific.cd8.clusters.hla[order(ra.specific.cd8.clusters.hla$hla_score),]

# Manual add of HLA allele

ra.specific.cd4.clusters.hla$hla_allele<-rep("NA",nrow(ra.specific.cd4.clusters.hla))
ra.specific.cd4.clusters.hla$hla_allele[1]<-"DRB1*15"
ra.specific.cd4.clusters.hla$hla_allele[2]<-"DQB1*02"
ra.specific.cd4.clusters.hla$hla_allele[3]<-"DRB1*04 and DQA1*03"
ra.specific.cd4.clusters.hla$hla_allele[4]<-"DRB1*11"

ra.specific.cd8.clusters.hla$hla_allele<-rep("NA",nrow(ra.specific.cd8.clusters.hla))
ra.specific.cd8.clusters.hla$hla_allele[1]<-"DRB1*04 and DQA1*03"
ra.specific.cd8.clusters.hla$hla_allele[2]<-"DQB1*03"
ra.specific.cd8.clusters.hla$hla_allele[3]<-"DQB1*06"
ra.specific.cd8.clusters.hla$hla_allele[4]<-"DQB1*03"

saveRDS(ra.specific.cd4.clusters.hla,"Output_Data/Sig_HLA_RA_CD4.rds")
saveRDS(ra.specific.cd8.clusters.hla,"Output_Data/Sig_HLA_RA_CD8.rds")

ractrl.data<-readRDS("Output_Data/Kmer_IMID_TRB_statistics.rds") # This file must be created by running the KmerAnalysis.R script
RespGM.data<-readRDS("Data_Output/Kmer_RespGM_TRB_statistics.rds") # This file must be created by running the KmerAnalysis.R script
ActB.data<-readRDS("Data_Output/Kmer_ActB_TRB_statistics.rds") # This file must be created by running the KmerAnalysis.R script
ACPA.data<-readRDS("Data_Output/Kmer_ACPA_TRB_statistics.rds") # This file must be created by running the KmerAnalysis.R script
RF.data<-readRDS("Data_Output/Kmer_RF_TRB_statistics.rds") # This file must be created by running the KmerAnalysis.R script

ra.specific.cd4.clusters.hla$IMID<-rep(NA,nrow(ra.specific.cd4.clusters.hla))
ra.specific.cd4.clusters.hla$RespGM<-rep(NA,nrow(ra.specific.cd4.clusters.hla))
ra.specific.cd4.clusters.hla$ActB<-rep(NA,nrow(ra.specific.cd4.clusters.hla))
ra.specific.cd4.clusters.hla$ACPA<-rep(NA,nrow(ra.specific.cd4.clusters.hla))
ra.specific.cd4.clusters.hla$RF<-rep(NA,nrow(ra.specific.cd4.clusters.hla))

phenos<-c("IMID","RespGM","ActB","ACPA","RF")

for (i in 1:nrow(ra.specific.cd4.clusters.hla)) {
  kmer<-as.character(ra.specific.cd4.clusters.hla[i,"pattern"])
  count=3
  for (ph in phenos) {
    count=count+1
    if (ph=="IMID") {data<-ractrl.data}
    if (ph=="RespGM") {data<-RespGM.data}
    if (ph=="ActB") {data<-ActB.data}
    if (ph=="ACPA") {data<-ACPA.data}
    if (ph=="RF") {data<-RF.data}
    pval<-data[data$Kmer==kmer,"Pval.hurdle"]
    ra.specific.cd4.clusters.hla[i,count]<-pval
  }
}

ra.specific.cd8.clusters.hla$IMID<-rep(NA,nrow(ra.specific.cd8.clusters.hla))
ra.specific.cd8.clusters.hla$RespGM<-rep(NA,nrow(ra.specific.cd8.clusters.hla))
ra.specific.cd8.clusters.hla$ActB<-rep(NA,nrow(ra.specific.cd8.clusters.hla))
ra.specific.cd8.clusters.hla$ACPA<-rep(NA,nrow(ra.specific.cd8.clusters.hla))
ra.specific.cd8.clusters.hla$RF<-rep(NA,nrow(ra.specific.cd8.clusters.hla))

phenos<-c("IMID","RespGM","ActB","ACPA","RF")

for (i in 1:nrow(ra.specific.cd8.clusters.hla)) {
  kmer<-as.character(ra.specific.cd8.clusters.hla[i,"pattern"])
  count=3
  for (ph in phenos) {
    count=count+1
    if (ph=="IMID") {data<-ractrl.data}
    if (ph=="RespGM") {data<-RespGM.data}
    if (ph=="ActB") {data<-ActB.data}
    if (ph=="ACPA") {data<-ACPA.data}
    if (ph=="RF") {data<-RF.data}
    pval<-data[data$Kmer==kmer,"Pval.hurdle"]
    if (length(pval)==0) {pval=NA}
    ra.specific.cd8.clusters.hla[i,count]<-pval
  }
}

saveRDS(ra.specific.cd4.clusters.hla, "Output_Data/Sig_HLA_RA_CD4_Assoc.rds")
saveRDS(ra.specific.cd8.clusters.hla, "Output_Data/Sig_HLA_RA_CD8_Assoc.rds")




