################  Loading libraries ########

suppressPackageStartupMessages({
  library(MAST);  library(immunarch);  library(plyr); library(data.table);  library(reshape);  library(scater); library(ComplexHeatmap);
  library(edgeR);  library(ROCR);  library(scRNA.seq.funcs);  library(gridExtra);  library(NMF);  library(RColorBrewer);  library(pheatmap);
  library(stringdist); library(Peptides);  library(ggpubr);   library(igraph);  library(stringdist);
  library(tcR);  library(purrr);  library(tidyr);  library(dendextend);  library(doParallel);  library(foreach);  library(motifStack)
  library(protr);  library(msa);  library(Biostrings);  library(ggseqlogo);  library(HDMD);  library(viridis);   library(readr)
})

################  Creating NoisET Data Format  ################

# Grouping by CDR3nt

dir.create("Data2analyze/NoisET_Analysis/")
names=c("Samples","Replicates")
chains<-c("IGK","IGH","IGL","TRB","TRA","TRG","TRD")

for (name in names) {
  for (c in chains){
    print(paste0(name," - ",c))
    inds<-grep("_R", list.files(paste0("Data2analyze/Immunarch_Format_BestVDJ/",c,"/")), value=T, invert=T)
    for (i in inds){
      if (name=="Samples") {
        data<-readRDS(paste0("Data2analyze/Immunarch_Format_BestVDJ/",c,"/",i))
        data2<-data %>% dplyr::group_by(CDR3.nt) %>% dplyr::summarise(Clones=sum(Clones), CDR3.aa=unique(CDR3.aa), V.name=paste(unlist(unique(V.name)), collapse=','), D.name=paste(unlist(unique(D.name)), collapse=','), J.name=paste(unlist(unique(J.name)), collapse=','), V.end=paste(unlist(unique(V.end)), collapse=','), D.start=paste(unlist(unique(D.start)), collapse=','), D.end=paste(unlist(unique(D.end)), collapse=','), J.start=paste(unlist(unique(J.start)), collapse=','), VJ.ins=paste(unlist(unique(VJ.ins)), collapse=','), VD.ins=paste(unlist(unique(VD.ins)), collapse=','), DJ.ins=paste(unlist(unique(DJ.ins)), collapse=','), Sequence=paste(unlist(unique(Sequence)), collapse=','))
        data2$Proportion<-data2$Clones/sum(data2$Clones, na.rm=T)
        ind.out<-gsub(".rds","_S",i)
        write.table(data2, file=paste0("Data2analyze/NoisET_Analysis/",c,"_",ind.out,".tsv"), quote=F, row.names=F, col.names=T, sep="\t")
      } else {
        i2read<-gsub(".rds","_R.rds",i)
        data<-readRDS(paste0("Data2analyze/Immunarch_Format_BestVDJ/",c,"/",i2read))
        data2<-data %>% dplyr::group_by(CDR3.nt) %>% dplyr::summarise(Clones=sum(Clones), CDR3.aa=unique(CDR3.aa), V.name=paste(unlist(unique(V.name)), collapse=','), D.name=paste(unlist(unique(D.name)), collapse=','), J.name=paste(unlist(unique(J.name)), collapse=','), V.end=paste(unlist(unique(V.end)), collapse=','), D.start=paste(unlist(unique(D.start)), collapse=','), D.end=paste(unlist(unique(D.end)), collapse=','), J.start=paste(unlist(unique(J.start)), collapse=','), VJ.ins=paste(unlist(unique(VJ.ins)), collapse=','), VD.ins=paste(unlist(unique(VD.ins)), collapse=','), DJ.ins=paste(unlist(unique(DJ.ins)), collapse=','), Sequence=paste(unlist(unique(Sequence)), collapse=','))
        data2$Proportion<-data2$Clones/sum(data2$Clones, na.rm=T)
        ind.out<-gsub(".rds","_R",i)
        write.table(data2, file=paste0("Data2analyze/NoisET_Analysis/",c,"_",ind.out,".tsv"), quote=F, row.names=F, col.names=T, sep="\t")
      }
    }
  }
}

################  Longitudinal Analysis Expansion  ################

setwd("Data2analyze/NoisET_Analysis/")
df.donor<-read.table("../Donor2RNA_Codes.txt", header=T)

chains<-c("TRB","TRA","TRG","TRD","IGK","IGH","IGL")

for (i in 1:nrow(df.donor)) {
  donor<-df.donor$Donor[i]
  wk0<-df.donor$RNA1[i]
  wk12<-df.donor$RNA2[i]
  for (c in chains) {
    
    # Step1
    comm1.noiset.wk0<-paste0("noiset-noise --path '' --f1 '",c,"_",wk0,"_S.tsv","' --f2 '",c,"_",wk0,"_R.tsv","' --specify --freq 'Proportion' --counts 'Clones' --ntCDR3 'CDR3.nt' --AACDR3 'CDR3.aa' --NB")
    system(comm1.noiset.wk0)
    system(paste0("mv nullpara1.txt NoisET_nullpara1_",c,"_",donor,".txt"))
    comm1.noiset.wk12<-paste0("noiset-noise --path '' --f1 '",c,"_",wk12,"_S.tsv","' --f2 '",c,"_",wk12,"_R.tsv","' --specify --freq 'Proportion' --counts 'Clones' --ntCDR3 'CDR3.nt' --AACDR3 'CDR3.aa' --NB")
    system(comm1.noiset.wk12)
    system(paste0("mv nullpara1.txt NoisET_nullpara2_",c,"_",donor,".txt"))
    
    # Step2
    wk0.sample.data<-read_tsv(paste0(c,"_",wk0,"_S.tsv"))
    wk0.replicate.data<-read_tsv(paste0(c,"_",wk0,"_R.tsv"))
    wk12.sample.data<-read_tsv(paste0(c,"_",wk12,"_S.tsv"))
    wk12.replicate.data<-read_tsv(paste0(c,"_",wk12,"_R.tsv"))
    wk0.nreads1<-sum(wk0.sample.data$Clones, na.rm=T)
    wk0.nreads2<-sum(wk0.replicate.data$Clones, na.rm=T)
    wk0.nclones<-length(unique(c(wk0.sample.data$CDR3.nt,wk0.replicate.data$CDR3.nt)))
    wk12.nreads1<-sum(wk12.sample.data$Clones, na.rm=T)
    wk12.nreads2<-sum(wk12.replicate.data$Clones, na.rm=T)
    wk12.nclones<-length(unique(c(wk12.sample.data$CDR3.nt,wk12.replicate.data$CDR3.nt)))    
    comm2.noiset.wk0<-paste0("noiset-nullgenerator --NB --nullpara 'NoisET_nullpara1_",c,"_",donor,".txt' --NreadsI ",wk0.nreads1," --NreadsII ",wk0.nreads2," --Nclones ",wk0.nclones," --output NoisET_",c,"_",donor,"_wk0")
    system(comm2.noiset.wk0)
    comm2.noiset.wk12<-paste0("noiset-nullgenerator --NB --nullpara 'NoisET_nullpara2_",c,"_",donor,".txt' --NreadsI ",wk12.nreads1," --NreadsII ",wk12.nreads2," --Nclones ",wk12.nclones," --output NoisET_",c,"_",donor,"_wk12")
    system(comm2.noiset.wk12)
    
    # Step3
    comm3.noiset.long<-paste0("noiset-detection --NB --nullpara1 'NoisET_nullpara1_",c,"_",donor,".txt' --nullpara2 'NoisET_nullpara2_",c,"_",donor,".txt' --path '' --f1 '",c,"_",wk0,"_S.tsv","' --f2 '",c,"_",wk12,"_S.tsv","' --specify --freq 'Proportion' --counts 'Clones' --ntCDR3 'CDR3.nt' --AACDR3 'CDR3.aa' --pval 0.5 --smedthresh 0 --output Detection_")
    system(comm3.noiset.long)
    system(paste0("mv Detection_",c,"_",wk0,"_S.tsv",c,"_",wk12,"_S.tsvtop_expanded.csv NoisET_Detection_",c,"_",donor,".csv"))
  }
}

################ Longitudinal Analysis Contraction  ################

df.donor<-read.table("../Donor2RNA_Codes.txt", header=T)
chains<-c("TRB","TRA","TRG","TRD","IGK","IGH","IGL")

dir.name<-"P05_Contracted/"
dir.create(dir.name)

for (i in 1:nrow(df.donor)) {
  donor<-df.donor$Donor[i]
  wk0<-df.donor$RNA2[i]
  wk12<-df.donor$RNA1[i]
  for (c in chains) {
    
    # Step1
    comm1.noiset.wk0<-paste0("noiset-noise --path '' --f1 '",c,"_",wk0,"_S.tsv","' --f2 '",c,"_",wk0,"_R.tsv","' --specify --freq 'Proportion' --counts 'Clones' --ntCDR3 'CDR3.nt' --AACDR3 'CDR3.aa' --NB")
    system(comm1.noiset.wk0)
    system(paste0("mv nullpara1.txt NoisET_nullpara1_",c,"_",donor,".txt"))
    comm1.noiset.wk12<-paste0("noiset-noise --path '' --f1 '",c,"_",wk12,"_S.tsv","' --f2 '",c,"_",wk12,"_R.tsv","' --specify --freq 'Proportion' --counts 'Clones' --ntCDR3 'CDR3.nt' --AACDR3 'CDR3.aa' --NB")
    system(comm1.noiset.wk12)
    system(paste0("mv nullpara1.txt NoisET_nullpara2_",c,"_",donor,".txt"))
    
    # Step2
    wk0.sample.data<-read_tsv(paste0(c,"_",wk0,"_S.tsv"))
    wk0.replicate.data<-read_tsv(paste0(c,"_",wk0,"_R.tsv"))
    wk12.sample.data<-read_tsv(paste0(c,"_",wk12,"_S.tsv"))
    wk12.replicate.data<-read_tsv(paste0(c,"_",wk12,"_R.tsv"))
    wk0.nreads1<-sum(wk0.sample.data$Clones, na.rm=T)
    wk0.nreads2<-sum(wk0.replicate.data$Clones, na.rm=T)
    wk0.nclones<-length(unique(c(wk0.sample.data$CDR3.nt,wk0.replicate.data$CDR3.nt)))
    wk12.nreads1<-sum(wk12.sample.data$Clones, na.rm=T)
    wk12.nreads2<-sum(wk12.replicate.data$Clones, na.rm=T)
    wk12.nclones<-length(unique(c(wk12.sample.data$CDR3.nt,wk12.replicate.data$CDR3.nt)))    
    comm2.noiset.wk0<-paste0("noiset-nullgenerator --NB --nullpara 'NoisET_nullpara1_",c,"_",donor,".txt' --NreadsI ",wk0.nreads1," --NreadsII ",wk0.nreads2," --Nclones ",wk0.nclones," --output NoisET_",c,"_",donor,"_wk0")
    system(comm2.noiset.wk0)
    comm2.noiset.wk12<-paste0("noiset-nullgenerator --NB --nullpara 'NoisET_nullpara2_",c,"_",donor,".txt' --NreadsI ",wk12.nreads1," --NreadsII ",wk12.nreads2," --Nclones ",wk12.nclones," --output NoisET_",c,"_",donor,"_wk12")
    system(comm2.noiset.wk12)
    
    # Step3
    comm3.noiset.long<-paste0("noiset-detection --NB --nullpara1 'NoisET_nullpara1_",c,"_",donor,".txt' --nullpara2 'NoisET_nullpara2_",c,"_",donor,".txt' --path '' --f1 '",c,"_",wk0,"_S.tsv","' --f2 '",c,"_",wk12,"_S.tsv","' --specify --freq 'Proportion' --counts 'Clones' --ntCDR3 'CDR3.nt' --AACDR3 'CDR3.aa' --pval 0.5 --smedthresh 0 --output Detection_")
    system(comm3.noiset.long)
    system(paste0("mv Detection_",c,"_",wk0,"_S.tsv",c,"_",wk12,"_S.tsvtop_expanded.csv NoisET_Detection_",c,"_",donor,".csv"))
  }
}

system(paste0("mv NoisET_* ",dir.name))
system(paste0("mv synthetic_* ",dir.name))

################  QC &  Redoing the analysis by averaging the parameters ########
#### Contracted ####

dir="P05_Contracted/"
pval="0.05"
type="contracted"

df.donor<-read.table("../Donor2RNA_Codes.txt", header=T)
metadata<-readRDS("../Metadata_IGH.rds")
colnames(metadata)[1]<-"RNA1"
df.pheno<-merge(df.donor, metadata, by="RNA1")

df.qc<-data.frame("Donor"="test", "Week0"="test", "Week12"="test", "Chain"="test", "File"="0", stringsAsFactors=F)
files.out<-system(paste0("ls ",dir," | grep NoisET_Detection_"), intern=T)
chains<-c("TRB","TRA","TRG","TRD","IGK","IGH","IGL")
count=0
for (i in 1:nrow(df.donor)) {
  donor<-as.character(df.donor$Donor[i])
  wk0<-as.character(df.donor$RNA1[i])
  wk12<-as.character(df.donor$RNA2[i])
  for (c in chains) {
    value<-length(grep(paste0("NoisET_Detection_",c,"_",donor,".csv"), files.out))
    count=count+1
    df.qc[count,]<-c(donor,wk0,wk12,c,value)
  }
}

for (c in chains) {
  # wk0
  avg.rho.wk0<-mean(as.numeric(system(paste0("grep alph_rho ",dir,"/NoisET_nullpara1_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.beta.wk0<-mean(as.numeric(system(paste0("grep beta ",dir,"/NoisET_nullpara1_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.alpha.wk0<-mean(as.numeric(system(paste0("grep alpha ",dir,"/NoisET_nullpara1_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.fmin.wk0<-mean(as.numeric(system(paste0("grep fmin ",dir,"/NoisET_nullpara1_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  df.wk0<-data.frame("label"=c("alph_rho","beta","alpha","fmin"), "value"=c(avg.rho.wk0,avg.beta.wk0,avg.alpha.wk0,avg.fmin.wk0))
  write.table(df.wk0,paste0("NoisET_nullpara1_",c,"_average.txt"), quote=F, sep="\t")
  system(paste0("sed 's/label/\tlabel/g' NoisET_nullpara1_",c,"_average.txt > NoisET_nullpara1_",c,"_Average.txt"))
  system(paste0("rm NoisET_nullpara1_",c,"_average.txt"))
  # wk12
  avg.rho.wk12<-mean(as.numeric(system(paste0("grep alph_rho ",dir,"/NoisET_nullpara2_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.beta.wk12<-mean(as.numeric(system(paste0("grep beta ",dir,"/NoisET_nullpara2_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.alpha.wk12<-mean(as.numeric(system(paste0("grep alpha ",dir,"/NoisET_nullpara2_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.fmin.wk12<-mean(as.numeric(system(paste0("grep fmin ",dir,"/NoisET_nullpara2_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  df.wk12<-data.frame("label"=c("alph_rho","beta","alpha","fmin"), "value"=c(avg.rho.wk12,avg.beta.wk12,avg.alpha.wk12,avg.fmin.wk12))
  write.table(df.wk12,paste0("NoisET_nullpara2_",c,"_average.txt"), quote=F, sep="\t")
  system(paste0("sed 's/label/\tlabel/g' NoisET_nullpara2_",c,"_average.txt > NoisET_nullpara2_",c,"_Average.txt"))
  system(paste0("rm NoisET_nullpara2_",c,"_average.txt"))
}

df.fail<-df.qc[df.qc$File==0,]
for (i in 1:nrow(df.fail)) {
  donor<-as.character(df.fail$Donor[i])
  for (c in chains) {
    if (type=="expanded") {wk0<-as.character(df.fail$Week0[i]); wk12<-as.character(df.fail$Week12[i]);  nullpara1<-paste0("NoisET_nullpara1_",c,"_Average.txt");  nullpara2<-paste0("NoisET_nullpara2_",c,"_Average.txt")}
    if (type=="contracted") {wk0<-as.character(df.fail$Week12[i]); wk12<-as.character(df.fail$Week0[i]);  nullpara1<-paste0("NoisET_nullpara2_",c,"_Average.txt");  nullpara2<-paste0("NoisET_nullpara1_",c,"_Average.txt")}
    
    # Step2
    wk0.sample.data<-read_tsv(paste0(c,"_",wk0,"_S.tsv"))
    wk0.replicate.data<-read_tsv(paste0(c,"_",wk0,"_R.tsv"))
    wk12.sample.data<-read_tsv(paste0(c,"_",wk12,"_S.tsv"))
    wk12.replicate.data<-read_tsv(paste0(c,"_",wk12,"_R.tsv"))
    wk0.nreads1<-sum(wk0.sample.data$Clones, na.rm=T)
    wk0.nreads2<-sum(wk0.replicate.data$Clones, na.rm=T)
    wk0.nclones<-length(unique(c(wk0.sample.data$CDR3.nt,wk0.replicate.data$CDR3.nt)))
    wk12.nreads1<-sum(wk12.sample.data$Clones, na.rm=T)
    wk12.nreads2<-sum(wk12.replicate.data$Clones, na.rm=T)
    wk12.nclones<-length(unique(c(wk12.sample.data$CDR3.nt,wk12.replicate.data$CDR3.nt)))    
    comm2.noiset.wk0<-paste0("noiset-nullgenerator --NB --nullpara '",nullpara1,"' --NreadsI ",wk0.nreads1," --NreadsII ",wk0.nreads2," --Nclones ",wk0.nclones," --output NoisET_",c,"_",donor,"_wk0")
    system(comm2.noiset.wk0)
    comm2.noiset.wk12<-paste0("noiset-nullgenerator --NB --nullpara '",nullpara2,"' --NreadsI ",wk12.nreads1," --NreadsII ",wk12.nreads2," --Nclones ",wk12.nclones," --output NoisET_",c,"_",donor,"_wk12")
    system(comm2.noiset.wk12)
    
    # Step3
    comm3.noiset.long<-paste0("noiset-detection --NB --nullpara1 '",nullpara1,"' --nullpara2 '",nullpara2,"' --path '' --f1 '",c,"_",wk0,"_S.tsv","' --f2 '",c,"_",wk12,"_S.tsv","' --specify --freq 'Proportion' --counts 'Clones' --ntCDR3 'CDR3.nt' --AACDR3 'CDR3.aa' --pval ",pval," --smedthresh 0 --output Detection_")
    system(comm3.noiset.long)
    system(paste0("mv Detection_",c,"_",wk0,"_S.tsv",c,"_",wk12,"_S.tsvtop_expanded.csv NoisET_Detection_",c,"_",donor,".csv"))
    system(paste0("mv NoisET_Detection_",c,"_",donor,".csv ",dir))
  }
}

#### Expanded ####

dir="P05_Expanded/"
pval="0.05"
type="expanded"

df.donor<-read.table("../Donor2RNA_Codes.txt", header=T)
metadata<-readRDS("../Metadata_IGH.rds")
colnames(metadata)[1]<-"RNA1"
df.pheno<-merge(df.donor, metadata, by="RNA1")

df.qc<-data.frame("Donor"="test", "Week0"="test", "Week12"="test", "Chain"="test", "File"="0", stringsAsFactors=F)
files.out<-system(paste0("ls ",dir," | grep NoisET_Detection_"), intern=T)
chains<-c("TRB","TRA","TRG","TRD","IGK","IGH","IGL")
count=0
for (i in 1:nrow(df.donor)) {
  donor<-as.character(df.donor$Donor[i])
  wk0<-as.character(df.donor$RNA1[i])
  wk12<-as.character(df.donor$RNA2[i])
  for (c in chains) {
    value<-length(grep(paste0("NoisET_Detection_",c,"_",donor,".csv"), files.out))
    count=count+1
    df.qc[count,]<-c(donor,wk0,wk12,c,value)
  }
}

for (c in chains) {
  # wk0
  avg.rho.wk0<-mean(as.numeric(system(paste0("grep alph_rho ",dir,"/NoisET_nullpara1_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.beta.wk0<-mean(as.numeric(system(paste0("grep beta ",dir,"/NoisET_nullpara1_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.alpha.wk0<-mean(as.numeric(system(paste0("grep alpha ",dir,"/NoisET_nullpara1_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.fmin.wk0<-mean(as.numeric(system(paste0("grep fmin ",dir,"/NoisET_nullpara1_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  df.wk0<-data.frame("label"=c("alph_rho","beta","alpha","fmin"), "value"=c(avg.rho.wk0,avg.beta.wk0,avg.alpha.wk0,avg.fmin.wk0))
  write.table(df.wk0,paste0("NoisET_nullpara1_",c,"_average.txt"), quote=F, sep="\t")
  system(paste0("sed 's/label/\tlabel/g' NoisET_nullpara1_",c,"_average.txt > NoisET_nullpara1_",c,"_Average.txt"))
  system(paste0("rm NoisET_nullpara1_",c,"_average.txt"))
  # wk12
  avg.rho.wk12<-mean(as.numeric(system(paste0("grep alph_rho ",dir,"/NoisET_nullpara2_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.beta.wk12<-mean(as.numeric(system(paste0("grep beta ",dir,"/NoisET_nullpara2_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.alpha.wk12<-mean(as.numeric(system(paste0("grep alpha ",dir,"/NoisET_nullpara2_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  avg.fmin.wk12<-mean(as.numeric(system(paste0("grep fmin ",dir,"/NoisET_nullpara2_",c,"_*.txt | awk '{print $3}'"), intern=T)), na.rm=T)
  df.wk12<-data.frame("label"=c("alph_rho","beta","alpha","fmin"), "value"=c(avg.rho.wk12,avg.beta.wk12,avg.alpha.wk12,avg.fmin.wk12))
  write.table(df.wk12,paste0("NoisET_nullpara2_",c,"_average.txt"), quote=F, sep="\t")
  system(paste0("sed 's/label/\tlabel/g' NoisET_nullpara2_",c,"_average.txt > NoisET_nullpara2_",c,"_Average.txt"))
  system(paste0("rm NoisET_nullpara2_",c,"_average.txt"))
}

df.fail<-df.qc[df.qc$File==0,]
for (i in 1:nrow(df.fail)) {
  donor<-as.character(df.fail$Donor[i])
  for (c in chains) {
    if (type=="expanded") {wk0<-as.character(df.fail$Week0[i]); wk12<-as.character(df.fail$Week12[i]);  nullpara1<-paste0("NoisET_nullpara1_",c,"_Average.txt");  nullpara2<-paste0("NoisET_nullpara2_",c,"_Average.txt")}
    if (type=="contracted") {wk0<-as.character(df.fail$Week12[i]); wk12<-as.character(df.fail$Week0[i]);  nullpara1<-paste0("NoisET_nullpara2_",c,"_Average.txt");  nullpara2<-paste0("NoisET_nullpara1_",c,"_Average.txt")}
    
    # Step2
    wk0.sample.data<-read_tsv(paste0(c,"_",wk0,"_S.tsv"))
    wk0.replicate.data<-read_tsv(paste0(c,"_",wk0,"_R.tsv"))
    wk12.sample.data<-read_tsv(paste0(c,"_",wk12,"_S.tsv"))
    wk12.replicate.data<-read_tsv(paste0(c,"_",wk12,"_R.tsv"))
    wk0.nreads1<-sum(wk0.sample.data$Clones, na.rm=T)
    wk0.nreads2<-sum(wk0.replicate.data$Clones, na.rm=T)
    wk0.nclones<-length(unique(c(wk0.sample.data$CDR3.nt,wk0.replicate.data$CDR3.nt)))
    wk12.nreads1<-sum(wk12.sample.data$Clones, na.rm=T)
    wk12.nreads2<-sum(wk12.replicate.data$Clones, na.rm=T)
    wk12.nclones<-length(unique(c(wk12.sample.data$CDR3.nt,wk12.replicate.data$CDR3.nt)))    
    comm2.noiset.wk0<-paste0("noiset-nullgenerator --NB --nullpara '",nullpara1,"' --NreadsI ",wk0.nreads1," --NreadsII ",wk0.nreads2," --Nclones ",wk0.nclones," --output NoisET_",c,"_",donor,"_wk0")
    system(comm2.noiset.wk0)
    comm2.noiset.wk12<-paste0("noiset-nullgenerator --NB --nullpara '",nullpara2,"' --NreadsI ",wk12.nreads1," --NreadsII ",wk12.nreads2," --Nclones ",wk12.nclones," --output NoisET_",c,"_",donor,"_wk12")
    system(comm2.noiset.wk12)
    
    # Step3
    comm3.noiset.long<-paste0("noiset-detection --NB --nullpara1 '",nullpara1,"' --nullpara2 '",nullpara2,"' --path '' --f1 '",c,"_",wk0,"_S.tsv","' --f2 '",c,"_",wk12,"_S.tsv","' --specify --freq 'Proportion' --counts 'Clones' --ntCDR3 'CDR3.nt' --AACDR3 'CDR3.aa' --pval ",pval," --smedthresh 0 --output Detection_")
    system(comm3.noiset.long)
    system(paste0("mv Detection_",c,"_",wk0,"_S.tsv",c,"_",wk12,"_S.tsvtop_expanded.csv NoisET_Detection_",c,"_",donor,".csv"))
    system(paste0("mv NoisET_Detection_",c,"_",donor,".csv ",dir))
  }
}

########  Summarizing Results  ####

df.donor<-read.table("../Donor2RNA_Codes.txt", header=T)
metadata<-readRDS("../Metadata_IGH.rds")
colnames(metadata)[1]<-"RNA1"
df.pheno<-merge(df.donor, metadata, by="RNA1")

df.res<-data.frame("Chain"="test", "Donor"="test", "Week0"="test", "Week12"="test", "All.Sig"=0, "Exp.Sig"=0, "Cont.Sig"=0, "Exp.Sig.Clones"="test", "Cont.Sig.Clones"="test", stringsAsFactors=F)
chains<-c("TRB","TRA","TRG","TRD","IGK","IGH","IGL")
count=0
for (i in 1:nrow(df.pheno)) {
  donor<-as.character(df.pheno$Donor[i])
  wk0<-as.character(df.pheno$RNA1[i])
  wk12<-as.character(df.pheno$RNA2[i])
  print(donor)
  for (c in chains) {
    data.wk0<-read_tsv(paste0(c,"_",wk0,"_S.tsv")) 
    data.wk12<-read_tsv(paste0(c,"_",wk12,"_S.tsv"))
    df.merged<-merge(data.wk0, data.wk12, by="CDR3.nt", all=T)
    df.merged.filt<-df.merged[,c("CDR3.nt","CDR3.aa.x","CDR3.aa.y","Clones.x","Clones.y","Proportion.x","Proportion.y")]
    df.merged.filt$CDR3.aa.x[is.na(df.merged.filt$CDR3.aa.x)]<-""
    df.merged.filt$CDR3.aa.y[is.na(df.merged.filt$CDR3.aa.y)]<-""
    new.vect<-list()
    for (n in 1:nrow(df.merged.filt)) {
      value1<-df.merged.filt$CDR3.aa.x[n]
      value2<-df.merged.filt$CDR3.aa.y[n]
      value3<-c(value1,value2)
      value3<-unique(value3[!(value3 %in% "")])
      new.vect[[n]]<-value3
    }
    df.merged.filt$CDR3.aa<-unlist(new.vect)
    df.merged.filt$Proportion.x[is.na(df.merged.filt$Proportion.x)]<-0
    df.merged.filt$Proportion.y[is.na(df.merged.filt$Proportion.y)]<-0
    df.merged.filt<-df.merged.filt[,c(1,8,4,5,6,7)]
    df.pval.exp<-read.csv(paste0("P05_Expanded/NoisET_Detection_",c,"_",donor,".csv"), sep="\t")[,c(10,12)]
    colnames(df.pval.exp)[1]<-"CDR3.nt"
    df.merged.final<-merge(df.merged.filt, df.pval.exp, by="CDR3.nt", all=T)
    df.merged.final$X.1.P.s.0..[is.na(df.merged.final$X.1.P.s.0..)]<-1
    colnames(df.merged.final)<-c("CDR3.nt","CDR3.aa","Clones.wk0","Clones.wk12","Proportion.wk0","Proportion.wk12","Pval.Exp")
    df.pval.cont<-read.csv(paste0("P05_Contracted/NoisET_Detection_",c,"_",donor,".csv"), sep="\t")[,c(10,12)]
    colnames(df.pval.cont)[1]<-"CDR3.nt"
    df.merged.final2<-merge(df.merged.final, df.pval.cont, by="CDR3.nt", all=T)
    df.merged.final2$X.1.P.s.0..[is.na(df.merged.final2$X.1.P.s.0..)]<-1
    colnames(df.merged.final2)[8]<-"Pval.Cont"
    df.merged.final2$Type<-rep("NoSig",nrow(df.merged.final2))
    df.merged.final2$Type[df.merged.final2$Pval.Exp<0.05]<-"Sig.Expanded"
    df.merged.final2$Type[df.merged.final2$Pval.Cont<0.05]<-"Sig.Contracted"
    all.sig<-nrow(df.merged.final2[(df.merged.final2$Type=="Sig.Expanded" | df.merged.final2$Type=="Sig.Contracted"),])
    exp.sig<-nrow(df.merged.final2[df.merged.final2$Type=="Sig.Expanded" & df.merged.final2$Proportion.wk0<df.merged.final2$Proportion.wk12,])
    cont.sig<-nrow(df.merged.final2[df.merged.final2$Type=="Sig.Contracted" & df.merged.final2$Proportion.wk0>df.merged.final2$Proportion.wk12,])
    exp.sig.clones<-do.call(paste, c(as.list(unique(df.merged.final2[df.merged.final2$Type=="Sig.Expanded" & df.merged.final2$Proportion.wk0<df.merged.final2$Proportion.wk12,"CDR3.aa"])), sep=";"))
    cont.sig.clones<-do.call(paste, c(as.list(unique(df.merged.final[df.merged.final2$Type=="Sig.Contracted" & df.merged.final2$Proportion.wk0>df.merged.final2$Proportion.wk12,"CDR3.aa"])), sep=";"))
    if (length(exp.sig.clones)==0) {exp.sig.clones<-NA}
    if (length(cont.sig.clones)==0) {cont.sig.clones<-NA}
    count=count+1
    df.res[count,]<-c(c,donor,wk0,wk12,all.sig,exp.sig,cont.sig,exp.sig.clones,cont.sig.clones)
  }
}
saveRDS(df.res,"../Output_Data/NoisET_Summary_ExpCont.rds")



########  Visualizing results  ####

df.res<-readRDS("../Output_Data/NoisET_Summary_ExpCont.rds")

#### Overview ####

df.res2plot<-df.res[,c(1:7)]
df.res2plot<-gather(df.res2plot, condition, measurement, All.Sig:Cont.Sig, factor_key=T)
df.res2plot$Chain<-as.factor(df.res2plot$Chain)
df.res2plot$measurement<-as.numeric(df.res2plot$measurement)
p.sample<-ggplot(df.res2plot, aes(x=Chain, y=measurement, fill=condition)) + geom_bar(stat="identity", position=position_dodge()) +
  labs(x=paste0("\n"), y="\nSignificant Clones (N, P<0.05)\n") + facet_wrap(~Donor) + theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("azure4","chartreuse4","brown3"))
jpeg("../Output_Data/NoisET_Summary_BySample.jpeg", res=300, width=7300, height=4800)
plot(p.sample)
dev.off()

p1<-ggboxplot(df.res2plot, x="Chain", y="measurement", color="condition", add="jitter", fill="condition") +
  labs(title="\n\n", x="", y=paste0("\nSignificant Clones (N, P<0.05)\n")) +
  scale_fill_manual(values=c("azure4","chartreuse4","brown3")) + scale_color_manual(values=c("black","black","black","black","black","black","black")) +
  theme(plot.title=element_text(hjust=0.5, size=16, face="bold"), legend.position="right", axis.text.x = element_text(angle=45, hjust=1, size=14))

df.chain2plot<-data.frame("Chain"=c("IGH","IGK","IGL","TRA","TRB","TRD","TRG"))
df.chain2plot$All.Sig[1]<-sum(df.res2plot[(df.res2plot$Chain=="IGH" & df.res2plot$condition=="All.Sig"),"measurement"], na.rm=T)
df.chain2plot$Exp.Sig[1]<-sum(df.res2plot[(df.res2plot$Chain=="IGH" & df.res2plot$condition=="Exp.Sig"),"measurement"], na.rm=T)
df.chain2plot$Cont.Sig[1]<-sum(df.res2plot[(df.res2plot$Chain=="IGH" & df.res2plot$condition=="Cont.Sig"),"measurement"], na.rm=T)
df.chain2plot$All.Sig[2]<-sum(df.res2plot[(df.res2plot$Chain=="IGK" & df.res2plot$condition=="All.Sig"),"measurement"], na.rm=T)
df.chain2plot$Exp.Sig[2]<-sum(df.res2plot[(df.res2plot$Chain=="IGK" & df.res2plot$condition=="Exp.Sig"),"measurement"], na.rm=T)
df.chain2plot$Cont.Sig[2]<-sum(df.res2plot[(df.res2plot$Chain=="IGK" & df.res2plot$condition=="Cont.Sig"),"measurement"], na.rm=T)
df.chain2plot$All.Sig[3]<-sum(df.res2plot[(df.res2plot$Chain=="IGL" & df.res2plot$condition=="All.Sig"),"measurement"], na.rm=T)
df.chain2plot$Exp.Sig[3]<-sum(df.res2plot[(df.res2plot$Chain=="IGL" & df.res2plot$condition=="Exp.Sig"),"measurement"], na.rm=T)
df.chain2plot$Cont.Sig[3]<-sum(df.res2plot[(df.res2plot$Chain=="IGL" & df.res2plot$condition=="Cont.Sig"),"measurement"], na.rm=T)
df.chain2plot$All.Sig[4]<-sum(df.res2plot[(df.res2plot$Chain=="TRA" & df.res2plot$condition=="All.Sig"),"measurement"], na.rm=T)
df.chain2plot$Exp.Sig[4]<-sum(df.res2plot[(df.res2plot$Chain=="TRA" & df.res2plot$condition=="Exp.Sig"),"measurement"], na.rm=T)
df.chain2plot$Cont.Sig[4]<-sum(df.res2plot[(df.res2plot$Chain=="TRA" & df.res2plot$condition=="Cont.Sig"),"measurement"], na.rm=T)
df.chain2plot$All.Sig[5]<-sum(df.res2plot[(df.res2plot$Chain=="TRB" & df.res2plot$condition=="All.Sig"),"measurement"], na.rm=T)
df.chain2plot$Exp.Sig[5]<-sum(df.res2plot[(df.res2plot$Chain=="TRB" & df.res2plot$condition=="Exp.Sig"),"measurement"], na.rm=T)
df.chain2plot$Cont.Sig[5]<-sum(df.res2plot[(df.res2plot$Chain=="TRB" & df.res2plot$condition=="Cont.Sig"),"measurement"], na.rm=T)
df.chain2plot$All.Sig[6]<-sum(df.res2plot[(df.res2plot$Chain=="TRD" & df.res2plot$condition=="All.Sig"),"measurement"], na.rm=T)
df.chain2plot$Exp.Sig[6]<-sum(df.res2plot[(df.res2plot$Chain=="TRD" & df.res2plot$condition=="Exp.Sig"),"measurement"], na.rm=T)
df.chain2plot$Cont.Sig[6]<-sum(df.res2plot[(df.res2plot$Chain=="TRD" & df.res2plot$condition=="Cont.Sig"),"measurement"], na.rm=T)
df.chain2plot$All.Sig[7]<-sum(df.res2plot[(df.res2plot$Chain=="TRG" & df.res2plot$condition=="All.Sig"),"measurement"], na.rm=T)
df.chain2plot$Exp.Sig[7]<-sum(df.res2plot[(df.res2plot$Chain=="TRG" & df.res2plot$condition=="Exp.Sig"),"measurement"], na.rm=T)
df.chain2plot$Cont.Sig[7]<-sum(df.res2plot[(df.res2plot$Chain=="TRG" & df.res2plot$condition=="Cont.Sig"),"measurement"], na.rm=T)

df.chain2plot<-gather(df.chain2plot, condition, measurement, All.Sig:Cont.Sig, factor_key=T)
p2<-ggplot(df.chain2plot, aes(x=Chain, y=measurement, fill=condition)) + geom_bar(stat="identity", colour="black", position=position_dodge()) +
  labs(x=paste0("\n"), y="\nSignificant Clones (N, P<0.05)\n") + theme_classic() + theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("azure4","chartreuse4","brown3"))

jpeg("../Output_Data/NoisET_Summary_BySampleChain.jpeg", res=300, width=2900, height=4000)
grid.arrange(p1,p2)
dev.off()


#### Single Sample ####

df.donor<-read.table("../Donor2RNA_Codes.txt", header=T)
metadata<-readRDS("../Metadata_IGH.rds")
colnames(metadata)[1]<-"RNA1"
df.pheno<-merge(df.donor, metadata, by="RNA1")

df.res<-data.frame("Chain"="test", "Donor"="test", "Week0"="test", "Week12"="test", "All.Sig"=0, "Exp.Sig"=0, "Cont.Sig"=0, "Exp.Sig.Clones"="test", "Cont.Sig.Clones"="test", stringsAsFactors=F)
chains<-c("TRB","TRA","TRG","TRD","IGK","IGH","IGL")
count=0
for (i in 1:nrow(df.pheno)) {
  donor<-as.character(df.pheno$Donor[i])
  wk0<-as.character(df.pheno$RNA1[i])
  wk12<-as.character(df.pheno$RNA2[i])
  print(donor)
  for (c in chains) {
    data.wk0<-read_tsv(paste0(c,"_",wk0,"_S.tsv")) 
    data.wk12<-read_tsv(paste0(c,"_",wk12,"_S.tsv"))
    df.merged<-merge(data.wk0, data.wk12, by="CDR3.nt", all=T)
    df.merged.filt<-df.merged[,c("CDR3.nt","CDR3.aa.x","CDR3.aa.y","Clones.x","Clones.y","Proportion.x","Proportion.y")]
    df.merged.filt$CDR3.aa.x[is.na(df.merged.filt$CDR3.aa.x)]<-""
    df.merged.filt$CDR3.aa.y[is.na(df.merged.filt$CDR3.aa.y)]<-""
    df.merged.filt$CDR3.aa<-paste0(df.merged.filt$CDR3.aa.x,df.merged.filt$CDR3.aa.y)
    df.merged.filt$Proportion.x[is.na(df.merged.filt$Proportion.x)]<-0
    df.merged.filt$Proportion.y[is.na(df.merged.filt$Proportion.y)]<-0
    df.merged.filt<-df.merged.filt[,c(1,8,4,5,6,7)]
    df.pval.exp<-read.csv(paste0("P05_Expanded/NoisET_Detection_",c,"_",donor,".csv"), sep="\t")[,c(10,12)]
    colnames(df.pval.exp)[1]<-"CDR3.nt"
    df.merged.final<-merge(df.merged.filt, df.pval.exp, by="CDR3.nt", all=T)
    df.merged.final$X.1.P.s.0..[is.na(df.merged.final$X.1.P.s.0..)]<-1
    colnames(df.merged.final)<-c("CDR3.nt","CDR3.aa","Clones.wk0","Clones.wk12","Proportion.wk0","Proportion.wk12","Pval.Exp")
    df.pval.cont<-read.csv(paste0("P05_Contracted/NoisET_Detection_",c,"_",donor,".csv"), sep="\t")[,c(10,12)]
    colnames(df.pval.cont)[1]<-"CDR3.nt"
    df.merged.final2<-merge(df.merged.final, df.pval.cont, by="CDR3.nt", all=T)
    df.merged.final2$X.1.P.s.0..[is.na(df.merged.final2$X.1.P.s.0..)]<-1
    colnames(df.merged.final2)[8]<-"Pval.Cont"
    df.merged.final2$Type<-rep("NoSig",nrow(df.merged.final2))
    df.merged.final2$Type[df.merged.final2$Pval.Exp<0.05]<-"Sig.Expanded"
    df.merged.final2$Type[df.merged.final2$Pval.Cont<0.05]<-"Sig.Contracted"
    
    match.nosig<-0;
    match.exp<-0;
    match.cont<-0;
    if ("NoSig" %in% unique(df.merged.final2$Type)) {match.nosig<-1}
    if ("Sig.Expanded" %in% unique(df.merged.final2$Type)) {match.exp<-1}
    if ("Sig.Contracted" %in% unique(df.merged.final2$Type)) {match.cont<-1}
    
    if (match.nosig==1 & match.exp==1 & match.cont==1) {colors<-c("azure4","brown3","chartreuse4"); shapes<-c(1,19,19)}
    if (match.nosig==1 & match.exp==0 & match.cont==0) {colors<-c("azure4"); shapes<-c(1)}
    if (match.nosig==1 & match.exp==1 & match.cont==0) {colors<-c("azure4","chartreuse4"); shapes<-c(1,19)}
    if (match.nosig==1 & match.exp==0 & match.cont==1) {colors<-c("azure4","brown3"); shapes<-c(1,19)}
    if (match.nosig==0 & match.exp==1 & match.cont==1) {colors<-c("brown3","chartreuse4"); shapes<-c(19,19)}
    
    p.scat<-ggscatter(df.merged.final2, x="Proportion.wk0", y="Proportion.wk12", color="Type", shape="Type") +
      scale_color_manual(values =colors) + scale_shape_manual(values=shapes) +
      ylim(min(c(df.merged.final$Proportion.wk0,df.merged.final$Proportion.wk12)),max(c(df.merged.final$Proportion.wk0,df.merged.final$Proportion.wk12))) +
      xlim(min(c(df.merged.final$Proportion.wk0,df.merged.final$Proportion.wk12)),max(c(df.merged.final$Proportion.wk0,df.merged.final$Proportion.wk12))) +
      labs(x="\nClone Proportion Week0\n", y="\nClone Proportion Week12\n", title=paste0("\n",donor," - ",c,"\n")) +
      geom_abline(intercept=0, slope=1, color="purple", linetype="dashed", size=0.1) +
      theme(plot.title=element_text(hjust=0.5, size=16, face="bold"), legend.position="right")
    
    jpeg(paste0("../Output_Data/NoisET_ExpCont_BySample/",donor,"_",c,".jpeg"), res=300, width=2800, height=2600)
    grid.arrange(p.scat)
    dev.off()
  }
}


########  Clone characterization ########

df.donor<-read.table("../Donor2RNA_Codes.txt", header=T)
metadata<-readRDS("../Metadata_IGH.rds")
colnames(metadata)[1]<-"RNA1"
df.pheno<-merge(df.donor, metadata, by="RNA1")
       
#### Expanded   ####

type="Expanded"
dir.create(paste0("../../Output_Data/NoisET_Characterization_",type,"/"))

get.lev.dc.pval<-function(chain.clones.cont.unique, chain.clones.all, c, type, nperm) {
  
  graph<-mutation.network(chain.clones.cont.unique, .method="lev", .max.errors=1)
  deg<-degree(graph, mode="all")
  V(graph)$size<-deg*0.6
  V(graph)$frame.color<-"white"
  V(graph)$color<-"orange"
  V(graph)$label.cex<-0.3
  E(graph)$arrow.mode<-0
  
  jpeg(paste0("../../Output_Data/NoisET_Characterization_",type,"/SimilarityNetwork",c,".jpeg"), res=300, height=3000, width=3000)
  plot(graph, layout=layout_with_kk, main=paste0("\n",type," Cloness (DC=", round(mean(degree(graph)),2),")"))
  dev.off()
  
  # Leveinshtein & Degree Centrality
  df.perm.lev.random.chain<-data.frame("Chain"="test", "Perm"=0, "Mean.LV"=0, stringsAsFactors=F)
  df.perm.dc.random.chain<-data.frame("Chain"="test", "Perm"=0, "Mean.DC"=0, stringsAsFactors=F)
  set.size<-length(chain.clones.cont.unique)
  for (i in 1:nperm) {
    print(paste0(c," - ",i))
    random.clones<-sample(chain.clones.all, set.size, replace=F)
    lv.mean<-mean(stringdistmatrix(random.clones, method="lv"))
    df.perm.lev.random.chain[i,]<-c(c, i, lv.mean)
    dc.mean<-mean(degree(mutation.network(random.clones, .method="lev", .max.errors=1)))
    df.perm.dc.random.chain[i,]<-c(c, i, dc.mean)
  }
  obs.lev<-mean(stringdistmatrix(chain.clones.cont.unique, method="lv"))
  obs.dc<-mean(degree(graph))
  
  random.lev.mean<-mean(as.numeric(df.perm.lev.random.chain$Mean.LV))
  random.dc.mean<-mean(as.numeric(df.perm.dc.random.chain$Mean.DC))
  
  pval.lev<-(nrow(df.perm.lev.random.chain[(as.numeric(df.perm.lev.random.chain$Mean.LV)<obs.lev),])+1)/(nperm+1)
  pval.dc<-(nrow(df.perm.dc.random.chain[(as.numeric(df.perm.dc.random.chain$Mean.DC)>obs.dc),])+1)/(nperm+1)
  
  df.out<-data.frame("Chain"=c(c,c),
                     "Measure"=c("Leveinshtein Distance","Degree Centrality"),
                     "N.Perm"=c(nperm,nperm),
                     "Measure.Obs"=c(obs.lev, obs.dc),
                     "Measure.Perm"=c(mean(as.numeric(df.perm.lev.random.chain$Mean.LV)),
                                      mean(as.numeric(df.perm.dc.random.chain$Mean.DC))),
                     "Perm.Lower.Obs"=c(nrow(df.perm.lev.random.chain[(as.numeric(df.perm.lev.random.chain$Mean.LV)<obs.lev),]),
                                        nrow(df.perm.dc.random.chain[(as.numeric(df.perm.dc.random.chain$Mean.DC)>obs.dc),])),
                     "Pval"=c(pval.lev,pval.dc))
  return(df.out)
}
get.aa.features<-function(cl) {
  mw.res<-mw(cl)
  tiny.res<-aaComp(cl)[[1]]["Tiny","Number"]
  small.res<-aaComp(cl)[[1]]["Small","Number"]
  aliphatic.res<-aaComp(cl)[[1]]["Aliphatic","Number"]
  aromatic.res<-aaComp(cl)[[1]]["Aromatic","Number"]
  nonpolar.res<-aaComp(cl)[[1]]["NonPolar","Number"]
  polar.res<-aaComp(cl)[[1]]["Polar","Number"]
  charged.res<-aaComp(cl)[[1]]["Charged","Number"]
  basic.res<-aaComp(cl)[[1]]["Basic","Number"]
  acidic.res<-aaComp(cl)[[1]]["Acidic","Number"]
  netcharge.res<-charge(cl, pH=7, pKscale="EMBOSS")
  pI.res<-pI(cl, pKscale="EMBOSS")
  aliphI.res<-aIndex(cl)
  insta.res<-instaIndex(cl)
  boman.res<-boman(cl)
  hydroI.res<-hydrophobicity(cl)
  hmom.res<-hmoment(cl)
  vect.out<-c(mw.res,tiny.res,small.res,aliphatic.res,aromatic.res,nonpolar.res,polar.res,charged.res,basic.res,acidic.res,netcharge.res,pI.res,aliphI.res,insta.res,boman.res,hydroI.res,hmom.res)
  names(vect.out)<-c("Molecular.Weight","Tiny","Small","Aliphatic","Aromatic","Non.Polar","Polar","Charged","Basic","Acidic","Net.Charge","Isoelectric.Point","Aliphatic.Index","Instability.Index","PPI.Index","Hydrophobicity.Index","HydrophobicMoment.Index")
  return(vect.out)
}
get.perm.features<-function(set.n, chain.clones.all, c, get.aa.features, nperm) {
  count.low<-0;  count.hi<-0;
  df.perm.all<-data.frame("Chain"="TEST", "Perm"=0,
                          "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                          "Aliphatic"=0, "Aromatic"=0,
                          "Non.Polar"=0, "Polar"=0, "Charged"=0,
                          "Basic"=0, "Acidic"=0,
                          "Net.Charge"=0, "Isoelectric.Point"=0,
                          "Aliphatic.Index"=0, "Instability.Index"=0,
                          "PPI.Index"=0, "Hydrophobicity.Index"=0,
                          "HydrophobicMoment.Index"=0, stringsAsFactors=F)
  for (i in 1:nperm) {
    print(paste0(c,"-",i))
    perm.clones<-sample(chain.clones.all, size=set.n, replace=F)
    df.perm.clones<-data.frame("Chain"="TEST", "Perm"=0,
                               "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                               "Aliphatic"=0, "Aromatic"=0,
                               "Non.Polar"=0, "Polar"=0, "Charged"=0,
                               "Basic"=0, "Acidic"=0,
                               "Net.Charge"=0, "Isoelectric.Point"=0,
                               "Aliphatic.Index"=0, "Instability.Index"=0,
                               "PPI.Index"=0, "Hydrophobicity.Index"=0,
                               "HydrophobicMoment.Index"=0, stringsAsFactors=F)
    count.perm=0
    for (cl in perm.clones) {
      count.perm=count.perm+1
      cl.perm.prop<-get.aa.features(cl)
      df.perm.clones[count.perm,]<-c(c,i,cl.perm.prop)
    }
    cols.interest.perm<-sapply(df.perm.clones[,c(3:ncol(df.perm.clones))], as.numeric)
    perm.prop<-colMeans(cols.interest.perm)
    df.perm.all[i,]<-c(c,i,perm.prop)
  }
  perm.prop<-df.perm.all
  return(perm.prop)
}
get.final.out<-function(bc.obs, bc.random, nperm){
  aa.propierties<-names(bc.obs)[3:length(bc.obs)]
  df.out<-data.frame("Biochemical.Feature"="TEST","Obs.Value"="TEST","Perm.Value"="TEST","PermValue.LT.ObsValue"="TEST","PermValue.GT.ObsValue"="TEST","Pval.PermValue.LT.Obs.Value"="TEST","Pval.PermValue.GT.ObsValue"="TEST", stringsAsFactors=F)
  count=0
  for (aap in aa.propierties) {
    obs.value<-as.numeric(bc.obs[aap])
    perm.value<-as.numeric(bc.random[,aap])
    pval1<-(length(which(perm.value<obs.value))+1)/(nperm+1)
    pval2<-(length(which(perm.value>=obs.value))+1)/(nperm+1)
    count=count+1
    df.out[count,]<-c(aap, unname(obs.value), mean(perm.value), length(which(perm.value<obs.value)), length(which(perm.value>obs.value)), pval1, pval2)
  }
  return(df.out)
}

df.res<-readRDS("../../Output_Data/NoisET_Summary_ExpCont.rds")
chains<-c("TRA","TRD","TRG","TRB","IGH","IGL","IGK")

samples2keep<-c(as.character(df.donor[,2]), as.character(df.donor[,3]))

for (c in chains) {
  df.res.chain<-df.res[df.res$Chain==c,]
  clones.cont<-unlist(strsplit(paste(df.res.chain[,"Exp.Sig.Clones"], collapse=';'),";"))
  clones.cont<-clones.cont[!(clones.cont %in% "NA")]
  chain.clones.cont.n<-length(unique(clones.cont))
  chain.clones.cont.max<-max(table(clones.cont))
  chain.clones.cont.unique<-unique(clones.cont)[!(unique(clones.cont) %in% "NA")]
  tbl.freq<-table(clones.cont)
  df.freq.sorted<-as.data.frame(tbl.freq[order(tbl.freq,decreasing = T)])
  df.freq.sorted<-df.freq.sorted[df.freq.sorted$clones.cont!="NA",]
  
  clones.total<-nrow(df.freq.sorted)
  clones.in50perc<-nrow(df.freq.sorted[df.freq.sorted$Freq>77*0.5,])
  clones.in40perc<-nrow(df.freq.sorted[df.freq.sorted$Freq>77*0.4,])
  clones.in30perc<-nrow(df.freq.sorted[df.freq.sorted$Freq>77*0.3,])
  clones.in20perc<-nrow(df.freq.sorted[df.freq.sorted$Freq>77*0.2,])
  
  p.freq<-ggplot(df.freq.sorted, aes(x=clones.cont, y=Freq, fill=clones.cont)) + geom_bar(stat="identity", position=position_dodge()) +
    labs(x=paste0("\n"), y="\nSamples (N)\n") + theme_classic() + ylim(0,max(df.freq.sorted$Freq)*1.5) +
    theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1), legend.position="none") +
    annotate("text", x=nrow(df.freq.sorted)*0.7, y=max(df.freq.sorted$Freq)*1.4, label=paste0("Total Clones (n) = ",clones.total)) +
    annotate("text", x=nrow(df.freq.sorted)*0.7, y=max(df.freq.sorted$Freq)*1.36, label=paste0("Clones in >50% samples (n) = ",clones.in50perc)) +
    annotate("text", x=nrow(df.freq.sorted)*0.7, y=max(df.freq.sorted$Freq)*1.32, label=paste0("Clones in >40% samples (n) = ",clones.in40perc)) +
    annotate("text", x=nrow(df.freq.sorted)*0.7, y=max(df.freq.sorted$Freq)*1.28, label=paste0("Clones in >30% samples (n) = ",clones.in30perc)) +
    annotate("text", x=nrow(df.freq.sorted)*0.7, y=max(df.freq.sorted$Freq)*1.24, label=paste0("Clones in >20% samples (n) = ",clones.in20perc))
  
  jpeg(paste0("../../Output_Data/NoisET_Characterization_",type,"/SampleDist_",c,".jpeg"), res=300, height=3500, width=4500)
  plot(p.freq)
  dev.off()
  
  irep.chain<-readRDS(paste0("../Clones_By_CDR3aa/",c,".rds"))
  meta_col<-dplyr::filter(irep.chain$meta, IMID=="RA");
  meta_col1<-meta_col[(meta_col$Sample %in% samples2keep),]
  samples_col<-meta_col1$Sample
  data_col<-irep.chain$data[match(samples_col,names(irep.chain$data))];
  irep_chain<-list(data_col,meta_col1);  names(irep_chain)<-c("data","meta")
  
  chain.clones.all<-as.data.frame(immunarch::pubRep(irep_chain$data, "aa", .verbose=F))[,"CDR3.aa"]
  
  # Leveinshtein & Degree Centrality
  print(paste0(c," - Leveinshtein & DC"))
  nperm<-10000
  dc.lev.out<-get.lev.dc.pval(chain.clones.cont.unique, chain.clones.all, c, type, nperm)
  saveRDS(dc.lev.out, paste0("../../Output_Data/NoisET_Characterization_",type,"/DC_Lev_",c,".rds"))
  
  # Biochemical Properties
  print(paste0(c," - Biochemical Prop"))
  nperm<-1000
  bc.random<-get.perm.features(length(chain.clones.cont.unique), chain.clones.all, c, get.aa.features, nperm)
  bc.obs<-data.frame("Chain"="TEST", "Perm"=0,
                     "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                     "Aliphatic"=0, "Aromatic"=0,
                     "Non.Polar"=0, "Polar"=0, "Charged"=0,
                     "Basic"=0, "Acidic"=0,
                     "Net.Charge"=0, "Isoelectric.Point"=0,
                     "Aliphatic.Index"=0, "Instability.Index"=0,
                     "PPI.Index"=0, "Hydrophobicity.Index"=0,
                     "HydrophobicMoment.Index"=0, stringsAsFactors=F)
  count=0
  for (cl in chain.clones.cont.unique) {
    count=count+1
    cl.prop<-get.aa.features(cl)
    bc.obs[count,]<-c(c,cl,cl.prop)
  }
  cols.interest<-sapply(bc.obs[,c(3:ncol(bc.obs))], as.numeric)
  bc.obs<-c(c,"Obs",colMeans(cols.interest))
  names(bc.obs)[1:2]<-c("Chain","Perm")
  bc.final.out<-get.final.out(bc.obs, bc.random, nperm)
  saveRDS(bc.final.out, paste0("../../Output_Data/NoisET_Characterization_",type,"/BiochemicalProp_",c,".rds"))
  
  # Amino Acid Motif
  print(paste0(c," - Amino Acid Motif"))
  count=0
  count.name=0
  out.data<-c()
  for (clone in chain.clones.cont.unique) {
    count=count+1
    count.name=count.name+1
    out.data[count]<-paste0(">Clone",count.name,"|Clone",count.name,"|Clone",count.name)
    out.data[count+1]<-clone
    count=count+1
  }
  df.out.data<-as.data.frame(out.data)
  write.table(df.out.data, paste0("../../Output_Data/NoisET_Characterization_",type,"/Seq_",c,".txt"), quote=F, col.names=F, row.names=F)
  seqs<-readAAStringSet(paste0("../../Output_Data/NoisET_Characterization_",type,"/Seq_",c,".txt"))
  seqs.algn<-msa(seqs, "ClustalW")
  conMat<-consensusMatrix(seqs.algn)
  jpeg(paste0("../../Output_Data/NoisET_Characterization_",type,"/AAMotif_",c,".jpeg"), res=300, height=2500, width=3500)
  plot(ggseqlogo(conMat, method='prob'))
  dev.off()
}

#### Notes ####

# the same code was used to characterize: 1) contracted clones, 2) expanded clones detected in responders, 3) expanded clones detected in non-responders, 4) contracted clones detected in responders, and 5) contracted clones detected in non-responders.

########  Getting Sequence Similarity Network with Individual Information  ######## 
#### Creating Clone-Individual Dataframe ####

df.donor<-read.table("../Donor2RNA_Codes.txt", header=T)

chains<-c("TRA","TRB","TRD","TRG","IGH","IGL","IGK")

for (c in chains) {
  df.out.chain<-list()
  for (i in 1:nrow(df.donor)) {
    donor<-as.character(df.donor$Donor[i])
    wk0<-as.character(df.donor$RNA1[i])
    wk12<-as.character(df.donor$RNA2[i])
    print(paste0(c," - ",donor))
    data.wk0<-read_tsv(paste0(c,"_",wk0,"_S.tsv"))
    data.wk12<-read_tsv(paste0(c,"_",wk12,"_S.tsv"))
    clones<-unique(c(data.wk0$CDR3.aa, data.wk12$CDR3.aa))
    df.ind<-data.frame("Donor"=rep(donor,length(clones)), "Clone"=clones, stringsAsFactors=F)
    df.out.chain[[i]]<-df.ind
  }
  names(df.out.chain)<-as.character(df.donor$Donor)
  df.out.chain<-do.call(rbind, df.out.chain)
  saveRDS(df.out.chain, paste0("../../Output_Data/NoisET_Donor_Clone_Dict_",c,".rds"))
}

#### Expanded ####

df.res<-readRDS("../../Output_Data/NoisET_Summary_ExpCont.rds")

type="Expanded"
samples2keep<-c(as.character(df.donor[,2]), as.character(df.donor[,3]))

chains<-c("TRA","TRB","TRD","TRG","IGH","IGL","IGK")
for (c in chains) {
  print(c)
  df.res.chain<-df.res[df.res$Chain==c,]
  clones.cont<-unlist(strsplit(paste(df.res.chain[,"Exp.Sig.Clones"], collapse=';'),";"))
  clones.cont<-unique(clones.cont[!(clones.cont %in% "NA")])
  
  if (length(clones.cont)>0) {
    graph<-mutation.network(clones.cont, .method="lev", .max.errors=1)
    deg<-degree(graph, mode="all")
    
    if (c=="IGH" | c=="IGL" | c=="IGK") {fact<-1.5} else {fact<-5}
    
    V(graph)$size<-deg*fact;    V(graph)$frame.color<-"white";   V(graph)$color<-"orange";  V(graph)$label.cex<-0.3;   E(graph)$arrow.mode<-0;
    
    dict.data<-readRDS(paste0("../../Output_Data/NoisET_Donor_Clone_Dict_",c,".rds"))
    values<-list()
    count=0
    for (lab in V(graph)$label) {
      ind<-as.numeric(unique(dict.data[dict.data$Clone==lab,"Donor.Num"]))
      count=count+1
      values[[count]]<-ind
    }
    names(values)<-V(graph)$label
    
    jpeg(paste0("../../Output_Data/NoisET_Characterization_",type,"/SimilarityNetwork",c,"_ByInd.jpeg"), res=300, height=3000, width=3000)
    plot(graph, layout=layout_with_kk, vertex.frame.color="darkblue", vertex.label=NA, vertex.shape="pie", vertex.pie=values,
         main=paste0("\n",type," Cloness (DC=", round(mean(degree(graph)),2),")"))
    dev.off()
  }

}

#### Notes ####

# the same code was used to draw the sequence similarity network for: 1) contracted clones, 2) expanded clones detected in responders, 3) expanded clones detected in non-responders, 4) contracted clones detected in responders, and 5) contracted clones detected in non-responders.





########  ALL: Expanded vs Contracted (weighted by abundance)  ########

samples2keep<-c(as.character(df.donor[,2]), as.character(df.donor[,3]))
donor.pheno<-df.pheno$Donor

out.expVScont.all<-data.frame("Chain"="NA","Biochemical.Prop"="NA","Beta"=0,"Pval"=0, stringsAsFactors=F)
chains<-c("TRG","TRB","TRD","TRA","IGH","IGL","IGK")
counter=0
for (chain in chains) {
  df.res.chain<-df.res[df.res$Chain==chain,]
  clones.exp<-unlist(strsplit(paste(df.res.chain[(df.res.chain$Donor %in% donor.pheno),"Exp.Sig.Clones"], collapse=';'),";"))
  clones.exp<-clones.exp[!(clones.exp %in% "NA")]
  clones.cont<-unlist(strsplit(paste(df.res.chain[(df.res.chain$Donor %in% donor.pheno),"Cont.Sig.Clones"], collapse=';'),";"))
  clones.cont<-clones.cont[!(clones.cont %in% "NA")]

  df.exp<-list(); df.cont<-list()
  exp.fs=grep(paste0("NoisET_Detection_",chain), list.files("P05_Expanded/"), value=T)
  cont.fs=grep(paste0("NoisET_Detection_",chain), list.files("P05_Contracted/"), value=T)
  count=0;  for (j in 1:length(exp.fs)) {f1<-read.csv(paste0("P05_Expanded/",exp.fs[j]),sep="\t"); count=count+1; df.exp[[count]]<-f1}
  count=0;  for (j in 1:length(cont.fs)) {f1<-read.csv(paste0("P05_Contracted/",cont.fs[j]),sep="\t"); count=count+1; df.cont[[count]]<-f1}
  chain.cont<-do.call(rbind,df.cont)
  chain.exp<-do.call(rbind,df.exp)
  chain.cont.sig<-chain.cont[chain.cont$X.1.P.s.0..<0.05,]
  chain.exp.sig<-chain.exp[chain.exp$X.1.P.s.0..<0.05,]
  chain.exp.sig$Change<-abs(chain.exp.sig$X.n_2.-chain.exp.sig$X.n_1.)
  chain.cont.sig$Change<-abs(chain.cont.sig$X.n_2.-chain.cont.sig$X.n_1.)
  chain.exp.sig.filt<-chain.exp.sig[,c("CDR3_AA","Change")]
  chain.cont.sig.filt<-chain.cont.sig[,c("CDR3_AA","Change")]
  chain.exp.sig.filt.group<-as.data.frame(chain.exp.sig.filt %>% group_by(CDR3_AA) %>% dplyr::summarise(Change=mean(Change)))
  chain.cont.sig.filt.group<-as.data.frame(chain.cont.sig.filt %>% group_by(CDR3_AA) %>% dplyr::summarise(Change=mean(Change)))
  
  df.clones.cont<-data.frame("Clones"=clones.cont,"CDR3_AA"=clones.cont,stringsAsFactors=F)
  df.clones.exp<-data.frame("Clones"=clones.exp,"CDR3_AA"=clones.exp,stringsAsFactors=F)
  
  dfm.exp<-merge(chain.exp.sig.filt.group,df.clones.exp,by="CDR3_AA")[,c(3,2)]
  dfm.cont<-merge(chain.cont.sig.filt.group,df.clones.cont,by="CDR3_AA")[,c(3,2)]
  
  clones.exp.weighted<-data.frame("Clone"=rep(dfm.exp$Clones,dfm.exp$Change), "Type"=rep("Expanded",length(rep(dfm.exp$Clones,dfm.exp$Change))))
  clones.cont.weighted<-data.frame("Clone"=rep(dfm.cont$Clones,dfm.cont$Change), "Type"=rep("Contracted",length(rep(dfm.cont$Clones,dfm.cont$Change))))
  
  df.all<-rbind(clones.exp.weighted,clones.cont.weighted)
  df.all$PPI<-0; df.all$Polar<-0; df.all$MW<-0; df.all$IsoP<-0; df.all$Insta<-0; df.all$Hydro<-0; df.all$Charged<-0; df.all$Basic<-0; df.all$Aromatic<-0; df.all$Aliphatic<-0; df.all$Acidic<-0; df.all$Length<-0;
  if (length(clones.cont)>0 & length(clones.exp)>0) {
    for (i in 1:nrow(df.all)) {
      print(paste0(chain," ",i,"/",nrow(df.all)))
      cl<-as.character(df.all$Clone[i])
      type<-as.character(df.all$Type[i])
      molw<-mw(cl)
      polar<-aaComp(cl)[[1]]["Polar","Number"]
      aliphatic<-aaComp(cl)[[1]]["Aliphatic","Number"]
      aromatic<-aaComp(cl)[[1]]["Aromatic","Number"]
      charged<-aaComp(cl)[[1]]["Charged","Number"]
      basic<-aaComp(cl)[[1]]["Basic","Number"]
      acidic<-aaComp(cl)[[1]]["Acidic","Number"]
      hydro<-hydrophobicity(cl)
      iso<-pI(cl, pKscale="EMBOSS")
      ppi<-boman(cl)
      insta<-instaIndex(cl)
      len<-nchar(cl)
      df.all[i,]<-c(cl,type,ppi,polar,molw,iso,insta,hydro,charged,basic,aromatic,aliphatic,acidic,len)
    }
    bp<-colnames(df.all)[3:13]
    for (b in bp) {
      print(paste0(chain,"-",b))
      counter=counter+1
      data2use<-df.all[,c("Clone","Type","Length",b)]
      data2use$Type<-as.character(data2use$Type)
      data2use$Type[data2use$Type=="Expanded"]<-0
      data2use$Type[data2use$Type=="Contracted"]<-1
      data2use$Type<-as.numeric(data2use$Type)
      colnames(data2use)[4]<-"BP"
      data2use$BP<-as.numeric(data2use$BP)
      data2use$Length<-as.numeric(data2use$Length)
      model<-glm(Type~BP+Length, data=data2use, family="binomial")
      pval<-coefficients(summary(model))["BP","Pr(>|z|)"]
      beta<-coefficients(summary(model))["BP","Estimate"] # beta<0 means lower in contracted
      out.expVScont.all[counter,]<-c(chain,b,beta,pval)
    }
  } else {counter=counter+1; out.expVScont.all[counter,]<-c(chain,NA,NA,NA)}
}
saveRDS(out.expVScont.all,"../../Output_Data/NoisET_Characterization_ALL_ExpVSCont.rds")

#### Notes ####
# The same code was used for the following comparisons: 1) Expanded vs Contracted in responders, and 2) Expanded vs Contracted in non-responders.
########  Summarizing Results V2 for graphical representation  ####

# Loading NoisET output
df.res<-readRDS("../../Output_Data/NoisET_Summary_ExpCont.rds")
df.res.bcr<-df.res[(df.res$Chain=="IGH" | df.res$Chain=="IGL" | df.res$Chain=="IGK"),]
length(unique(unlist(strsplit(df.res.bcr$Cont.Sig.Clones,";"))))
df.res.tcr<-df.res[(df.res$Chain!="IGH" & df.res$Chain!="IGL" & df.res$Chain!="IGK"),]
length(unique(unlist(strsplit(df.res.tcr$Cont.Sig.Clones,";"))))

chains<-c("TRA","TRB","TRD","TRG","IGH","IGL","IGK")
df.clones<-data.frame("Chain"="test","Clones"=0, stringsAsFactors=F)
count=0
for (c in chains) {
  print(c)
  ra.data<-readRDS(paste0("../Clones_By_CDR3aa/",c,".rds"))
  ra.ltgd<-ra.data$data[names(ra.data$data) %in% ra.all]
  df.lgtd<-do.call(rbind,ra.ltgd)
  n.clones<-length(unique(df.lgtd$CDR3.aa,","))
  count=count+1
  df.clones[count,]<-c(c,n.clones)
}

# Plots
df.chain2plot<-data.frame("Chain"=c("IGH","IGK","IGL","TRA","TRB","TRD","TRG"))
df.chain2plot$Cont.Sig[1]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="IGH","Cont.Sig.Clones"],";"))))
df.chain2plot$Cont.Sig[2]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="IGK","Cont.Sig.Clones"],";"))))  
df.chain2plot$Cont.Sig[3]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="IGL","Cont.Sig.Clones"],";"))))  
df.chain2plot$Cont.Sig[4]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="TRA","Cont.Sig.Clones"],";"))))
df.chain2plot$Cont.Sig[5]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="TRB","Cont.Sig.Clones"],";"))))
df.chain2plot$Cont.Sig[6]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="TRD","Cont.Sig.Clones"],";"))))
df.chain2plot$Cont.Sig[7]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="TRG","Cont.Sig.Clones"],";"))))
df.chain2plot$Receptor<-c("BCR","BCR","BCR","TCR","TCR","TCR","TCR")
df.chain2plot$Cont.Sig<-as.numeric(df.chain2plot$Cont.Sig)
df.chain2plot$Chain<-as.factor(df.chain2plot$Chain)
df.chain2plot$Receptor<-as.factor(df.chain2plot$Receptor)

df.chain2plot$Exp.Sig[1]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="IGH","Exp.Sig.Clones"],";"))))
df.chain2plot$Exp.Sig[2]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="IGK","Exp.Sig.Clones"],";"))))  
df.chain2plot$Exp.Sig[3]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="IGL","Exp.Sig.Clones"],";"))))  
df.chain2plot$Exp.Sig[4]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="TRA","Exp.Sig.Clones"],";"))))
df.chain2plot$Exp.Sig[5]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="TRB","Exp.Sig.Clones"],";"))))
df.chain2plot$Exp.Sig[6]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="TRD","Exp.Sig.Clones"],";"))))
df.chain2plot$Exp.Sig[7]<-length(unique(unlist(strsplit(df.res[df.res$Chain=="TRG","Exp.Sig.Clones"],";"))))

write.csv(df.chain2plot[,c(3,1,2,4)], "../../Output_Data/Table_ExpCont_Summary_part1.csv", row.names=F, quote=F)

df.chain2plot$Mean.Delta.Cont<-0
df.chain2plot$Mean.Delta.Exp<-0
df.chain2plot$Disappear.N<-0
df.chain2plot$Appear.N<-0

for (i in 1:nrow(df.chain2plot)) {
  print(i)
  chain<-as.character(df.chain2plot$Chain[i])
  df.exp<-list()
  df.cont<-list()
  exp.fs=grep(paste0("NoisET_Detection_",chain), list.files("P05_Expanded/"), value=T)
  cont.fs=grep(paste0("NoisET_Detection_",chain), list.files("P05_Contracted/"), value=T)
  count=0
  for (j in 1:length(exp.fs)) {f1<-read.csv(paste0("P05_Expanded/",exp.fs[j]),sep="\t"); count=count+1; df.exp[[count]]<-f1}
  count=0
  for (j in 1:length(cont.fs)) {f1<-read.csv(paste0("P05_Contracted/",cont.fs[j]),sep="\t"); count=count+1; df.cont[[count]]<-f1}
  chain.cont<-do.call(rbind,df.cont)
  chain.exp<-do.call(rbind,df.exp)
  chain.cont.sig<-chain.cont[chain.cont$X.1.P.s.0..<0.05,]
  chain.exp.sig<-chain.exp[chain.exp$X.1.P.s.0..<0.05,]
  chain.exp.sig$Change<-chain.exp.sig$X.f_2.-chain.exp.sig$X.f_1.
  chain.cont.sig$Change<-chain.cont.sig$X.f_2.-chain.cont.sig$X.f_1.
  df.chain2plot$Mean.Delta.Cont[i]<-mean(chain.cont.sig$Change)
  df.chain2plot$Mean.Delta.Exp[i]<-mean(chain.exp.sig$Change)
  n.appear.exp<-length(which(chain.exp.sig$X.f_1.==0))
  n.disappear.cont<-length(which(chain.cont.sig$X.f_1.==0))
  df.chain2plot$Disappear.N[i]<-n.disappear.cont
  df.chain2plot$Appear.N[i]<-n.appear.exp
}
df.chain2plot$Delta.Magnitude<-round(df.chain2plot$Mean.Delta.Cont/df.chain2plot$Mean.Delta.Exp,2)
df.chain.final<-merge(df.chain2plot,df.clones,by="Chain")
df.chain.final<-df.chain.final[,c(3,1,10,2,4,5,6,9)]
colnames(df.chain.final)<-c("Receptor","Chain","Total.Clones","Contracted.Clones.Sig","Expanded.Clones.Sig","Delta.Contraction.Mean","Delta.Expansion.Mean","Delta.Magnitude.Contracted/Expanded")
write.csv(df.chain.final, "../../Output_Data/Table_ExpCont_Summary.csv", row.names=F, quote=F)

df.chain2plot2use<-df.chain2plot[,c(1,3,2,4)]
colnames(df.chain2plot2use)[3:4]<-c("Contraction","Expansion")
df.chain2plot.long<-gather(df.chain2plot2use, condition, value, Contraction:Expansion, factor_key=T)

p1<-ggplot(df.chain2plot.long, aes(fill=Chain, y=value, x=condition)) + coord_flip() +
  geom_bar(position="stack", stat="identity") + theme_classic() + theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
  scale_fill_manual(values=c("honeydew","honeydew2","honeydew4",brewer.pal(4, "Reds"))) + ggtitle("Studying 4 species..") +
  labs(y="\nClones modulated by TNF inhibition (N, P<0.05)\n", x="\n", title="\n") +
  annotate(geom="text", x=2.55, y=58900, size=4, label="P<2.2e-16", color="black")



corresp<-unique(df.res$Donor)
names(corresp)<-paste("RA_sample_",seq(1:length(corresp)))
corresp<-data.frame("Donor"=as.character(corresp), "CODE"=names(corresp))
df.res.cont<-merge(corresp,df.res[,c(1,2,7,9)],by="Donor")
write.csv(df.res.cont[,c(2,3,4,5)], "Manuscript_Data/Table_Contraction.csv", row.names=F)

df.res.exp<-merge(corresp,df.res[,c(1,2,6,8)],by="Donor")
write.csv(df.res.exp[,c(2,3,4,5)], "Manuscript_Data/Table_Expansion.csv", row.names=F)

df2plot<-data.frame("Chain"=c("IGH","IGK","IGL","TRA","TRB","TRD","TRG"),"Ratio"=c(1.35,1.32,1.17,1.30,1.69,0.34,1.33))
p3<-ggplot(df2plot, aes(fill=Chain, y=Ratio, x=Chain)) +
  geom_bar(position="stack", stat="identity") + theme_classic() + theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)) +
  scale_fill_manual(values=c("honeydew","honeydew2","honeydew4",brewer.pal(4, "Reds"))) + ggtitle("Studying 4 species..") +
  geom_hline(yintercept=1, linetype="dashed", color = "black") +
  labs(y="\nContraction/Expansion Ratio\n", x="\n", title="\n")

jpeg(paste0("../../Output_Data/Summary_Ratio.jpeg"), width=2000, height=2600, res=300)
plot(p3)
dev.off()

# Networks

igk.cont<-unique(unlist(strsplit(df.res[df.res$Chain=="IGK","Cont.Sig.Clones"],";")))
graph<-mutation.network(igk.cont, .method="lev", .max.errors=1)
deg<-degree(graph, mode="all")
V(graph)$size<-deg*0.6
V(graph)$frame.color<-"white"
V(graph)$color<-"honeydew2"
V(graph)$label.cex<-0.12
V(graph)$label.color<-"black"
E(graph)$arrow.mode<-0
E(graph)$arrow.width<-0.5
jpeg("../../Output_Data/AdditionalFigure1.jpeg", res=350, height=4000, width=4000)
plot(graph, layout=layout_with_kk)
dev.off()

igk.cont<-unique(unlist(strsplit(df.res[df.res$Chain=="IGL","Cont.Sig.Clones"],";")))
graph<-mutation.network(igk.cont, .method="lev", .max.errors=1)
deg<-degree(graph, mode="all")
V(graph)$size<-deg*0.6
V(graph)$frame.color<-"white"
V(graph)$color<-"honeydew4"
V(graph)$label.cex<-0.12
V(graph)$label.color<-"black"
E(graph)$arrow.mode<-0
E(graph)$arrow.width<-0.5
jpeg("../../Output_Data/AdditionalFigure2.jpeg", res=350, height=4000, width=4000)
plot(graph, layout=layout_with_kk)
dev.off()

igk.cont<-unique(unlist(strsplit(df.res[df.res$Chain=="IGH","Cont.Sig.Clones"],";")))
graph<-mutation.network(igk.cont, .method="lev", .max.errors=1)
deg<-degree(graph, mode="all")
V(graph)$size<-deg*0.6
V(graph)$frame.color<-"white"
V(graph)$color<-"honeydew"
V(graph)$label.cex<-0.12
V(graph)$label.color<-"black"
E(graph)$arrow.mode<-0
E(graph)$arrow.width<-0.5
jpeg("../../Output_Data/AdditionalFigure3.jpeg", res=300, height=4000, width=4000)
plot(graph, layout=layout_with_kk)
dev.off()

out.sim.exp<-list()
out.sim.cont<-list()
count=0
for (ch in chains) {
  count=count+1
  exp.data<-readRDS(paste0("../../Output_Data/NoisET_Characterization_Expanded/DC_Lev_",ch,".rds"))
  out.sim.exp[[count]]<-exp.data
  cont.data<-readRDS(paste0("../../Output_Data/NoisET_Characterization_Contracted/DC_Lev_",ch,".rds"))
  out.sim.cont[[count]]<-cont.data
}
out.final.exp<-do.call(rbind,out.sim.exp)
out.final.cont<-do.call(rbind,out.sim.cont)
write.csv(out.final.exp, "../../Output_Data/NoisET_AdditionalTable3.csv", row.names=F, quote=F)
write.csv(out.final.cont, "../../Output_Data/NoisET_AdditionalTable4.csv", row.names=F, quote=F)


out.sim.exp.resp<-list()
out.sim.cont.resp<-list()
out.sim.exp.noresp<-list()
out.sim.cont.noresp<-list()
count=0
for (ch in chains) {
  count=count+1
  exp.data.resp<-readRDS(paste0("../../Output_Data/NoisET_Characterization_Expanded_RespGM/DC_Lev_",ch,".rds"))
  out.sim.exp.resp[[count]]<-exp.data.resp
  cont.data.resp<-readRDS(paste0("../../Output_Data/NoisET_Characterization_Contracted_RespGM/DC_Lev_",ch,".rds"))
  out.sim.cont.resp[[count]]<-cont.data.resp
  exp.data.noresp<-readRDS(paste0("../../Output_Data/NoisET_Characterization_Expanded_NoResp/DC_Lev_",ch,".rds"))
  out.sim.exp.noresp[[count]]<-exp.data.noresp
  cont.data.noresp<-readRDS(paste0("../../Output_Data/NoisET_Characterization_Contracted_NoResp/DC_Lev_",ch,".rds"))
  out.sim.cont.noresp[[count]]<-cont.data.noresp
}
out.final.exp.resp<-do.call(rbind,out.sim.exp.resp)
out.final.cont.resp<-do.call(rbind,out.sim.cont.resp)
out.final.exp.noresp<-do.call(rbind,out.sim.exp.noresp)
out.final.cont.noresp<-do.call(rbind,out.sim.cont.noresp)
write.csv(out.final.exp.resp, "../../Output_Data/NoisET_AdditionalTable5.csv", row.names=F, quote=F)
write.csv(out.final.cont.resp, "../../Output_Data/NoisET_AdditionalTable6.csv", row.names=F, quote=F)
write.csv(out.final.exp.noresp, "../../Output_Data/NoisET_AdditionalTable7.csv", row.names=F, quote=F)
write.csv(out.final.cont.noresp, "../../Output_Data/NoisET_AdditionalTable8.csv", row.names=F, quote=F)


# Graphical representation biochemical properties

type<-c("Expanded","Contracted")
chains<-c("TRA","TRB","TRD","TRG","IGH","IGL","IGK")
df.all<-list()
count=0
for (t in type) {
  for (c in chains) {
    count=count+1
    file.name<-paste0("../../Output_Data/NoisET_Characterization_",t,"/BiochemicalProp_",c,".rds")
    res.data<-readRDS(file.name)
    print(paste0(t," - ",c," - ",max(as.numeric(c(res.data[,4],res.data[,5])))))
    res.data$Chain<-rep(c,nrow(res.data))
    res.data$TNFeff<-rep(t,nrow(res.data))
    res.data.depl<-res.data[,c(1,6,8,9)]
    res.data.enrich<-res.data[,c(1,7,8,9)]
    res.data.depl$BiochEff<-rep("Depletion",nrow(res.data.depl))
    res.data.enrich$BiochEff<-rep("Enrichment",nrow(res.data.enrich))
    colnames(res.data.depl)[2]<-"Pval"
    colnames(res.data.enrich)[2]<-"Pval"
    res.data.depl$Pval<-(log10(as.numeric(res.data.depl$Pval)))
    res.data.enrich$Pval<-(-log10(as.numeric(res.data.enrich$Pval)))
    res.data.enrich$key<-paste(res.data.enrich$Biochemical.Feature, res.data.enrich$Chain, res.data.enrich$TNFeff)
    res.data.depl$key<-paste(res.data.depl$Biochemical.Feature, res.data.depl$Chain, res.data.depl$TNFeff)
    res2use<-merge(res.data.depl,res.data.enrich, by="key")
    res2use$enrich.score<-res2use$Pval.y+res2use$Pval.x
    res2use<-res2use[,c(2,4,5,12)]
    colnames(res2use)<-c("Biochemical.Feature","Chain","TNFeff","Biochemical.Enrichment.Score")
    df.all[[count]]<-res2use
  }
}
data2plot<-do.call(rbind, df.all)
bp2remove<-c("Tiny","Aliphatic.Index","Non.Polar","HydrophobicMoment.Index","Small","Net.Charge")
data2plot<-data2plot[!(data2plot$Biochemical.Feature %in% bp2remove),]
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="Hydrophobicity.Index"]<-"Hydrophobic"
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="Instability.Index"]<-"Instability"
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="PPI.Index"]<-"Protein-Protein Interaction Index"
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="Molecular.Weight"]<-"Molecular weight"
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="Isoelectric.Point"]<-"Isoelectric point"

saveRDS(data2plot[data2plot$TNFeff=="Contracted",],"../../Output_Data/NoisET_Table_BiochProf_BES_ContVSRandom.rds")
saveRDS(data2plot[data2plot$TNFeff=="Expanded",],"../../Output_Data/NoisET_Table_BiochProf_BES_ExpVSRandom.rds")

# only plot if BEP (Biochemical Enrichment Score is > -log10(0.05))

data2plot$Biochemical.Enrichment.Score[(data2plot$Biochemical.Enrichment.Score<1.3 & data2plot$Biochemical.Enrichment.Score>-1.3)]<-NA
colnames(data2plot)[4]<-"BES"

p.all.chain.comp<-ggplot(data=data2plot, aes(x=factor(Chain), y=factor(Biochemical.Feature), size=abs(BES), fill=BES)) +
  geom_point(alpha=1, shape=21, color="black") +
  scale_size(range=c(.1, 8), name="Significance") +
  facet_wrap(~TNFeff) + labs(y="", x="") +
  scale_fill_gradientn(colours = rev(colorspace::diverge_hcl(15))) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face="bold", size=16), axis.text.x=element_text(vjust=0.5, size=13), strip.text=element_text(size=14, face="bold"))

p.all.expcont.comp<-ggplot(data=data2plot, aes(x=factor(TNFeff), y=factor(Biochemical.Feature), size=abs(BES), fill=BES)) +
  geom_point(alpha=1, shape=21, color="black") +
  scale_size(range=c(.1, 8), name="Significance") +
  facet_wrap(~Chain) + labs(y="", x="") +
  scale_fill_gradientn(colours = rev(colorspace::diverge_hcl(15))) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face="bold", size=16), axis.text.x=element_text(vjust=0.5, size=13), strip.text=element_text(size=14, face="bold"))

data2plot.v2<-data2plot[(data2plot$Chain %in% c("IGH","IGL","IGK") & data2plot$TNFeff=="Contracted"),]
p.cont.int<-ggplot(data=data2plot.v2, aes(x=factor(Chain), y=factor(Biochemical.Feature), size=abs(BES), fill=BES)) +
  geom_point(alpha=1, shape=21, color="black") +
  scale_size(range=c(.1, 10), name="Significance") +
  facet_wrap(~TNFeff) + labs(y="", x="") +
  scale_fill_gradientn(colours = rev(colorspace::diverge_hcl(15))) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face="bold", size=16), axis.text.x=element_text(vjust=0.5, size=13), strip.text=element_text(size=14, face="bold"))

jpeg("../../Output_Data/NoisET_BP1.jpeg", res=300, height=4000, width=8300)
grid.arrange(p.cont.int, p.all.chain.comp, p.all.expcont.comp, ncol=3)
dev.off()





res.data<-readRDS("../../Output_Data/NoisET_Characterization_Contracted/BiochemicalProp_IGK.rds")
res2plot.v1<-res.data[,c(1,6)]
res2plot.v1$Type<-rep("Depletion",nrow(res2plot.v1))
res2plot.v2<-res.data[,c(1,7)]
res2plot.v2$Type<-rep("Enrichment",nrow(res2plot.v1))
colnames(res2plot.v1)[2]<-"Pval"
colnames(res2plot.v2)[2]<-"Pval"
res2plot<-rbind(res2plot.v1,res2plot.v2)
res2plot$Pval<-as.numeric(res2plot$Pval)
res2plot<-res2plot[order(res2plot$Pval),]
res2plot$FDR<-p.adjust(res2plot$Pval, method="fdr")
res2plot$Pval2plot<--log10(res2plot$Pval)
res2plot$Type<-as.factor(res2plot$Type)
res2plot.igk<-res2plot

res.data<-readRDS("../../Output_Data/NoisET_Characterization_Contracted/BiochemicalProp_IGL.rds")
res2plot.v1<-res.data[,c(1,6)]
res2plot.v1$Type<-rep("Depletion",nrow(res2plot.v1))
res2plot.v2<-res.data[,c(1,7)]
res2plot.v2$Type<-rep("Enrichment",nrow(res2plot.v1))
colnames(res2plot.v1)[2]<-"Pval"
colnames(res2plot.v2)[2]<-"Pval"
res2plot<-rbind(res2plot.v1,res2plot.v2)
res2plot$Pval<-as.numeric(res2plot$Pval)
res2plot<-res2plot[order(res2plot$Pval),]
res2plot$FDR<-p.adjust(res2plot$Pval, method="fdr")
res2plot$Pval2plot<--log10(res2plot$Pval)
res2plot$Type<-as.factor(res2plot$Type)
res2plot.igl<-res2plot

res.data<-readRDS("../../Output_Data/NoisET_NoisET/Characterization_Contracted/BiochemicalProp_IGH.rds")
res2plot.v1<-res.data[,c(1,6)]
res2plot.v1$Type<-rep("Depletion",nrow(res2plot.v1))
res2plot.v2<-res.data[,c(1,7)]
res2plot.v2$Type<-rep("Enrichment",nrow(res2plot.v1))
colnames(res2plot.v1)[2]<-"Pval"
colnames(res2plot.v2)[2]<-"Pval"
res2plot<-rbind(res2plot.v1,res2plot.v2)
res2plot$Pval<-as.numeric(res2plot$Pval)
res2plot<-res2plot[order(res2plot$Pval),]
res2plot$FDR<-p.adjust(res2plot$Pval, method="fdr")
res2plot$Pval2plot<--log10(res2plot$Pval)
res2plot$Type<-as.factor(res2plot$Type)
res2plot.igh<-res2plot

res2plot.igl$Chain<-rep("IGL",nrow(res2plot.igl))
res2plot.igk$Chain<-rep("IGK",nrow(res2plot.igk))
res2plot.igh$Chain<-rep("IGH",nrow(res2plot.igh))


res2plot<-rbind(res2plot.igl,res2plot.igk)
res2plot<-rbind(res2plot,res2plot.igh)

res2plot2use<-res2plot[,c(1,3,5,6)]
res2plot2use$TNFeff<-rep("Contracted",nrow(res2plot2use))
colnames(res2plot2use)[3]<-"Pval"
res2plot2use$Pheno<-rep("Pheno1",nrow(res2plot2use))

res2plot2use.exp<-res2plot2use
res2plot2use.exp$TNFeff<-"Expanded"

res2plot.p1<-rbind(res2plot2use, res2plot2use.exp)
res2plot.p2<-res2plot.p1
res2plot.p2$Pheno<-"Pheno2"

test.df<-rbind(res2plot.p1,res2plot.p2)
test.df$key<-paste0(test.df$Pheno,"-",test.df$Chain)
saveRDS(test.df, "Test.rds")

ggplot(data = test.df,aes(x=factor(key), y=factor(Biochemical.Feature), size=Pval, fill=Type)) +
  geom_point(alpha=0.55, shape=21, color="black") +
  scale_size(range=c(.1, 8), name="-log10(pval)") +
  facet_wrap(~TNFeff) + labs(y="", x="") +
  scale_fill_manual(values=c("red3","blue3")) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face="bold", size=16), axis.text.x=element_text(vjust=0.5, size=13), strip.text=element_text(size=14, face="bold"))

#### Example of the code used for characterizing all contracted clones ####

get.aa.features<-function(cl) {
  mw.res<-mw(cl)
  tiny.res<-aaComp(cl)[[1]]["Tiny","Number"]
  small.res<-aaComp(cl)[[1]]["Small","Number"]
  aliphatic.res<-aaComp(cl)[[1]]["Aliphatic","Number"]
  aromatic.res<-aaComp(cl)[[1]]["Aromatic","Number"]
  nonpolar.res<-aaComp(cl)[[1]]["NonPolar","Number"]
  polar.res<-aaComp(cl)[[1]]["Polar","Number"]
  charged.res<-aaComp(cl)[[1]]["Charged","Number"]
  basic.res<-aaComp(cl)[[1]]["Basic","Number"]
  acidic.res<-aaComp(cl)[[1]]["Acidic","Number"]
  netcharge.res<-charge(cl, pH=7, pKscale="EMBOSS")
  pI.res<-pI(cl, pKscale="EMBOSS")
  aliphI.res<-aIndex(cl)
  insta.res<-instaIndex(cl)
  boman.res<-boman(cl)
  hydroI.res<-hydrophobicity(cl)
  hmom.res<-hmoment(cl)
  vect.out<-c(mw.res,tiny.res,small.res,aliphatic.res,aromatic.res,nonpolar.res,polar.res,charged.res,basic.res,acidic.res,netcharge.res,pI.res,aliphI.res,insta.res,boman.res,hydroI.res,hmom.res)
  names(vect.out)<-c("Molecular.Weight","Tiny","Small","Aliphatic","Aromatic","Non.Polar","Polar","Charged","Basic","Acidic","Net.Charge","Isoelectric.Point","Aliphatic.Index","Instability.Index","PPI.Index","Hydrophobicity.Index","HydrophobicMoment.Index")
  return(vect.out)
}
get.perm.features<-function(set.n, chain.clones.all, c, get.aa.features, nperm) {
  count.low<-0;  count.hi<-0;
  df.perm.all<-data.frame("Chain"="TEST", "Perm"=0,
                          "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                          "Aliphatic"=0, "Aromatic"=0,
                          "Non.Polar"=0, "Polar"=0, "Charged"=0,
                          "Basic"=0, "Acidic"=0,
                          "Net.Charge"=0, "Isoelectric.Point"=0,
                          "Aliphatic.Index"=0, "Instability.Index"=0,
                          "PPI.Index"=0, "Hydrophobicity.Index"=0,
                          "HydrophobicMoment.Index"=0, stringsAsFactors=F)
  for (i in 1:nperm) {
    print(paste0(c,"-",i))
    perm.clones<-sample(chain.clones.all, size=set.n, replace=F)
    df.perm.clones<-data.frame("Chain"="TEST", "Perm"=0,
                               "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                               "Aliphatic"=0, "Aromatic"=0,
                               "Non.Polar"=0, "Polar"=0, "Charged"=0,
                               "Basic"=0, "Acidic"=0,
                               "Net.Charge"=0, "Isoelectric.Point"=0,
                               "Aliphatic.Index"=0, "Instability.Index"=0,
                               "PPI.Index"=0, "Hydrophobicity.Index"=0,
                               "HydrophobicMoment.Index"=0, stringsAsFactors=F)
    count.perm=0
    for (cl in perm.clones) {
      count.perm=count.perm+1
      cl.perm.prop<-get.aa.features(cl)
      df.perm.clones[count.perm,]<-c(c,i,cl.perm.prop)
    }
    cols.interest.perm<-sapply(df.perm.clones[,c(3:ncol(df.perm.clones))], as.numeric)
    perm.prop<-colMeans(cols.interest.perm)
    df.perm.all[i,]<-c(c,i,perm.prop)
  }
  perm.prop<-df.perm.all
  return(perm.prop)
}
get.final.out<-function(bc.obs, bc.random, nperm){
  aa.propierties<-names(bc.obs)[3:length(bc.obs)]
  df.out<-data.frame("Biochemical.Feature"="TEST","Obs.Value"="TEST","Perm.Value"="TEST","PermValue.LT.ObsValue"="TEST","PermValue.GT.ObsValue"="TEST","Pval.PermValue.LT.Obs.Value"="TEST","Pval.PermValue.GT.ObsValue"="TEST", stringsAsFactors=F)
  count=0
  for (aap in aa.propierties) {
    obs.value<-as.numeric(bc.obs[aap])
    perm.value<-as.numeric(bc.random[,aap])
    pval1<-(length(which(perm.value<obs.value))+1)/(nperm+1)
    pval2<-(length(which(perm.value>=obs.value))+1)/(nperm+1)
    count=count+1
    df.out[count,]<-c(aap, unname(obs.value), mean(perm.value), length(which(perm.value<obs.value)), length(which(perm.value>obs.value)), pval1, pval2)
  }
  return(df.out)
}

type="Contracted"
dir.create(paste0("../../Output_Data/NoisET_Characterization_",type,"_ALL/"))
nperm<-10000

df.res<-readRDS("../../Output_Data/NoisET_Summary_ExpCont.rds")
chains<-c("IGH","IGL","IGK")
samples2keep<-c(as.character(df.pheno[,1]), as.character(df.pheno[,3]))
for (c in chains) {
  df.res.chain<-df.res[df.res$Chain==c,]
  clones.cont<-unlist(strsplit(paste(df.res.chain[(df.res.chain$Donor %in% donor.pheno),"Cont.Sig.Clones"], collapse=';'),";"))
  clones.cont<-clones.cont[!(clones.cont %in% "NA")]
  chain.clones.cont.n<-length(unique(clones.cont))
  chain.clones.cont.max<-max(table(clones.cont))
  chain.clones.cont.unique<-unique(clones.cont)[!(unique(clones.cont) %in% "NA")]
  tbl.freq<-table(clones.cont)
  df.freq.sorted<-as.data.frame(tbl.freq[order(tbl.freq,decreasing = T)])
  df.freq.sorted<-df.freq.sorted[df.freq.sorted$clones.cont!="NA",]
  
  clones.total<-nrow(df.freq.sorted)
  
  irep.chain<-readRDS(paste0("../Data2analyzeClones_By_CDR3aa/",c,".rds"))
  meta_col<-dplyr::filter(irep.chain$meta, IMID=="RA");
  meta_col1<-meta_col[(meta_col$Sample %in% samples2keep),]
  samples_col<-meta_col1$Sample
  data_col<-irep.chain$data[match(samples_col,names(irep.chain$data))];
  irep_chain<-list(data_col,meta_col1);  names(irep_chain)<-c("data","meta")
  
  chain.clones.all<-as.data.frame(immunarch::pubRep(irep_chain$data, "aa", .verbose=F))[,"CDR3.aa"]
  
  # Biochemical Properties
  print(paste0(c," - Biochemical Prop"))
  bc.random<-get.perm.features(length(chain.clones.cont.unique), chain.clones.all, c, get.aa.features, nperm)
  bc.obs<-data.frame("Chain"="TEST", "Perm"=0,
                     "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                     "Aliphatic"=0, "Aromatic"=0,
                     "Non.Polar"=0, "Polar"=0, "Charged"=0,
                     "Basic"=0, "Acidic"=0,
                     "Net.Charge"=0, "Isoelectric.Point"=0,
                     "Aliphatic.Index"=0, "Instability.Index"=0,
                     "PPI.Index"=0, "Hydrophobicity.Index"=0,
                     "HydrophobicMoment.Index"=0, stringsAsFactors=F)
  count=0
  for (cl in chain.clones.cont.unique) {
    count=count+1
    cl.prop<-get.aa.features(cl)
    bc.obs[count,]<-c(c,cl,cl.prop)
  }
  cols.interest<-sapply(bc.obs[,c(3:ncol(bc.obs))], as.numeric)
  bc.obs<-c(c,"Obs",colMeans(cols.interest))
  names(bc.obs)[1:2]<-c("Chain","Perm")
  bc.final.out<-get.final.out(bc.obs, bc.random, nperm)
  saveRDS(bc.final.out, paste0("../../Output_Data/NoisET_Characterization_",type,"_ALL/BiochemicalProp_",c,".rds"))
}




















