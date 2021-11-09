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

# RA vs Controls

df.div.chain<-readRDS("Data2analyze/DivBaselineControls.rds")
chains<-c("TRA","TRB","TRD","TRG","IGH","IGL","IGK")
indexs<-c("D20","D50","Gini.Simp.Logit","Inv.Simp","Shannon")
df.pvals<-data.frame("Chain"=rep("NA",7*5), "Div.Measure"=rep("NA",7*5), "Beta"=rep("NA",7*5), "Pval"=rep("NA",7*5), stringsAsFactors=F)
count=0
for (i in indexs) {
  print(i)
  for (c in chains) {
    count=count+1
    df.chain<-df.div.chain[df.div.chain$Chain==c,]
    df.chain$IMID.num<-as.numeric(df.chain$IMID)
    df.chain$IMID.num[df.chain$IMID.num==1]<-0
    df.chain$IMID.num[df.chain$IMID.num==2]<-1
    df.chain$Gender.num<-as.numeric(df.chain$Gender)
    df.chain$Age.num<-as.numeric(df.chain$Age)
    df.chain$Reads.Chain<-as.numeric(df.chain$Reads.Chain)
    if (i=="D20" | i=="D50" | i=="Inv.Simp") {
      data<-df.chain[,c("IMID.num","Age.num","Gender.num","Reads.Chain",i)]   # data<-df.chain[,c("IMID.num","Age.num","Gender.num","Reads.Chain","Plate","RNA.Integrity",i)]
      colnames(data)[5]<-"Div"
      data$Div<-as.numeric(data$Div)
      model<-glm.nb(Div~IMID.num+Age.num+Gender.num+Reads.Chain, data=data)
      pval<-coef(summary(model))["IMID.num",'Pr(>|z|)']
      beta<-coef(summary(model))["IMID.num",'Estimate']
    }
    if (i=="Gini" | i=="Gini.Simp.Logit" | i=="Shannon") {
      data<-df.chain[,c("IMID.num","Age.num","Gender.num","Reads.Chain",i)] #data<-df.chain[,c("IMID.num","Age.num","Gender.num","Reads.Chain","Plate","RNA.Integrity",i)]
      colnames(data)[5]<-"Div"
      data$Div<-as.numeric(data$Div)
      model<-glm(IMID.num~Div+Age.num+Gender.num+Reads.Chain, data=data, family="binomial")
      pval<-coef(summary(model))["Div",'Pr(>|z|)']
      beta<-coef(summary(model))["Div",'Estimate']
    }
    df.pvals[count,]<-c(c,i,beta,pval)
  }
}
df.pvals$Div.Measure[df.pvals$Div.Measure=="Shannon"]<-"Shannon Entropy"
df.pvals$Div.Measure[df.pvals$Div.Measure=="Inv.Simp"]<-"Inverse-Simpson"
df.pvals$Div.Measure[df.pvals$Div.Measure=="Gini.Simp.Logit"]<-"Gini-Simpson"

write.csv(df.pvals, "Output_Data/Table_Div_RACTRL.csv", row.names=F)

data.summary.pvals<-list()
data.summary.beta<-list()
count=0
for (c in chains) {
  chain.summary.beta<-as.numeric(df.pvals[df.pvals$Chain==c,"Beta"])
  names(chain.summary.beta)<-df.pvals[df.pvals$Chain==c,"Div.Measure"]
  chain.summary.pval<-as.numeric(df.pvals[df.pvals$Chain==c,"Pval"])
  names(chain.summary.pval)<-df.pvals[df.pvals$Chain==c,"Div.Measure"]
  count=count+1
  data.summary.pvals[[count]]<-chain.summary.pval
  data.summary.beta[[count]]<-chain.summary.beta
}
df.p<-do.call("rbind", data.summary.pvals);  rownames(df.p)<-chains;  df.pval.mat<-as.matrix(df.p)
df.b<-do.call("rbind", data.summary.beta);  rownames(df.b)<-chains;  df.beta.mat<-as.matrix(df.b)

lims<-c(-round(abs(min(df.beta.mat)),2)-0.2, round(abs(min(df.beta.mat)),2)+0.2)
cols<-colorRampPalette(c("firebrick2","white","forestgreen"))(30)
sigs<-c(.000005,.00005,.0005,.005,.05)
jpeg("Output_Data/Div_RACTRL.jpeg", width=2200, height=2200, res=300)
corrplot(df.beta.mat, p.mat=df.pval.mat, is.corr=F, method="color", cl.lim=lims, col=cols, pch.cex=.9, tl.cex=1.2, tl.col="black",insig="label_sig",sig.level=sigs, mar=c(0,0,2,0))
dev.off()


# Longitudinal
  
  div.wk0<-readRDS("Data2analyze/DivBaselineControls.rds")
  div.wk12<-readRDS("Data2analyze/DivWeek12Controls.rds")
  div.wk0wk12<-rbind(div.wk0, div.wk12[div.wk12$IMID!="CTRL",])
  div.lgtd<-div.wk0wk12[(as.character(div.wk0wk12$Sample) %in% c(ra.wk0,ra.wk12,ctrl)),]
  df.lgtd<-gather(div.lgtd, condition, measurement, D50:Shannon, factor_key=T)
  df.lgtd$Week<-as.character(df.lgtd$Week)
  df.lgtd$Week[df.lgtd$Week=="wk0"]<-"0.Baseline"
  df.lgtd$Week[df.lgtd$Week=="wk12"]<-"1.Treated"
  df.lgtd$Week[df.lgtd$Week=="CTRL"]<-"2.Control"
  donor.data<-read.table("Data2analyze/Donor2RNA_Codes.txt", header=T)
  donor.data.wk0<-donor.data[,c(1,2)];    colnames(donor.data.wk0)<-c("Donor","Sample")
  donor.data.wk12<-donor.data[,c(1,3)];   colnames(donor.data.wk12)<-c("Donor","Sample")
  donor.data2<-rbind(donor.data.wk0, donor.data.wk12)
  donor.data2<-donor.data2[donor.data2$Sample %in% c(ra.wk0,ra.wk12),]
  
  vars<-c("D20","D50","Gini.Simp.Logit","Inv.Simp","Shannon")
  chains<-c("TRA","TRB","TRD","TRG","IGH","IGL","IGK")
  pvals.df<-data.frame("Chain"=rep("TEST",5*7), "Variable"=rep("TEST",5*7), "Beta"=rep(0,5*7), "Pval"=rep(0,5*7), stringsAsFactors=F)
  count=0
  for (v in vars) {
    for (c in chains) {
      count=count+1
      data2use<-df.lgtd[(df.lgtd$Week!="2.Control" & df.lgtd$Chain==c & df.lgtd$condition==v),]
      data2analyze<-merge(data2use, donor.data2, by="Sample")
      if (v=="D20" | v=="D50" | v=="Inv.Simp") {
        data<-data2analyze[,c("Donor","Sample","Week","Age","Gender","Reads.Chain","condition","measurement")]
        data$Week[data$Week=="0.Baseline"]<-0
        data$Week[data$Week=="1.Treated"]<-1
        data$Week<-as.numeric(data$Week)
        data$Gender<-as.numeric(data$Gender)
        model1<-glmer.nb(measurement~Week+Age+Gender+Reads.Chain+(1|Donor), nAGQ=0, data=data)
        pval<-summary(model1)$coefficients["Week","Pr(>|z|)"]
        beta<-summary(model1)$coefficients["Week","Estimate"]
      }
      if (v=="Gini" | v=="Gini.Simp.Logit" | v=="Shannon") {
        data<-data2analyze[,c("Donor","Sample","Week","Age","Gender","Reads.Chain","condition","measurement")]
        data$Week[data$Week=="0.Baseline"]<-0
        data$Week[data$Week=="1.Treated"]<-1
        data$Week<-as.numeric(data$Week)
        data$Gender<-as.numeric(data$Gender)
        model1<-glmer(Week~measurement+Age+Gender+Reads.Chain+(1|Donor),  family=binomial, data=data)
        pval<-summary(model1)$coefficients["measurement","Pr(>|z|)"]
        beta<-summary(model1)$coefficients["measurement","Estimate"]
      }
      pvals.df[count,]<-c(c,v,beta,pval)
    }
  }
  pvals.df$Variable[pvals.df$Variable=="Shannon"]<-"Shannon Entropy"
  pvals.df$Variable[pvals.df$Variable=="Gini.Simp.Logit"]<-"Gini-Simpson"
  pvals.df$Variable[pvals.df$Variable=="Inv.Simp"]<-"Inverse-Simpson"
  write.csv(pvals.df, "Output_Data/Table_Div_wk0wk12.csv", row.names=F)

data.summary.pvals.all<-list()
data.summary.beta.all<-list()
count.all=0
vars<-unique(pvals.df$Variable)
for (v in vars) {
  pvals.df.int<-pvals.df[pvals.df$Variable==v,]
  count.var=0
  data.summary.pvals<-list()
  data.summary.beta<-list()
  for (c in chains) {
    chain.summary.beta<-as.numeric(pvals.df.int[pvals.df.int$Chain==c,"Beta"])
    names(chain.summary.beta)<-pvals.df.int[pvals.df.int$Chain==c,"Variable"]
    chain.summary.pval<-as.numeric(pvals.df.int[pvals.df.int$Chain==c,"Pval"])
    names(chain.summary.pval)<-pvals.df.int[pvals.df.int$Chain==c,"Variable"]
    count.var=count.var+1
    data.summary.pvals[[count.var]]<-chain.summary.pval
    data.summary.beta[[count.var]]<-chain.summary.beta
  }
  df.p<-do.call("rbind", data.summary.pvals);  rownames(df.p)<-chains;  df.pval.mat<-as.matrix(df.p)
  df.b<-do.call("rbind", data.summary.beta);  rownames(df.b)<-chains;  df.beta.mat<-as.matrix(df.b)
  count.all=count.all+1
  data.summary.pvals.all[[count.all]]<-df.pval.mat
  data.summary.beta.all[[count.all]]<-df.beta.mat
}
df.p.all<-do.call("cbind", data.summary.pvals.all);  df.pval.mat.all<-as.matrix(df.p.all)
df.b.all<-do.call("cbind", data.summary.beta.all);   df.beta.mat.all<-as.matrix(df.b.all)

lims<-c(-max(abs(df.beta.mat.all)), max(abs(df.beta.mat.all)))
cols<-colorRampPalette(c("firebrick2","white","forestgreen"))(30)
sigs<-c(.000005,.00005,.0005,.005,.05)
jpeg("Output_Data/Div_Long.jpeg", width=2200, height=2200, res=300)
corrplot(df.beta.mat.all, p.mat=df.pval.mat.all, is.corr=F, method="color", cl.lim=lims, col=cols, pch.cex=.9, tl.cex=1.2, tl.col="black",insig="label_sig",sig.level=sigs, mar=c(0,0,2,0))
dev.off()


# Boxplots

data2use<-df.lgtd[(df.lgtd$Week!="2.Control" & df.lgtd$Chain=="IGK" & df.lgtd$condition=="D50"),]
data2analyze<-merge(data2use, donor.data2, by="Sample")
data2analyze$Week[data2analyze$Week=="0.Baseline"]<-"Baseline"
data2analyze$Week[data2analyze$Week=="1.Treated"]<-"TNFi therapy"
igk<-ggplot(data2analyze, aes(x=Week, y=measurement, fill=Week)) + 
  geom_boxplot(aes(group=Week)) + theme_classic() +
  scale_fill_manual(values=c("#F69191", "#D8EBD8")) +
  geom_point(aes(color=Week), size=1) +
  geom_line(aes(group=Donor), color="black", size=0.05) +
  scale_color_manual(values=c("#F69191", "#D8EBD8")) +
  labs(title="", x="", y="\nIGK Diversity (D50) \n") +
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.position="none") +
  annotate("text", x=1.5, y=820, label="P=2.72e-07\n")
data2use<-df.lgtd[(df.lgtd$Week!="2.Control" & df.lgtd$Chain=="IGL" & df.lgtd$condition=="D50"),]
data2analyze<-merge(data2use, donor.data2, by="Sample")
data2analyze$Week[data2analyze$Week=="0.Baseline"]<-"Baseline"
data2analyze$Week[data2analyze$Week=="1.Treated"]<-"TNFi therapy"
igl<-ggplot(data2analyze, aes(x=Week, y=measurement, fill=Week)) + 
  geom_boxplot(aes(group=Week)) + theme_classic() +
  scale_fill_manual(values=c("#F69191", "#D8EBD8")) +
  geom_point(aes(color=Week), size=1) +
  geom_line(aes(group=Donor), color="black", size=0.05) +
  scale_color_manual(values=c("#F69191", "#D8EBD8")) +
  labs(title="", x="", y="\nIGL Diversity (D50) \n") +
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.position="none") +
  annotate("text", x=1.5, y=750, label="P=2.16e-03\n")
data2use<-df.lgtd[(df.lgtd$Week!="2.Control" & df.lgtd$Chain=="IGH" & df.lgtd$condition=="D50"),]
data2analyze<-merge(data2use, donor.data2, by="Sample")
data2analyze$Week[data2analyze$Week=="0.Baseline"]<-"Baseline"
data2analyze$Week[data2analyze$Week=="1.Treated"]<-"TNFi therapy"
igh<-ggplot(data2analyze, aes(x=Week, y=measurement, fill=Week)) + 
  geom_boxplot(aes(group=Week)) + theme_classic() +
  scale_fill_manual(values=c("#F69191", "#D8EBD8")) +
  geom_point(aes(color=Week), size=1) +
  geom_line(aes(group=Donor), color="black", size=0.05) +
  scale_color_manual(values=c("#F69191", "#D8EBD8")) +
  labs(title="", x="", y="\nIGH Diversity (D50) \n") +
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.position="none") +
  annotate("text", x=1.5, y=7700, label="P=3.63e-03\n") + ylim(0,8000) +
  jpeg("Output_Data/Boxplot_1.jpeg", width=1000, height=3300, res=300)
grid.arrange(igk,igl,igh, ncol=1)
dev.off()

data2analyze<-df.lgtd[(df.lgtd$Week=="2.Control" & df.lgtd$Chain=="IGK" & df.lgtd$condition=="D50"),]
data2analyze$Week<-data2analyze$Week[data2analyze$Week=="2.Control"]<-"Controls"
igk<-ggplot(data2analyze, aes(x=Week, y=measurement, fill=Week)) + 
  geom_boxplot(aes(group=Week)) + theme_classic() +
  scale_fill_manual(values=c("#7DBB7D")) +
  geom_point(aes(color=Week), size=1) +
  scale_color_manual(values=c("#7DBB7D")) +
  labs(title="", x="", y="\nIGK Diversity (D50) \n") + ylim(0,800) +
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.position="none")
data2analyze<-df.lgtd[(df.lgtd$Week=="2.Control" & df.lgtd$Chain=="IGL" & df.lgtd$condition=="D50"),]
data2analyze$Week<-data2analyze$Week[data2analyze$Week=="2.Control"]<-"Controls"
igl<-ggplot(data2analyze, aes(x=Week, y=measurement, fill=Week)) + 
  geom_boxplot(aes(group=Week)) + theme_classic() +
  scale_fill_manual(values=c("#7DBB7D")) +
  geom_point(aes(color=Week), size=1) +
  scale_color_manual(values=c("#7DBB7D")) +
  labs(title="", x="", y="\nIGL Diversity (D50) \n") + ylim(0,800) +
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.position="none")
data2analyze<-df.lgtd[(df.lgtd$Week=="2.Control" & df.lgtd$Chain=="IGH" & df.lgtd$condition=="D50"),]
data2analyze$Week<-data2analyze$Week[data2analyze$Week=="2.Control"]<-"Controls"
igh<-ggplot(data2analyze, aes(x=Week, y=measurement, fill=Week)) + 
  geom_boxplot(aes(group=Week)) + theme_classic() +
  scale_fill_manual(values=c("#7DBB7D")) +
  geom_point(aes(color=Week), size=1) +
  scale_color_manual(values=c("#7DBB7D")) +
  labs(title="", x="", y="\nIGH Diversity (D50) \n") +
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"), legend.position="none")
jpeg("Output_Data/Boxplots_2.jpeg", width=700, height=3300, res=300)
grid.arrange(igk,igl,igh, ncol=1)
dev.off()


# Longitudinal stratified

phenos<-c("RespGM")
vars<-c("D50")
chains<-c("IGK","IGL","IGH")
pvals.pheno<-data.frame("Variable"=rep("TEST", 3), "Chain"=rep("TEST", 3),
                        "Pheno"=rep("TEST", 3), "Subpheno"=rep("TEST", 3),
                        "Beta"=rep(0, 3), "Pval"=rep(0, 3), stringsAsFactors=F)
count=0
for (v in vars) {
  for (c in chains) {
    for (ph in phenos) {
      print(paste0(v," - ",c," - ",ph))
      data2use<-df.lgtd[(df.lgtd$Week!="2.Control" & df.lgtd$Chain==c & df.lgtd$condition==v),]
      data2analyze<-merge(data2use, donor.data2, by="Sample")
      data2analyze$RespGM<-data2analyze$Resp
      data2analyze$RespGM[data2analyze$RespGM=="MOD"]<-"GOOD"
      data2analyze$Resp[data2analyze$Resp=="MOD"]<-NA
      data2test<-data2analyze[,c("Donor","Sample","Week","Age","Gender","Reads.Chain","condition","measurement",ph)]
      colnames(data2test)[9]<-"Pheno"
      data2test<-data2test[!is.na(data2test$Pheno),]
      data2test$Pheno<-droplevels(data2test$Pheno)
      subphenos<-levels(data2test$Pheno)
      if (ph=="ActB") {
        df.baseline<-data2test[data2test$Week=="0.Baseline",]
        df.baseline$Donor<-droplevels(df.baseline$Donor)
        donor.act0<-table(df.baseline$Donor, df.baseline$Pheno)[,"High"];   donor.act0<-gsub(0,"Low",donor.act0);   donor.act0<-gsub(1,"High",donor.act0);
        donor.act<-data.frame("Donor"=names(table(df.baseline$Donor, df.baseline$Pheno)[,"High"]), "ActB"=donor.act0)
        data2test.expl<-merge(data2test, donor.act, by="Donor")
        data2test.expl$Pheno<-data2test.expl$ActB
        data2test<-data2test.expl
      }
      for (sp in subphenos) {
        data2test.sp<-data2test[data2test$Pheno==sp,]
        if (v=="D20" | v=="D50" | v=="Inv.Simp") {
          data<-data2test.sp[,c("Donor","Sample","Week","Age","Gender","Reads.Chain","condition","measurement")]
          data$Week[data$Week=="0.Baseline"]<-0
          data$Week[data$Week=="1.Treated"]<-1
          data$Week<-as.numeric(data$Week)
          data$Gender<-as.numeric(data$Gender)
          model1<-glmer.nb(measurement~Week+Age+Gender+Reads.Chain+(1|Donor), nAGQ=0, data=data)
          pval<-summary(model1)$coefficients["Week","Pr(>|z|)"]
          beta<-summary(model1)$coefficients["Week","Estimate"]
        }
        if (v=="Gini" | v=="Gini.Simp.Logit" | v=="Shannon") {
          data<-data2test.sp[,c("Donor","Sample","Week","Age","Gender","Reads.Chain","condition","measurement")]
          data$Week[data$Week=="0.Baseline"]<-0
          data$Week[data$Week=="1.Treated"]<-1
          data$Week<-as.numeric(data$Week)
          data$Gender<-as.numeric(data$Gender)
          model1<-glmer(Week~measurement+Age+Gender+Reads.Chain+(1|Donor),  family=binomial, data=data)
          pval<-summary(model1)$coefficients["measurement","Pr(>|z|)"]
          beta<-summary(model1)$coefficients["measurement","Estimate"]
        }
        count=count+1
        pvals.pheno[count,]<-c(v,c,ph,sp,beta,pval)
      }
    }
  }
}
pvals.pheno$Pval<-as.numeric(pvals.pheno$Pval)
pvals.pheno$Beta<-as.numeric(pvals.pheno$Beta)

sig.data<-pvals.pheno[pvals.pheno$Pval<0.05 & pvals.pheno$Chain!="TRD",]
sig.phenos<-unique(sig.data$Pheno)
sig.data.pheno<-sig.data[sig.data$Pheno=="RespGM",]
subph.sig<-unique(sig.data.pheno$Subpheno)
sig.pheno.vars<-unique(sig.data.pheno$Variable)
sig.data.pheno.vars<-sig.data.pheno[sig.data.pheno$Variable=="D50",]
chains<-c("IGK","IGL","IGH")
spv<-unique(sig.data$Variable)
for (chain in chains) {
  pval.sig<-sig.data.pheno.vars[sig.data.pheno.vars$Chain==chain,"Pval"]
  subph.nosig<-"NON_RESP"
  pval.nosig<-pvals.pheno[pvals.pheno$Pheno=="RespGM" & pvals.pheno$Subpheno==subph.nosig & pvals.pheno$Chain==chain & pvals.pheno$Variable=="D50","Pval"]
  data0<-df.lgtd[(df.lgtd$Week!="2.Control" & df.lgtd$Chain==chain & df.lgtd$condition==spv),]
  data00<-merge(data0, donor.data2, by="Sample")
  data00$RespGM<-data00$Resp
  data00$RespGM[data00$RespGM=="MOD"]<-"GOOD"
  data00$Resp[data00$Resp=="MOD"]<-NA
  data2test<-data00[,c("Donor","Sample","Week","Chain","Age","Gender","Reads.Chain","condition","measurement",ph)]
  colnames(data2test)[10]<-"Pheno"
  data2test<-data2test[!is.na(data2test$Pheno),]
  data2test$Pheno<-droplevels(data2test$Pheno)
  data2use<-data2test[(data2test$Chain==chain & data2test$condition==spv), c("Donor","Sample","Chain","Week","condition","measurement","Pheno")]
  data2use.ra<-data2use[!is.na(data2use$Pheno),]
  data2use.ctrls<-df.lgtd[(df.lgtd$Week=="2.Control" & df.lgtd$Chain==chain & df.lgtd$condition==spv), c("Sample","Sample","Chain","Week","condition","measurement","Week")]
  colnames(data2use.ctrls)<-c("Donor","Sample","Chain","Week","condition","measurement","Pheno")
  data2use.all<-rbind(data2use.ra,data2use.ctrls)
  data2use.all$Pheno<-as.character(data2use.all$Pheno)
  data2use.all$Pheno[data2use.all$Pheno=="2.Control"]<-"Controls"
  pvals2use<-c(textPval<-paste0("P=",formatC(pval.nosig,format="e", digits=2)), textPval<-paste0("P=",formatC(pval.sig,format="e", digits=2)), "-")
  data2use.all$Pheno[data2use.all$Pheno=="GOOD"]<-"Responders"
  data2use.all$Pheno[data2use.all$Pheno=="NON_RESP"]<-"Non Responders"
  names(pvals2use)<-c("Non Responders", "Responders", "Controls")
  data2use.all$Pheno<-factor(data2use.all$Pheno, levels=names(pvals2use))
  plot.sig<-ggboxplot(data2use.all, x="Pheno", y="measurement", color="Week", add="jitter", fill="Week") +
    labs(title="", x="", y=paste0("\n",chain," Diversity (",spv,")\n")) + scale_color_manual(values=c("black","black","black")) +
    scale_fill_manual(values=c("slategray1", "plum3", "#7DBB7D")) + ylim(0,max(data2use.all$measurement)*1.1) +
    theme(plot.title=element_text(hjust=0.5, size=16, face="bold"), legend.position="right") +
    annotate(geom="text", x=1, y=max(data2use.all$measurement)*1.1, label=pvals2use[1], color="black", size=4) +
    annotate(geom="text", x=2, y=max(data2use.all$measurement)*1.1, label=pvals2use[2], color="black", size=4)
  jpeg(paste0("Output_Data/Long_strat_",spv,"_",chain,".jpeg"), width=2100, height=2000, res=300)
  plot(plot.sig)
  dev.off()
}


# RA pheno

df.div.chain<-readRDS("Data2analyze/DivBaselineControls.rds")
phenos<-c("RespGM","ActB","ACPA","RF")
chains<-c("TRA","TRB","TRD","TRG","IGH","IGL","IGK")
indexs<-c("D20","D50","Gini.Simp.Logit","Inv.Simp","Shannon")
df.div.chain$RespGM<-df.div.chain$Resp
df.div.chain$RespGM[df.div.chain$RespGM=="MOD"]<-"GOOD"
df.div.chain$Resp[df.div.chain$Resp=="MOD"]<-NA

df.pvals.pheno<-data.frame("Pheno"=rep("NA",7*5*4),"Chain"=rep("NA",7*5*4), "Div.Measure"=rep("NA",7*5*4), "Beta"=rep("NA",7*5*4), "Pval"=rep("NA",7*5*4), stringsAsFactors=F)
count=0
for (ph in phenos) {
  for (i in indexs) {
    print(paste0(ph,"-",i))
    for (c in chains) {
      count=count+1
      df.chain<-df.div.chain[(df.div.chain$Chain==c & df.div.chain$IMID=="RA"),]
      df.chain$IMID.num<-as.numeric(df.chain$IMID)
      df.chain$IMID.num[df.chain$IMID.num==1]<-0
      df.chain$IMID.num[df.chain$IMID.num==2]<-1
      df.chain$Gender.num<-as.numeric(df.chain$Gender)
      df.chain$Age.num<-as.numeric(df.chain$Age)
      df.chain$Reads.Chain<-as.numeric(df.chain$Reads.Chain)
      if (i=="D20" | i=="D50" | i=="Inv.Simp") {
        data<-df.chain[,c(ph,"Age.num","Gender.num","Reads.Chain",i)]
        colnames(data)[1]<-"Pheno"; colnames(data)[5]<-"Div"
        data$Div<-as.numeric(data$Div)
        data2use<-data[!is.na(data$Pheno),]; data2use$Pheno<-droplevels(data2use$Pheno);  data2use$Pheno<-as.numeric(data2use$Pheno)
        model<-glm.nb(Div~Pheno+Age.num+Gender.num+Reads.Chain, data=data2use)
        pval<-coef(summary(model))["Pheno",'Pr(>|z|)']
        beta<-coef(summary(model))["Pheno",'Estimate']
      }
      if (i=="Gini" | i=="Gini.Simp.Logit" | i=="Shannon") {
        data<-df.chain[,c(ph,"Age.num","Gender.num","Reads.Chain",i)]
        colnames(data)[1]<-"Pheno"; colnames(data)[5]<-"Div"
        data$Div<-as.numeric(data$Div)
        data2use<-data[!is.na(data$Pheno),]; data2use$Pheno<-droplevels(data2use$Pheno);  data2use$Pheno<-as.numeric(data2use$Pheno)
        data2use$Pheno[data2use$Pheno==1]<-0;  data2use$Pheno[data2use$Pheno==2]<-1
        model<-glm(Pheno~Div+Age.num+Gender.num+Reads.Chain, data=data2use, family="binomial")
        pval<-coef(summary(model))["Div",'Pr(>|z|)']
        beta<-coef(summary(model))["Div",'Estimate']
      }
      df.pvals.pheno[count,]<-c(ph,c,i,beta,pval)
    }
  }
}
df.pvals.pheno$Div.Measure[df.pvals.pheno$Div.Measure=="Gini.Simp.Logit"]<-"Gini-Simpson"
df.pvals.pheno$Div.Measure[df.pvals.pheno$Div.Measure=="Inv.Simp"]<-"Inverse-Simpson"
df.pvals.pheno$Div.Measure[df.pvals.pheno$Div.Measure=="Shannon"]<-"Shannon Entropy"

write.csv(df.pvals.pheno, "Output_Data/Table_Div_RApheno.csv", row.names=F)

pdf("Output_Data/RApheno_Figure_Div.pdf")
for (ph in phenos) {
  df.pvals.pheno.int<-df.pvals.pheno[df.pvals.pheno$Pheno==ph,]
  data.summary.pvals<-list()
  data.summary.beta<-list()
  count=0
  for (c in chains) {
    chain.summary.beta<-as.numeric(df.pvals.pheno.int[df.pvals.pheno.int$Chain==c,"Beta"])
    names(chain.summary.beta)<-df.pvals.pheno.int[df.pvals.pheno.int$Chain==c,"Div.Measure"]
    chain.summary.pval<-as.numeric(df.pvals.pheno.int[df.pvals.pheno.int$Chain==c,"Pval"])
    names(chain.summary.pval)<-df.pvals.pheno.int[df.pvals.pheno.int$Chain==c,"Div.Measure"]
    count=count+1
    data.summary.pvals[[count]]<-chain.summary.pval
    data.summary.beta[[count]]<-chain.summary.beta
  }
  df.p<-do.call("rbind", data.summary.pvals);  rownames(df.p)<-chains;  df.pval.mat<-as.matrix(df.p)
  df.b<-do.call("rbind", data.summary.beta);  rownames(df.b)<-chains;  df.beta.mat<-as.matrix(df.b)
  
  lims<-c(-max(abs(df.beta.mat)), max(abs(df.beta.mat)))
  cols<-colorRampPalette(c("firebrick2","white","forestgreen"))(30)
  sigs<-c(.000005,.00005,.0005,.005,.05)
  if (ph=="RespGM") {pheno.name<-"Response to TNFi therapy"}
  if (ph=="ActB") {pheno.name<-"Disease activity"}
  if (ph=="ACPA") {pheno.name<-"ACPA"}
  if (ph=="RF") {pheno.name<-"RF"}
  corrplot(df.beta.mat, p.mat=df.pval.mat, is.corr=F, method="color", cl.lim=lims, col=cols, pch.cex=.9, tl.cex=1.2, tl.col="black", insig="label_sig",sig.level=sigs, mar=c(0,0,2,0), title=paste0("\n\n",pheno.name,"\n\n"))
}
dev.off()

    
    phenos<-c("RespGM")
    vars<-c("D20","D50","Gini.Simp.Logit","Inv.Simp","Shannon")
    chains<-c("TRA","TRB","TRD","TRG","IGH","IGK","IGL")
    pvals.pheno<-data.frame("Variable"=rep("TEST", 420), "Chain"=rep("TEST",420),
                            "Pheno"=rep("TEST", 420), "Subpheno"=rep("TEST", 420),
                            "Beta"=rep(0, 420), "Pval"=rep(0,420), stringsAsFactors=F)
    count=0
    for (v in vars) {
      for (c in chains) {
        for (ph in phenos) {
          print(paste0(v," - ",c," - ",ph))
          data2use<-df.lgtd[(df.lgtd$Week!="2.Control" & df.lgtd$Chain==c & df.lgtd$condition==v),]
          data2analyze<-merge(data2use, donor.data2, by="Sample")
          data2analyze$RespGM<-data2analyze$Resp
          data2analyze$RespGM[data2analyze$RespGM=="MOD"]<-"GOOD"
          data2analyze$Resp[data2analyze$Resp=="MOD"]<-NA
          data2test<-data2analyze[,c("Donor","Sample","Week","Age","Gender","Reads.Chain","condition","measurement",ph)]
          colnames(data2test)[9]<-"Pheno"
          data2test<-data2test[!is.na(data2test$Pheno),]
          data2test$Pheno<-droplevels(data2test$Pheno)
          subphenos<-levels(data2test$Pheno)
          if (ph=="ActB") {
            df.baseline<-data2test[data2test$Week=="0.Baseline",]
            df.baseline$Donor<-droplevels(df.baseline$Donor)
            donor.act0<-table(df.baseline$Donor, df.baseline$Pheno)[,"High"];   donor.act0<-gsub(0,"Low",donor.act0);   donor.act0<-gsub(1,"High",donor.act0);
            donor.act<-data.frame("Donor"=names(table(df.baseline$Donor, df.baseline$Pheno)[,"High"]), "ActB"=donor.act0)
            data2test.expl<-merge(data2test, donor.act, by="Donor")
            data2test.expl$Pheno<-data2test.expl$ActB
            data2test<-data2test.expl
          }
          for (sp in subphenos) {
            data2test.sp<-data2test[data2test$Pheno==sp,]
            if (v=="D20" | v=="D50" | v=="Inv.Simp") {
              data<-data2test.sp[,c("Donor","Sample","Week","Age","Gender","Reads.Chain","condition","measurement")]
              data$Week[data$Week=="0.Baseline"]<-0
              data$Week[data$Week=="1.Treated"]<-1
              data$Week<-as.numeric(data$Week)
              data$Gender<-as.numeric(data$Gender)
              model1<-glmer.nb(measurement~Week+Age+Gender+Reads.Chain+(1|Donor), nAGQ=0, data=data)
              pval<-summary(model1)$coefficients["Week","Pr(>|z|)"]
              beta<-summary(model1)$coefficients["Week","Estimate"]
            }
            if (v=="Gini" | v=="Gini.Simp.Logit" | v=="Shannon") {
              data<-data2test.sp[,c("Donor","Sample","Week","Age","Gender","Reads.Chain","condition","measurement")]
              data$Week[data$Week=="0.Baseline"]<-0
              data$Week[data$Week=="1.Treated"]<-1
              data$Week<-as.numeric(data$Week)
              data$Gender<-as.numeric(data$Gender)
              model1<-glmer(Week~measurement+Age+Gender+Reads.Chain+(1|Donor),  family=binomial, data=data)
              pval<-summary(model1)$coefficients["measurement","Pr(>|z|)"]
              beta<-summary(model1)$coefficients["measurement","Estimate"]
            }
            count=count+1
            pvals.pheno[count,]<-c(v,c,ph,sp,beta,pval)
          }
        }
      }
    }
    pvals.pheno$Pval<-as.numeric(pvals.pheno$Pval)
    pvals.pheno$Beta<-as.numeric(pvals.pheno$Beta)
    pvals.pheno.long<-pvals.pheno[pvals.pheno$Chain!="TEST",]
    
    write.csv(pvals.pheno.long, "Output_Data/Table_Div_Pheno_wk0wk12.csv", row.names=F)



    
    
    
