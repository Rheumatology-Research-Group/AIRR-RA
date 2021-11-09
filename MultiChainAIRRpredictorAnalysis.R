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

# Analysis
# The script provided in this section is related to the multi-chain AIRR predictor for RA diagnosis. It must be modified in order to evaluate the predictive potential of AIRR-seq data for RA clinical phenotypes

df.data<-readRDS("Data2analyze/SummaryAIRR.rds")
 
loocv_probs.bal<-function(x, y, cv, var_selection=F, filter_var=0.05, nb_cores, n_top, pheno2pred){
  
  cl=makeCluster(as.numeric(nb_cores), outfile="")
  registerDoParallel(cl) 
  
  results<-foreach(fold=cv, .packages="randomForest") %dopar% {
    
    trainPred=x[-c(fold),,drop=F]
    trainPheno=as.factor(as.character(y[-c(fold)]))
    
    testPred=x[fold,,drop=F]
    testPheno=y[fold]
    
    model=randomForest(as.data.frame(trainPred), trainPheno, importance=T, do.trace=F, probs=T, strata=trainPheno, sampsize=rep(min(table(trainPheno)),2))
    list(proba=predict(model, testPred, type="prob"))
  }
  
  stopCluster(cl)
  
  prob=foreach(fold.result=results, fold.num=icount(), .combine=rbind) %do%{as.data.frame(fold.result$proba)}
  return(prob[,pheno2pred])
}

get.pheno<-function(x, y, prob.pred) {
  df.out.pred<-data.frame("Sample"=rownames(x), "Obs.Status"=y, "RF.prob.bal"=prob.pred, "RF.Pred.Status.bal"=rep(NA,nrow(x)))
  df.out.pred$RF.Pred.Status.bal[df.out.pred$RF.prob.bal<0.5]<-phenoOther
  df.out.pred$RF.Pred.Status.bal[df.out.pred$RF.prob.bal>=0.5]<-pheno2pred
  return(df.out.pred)
}
get.metrics<-function(conf.matrix) {
  acc<-conf.matrix$overall["Accuracy"][[1]]
  sens<-conf.matrix$byClass["Sensitivity"][[1]]
  spec<-conf.matrix$byClass["Specificity"][[1]]
  prec<-conf.matrix$byClass["Precision"][[1]]
  fpr<-(1-sens)
  metrics<-c(acc, sens, spec, prec, fpr)
  names(metrics)<-c("Accuracy", "Sensitivity", "Specificity", "Precision", "FPR")
  return(metrics)
}

pheno.name="IMID";  pheno2pred="RA";  phenoOther="CTRL";  covars="";  week2keep="wk0";  week2discard="wk12"  # Phenotype information should be adapted in accordance to the predictor of interest.
  
pheno.info<-read.csv("Data2analyze/Metadata.csv", head=T)
pheno.info2<-pheno.info[(pheno.info$Week!=week2discard), c("Sample",pheno.name)];  colnames(pheno.info2)<-c("Sample","Pheno")
df.pheno<-pheno.info2[(pheno.info2$Pheno==pheno2pred | pheno.info2$Pheno==phenoOther),]
df.pheno<-df.pheno[!(is.na(df.pheno$Pheno)),]
  
df.data<-merge(df.data, df.pheno, by="Sample")
  
preds.list<-list(c("Reads.All"),
                 colnames(epi.wide)[-1],
                 colnames(df.is.perc.wide)[-1],
                 colnames(df.csr.wide)[-1],
                 colnames(df.shm.wide)[-1],
                 colnames(df.IgDMmut.perc.wide)[-1],
                 colnames(df.ucn.wide)[-1],
                 colnames(cu.wide)[-1],
                 colnames(div.D50.wide)[-1],
                 colnames(div.D20.wide)[-1],
                 colnames(div.Gini.wide)[-1],
                 colnames(div.GiniSimp.wide)[-1],
                 colnames(div.InvSimp.wide)[-1],
                 colnames(div.Shannon.wide)[-1],
                 colnames(df.gu.V.all)[-ncol(df.gu.V.all)],
                 colnames(df.gu.J.all)[-ncol(df.gu.J.all)],
                 colnames(df.gu.VJ.all)[-ncol(df.gu.VJ.all)],
                 colnames(geno.wide)[grep("^HLA_",colnames(geno.wide))],
                 colnames(geno.wide)[grep("^AA_",colnames(geno.wide))],  
                 colnames(geno.wide)[grep("^SNP_",colnames(geno.wide))],
                 colnames(geno.wide)[grep("^rs",colnames(geno.wide))])

names(preds.list)<-c("Reads_Technical",
                     "AgeGender_Epidemiological",
                     "IsotypePerc_iRepBCR",
                     "CSRIndex_iRepBCR",
                     "SHMPerc_iRepBCR",
                     "MutIgDMPerc_iRepBCR",
                     "UniqueCloneNumber_iRepUCN",
                     "ChainUsage_iRepCU",
                     "D50_iRepDiv",
                     "D20_iRepDiv",
                     "Gini_iRepDiv",
                     "GiniSimp_iRepDiv",
                     "InvSimp_iRepDiv",
                     "Shannon_iRepDiv",
                     "GeneUsageV_iRepGU",
                     "GeneUsageJ_iRepGU",
                     "GeneUsageVJ_iRepGU",
                     "ClassicHLA_GenoHLA",
                     "AminoAcidsHLA_GenoHLA",
                     "IndelsHLA_GenoHLA",
                     "SNPsHLA_GenoHLA")

df.perf<-data.frame("Pred.Name"=rep("TEST",length(preds.list)),
                    "Pred.Type"=rep("TEST",length(preds.list)),
                    "Accuracy"=rep("TEST",length(preds.list)),
                    "Sensitivity"=rep("TEST",length(preds.list)),
                    "Specificity"=rep("TEST",length(preds.list)),
                    "Precision"=rep("TEST",length(preds.list)),
                    "FPR"=rep("TEST",length(preds.list)), stringsAsFactors=F)
count=0
for (n in names(preds.list)) {
  print(n)
  pred<-preds.list[[n]]
  x<-as.data.frame(df.data[,pred], row.names<-as.character(df.data$Sample))
  if ("Gender" %in% colnames(x)) {x$Gender<-as.numeric(x$Gender)}
  y<-factor(df.data[,"Pheno"], levels=c(pheno2pred,phenoOther))
  cv<-as.list(1:NROW(x));  nb_cores=7;  n_top=10
  prob.bal.pred<-loocv_probs.bal(x, y, cv, var_selection=F, filter_var=0.05, nb_cores, n_top, pheno2pred)
  df.out.pred<-get.pheno(x, y, prob.bal.pred)
  conf.matrix.bal.pred<<-confusionMatrix(factor(df.out.pred$RF.Pred.Status.bal, levels=c(pheno2pred,phenoOther)), factor(df.out.pred$Obs.Status, levels=c(pheno2pred,phenoOther)))
  metrics.pred<-c(strsplit(n,"_")[[1]][1], strsplit(n,"_")[[1]][2], get.metrics(conf.matrix.bal.pred))
  count=count+1
  df.perf[count,]<-metrics.pred
}

acc.thrs<-c(0.9,0.8,0.7,0.6,0.5)
count=nrow(df.perf)
for (th in acc.thrs) {
  th.name<-gsub('\\.','',th);   print(th.name);
  df.perf.ord<-df.perf[order(df.perf$Accuracy, decreasing=T),]
  types<-unique(df.perf$Pred.Type)
  
  preds.sel<-c();   c=0;
  for (t in types) {
    data.pred<-df.perf.ord[(df.perf.ord$Accuracy>th & df.perf.ord$Pred.Type==t),c("Pred.Name","Pred.Type")][1,]
    pred2use<-paste0(data.pred[,1],"_",data.pred[,2])
    if (pred2use!="NA_NA") {c=c+1; preds.sel[c]<-pred2use}
  }
  
  preds.agg<-unlist(preds.list[names(preds.list) %in% preds.sel])
  names.preds.agg<-paste(names(preds.list[names(preds.list) %in% preds.sel]), collapse='/')
  
  x<-as.data.frame(df.data[,preds.agg], row.names<-as.character(df.data$Sample))
  if ("Gender" %in% colnames(x)) {x$Gender<-as.numeric(x$Gender)}
  y<-factor(df.data[,"Pheno"], levels=c(pheno2pred,phenoOther))
  cv<-as.list(1:NROW(x));  nb_cores=7;  n_top=10
  prob.bal.pred<-loocv_probs.bal(x, y, cv, var_selection=F, filter_var=0.05, nb_cores, n_top, pheno2pred)
  df.out.pred<-get.pheno(x, y, prob.bal.pred)
  conf.matrix.bal.pred<<-confusionMatrix(factor(df.out.pred$RF.Pred.Status.bal, levels=c(pheno2pred,phenoOther)), factor(df.out.pred$Obs.Status, levels=c(pheno2pred,phenoOther)))
  metrics.pred<-c(paste0("AggPred",th.name,":",names.preds.agg), "MultiRep.Predictor", get.metrics(conf.matrix.bal.pred))
  count=count+1
  df.perf[count,]<-metrics.pred
}
saveRDS(df.perf, paste0("Output_Data/",pheno.name,"_Perf_All.rds"))

make.plot<-function(var, df.perf.long) {
  df.pred.type<-df.perf.long[df.perf.long$Pred.Type==var,]
  n.names<-length(unique(df.pred.type$Pred.Name))
  cols<-brewer.pal(n=n.names, name="BuPu")
  if (n.names==1) {cols<-cols[2]}
  p.metrics.all<-ggplot(df.pred.type, aes(x=condition, y=measurement, fill=Pred.Name)) + theme_classic() +
    geom_text(aes(label=round(measurement,2), y=measurement+0.02), position=position_dodge(1), size=2) +
    geom_bar(stat="identity", color="black", position=position_dodge())+ylim(0,1) +
    scale_fill_manual(values=cols) +
    theme(legend.title=element_blank(), axis.text=element_text(size=12), plot.title=element_text(hjust=0.5, face="bold", size=15)) +
    labs(title=paste0("\n",unique(df.pred.type$Pred.Type)," variables \n"), x="", y="\nMetric value\n")
  return(p.metrics.all)
}

df.perf.long<-gather(df.perf, condition, measurement, Accuracy:FPR, factor_key=T)
df.perf.long$measurement<-as.numeric(df.perf.long$measurement)

vars<-unique(df.perf.long$Pred.Type)

p1<-make.plot(vars[1], df.perf.long)
p2<-make.plot(vars[2], df.perf.long)
p3<-make.plot(vars[3], df.perf.long)
p4<-make.plot(vars[4], df.perf.long)
p5<-make.plot(vars[5], df.perf.long)
p6<-make.plot(vars[6], df.perf.long)
p7<-make.plot(vars[7], df.perf.long)
p8<-make.plot(vars[8], df.perf.long)
p9<-make.plot(vars[9], df.perf.long)

jpeg(paste0("Output_Data/",pheno.name,"_Perf_All_Metrics.jpeg"), width=5300, height=6000, res=300)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=2, top=textGrob("\nComparative Predictors Performance\n", gp=gpar(fontsize=20,fontface="bold")))
dev.off()

jpeg(paste0("Output_Data/",pheno.name,"_Perf_AggPred_Metrics.jpeg"), width=5000, height=3000, res=300)
grid.arrange(p9, ncol=1, top=textGrob("\nComparative Predictors Performance\n", gp=gpar(fontsize=20,fontface="bold")))
dev.off()


# The following section is specific for characterizing the best multi-chain AIRR predictor for RA diagnosis (i.e. chain usage + V gene segment usage)

preds<-c(preds.list[["ChainUsage_iRepCU"]],preds.list[["GeneUsageV_iRepGU"]])
                                                  
df.perf.best<-data.frame("Pred.Name"="TEST", "Pred.Type"="TEST", "Accuracy"="TEST", "Sensitivity"="TEST", "Specificity"="TEST", "Precision"="TEST", "FPR"="TEST", stringsAsFactors=F)

x<-as.data.frame(df.data[,preds], row.names<-as.character(df.data$Sample))
if ("Gender" %in% colnames(x)) {x$Gender<-as.numeric(x$Gender)}
y<-factor(df.data[,"Pheno.x"], levels=c(pheno2pred,phenoOther))
cv<-as.list(1:NROW(x));  nb_cores=7;  n_top=10
prob.bal.pred<-loocv_probs.bal(x, y, cv, var_selection=F, filter_var=0.05, nb_cores, n_top, pheno2pred)
df.out.pred<-get.pheno(x, y, prob.bal.pred)
conf.matrix.bal.pred<<-confusionMatrix(factor(df.out.pred$RF.Pred.Status.bal, levels=c(pheno2pred,phenoOther)), factor(df.out.pred$Obs.Status, levels=c(pheno2pred,phenoOther)))
conf.matrix.bal.pred$table


metrics.pred<-c(pred.name.best, pred.type.best, get.metrics(conf.matrix.bal.pred))
df.perf.best[1,]<-metrics.pred

df.perf<-readRDS("Output_Data/IMID_Perf_All.rds")
df.perf<-df.perf[order(df.perf$Accuracy, decreasing=T),]
data2plot<-data.frame("Metric"=c("Accuracy","Sensitivity","Specificity","Precision","1-FPR"), "Data"=unlist(c(unname(df.perf[1,3:7]))), stringsAsFactors=F)
data2plot$Data<-as.numeric(data2plot$Data)
data2plot$Data[5]<-(1-data2plot$Data[5])
data2plot$Data=100*data2plot$Data
data2plot$Metric<-factor(data2plot$Metric, levels=c("Accuracy","Sensitivity","Specificity","Precision","1-FPR"))

p3<-ggplot(data2plot, aes(fill=Metric, y=Data, x=Metric)) +
    geom_text(aes(label=round(Data,2), y=Data+12), position=position_dodge(1), size=4) +
    geom_bar(position="stack", stat="identity") + theme_classic() + theme(axis.text=element_text(size=12), axis.title.x = element_text(size = 13)) +
    scale_fill_manual(values=rev(brewer.pal(n=5, name="PuBu"))) + ylim(0,120) +
    labs(y="\nMetric value\n", x="\n", title="\n") + geom_hline(yintercept=90, linetype="dashed", color = "black") + coord_flip()

jpeg(paste0("Output_Data/ConfMatrix_BestPred_Metrics.jpeg"), width=1700, height=2600, res=300)
plot(p3)
dev.off()

df.pred.best<-df.perf[order(df.perf$Accuracy, decreasing=T),][1,]
pred.name.best<-df.pred.best[,"Pred.Name"]
pred.type.best<-df.pred.best[,"Pred.Type"]
bestname2get<-paste0(pred.name.best,"_",pred.type.best)

pr.best.pred<-pr.curve(scores.class0=prob.bal.pred, weights.class0=ifelse(y==pheno2pred, 1, 0), curve=T)
pr.best.pred.df<-data.frame("Recall"=pr.best.pred$curve[,1], "Precision"=pr.best.pred$curve[,2], "Predictor"=rep(paste0("Pred: ",pred.name.best),length(pr.best.pred$curve[,2])))
pr.best.pred.df$Predictor<-as.character(pr.best.pred.df$Predictor)
pr.best.pred.df$Predictor<-paste0(pr.best.pred.df$Predictor,"\n(PRAUC=", round(pr.best.pred$auc.integral,2), ")")

jpeg(paste0("Output_Data/IMID_PRAUC.jpeg"), width=1700, height=1700, res=300)
ggplot(pr.best.pred.df, aes(x=Recall, y=Precision, group=Predictor)) + geom_line(size=0.7, alpha=1, aes(color=Predictor)) + theme_classic() +
  labs(title="\nPrecision-Recall Curve\n", x="\nRecall\n", y="\nPrecision\n") + ylim(0,1) + xlim(0,1) + scale_color_manual(values="maroon4") +
  theme(legend.title=element_blank(), legend.position = c(.5,.1), axis.text=element_text(size=12), plot.title=element_text(hjust=0.5, face="bold", size=15))
dev.off()

best.pred.rc<-df.out.pred$RF.prob.bal
best.pred.pred<-prediction(best.pred.rc, df.out.pred$Obs.Status)
best.pred.perf<-performance(best.pred.pred, "tpr", "fpr")
best.pred.rocauc<-performance(best.pred.pred, measure = "auc")@y.values[[1]]

df.roc<-data.frame("Pred"=paste0("Pred: ",pred.name.best,"\n(ROC=",round(best.pred.rocauc,2),")"),
                   "FPR"=best.pred.perf@x.values[[1]],
                   "TPR"=best.pred.perf@y.values[[1]])

jpeg(paste0("Output_Data/IMID_ROCAUC.jpeg"), width=1700, height=1700, res=300)
ggplot(df.roc, aes(FPR, TPR, color=Pred)) + geom_line(size=0.7, alpha=1) + theme_classic() + ylim(0,1) + xlim(0,1) + scale_color_manual(values="maroon4") +
  theme(legend.title=element_blank(), legend.position=c(.83,.1), axis.text=element_text(size=12), plot.title=element_text(hjust=0.5, face="bold", size=15)) +
  labs(title="\nROC\n", x="\nFalse Positive Rate (FPR)\n", y="\nTrue Positive Rate (TPR)\n")
dev.off()


