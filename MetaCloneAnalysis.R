# Loading libraries

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

pheno.name<-"IMID"
pheno1<-"RA"
pheno2<-"CTRL"
chain<-"TRD"

# In order to perform the case-case analysis, the previous variables and the code shown below need to be accordingly modified

message("\n\nLoading count input data ...")

grr<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',chain,'.rds'))
meta<-dplyr::filter(grr$meta, (IMID=="CTRL" | IMID=="RA"));
data2<-grr$data[match(meta$Sample, names(grr$data))]; 
grr<-list(data2, meta);   names(grr)<-c("data","meta")

ct.data<-as.data.frame(pubRep(grr$data, "aa", .verbose=F))
clones<-ct.data[,1]

# Computing network & Identifying clusters

message(paste0("\n\n",chain," - Computing network..."))

graph<-mutation.network(clones, .method="lev", .max.errors=1)
net.comp<-length(components(graph)$csize)
net.comp.size<-table(components(graph)$csize)

message(paste0("\n\n",chain," - Identifying clusters..."))
cfg<-cluster_fast_greedy(graph)
cfg.mod<-modularity(cfg)
vect.clust<-cfg$membership
names(vect.clust)<-clones

#  Cluster mapping

message(paste0("\n\n",chain," - Cluster mapping..."))

clust.list<-split(names(vect.clust), vect.clust)
names(clust.list)[which(lengths(clust.list)>1)]<-paste0("Clust_",which(lengths(clust.list)>1))
names(clust.list)[which(lengths(clust.list)==1)]<-unlist(clust.list[which(lengths(clust.list)==1)])

saveRDS(clust.list, paste0("Output_Data/MetaClone_",pheno.name,"_",chain,"_cluster_mapping.rds"))

# Cluster-based count matrix

message(paste0("\n\n",chain," - Creating cluster-baseed count matrix..."))

clust.data<-clust.list
rownames(ct.data)<-ct.data$CDR3.aa;  ct.data<-ct.data[,c(-1,-2)];  ct.data<-t(ct.data)

clust.size<-vector(mode="list", length=length(clust.data))
for (i in 1:length(clust.size)) {
  names(clust.size)[i]<-names(clust.data)[i]
  clust.size[[i]]<-length(clust.data[[i]])
}

lst2<-unlist(clust.size[grep("Clust",names(clust.size))], use.names=T)
df2density<-data.frame("Cluster"=names(lst2), "Clones"=unname(lst2))

barplot.clust.cont<-ggplot(df2density, aes(x=reorder(Cluster, -Clones), y=Clones, fill=reorder(Cluster, -Clones))) + 
                    geom_bar(stat="identity", show.legend=F) + theme_classic() + theme(axis.text.x=element_text(angle=90, hjust=1, size=2), plot.title=element_text(hjust=0.5, face="bold", size=15)) + 
                    labs(title=paste0("\nCluster Clone Content - ",chain,"\n"), x=paste0("\nClusters (N=",length(lst2),")\n"), y="\nClones (N)\n")

jpeg(paste0("Output_Data/MetaClone_",pheno.name,"_",chain,"_cluster_content.jpeg"), width=2800, height=1600, res=300)
grid.arrange(barplot.clust.cont)
dev.off()

# df1: clones
df1<-ct.data[,names(clust.size[clust.size==1])]
df1[is.na(df1)]<-0

# df2: clusters
clusts<-names(clust.size[clust.size>1])
ct.data2filter<-ct.data[,unname(unlist(clust.data[names(clust.size[clust.size>1])]))]
df2<-matrix(ncol=length(clusts), nrow=nrow(ct.data2filter)) # Treballem en la matriu i no en dataframe, és més ràpid
rownames(df2)<-rownames(ct.data2filter)
colnames(df2)<-clusts
system.time(
  for (i in 1:length(clusts)) {df2[,i]<-apply(ct.data2filter[, clust.data[[clusts[i]]]], 1, function(x) (sum(x, na.rm=T)))}
)[[1]]

# df3: merge
df.exp<-merge(df2,df1,by='row.names')
rownames(df.exp)<-df.exp$Row.names
df.exp<-df.exp[,-1]
saveRDS(df.exp, paste0(dir.out,"/",chain,"_cluster_exp.rds"))

#  Object creation for Hurdle test

message(paste0("\n\n",chain," - Creating object for Hurdle test..."))

#  Clone Info
clone.info<-rowSums(t(df.exp)>0)
clone.info<-data.frame("primerid"=names(clone.info), "Samples"=unname(clone.info))
colnames(clone.info)[1]<-"primerid"

#  Clone Reads
count.info<-as.matrix(df.exp);  dimnames(count.info)<-list(rownames(df.exp), colnames(df.exp))
reads.sample<-colSums(t(count.info), na.rm=T)
count.info.norm<-t(apply(t(count.info), 1, function(x) (x/reads.sample)*1000000))
count.info.cpm.plus<-count.info.norm+1
count.info.cpm.log<-log2(count.info.cpm.plus)

#  Sample Data --> cData
div.data<-readRDS("Data2analyze/DivStatistics.rds")[,c(2,1,6,5,3)]
div.data.chain<-div.data[div.data$Chain==chain,]
sample.info<-div.data.chain[(div.data.chain$Sample %in% grr$meta$Sample),]
colnames(sample.info)[2]<-"wellKey"

#  Creating Assay object
vbeta<-FromMatrix(count.info.cpm.log, cData=sample.info, fData=clone.info)
colData(vbeta)$CDR<-scale(colSums(assay(vbeta)>0, na.rm=T))[,1]

#  Creating vbeta by the phenotype to compare
pheno.values<-as.data.frame(grr$meta[match(colData(vbeta)$wellKey, as.character(grr$meta$Sample)),"IMID"])[,1]
colData(vbeta)$Pheno<-as.factor(pheno.values)

#  Filtering out clones detected in <1% of samples
filt.value<-round(ncol(assay(vbeta))*0.01)
n.clones.pheno<-length(which(rowSums(assay(vbeta)>0)!=0))
clones2keep<-rowData(vbeta)[which(rowData(vbeta)$Samples>filt.value),"primerid"]
vbeta<-vbeta[clones2keep,]

# Statistical analysis
message(paste0("\n\n",chain," - Statistical analysis..."))
options(warn=-1)
zlm.output<-zlm(~Pheno+Gender+Age+CDR, vbeta, method='bayesglm', ebayes=T)
contrasts<-unique(summary(zlm.output)$datatable$contrast)
comp<-contrasts[grep("Pheno",contrasts)]
summaryDt<-as.data.frame(summary(zlm.output, doLRT=comp)$datatable)
options(warn=0)
# hurdle results
pvals<-unlist(summaryDt[summaryDt$component=="H",4])
names(pvals)<-unlist(summaryDt[summaryDt$component=="H",1])
pvals.sorted<-sort(pvals)
fdrs<-p.adjust(pvals.sorted, method="fdr")
hurdle.df<-data.frame("Clone"=names(pvals.sorted), "Pval.hurdle"=unname(pvals.sorted), "FDR.hurdle"=fdrs)

# discrete results
pvals<-unlist(summaryDt[(summaryDt$component=="D" & summaryDt$contrast==comp),4])
names(pvals)<-unlist(summaryDt[(summaryDt$component=="D" & summaryDt$contrast==comp),1])
pvals.sorted<-sort(pvals)
fdrs<-p.adjust(pvals.sorted, method="fdr")
disc.df<-data.frame("Clone"=names(pvals.sorted), "Pval.disc"=unname(pvals.sorted), "FDR.disc"=fdrs)

# continuous results
pvals<-unlist(summaryDt[(summaryDt$component=="C" & summaryDt$contrast==comp),4])
names(pvals)<-unlist(summaryDt[(summaryDt$component=="C" & summaryDt$contrast==comp),1])
pvals.sorted<-sort(pvals)
fdrs<-p.adjust(pvals.sorted, method="fdr")
cont.df<-data.frame("Clone"=names(pvals.sorted), "Pval.cont"=unname(pvals.sorted), "FDR.cont"=fdrs)

# merging data.frames
dfm1<-merge(hurdle.df,disc.df,by="Clone")
dfm<-merge(dfm1,cont.df,by="Clone")
dfm<-dfm[order(dfm$Pval.hurdle),]
saveRDS(dfm, paste0("Output_Data/MetaClone_",chain,"_",pheno.name,"_statistics.rds"))

# Summary statistical analysis
features1<-c("Chain", "Clinical Phenotype", "Reads (n)", "Clones (n)", "Clusters (n)","Reads per sample (m)", "CDR per sample (m)", "Filtered clones (>1% samples, n)", "Filtered clusters (>1% samples, n)")
features2<-c("Hurdle", "Continuous", "Discrete")
total.chain.reads<-sum(df.exp[colData(vbeta)$wellKey,], na.rm=T)
total.chain.clones<-sum(unlist(clust.size))
total.chain.clusters<-length(grep("Clust", names(clust.size)))
reads.sample<-round(sum(df.exp[colData(vbeta)$wellKey,], na.rm=T)/nrow(colData(vbeta)),2)
cdr.sample<-round(mean(colSums(assay(vbeta)>0)),2)
chain.clones.filt<-length(grep("Clust_",names(pvals.sorted), invert=T))
chain.clust.filt<-length(grep("Clust_",names(pvals.sorted)))
h.sigP<-length(which(dfm$Pval.hurdle<0.05))
h.sigFDR<-length(which(dfm$FDR.hurdle<0.05))
d.sigP<-length(which(dfm$Pval.disc<0.05))
d.sigFDR<-length(which(dfm$FDR.disc<0.05))
c.sigP<-length(which(dfm$Pval.cont<0.05))
c.sigFDR<-length(which(dfm$FDR.cont<0.05))
h.sig<-paste0(h.sigP," (",h.sigFDR,")")
c.sig<-paste0(c.sigP," (",c.sigFDR,")")
d.sig<-paste0(d.sigP," (",d.sigFDR,")")

values1<-c(chain, paste0(pheno.name,": ",pheno1," vs ",pheno2), total.chain.reads, total.chain.clones, total.chain.clusters, reads.sample, cdr.sample, chain.clones.filt, chain.clust.filt)
values2<-c(h.sig, c.sig, d.sig)
summ.table.1<-data.frame("Feature iRep Analysis"=features1, "Value"=values1)
summ.table.2<-data.frame("Significant clones"=features2, "P<0.05 (FDR<0.05)"=values2)
colnames(summ.table.1)<-c("Feature iRep Analysis", "Value")
colnames(summ.table.2)<-c("Significant clones", "P<0.05 (FDR<0.05)")
st1<-tableGrob(summ.table.1, rows=NULL)
st2<-tableGrob(summ.table.2, rows=NULL)
st<-gtable_combine(st1, st2, along=2)

dfm.h.sorted<-dfm[order(dfm$Pval.hurdle),]
dfm.h.top<-dfm.h.sorted[c(1:12),c(1,2,4,6)]
dfm.melt<-melt(dfm.h.top)
p.top.hurdle<-ggplot(dfm.melt, aes(x=Clone, y=-log10(value))) + 
  geom_bar(stat='identity') + facet_wrap(~variable) + 
  coord_flip() + ggtitle("\nTop-10 Hurdle\n") + theme_bw() +
  xlab("") + ylab("\n-Log10(Pval)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold"), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

dfm.d.sorted<-dfm[order(dfm$Pval.disc),]
dfm.d.top<-dfm.d.sorted[c(1:12),c(1,2,4,6)]
dfm.melt<-melt(dfm.d.top)
p.top.disc<-ggplot(dfm.melt, aes(x=Clone, y=-log10(value))) + 
  geom_bar(stat='identity') + facet_wrap(~variable) + 
  coord_flip() + ggtitle("Top-10 Discrete\n") + theme_bw() +
  xlab("") + ylab("\n-Log10(Pval)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold"), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

dfm.c.sorted<-dfm[order(dfm$Pval.cont),]
dfm.c.top<-dfm.c.sorted[c(1:12),c(1,2,4,6)]
dfm.melt<-melt(dfm.c.top)
p.top.cont<-ggplot(dfm.melt, aes(x=Clone, y=-log10(value))) + 
  geom_bar(stat='identity') + facet_wrap(~variable) + 
  coord_flip() + ggtitle("Top-10 Continuous\n") + theme_bw() +
  xlab("") + ylab("\n-Log10(Pval)\n") +
  theme(plot.title=element_text(hjust=0.5, face="bold"), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

jpeg(paste0("Output_Data/MetaClone_",chain,"_",pheno.name,"_summary.jpeg"), res=300, width=3000, height=3000)
grid.arrange(st, p.top.hurdle, p.top.disc, p.top.cont, ncol=2)
dev.off()

# Graphical representation significant meta-clones in response to TNFi phenotype

clust.map<-readRDS("Output_Data/MetaClone_RespGM_",chain,"_cluster_mapping.rds")
clust.clones<-clust.map$Clust_9
graph<-mutation.network(clust.clones, .method="lev", .max.errors=1)
deg<-degree(graph, mode="all")
V(graph)$size<-deg*1.5
V(graph)$frame.color<-"white"
V(graph)$color<-"orange"
E(graph)$arrow.mode<-0
jpeg("Output_Data/MetaClone_Network_Sig_RespTNF.jpeg", res=300, height=3000, width=3000)
plot(graph, layout=layout_with_kk)
dev.off()

# Meta-clone characterization. The variables shown below should be modified accordint to the clinical phenotype of interest.

pheno.name<-"IMID"
pheno1<-"RA"
pheno2<-"CTRL"

dir.data="Output_Data"

files<-c(paste0("MetaClone_TRA_",pheno.name,"_statistics.rds"), paste0("MetaClone_TRB_",pheno.name,"_statistics.rds"), paste0("MetaClone_TRD_",pheno.name,"_statistics.rds"), paste0("MetaClone_TRG_",pheno.name,"_statistics.rds"), paste0("MetaClone_IGL_",pheno.name,"_statistics.rds"), paste0("MetaClone_IGH_",pheno.name,"_statistics.rds"), paste0("MetaClone_IGK_",pheno.name,"_statistics.rds"))

make.plot.cont<-function(x, m, chain, exp, meta, pheno.name){
  clusts<-as.character(x$Clone)
  for (cl in clusts) {
    clust.content<-chain.clust.map[cl]
    clust.name<-names(clust.content)
    clust.size<-length(clust.content[[clust.name]])
    clust.exp<-as.data.frame(exp[clust.name,]);  colnames(clust.exp)<-clust.name;
    df.exp<-merge(clust.exp,metadata,by='row.names')
    
    if (m=="Cont") {
      pval<-x[x$Clone==cl ,c("Pval.cont","FDR.cont")]
      pval2plot<-paste0("P=",formatC(pval$Pval.cont, format="e", digits=2), " ; FDR=", formatC(pval$FDR.cont, format="e", digits=2))
    }
    if (m=="Hurdle") {
      pval<-x[x$Clone==cl ,c("Pval.hurdle","FDR.hurdle")]
      pval2plot<-paste0("P=",formatC(pval$Pval.hurdle, format="e", digits=2), " ; FDR=", formatC(pval$FDR.hurdle, format="e", digits=2))
    }
    plot.irep<-ggplot(df.exp, aes(x=df.exp[,pheno.name], y=df.exp[,2], fill=df.exp[,pheno.name])) + theme_classic() + geom_jitter(shape=20, position=position_jitter(0)) + 
               theme(legend.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold")) + ylab("\nLog2(CPM)\n") + xlab("") + ylim(-10, max(df.exp[,2]*2)) +
               ggtitle(paste0("\nDifferentially Expressed MetaClone - ",clust.name,"\n")) + geom_violin(trim=F) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red", size=0.1) +
               scale_fill_manual(values=c("#E69F00", "#56B4E9")) + annotate(geom="text", x=1.5, y=max(df.exp[,2]*1.85), label=pval2plot, color="black", size=4)
    jpeg(paste0("Output_Data/MetaClone_",pheno.name,"_",chain,"_",model,"_",clust.name,".jpeg"), res=300, height=2400, width=2200)
    print(plot.irep)
    dev.off()
  }
}
make.plot.disc<-function(x, m, chain, exp, meta, pheno.name){
  clusts<-as.character(x$Clone)
  for (cl in clusts) {
    clust.content<-chain.clust.map[cl]
    clust.name<-names(clust.content)
    clust.size<-length(clust.content[[clust.name]])
    clust.exp<-as.data.frame(exp[clust.name,]);  colnames(clust.exp)<-clust.name;
    df.exp<-merge(clust.exp,metadata,by='row.names')
    df.exp$Detection<-df.exp[,2]
    df.exp$Detection[df.exp$Detection>0]<-"Detected"
    df.exp$Detection[df.exp$Detection==0]<-"Undetected"
    
    table.det<-table(df.exp$IMID, as.character(df.exp$Detection))
    df.det<-as.data.frame(table.det)
    colnames(df.det)<-c("IMID","Clust.Det","Freq")
    
    pval<-x[x$Clone==cl ,c("Pval.disc","FDR.disc")]
    pval2plot<-paste0("P=",formatC(pval$Pval.disc, format="e", digits=2), " ; FDR=", formatC(pval$FDR.disc, format="e", digits=2))
    plot.irep<-ggplot(data=df.det, aes(x=IMID, y=Freq, fill=Clust.Det)) + geom_bar(stat="identity", position=position_dodge(preserve='single')) +
      geom_text(aes(label=Freq, y=Freq+3), position=position_dodge(1), size=3.5) + scale_color_manual(values=c("black","black")) +
      scale_fill_manual(values=c("#E69F00", "#56B4E9")) + theme_classic() +
      labs(title=paste0("\nDifferentially Expressed MetaClone - ",clust.name,"\n"), x="", y="\nIndividuals (N)\n") + ylim(0,90) +
      theme(legend.position="top", legend.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold", size=16)) +
      annotate(geom="text", x=1.5, y=88, label=pval2plot, color="black")
    
    jpeg(paste0("Output_Data/MetaClone_",pheno.name,"_",chain,"_",model,"_",clust.name,".jpeg"), res=300, height=2800, width=2500)
    print(plot.irep)
    dev.off()
  }
}

# Violin Plots
for (i in 1:length(files)){
  x1<-readRDS(paste0(dir.data,"/",files[i]))
  chain<-unlist(strsplit(files[i],"_"))[1]
  
  print(chain)
  
  x.h<-x1[(x1$FDR.hurdle<0.05),]
  x.c<-x1[(x1$FDR.cont<0.05),]
  x.d<-x1[(x1$FDR.disc<0.05),]
  
  if(nrow(x.h)>0 | nrow(x.c)>0 | nrow(x.d)>0) {
    chain.clust.map<-readRDS(paste0(dir.data,"/",chain,"_cluster_mapping.rds"))
    metadata<-read.csv("Data2analyze/Metadata.csv", sep="\t")
    metadata$RespGM<-metadata$Resp
    metadata$RespGM[metadata$RespGM=="MOD"]<-"GOOD"
    metadata$Resp[metadata$Resp=="MOD"]<-NA

    count.info<-as.matrix(readRDS(paste0(dir.data,"/",chain,"_cluster_exp.rds")))
    reads.sample<-colSums(t(readRDS(paste0(dir.data,"/",chain,"_cluster_exp.rds"))), na.rm=T)
    count.info.norm<-t(apply(t(count.info), 1, function(x) (x/reads.sample)*1000000))
    count.info.cpm.plus<-count.info.norm+1
    count.info.cpm.log<-log2(count.info.cpm.plus)
    
    if (nrow(x.h)>0) {
      model<-"Hurdle"
      make.plot.cont(x.h, model, chain, count.info.cpm.log, metadata, pheno.name)
    }
    if (nrow(x.c)>0) {
      model<-"Cont"
      make.plot.cont(x.c, model, chain, count.info.cpm.log, metadata, pheno.name)
    }
    if (nrow(x.d)>0) {
      model<-"Disc"
      make.plot.disc(x.d, mode, chain, count.info.cpm.log, metadata, pheno.name)
    }
  } else {print(paste0("No significant cluster for ",chain))}
}

# CDR Plots
for (i in 1:length(files)){
  x1<-readRDS(paste0(dir.data,"/",files[i]))
  chain<-unlist(strsplit(files[i],"_"))[1]
  
  print(chain)
  
  x.h<-x1[(x1$FDR.hurdle<0.05),]
  x.c<-x1[(x1$FDR.cont<0.05),]
  x.d<-x1[(x1$FDR.disc<0.05),]
  
  metadata<-readRDS(paste0("Data2analyze/Metadata_",chain,".rds"), head=T, row.names=1)
  metadata$RespGM<-metadata$Resp
  metadata$RespGM[metadata$RespGM=="MOD"]<-"GOOD"
  metadata$Resp[metadata$Resp=="MOD"]<-NA
  
  if(nrow(x.h)>0 | nrow(x.c)>0 | nrow(x.d)>0) {
    sig.all<-rbind(x.h,x.c,x.d)
    chain.clust.map<-readRDS(paste0(dir.data,"/",chain,"_cluster_mapping.rds"))
    
    count.info<-as.matrix(readRDS(paste0(dir.data,"/",chain,"_cluster_exp.rds")))
    reads.sample<-colSums(t(readRDS(paste0(dir.data,"/",chain,"_cluster_exp.rds"))), na.rm=T)
    count.info.norm<-t(apply(t(count.info), 1, function(x) (x/reads.sample)*1000000))
    count.info.cpm.plus<-count.info.norm+1
    count.info.cpm.log<-log2(count.info.cpm.plus)
    
    div.data<-readRDS("Data2analyze/DivStatistics.rds")[,c(2,1,6,5,3)]
    div.data.chain<-div.data[div.data$Chain==chain,]
    sample.info<-div.data.chain[(div.data.chain$Sample %in% colnames(count.info.cpm.log)),]
    colnames(sample.info)[2]<-"wellKey"
    
    clone.info<-rowSums(count.info.cpm.log>0)
    clone.info<-data.frame("primerid"=names(clone.info), "Samples"=unname(clone.info))
    colnames(clone.info)[1]<-"primerid"
    
    vbeta<-FromMatrix(count.info.cpm.log, cData=sample.info, fData=clone.info)
    colData(vbeta)$CDR<-scale(colSums(assay(vbeta)>0, na.rm=T))[,1]
    
    pheno.values<-as.data.frame(metadata[match(colData(vbeta)$wellKey, as.character(rownames(metadata))),pheno.name])[,1]
    colData(vbeta)$Pheno<-pheno.values
    
    ##  Keeing significant meta-clones
    vbeta<-vbeta[as.character(sig.all$Clone),]
    
    ##  Plot
    flat.dat<-data.frame("Meta.Clone"=rep(rowData(vbeta)[,"primerid"], nrow(colData(vbeta))),
                         "Expr"=unname(assay(vbeta)[1,]),
                         "CDR"=colData(vbeta)$CDR,
                         "Pheno"=colData(vbeta)$Pheno)
    
    cdr.plot<-ggplot(flat.dat, aes(x=CDR, y=Expr, color=Pheno)) + facet_wrap(~Meta.Clone) + theme_bw() +
              ggtitle("\nClone abundance by detection rate\n") + theme(plot.title=element_text(hjust=0.5, face="bold")) +
              geom_jitter() + xlab('\nStandarized MetaClone Detection Rate (CDR)\n') + 
              ylab('\nLog2(CPM)\n') + scale_color_manual(values=c("#E69F00", "#56B4E9"))
    
    jpeg(paste0("Output_Data/MetaClone_",pheno.name,"_CDR_",chain,".jpeg"), res=300, height=2000, width=2000)
    plot(cdr.plot)
    dev.off()
  } else {print(paste0("No significant cluster for ",chain))}
}

####  MetaClone Content  ####

sig.clusts<-list()
count=0
for (i in 1:length(files)){
  x1<-readRDS(paste0(dir.data,"/",files[i]))
  chain<-unlist(strsplit(files[i],"_"))[1]
  x.h<-x1[(x1$FDR.hurdle<0.05),]
  x.c<-x1[(x1$FDR.cont<0.05),]
  x.d<-x1[(x1$FDR.disc<0.05),]
  sig.all<-rbind(x.h, x.c, x.d)
  if (nrow(sig.all)>0) {count=count+1;  sig.clusts[[count]]<-paste0(chain,"_",sig.all$Clone)}
}

sig.clusts<-unlist(sig.clusts)
clust.cont.out<-data.frame("Sig.Clusts"=sig.clusts, "Clones.N"=rep(2,length(sig.clusts)), stringsAsFactors=F)
count=0
for (cl in sig.clusts) {
  count=count+1
  cl.data<-unlist(strsplit(cl,"_"))
  chain<-cl.data[1];  cl.num<-cl.data[3];  cl.name<-paste0("Clust_",cl.num)
  cl.map.data<-readRDS(paste0(dir.data,"/",chain,"_cluster_mapping.rds"))
  cl.cont<-cl.map.data[which(names(cl.map.data)==cl.name)][[1]]
  cl.clones<-length(cl.cont)
  clust.cont.out[count,2]<-cl.clones
}
saveRDS(clust.cont.out, paste0("Output_Data/MetaClone_",pheno.name,"_Clones_per_Cluster.rds"))

# Sequence similarity networks

sig.clusts<-list()
count=0
for (i in 1:length(files)){
  x1<-readRDS(paste0(dir.data,"/",files[i]))
  chain<-unlist(strsplit(files[i],"_"))[1]
  x.h<-x1[(x1$FDR.hurdle<0.05),]
  x.c<-x1[(x1$FDR.cont<0.05),]
  x.d<-x1[(x1$FDR.disc<0.05),]
  sig.all<-rbind(x.h, x.c, x.d)
  if (nrow(sig.all)>0) {count=count+1;  sig.clusts[[count]]<-paste0(chain,"_",sig.all$Clone)}
}
sig.clusts<-unlist(sig.clusts)
for (cl in sig.clusts) {
  print(cl)
  cl.data<-unlist(strsplit(cl,"_"))
  chain<-cl.data[1];  cl.num<-cl.data[3];  cl.name<-paste0("Clust_",cl.num)
  cl.map.data<-readRDS(paste0(dir.data,"/",chain,"_cluster_mapping.rds"))
  cl.cont<-cl.map.data[which(names(cl.map.data)==cl.name)][[1]]
  if (length(cl.cont)<5000) {
    graph<-mutation.network(cl.cont, .method="lev", .max.errors=1)
    deg<-degree(graph, mode="all")
    V(graph)$size<-deg*1.5
    V(graph)$frame.color<-"white"
    V(graph)$color<-"orange"
    E(graph)$arrow.mode<-0
    jpeg(paste0("Output_Data/MetaClone_",pheno.name,"_Network_",chain,"_MetaClone_",cl.num,".jpeg"), res=300, height=3000, width=3000)
    plot(graph, layout=layout_with_kk)
    dev.off()
  }
}

####  Amino acid motifs  ####
sig.clusts<-list()
count=0
for (i in 1:length(files)){
  x1<-readRDS(paste0(dir.data,"/",files[i]))
  chain<-unlist(strsplit(files[i],"_"))[1]
  x.h<-x1[(x1$FDR.hurdle<0.05),]
  x.c<-x1[(x1$FDR.cont<0.05),]
  x.d<-x1[(x1$FDR.disc<0.05),]
  sig.all<-rbind(x.h, x.c, x.d)
  if (nrow(sig.all)>0) {count=count+1;  sig.clusts[[count]]<-paste0(chain,"_",sig.all$Clone)}
}
sig.clusts<-unlist(sig.clusts)
for (cl in sig.clusts) {
  print(cl)
  cl.data<-unlist(strsplit(cl,"_"))
  chain<-cl.data[1];  cl.num<-cl.data[3];  cl.name<-paste0("Clust_",cl.num)
  cl.map.data<-readRDS(paste0(dir.data,"/",chain,"_cluster_mapping.rds"))
  cl.cont<-cl.map.data[which(names(cl.map.data)==cl.name)][[1]]
  count=0
  count.name=0
  out.data<-c()
  for (clone in cl.cont) {
    count=count+1
    count.name=count.name+1
    out.data[count]<-paste0(">Clone",count.name,"|Clone",count.name,"|Clone",count.name)
    out.data[count+1]<-clone
    count=count+1
  }
  df.out.data<-as.data.frame(out.data)
  write.table(df.out.data, paste0("Output_Data/MetaClone_",pheno.name,"_Seq_",cl,".txt"), quote=F, col.names=F, row.names=F)
  seqs<-readAAStringSet(paste0("Output_Data/MetaClone_",pheno.name,"_Seq_",cl,".txt"))
  seqs.algn<-msa(seqs, "ClustalW")
  conMat<-consensusMatrix(seqs.algn)
  jpeg(paste0("Output_Data/MetaClone_",pheno.name,"_AAMotif_",cl,".jpeg"), res=300, height=2500, width=3500)
  plot(ggseqlogo(conMat, method='prob'))
  dev.off()
}

# Characterization of the meta-clone 9 from TRG that showed a significant association with clinical response to TNF inhibition
# In order to execute this code, the data that is read must be obtained after performing the meta-clone association analysis with clinical response to TNF inhibition

seqs<-readAAStringSet("Output_Data/MetaClone_RespGM_Seq_TRG_Clust_9")
seqs.algn<-msa(seqs, "ClustalW")
conMat<-consensusMatrix(seqs.algn)
jpeg("Output_Data/MetaClone_RespGM_AAMotif_TRG_Clust_9.jpeg", res=300, height=2500, width=3500)
plot(ggseqlogo(conMat, method='prob'))
dev.off()

chain.clust.map<-readRDS(paste0(dir.data,"/TRG_cluster_mapping.rds"))
metadata<-readRDS("Data2analyze/Metadata_TRG.rds")
metadata$RespGM<-metadata$Resp
metadata$RespGM[metadata$RespGM=="MOD"]<-"GOOD"
metadata$Resp[metadata$Resp=="MOD"]<-NA
count.info<-as.matrix(readRDS("Output_Data/TRG_cluster_exp.rds"))
reads.sample<-colSums(t(readRDS("Output_Data/TRG_cluster_exp.rds")), na.rm=T)
count.info.norm<-t(apply(t(count.info), 1, function(x) (x/reads.sample)*1000000))
count.info.cpm.plus<-count.info.norm+1
count.info.cpm.log<-log2(count.info.cpm.plus)
clust.content<-chain.clust.map$Clust_9
clust.name<-"Clust_9"
clust.size<-length(clust.content)
clust.exp<-as.data.frame(count.info.cpm.log[clust.name,]);  colnames(clust.exp)<-clust.name;
df.exp<-merge(clust.exp,metadata,by='row.names')
pval2plot<-"\nP=3.07e-04\n"
df.cont<-df.exp[df.exp$Clust_9>0,]
df.cont$RespGM<-as.character(df.cont$RespGM)
df.cont$RespGM[df.cont$RespGM=="GOOD"]<-"Responders"
df.cont$RespGM[df.cont$RespGM=="NON_RESP"]<-"Non responders"
p1<-ggboxplot(df.cont, x="RespGM", y="Clust_9", color="RespGM", fill="RespGM") + 
  labs(x="", y="\nExpression Meta-clone 9 (Log2(CPM))\n") + scale_color_manual(values=c("black","black")) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) + theme(legend.position="none") + ylim(7.5,14.5) +
  annotate("text", x=1.5, y=14.3, label=pval2plot, size=6)
jpeg("Output_Data/MetaClone9_RespTNF_Expr.jpeg", res=300, height=2300, width=2300)
plot(p1)
dev.off()


div.data<-readRDS("Data2analyze/DivStatistics.rds")[,c(2,1,6,5,3)]
div.data.chain<-div.data[div.data$Chain==chain,]
sample.info<-div.data.chain[(div.data.chain$Sample %in% colnames(count.info.cpm.log)),]
colnames(sample.info)[2]<-"wellKey"
clone.info<-rowSums(count.info.cpm.log>0)
clone.info<-data.frame("primerid"=names(clone.info), "Samples"=unname(clone.info))
colnames(clone.info)[1]<-"primerid"
vbeta<-FromMatrix(count.info.cpm.log, cData=sample.info, fData=clone.info)
colData(vbeta)$CDR<-scale(colSums(assay(vbeta)>0, na.rm=T))[,1]
pheno.values<-as.data.frame(metadata[match(colData(vbeta)$wellKey, as.character(rownames(metadata))),"RespGM"])[,1]
colData(vbeta)$Pheno<-pheno.values
vbeta<-vbeta["Clust_9",]
flat.dat<-data.frame("Meta.Clone"=rep(rowData(vbeta)[,"primerid"], nrow(colData(vbeta))),
                     "Expr"=unname(assay(vbeta)[1,]),
                     "CDR"=colData(vbeta)$CDR,
                     "Pheno"=colData(vbeta)$Pheno)
flat.dat<-flat.dat[flat.dat$Expr>0,]
flat.dat$Pheno<-as.character(flat.dat$Pheno)
flat.dat$Pheno[flat.dat$Pheno=="GOOD"]<-"Responders"
flat.dat$Pheno[flat.dat$Pheno=="NON_RESP"]<-"Non responders"
p2<-ggplot(flat.dat, aes(x=CDR, y=Expr, color=Pheno)) + theme_classic() +
  theme(legend.title=element_blank(), plot.title=element_text(hjust=0.5, face="bold")) +
  labs(x='\nStandarized Clone Detection Rate (CDR)', y='Expression Meta-clone 9 (Log2(CPM))\n') +
  geom_jitter() + scale_color_manual(values=c("#E69F00", "#56B4E9")) + ylim(8,14)
jpeg("Output_Data/MetaClone9_RespTNF_CDR.jpeg", res=300, height=2150, width=2150)
plot(p2)
dev.off()
