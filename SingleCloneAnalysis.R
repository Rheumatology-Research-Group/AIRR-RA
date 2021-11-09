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
week<-"wk0"

# In order to perform the case-case analysis, the previous variables and the code shown below need to be accordingly modified
              
df.all<-readRDS("Data2analyze/CloneProp.rds")
        
x1<-read.csv("Data2analyze/Metadata.csv", sep="\t")
x2<-read.table("Data2analyze/Donor2RNA_Codes.txt", header=T)

ra.wk0<-as.character(x1[x1$Week=="wk0", "Sample"]);
ra.wk12<-as.character(x2$RNA2);
ctrl<-as.character(x1[x1$IMID=="CTRL", "Sample"]);
samples2keep<-c(ra.wk0,ra.wk12,ctrl)
        
df.all<-df.all[(df.all$sample.vect %in% samples2keep),]
        
    
chains<-c("TRG","TRD","TRA","TRB","IGH","IGL","IGK")
          
for (chain in chains) {
            
            #  Data preprocessing
            
            message(paste0("\n",chain," - Data preprocessing..."))
            
            #  Loading Data  --> immunarch
            inp.name<-paste0('Data2analyze/Clones_By_CDR3aa/',chain,'.rds')
            grr<-readRDS(inp.name)
            meta_col<-grr$meta
            meta_col<-meta_col[(meta_col$Sample %in% samples2keep),]
            samples_col<-meta_col$Sample;  data_col<-grr$data[match(samples_col,names(grr$data))];
            grr<-list(data_col,meta_col);  names(grr)<-c("data","meta")
            
            #  Clone Info --> fData
            count.data<-pubRep(grr$data, "aa", .verbose=F)
            clone.info<-as.data.frame(count.data[,c(1,2)])
            colnames(clone.info)[1]<-"primerid"
          
            #  Clone Reads --> df
            count.info<-as.matrix(count.data[,c(-1,-2)])
            dimnames(count.info)<-list(count.data$CDR3.aa, colnames(count.info)) 
            count.info.cpm<-apply(count.info, 2, function(x) (x/sum(x, na.rm=T))*1000000)
            count.info.cpm.plus<-count.info.cpm+1
            count.info.cpm.plus[is.na(count.info.cpm.plus)]<-1
            count.info.cpm.log<-log2(count.info.cpm.plus)
          
            #  Sample Data --> cData
            df.chain<-df.all[(df.all$chain.vect==chain),]
            colnames(df.chain)<-c("Chain","Sample","Clone","IMID","Age","Gender","Prop","Reads","Total.Reads")
            df.chain<-df.chain[,c(-3,-5,-7,-8)]
            sample.info<-df.chain[!duplicated(df.chain), ]
            colnames(sample.info)[2]<-"wellKey"
          
            #  Creating SingleCellAssay object
            vbeta<-FromMatrix(count.info.cpm.log, cData=sample.info, fData=clone.info)
            colData(vbeta)$CDR<-scale(colSums(assay(vbeta)>0, na.rm=T))[,1]
          
            #  Creating vbeta by the phenotype to compare
            metadata<-read.csv("Data2analyze/Metadata.csv", sep="\t")
            pheno.values<-as.character(metadata[match(colData(vbeta)$wellKey, metadata$Sample),pheno.name])
            colData(vbeta)$Pheno<-pheno.values
            colData(vbeta)$Pheno[is.na(colData(vbeta)$Pheno)]<-"NA"
            command<-paste0('subset(vbeta, Pheno=="',pheno1,'" | Pheno=="',pheno2,'")')
            vbeta.subset<-eval(parse(text=command))
            vbeta<-vbeta.subset
            age.values<-as.numeric(metadata[match(colData(vbeta)$wellKey, metadata$Sample),"Age"])
            colData(vbeta)$Age<-age.values
            
            colData(vbeta)$Pheno<-as.factor(colData(vbeta)$Pheno)
            colData(vbeta)$Gender<-as.factor(colData(vbeta)$Gender)
            
            #  Filtering out clones detected in <1% of samples
            filt.value<-round(ncol(assay(vbeta))*0.01)
            n.clones.pheno<-length(which(rowSums(assay(vbeta)>0)!=0))
            clones2keep<-rowData(vbeta)[which(rowData(vbeta)$Samples>filt.value),"primerid"]
            vbeta<-vbeta[clones2keep,]
          
            #  Statistical analysis
            message(paste0("\n\n",chain," - Statistical analysis..."))
            
            zlm.output<-zlm(~Pheno+Age+Gender+CDR, vbeta, method='bayesglm', ebayes=T, parallel=T)
        
            contrasts<-unique(summary(zlm.output)$datatable$contrast)
            comp<-contrasts[grep("Pheno",contrasts)]
            summaryDt<-as.data.frame(summary(zlm.output, doLRT=comp)$datatable)
          
            #  Hurdle results
            pvals<-unlist(summaryDt[summaryDt$component=="H",4])
            names(pvals)<-unlist(summaryDt[summaryDt$component=="H",1])
            pvals.sorted<-sort(pvals)
            fdrs<-p.adjust(pvals.sorted, method="fdr")
            hurdle.df<-data.frame("Clone"=names(pvals.sorted), "Pval.hurdle"=unname(pvals.sorted), "FDR.hurdle"=fdrs)
          
            #  Discrete results
            pvals<-unlist(summaryDt[(summaryDt$component=="D" & summaryDt$contrast==comp),4])
            names(pvals)<-unlist(summaryDt[(summaryDt$component=="D" & summaryDt$contrast==comp),1])
            pvals.sorted<-sort(pvals)
            fdrs<-p.adjust(pvals.sorted, method="fdr")
            disc.df<-data.frame("Clone"=names(pvals.sorted), "Pval.disc"=unname(pvals.sorted), "FDR.disc"=fdrs)
          
            #  Continuous results
            pvals<-unlist(summaryDt[(summaryDt$component=="C" & summaryDt$contrast==comp),4])
            names(pvals)<-unlist(summaryDt[(summaryDt$component=="C" & summaryDt$contrast==comp),1])
            pvals.sorted<-sort(pvals)
            fdrs<-p.adjust(pvals.sorted, method="fdr")
            cont.df<-data.frame("Clone"=names(pvals.sorted), "Pval.cont"=unname(pvals.sorted), "FDR.cont"=fdrs)
          
            #  Merging results
            dfm1<-merge(hurdle.df,disc.df,by="Clone")
            dfm<-merge(dfm1,cont.df,by="Clone")
          
            saveRDS(dfm, paste0("Output_Data/SingleClone_",chain,"_",pheno.name,"_statistics.rds"))
}


# 
# Graphical representation 

div.data<-readRDS("Data2analyze/DivStatistics.rds")[,c(2,1,6,5,3)] # This file includes the diversity statistics that can be downloaded from the website.
chains<-c("IGL","IGK");
models<-c("Hurdle","Cont","Disc")

df.models<-vector("list",length(models)); names(df.models)<-models
df.sig<-vector("list",length(chains)); names(df.sig)<-chains

for (ch in names(df.sig)) {
  df.sig[[ch]]<-df.models
  df.results<-readRDS(paste0("Output_Data/SingleClone_",chain,"_",pheno.name,"_statistics.rds"))
  for (mod in models) {
    if (mod=="Hurdle") {mod.col.pval<-2; mod.col.fdr<-3}
    if (mod=="Disc") {mod.col.pval<-4; mod.col.fdr<-5}
    if (mod=="Cont") {mod.col.pval<-6; mod.col.fdr<-7}
    sig.data<-df.results[(df.results[,mod.col.fdr]<0.05), ]
    df.sig[[ch]][[mod]]<-sig.data
  }
}

grr.tra<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/TRA.rds'))
grr.trb<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/TRB.rds'))
grr.trd<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/TRD.rds'))
grr.trg<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/TRG.rds'))
grr.igh<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/IGH.rds'))
grr.igl<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/IGL.rds'))
grr.igk<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/IGK.rds'))

pdf("Output_Data/SingleClone_Plots.pdf", width=10)
for (ch in names(df.sig)) {
  print(ch)
  if (ch=="TRA") {grr<-grr.tra}
  if (ch=="TRB") {grr<-grr.trb}
  if (ch=="TRD") {grr<-grr.trd}
  if (ch=="TRG") {grr<-grr.trg}
  if (ch=="IGH") {grr<-grr.igh}
  if (ch=="IGL") {grr<-grr.igl}
  if (ch=="IGK") {grr<-grr.igk}
  
  count.data<-pubRep(grr$data, "aa", .verbose=F);
  clone.info<-as.data.frame(count.data[,c(1,2)]);  colnames(clone.info)[1]<-"primerid"
  
  count.info<-as.matrix(count.data[,c(-1,-2)]);  dimnames(count.info)<-list(count.data$CDR3.aa, colnames(count.info)) 
  count.info.cpm<-apply(count.info, 2, function(x) (x/sum(x, na.rm=T))*1000000)
  count.info.cpm.plus<-count.info.cpm+1
  count.info.cpm.plus[is.na(count.info.cpm.plus)]<-1
  count.info.cpm.log<-log2(count.info.cpm.plus)
  
  sample.info<-div.data[(div.data$Chain==ch),];
  colnames(sample.info)[2]<-"wellKey"
  
  count.info.cpm.log<-count.info.cpm.log[,colnames(count.info.cpm.log) %in% sample.info$wellKey]
  
  vbeta<-FromMatrix(count.info.cpm.log, cData=sample.info, fData=clone.info)
  colData(vbeta)$CDR<-scale(colSums(assay(vbeta)>0, na.rm=T))[,1]
  
  metadata<-readRDS(paste0("Data2analyze/Metadata_",ch,".rds"))
  metadata$RespGM<-metadata$Resp
  metadata$RespGM[metadata$RespGM=="MOD"]<-"GOOD"
  metadata$Resp[metadata$Resp=="MOD"]<-NA
  
  for (mod in models) {
    data<-df.sig[[ch]][[mod]]
    colnames(data)<-c("Clone","Pval.Hurdle","FDR.Hurdle","Pval.Disc","FDR.Disc","Pval.Cont","FDR.Cont")
    # Summary Plot
    clones<-rep(data$Clone,3*2)
    mods<-c(rep("Hurdle",2*nrow(data)), rep("Discrete",2*nrow(data)), rep("Continuous",2*nrow(data)))
    stats<-c(rep("Pval",nrow(data)), rep("FDR",nrow(data)), rep("Pval",nrow(data)), rep("FDR",nrow(data)), rep("Pval",nrow(data)), rep("FDR",nrow(data)))
    values<-c(data[,"Pval.Hurdle"], data[,"FDR.Hurdle"], data[,"Pval.Disc"], data[,"FDR.Disc"],  data[,"Pval.Cont"],data[,"FDR.Cont"])
    df2plot<-data.frame("Clone"=clones, "Model"=mods, "Significance"=stats, "Value"=values)

    for (i in 1:nrow(data)) {
      clone<-data[i,"Clone"];    pheno<-"IMID";    pval<-data[i,paste0("Pval.",mod)];
      
      col.pheno<-which(colnames(metadata)==pheno)
      meta2use<-metadata[,c(1,col.pheno)]
      
      pheno.values<-meta2use[match(colData(vbeta)$wellKey, meta2use$Sample),]
      colData(vbeta)$Pheno<-pheno.values[,pheno]
      
      dupl.sample.exc<-c("IXTCB01476","IXTCB01395","IXTCB01434")
      
      vbeta.s1<-subset(vbeta, wellKey!=dupl.sample.exc[1])
      vbeta.s2<-subset(vbeta.s1, wellKey!=dupl.sample.exc[2])
      vbeta.s3<-subset(vbeta.s2, wellKey!=dupl.sample.exc[3])
      vbeta.subset<-vbeta.s3[as.character(clone),]
      
      vbeta.pheno<-subset(vbeta.subset, (Pheno=="CTRL" | Pheno=="RA"))
      flat.dat.pheno<-as(vbeta.pheno[as.character(clone),], 'data.table')
      
      pval.pheno.text<-paste0("P=",formatC(as.numeric(pval), format="e", digits=2))
      
      if (mod=="Disc") {
        flat.dat.pheno.count<-flat.dat.pheno
        flat.dat.pheno.count$Detection<-as.character(flat.dat.pheno.count$Pheno)
        flat.dat.pheno.count$Detection[flat.dat.pheno.count$value>0]<-"Yes"
        flat.dat.pheno.count$Detection[flat.dat.pheno.count$value==0]<-"No"
        flat.dat.pheno.count$Pheno<-droplevels(flat.dat.pheno.count$Pheno)
        
        tab.ph<-as.data.frame(table(flat.dat.pheno.count$Pheno, flat.dat.pheno.count$Detection))
        
        exp.ph<-ggplot(data=tab.ph, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat="identity", position=position_dodge()) +
          ylab("\nIndividuals with clone detected (N)\n") + xlab("") + theme_classic() + labs (fill="Detection") +
          ggtitle("\nClone Detection\n") + theme(plot.title=element_text(hjust=0.5, face="bold")) +
          annotate(geom="text", x=1.5, y=max(tab.ph$Freq)*1.1, size=5, label=pval.pheno.text, color="black")
        
        cdr.ph<-ggplot(flat.dat.pheno, aes(x=CDR, y=value, color=Pheno)) + theme_classic() +
          ggtitle("\nAbundance by CDR\n") + theme(plot.title=element_text(hjust=0.5, face="bold"), legend.title=element_blank()) +
          geom_jitter() + xlab('\nStandardized Clone Detection Rate (CDR)\n') +
          ylab('\nLog2(CPM)\n') + scale_color_manual(name=pheno, values=c("#D8EBD8","#F69191"))
        
        grid.arrange(exp.ph, cdr.ph, ncol=2, top=textGrob(paste0("\n",clone," from ",ch," chain significant in ",mod," model\n"),gp=gpar(fontsize=18,fontface="bold")))
        
      } else {
        
        exp.ph<-ggplot(flat.dat.pheno, aes(x=Pheno, y=value, fill=Pheno)) + 
          theme_classic() + geom_jitter(shape=20, position=position_jitter(0)) + 
          theme(plot.title=element_text(hjust=0.5, face="bold"), legend.title=element_blank()) + ylab("\nLog2(CPM)\n") + xlab("") +
          ggtitle("\nClone Expression\n") + geom_violin(trim=F) +
          stat_summary(fun.data=mean_sdl, geom="pointrange", color="red", size=0.1) +
          scale_fill_manual(values=c("#D8EBD8","#F69191")) +
          annotate(geom="text", x=1.5, y=max(flat.dat.pheno$value)*1.1, size=5, label=pval.pheno.text, color="black")
        
        cdr.ph<-ggplot(flat.dat.pheno, aes(x=CDR, y=value, color=Pheno)) + theme_classic() +
          ggtitle("\nAbundance by CDR\n") + theme(plot.title=element_text(hjust=0.5, face="bold"), legend.title=element_blank()) +
          geom_jitter() + xlab('\nStandardized Clone Detection Rate (CDR)\n') + 
          ylab('\nLog2(CPM)\n') + scale_color_manual(name=pheno, values=c("#D8EBD8","#F69191"))
        
        grid.arrange(exp.ph, cdr.ph, ncol=2, top=textGrob(paste0("\n",clone," from ",ch," chain significant in ",mod," model\n"),gp=gpar(fontsize=18,fontface="bold")))
      }
    }
  }
}
dev.off()

# Sequence similarity networks - significant clones IGL, IGK     
  
chains.sig<-c("IGL","IGK")

# DC
nperm<-10000000;
for (chain in chains.sig) {
  print(chain)
  res.chain<-readRDS(paste0("Output_Data/SingleClone_",chain,"_",pheno.name,"_statistics.rds"))
  sig.clones<-unique(as.character(res.chain[(res.chain$FDR.hurdle<0.05 | res.chain$FDR.disc<0.05 | res.chain$FDR.cont<0.05),"Clone"]))
  grr<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',chain,'.rds'))
  samples2keep<-as.data.frame(grr$meta[grr$meta$IMID=="RA","Sample"])[,1]
  meta_col<-grr$meta
  meta_col<-meta_col[(meta_col$Sample %in% samples2keep),]
  samples_col<-meta_col$Sample;  data_col<-grr$data[match(samples_col,names(grr$data))];
  grr<-list(data_col,meta_col);  names(grr)<-c("data","meta")
  chain.clones<-as.data.frame(pubRep(grr$data, "aa", .verbose=F))[,1]
  obs.graph<-mutation.network(sig.clones, .method="lev", .max.errors=1)
  deg<-degree(obs.graph, mode="all")
  V(obs.graph)$size<-deg*1.5
  V(obs.graph)$frame.color<-"white"
  if (chain=="IGL") {col<-"honeydew4"} else {col<-"honeydew2"}
  V(obs.graph)$color<-col
  E(obs.graph)$arrow.mode<-0
  E(obs.graph)$arrow.width<-0.1
  V(obs.graph)$label.color<-"gray1"
  jpeg(paste0("Output_Data/RAassoc_singleclones_Similarity_Network_",chain,".jpeg"), res=300, height=4000, width=4000)
  plot(obs.graph, layout=layout_with_kk)
  dev.off()
  obs.deg<-mean(degree(obs.graph, mode="all"))
  set.n<-length(sig.clones)
  count.low=0;  count.hi=0; 
  perm.deg.all<-list()
  for (i in 1:nperm) {
    print(paste0(chain,"-",i))
    x<-sample(chain.clones, size=set.n, replace=F)
    perm.graph<-mutation.network(x, .method="lev", .max.errors=1)
    perm.deg<-mean(degree(perm.graph, mode="all"))
    perm.deg.all[[i]]<-perm.deg
    perm.deg<-unlist(perm.deg.all)
  }
  df.out<-data.frame("Similarity.Metric"="Degree Centrality", "Obs.Value"=obs.deg, "Perm.Value"=mean(perm.deg),
                       "PermValue.LT.ObsValue"=length(which(perm.deg<obs.deg)),"PermValue.GT.ObsValue"=length(which(perm.deg>obs.deg)),
                       "Pval"=(length(which(perm.deg>obs.deg))+1)/(nperm+1), stringsAsFactors=F)
  write.csv(df.out,paste0("Output_Data/RAassoc_singleclones_Similarity_DC_",chain,".csv"), row.names=F)
}

# AminoAcid Motifs

for (chain in chains.sig) {
  print(chain)
  res.chain<-readRDS(paste0("Output_Data/SingleClone_",chain,"_",pheno.name,"_statistics.rds"))
  sig.clones<-unique(as.character(res.chain[(res.chain$FDR.hurdle<0.05 | res.chain$FDR.disc<0.05 | res.chain$FDR.cont<0.05),"Clone"]))
  count=0
  count.name=0
  out.data<-c()
  for (clone in sig.clones) {
    count=count+1
    count.name=count.name+1
    out.data[count]<-paste0(">Clone",count.name,"|Clone",count.name,"|Clone",count.name)
    out.data[count+1]<-clone
    count=count+1
  }
  df.out.data<-as.data.frame(out.data)
  write.table(df.out.data, paste0("Output_Data/SingleClone_Seq_",chain,".txt"), quote=F, col.names=F, row.names=F)
  seqs<-readAAStringSet(paste0("Output_Data/SingleClone_Seq_",chain,".txt"))
  seqs.algn<-msa(seqs, "ClustalW")
  conMat<-consensusMatrix(seqs.algn)
  jpeg(paste0("Output_Data/SingleClone_AAMotif_",chain,".jpeg"), res=300, height=2500, width=3500)
  plot(ggseqlogo(conMat, method='prob'))
  dev.off()
}

# Biochemical properties

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
for (chain in chains.sig) {
  res.chain<-readRDS(paste0("Output_Data/SingleClone_",chain,"_",pheno.name,"_statistics.rds"))
  sig.clones<-unique(as.character(res.chain[(res.chain$FDR.hurdle<0.05 | res.chain$FDR.disc<0.05 | res.chain$FDR.cont<0.05),"Clone"]))
  df.sig.clones.prop<-data.frame("Chain"="TEST", "Clone"="TEST", "Clone.Type"="TEST",
                                 "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                                 "Aliphatic"=0, "Aromatic"=0,
                                 "Non.Polar"=0, "Polar"=0, "Charged"=0,
                                 "Basic"=0, "Acidic"=0,
                                 "Net.Charge"=0, "Isoelectric.Point"=0,
                                 "Aliphatic.Index"=0, "Instability.Index"=0,
                                 "PPI.Index"=0, "Hydrophobicity.Index"=0,
                                 "HydrophobicMoment.Index"=0, stringsAsFactors=F)
  count=0
  for (cl in sig.clones) {
    count=count+1
    cl.prop<-get.aa.features(cl)
    cl.prop.norm<-cl.prop/(lengthpep(cl))
    df.sig.clones.prop[count,]<-c(chain,cl,"Sig.Clone",cl.prop.norm)
  }
  cols.interest<-sapply(df.sig.clones.prop[,c(4:ncol(df.sig.clones.prop))], as.numeric)
  obs.prop<-colMeans(cols.interest)
  
  grr<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',chain,'.rds'))
  samples2keep<-as.data.frame(grr$meta[grr$meta$IMID=="RA","Sample"])[,1]
  meta_col<-grr$meta
  meta_col<-meta_col[(meta_col$Sample %in% samples2keep),]
  samples_col<-meta_col$Sample;  data_col<-grr$data[match(samples_col,names(grr$data))];
  grr<-list(data_col,meta_col);  names(grr)<-c("data","meta")
  chain.clones<-as.data.frame(pubRep(grr$data, "aa", .verbose=F))[,1]
  
  set.n<-length(sig.clones)
  count.low<-0;  count.hi<-0;  nperm<-10; # nperm can be modified
  df.perm.all<-data.frame("Chain"="TEST", "Clone.Type"="TEST", "Perm"=0,
                          "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                          "Aliphatic"=0, "Aromatic"=0,
                          "Non.Polar"=0, "Polar"=0, "Charged"=0,
                          "Basic"=0, "Acidic"=0,
                          "Net.Charge"=0, "Isoelectric.Point"=0,
                          "Aliphatic.Index"=0, "Instability.Index"=0,
                          "PPI.Index"=0, "Hydrophobicity.Index"=0,
                          "HydrophobicMoment.Index"=0, stringsAsFactors=F)
  for (i in 1:nperm) {
    print(paste0(chain,"-",i))
    perm.clones<-sample(chain.clones, size=set.n, replace=F)
    df.perm.clones.prop<-data.frame("Chain"="TEST", "Clone"="TEST", "Clone.Type"="TEST",
                                   "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                                   "Aliphatic"=0, "Aromatic"=0,
                                   "Non.Polar"=0, "Polar"=0, "Charged"=0,
                                   "Basic"=0, "Acidic"=0,
                                   "Net.Charge"=0, "Isoelectric.Point"=0,
                                   "Aliphatic.Index"=0, "Instability.Index"=0,
                                   "PPI.Index"=0, "Hydrophobicity.Index"=0,
                                   "HydrophobicMoment.Index"=0, stringsAsFactors=F)
    count.perm=0
    for (pc in perm.clones) {
      count.perm=count.perm+1
      cl.perm.prop<-get.aa.features(pc)
      cl.perm.prop.norm<-cl.perm.prop/(lengthpep(pc))
      df.perm.clones.prop[count.perm,]<-c(chain,cl,"Perm.Clone",cl.perm.prop.norm)
    }
    cols.interest.perm<-sapply(df.perm.clones.prop[,c(4:ncol(df.perm.clones.prop))], as.numeric)
    perm.prop<-colMeans(cols.interest.perm)
    df.perm.all[i,]<-c(chain,"Perm.Clone",i,perm.prop)
  }
  perm.prop<-df.perm.all
  
  aa.propierties<-names(obs.prop)
  
  df.out<-data.frame("Biochemical.Feature"="TEST","Obs.Value"="TEST","Perm.Value"="TEST","PermValue.LT.ObsValue"="TEST","PermValue.GT.ObsValue"="TEST","Pval.PermValue.LT.Obs.Value"="TEST","Pval.PermValue.GT.ObsValue"="TEST", stringsAsFactors=F)
  count.final=0
  for (aap in aa.propierties) {
    obs.value<-obs.prop[aap]
    perm.value<-as.numeric(perm.prop[,aap])
    pval1<-(length(which(perm.value<obs.value))+1)/(nperm+1)
    pval2<-(length(which(perm.value>obs.value))+1)/(nperm+1)
    count.final=count.final+1
    df.out[count.final,]<-c(aap, unname(obs.value), mean(perm.value), length(which(perm.value<obs.value)), length(which(perm.value>obs.value)), pval1, pval2)
  }
  saveRDS(df.out,paste0("Output_Data/SingleClone_AminoAcid_Prop_",chain,".rds"))
}

chains<-c("IGL","IGK")
df.all<-list()
count=0
for (c in chains) {
    count=count+1
    file.name<-paste0("Output_Data/SingleClone_AminoAcid_Prop_",c,".rds")
    res.data<-readRDS(file.name)
    res.data$Chain<-rep(c,nrow(res.data))
    res.data.depl<-res.data[,c(1,6,8)]
    res.data.enrich<-res.data[,c(1,7,8)]
    res.data.depl$BiochEff<-rep("Depletion",nrow(res.data.depl))
    res.data.enrich$BiochEff<-rep("Enrichment",nrow(res.data.enrich))
    colnames(res.data.depl)[2]<-"Pval"
    colnames(res.data.enrich)[2]<-"Pval"
    res.data.depl$Pval<-(log10(as.numeric(res.data.depl$Pval)))
    res.data.enrich$Pval<-(-log10(as.numeric(res.data.enrich$Pval)))
    res.data.enrich$key<-paste(res.data.enrich$Biochemical.Feature, res.data.enrich$Chain)
    res.data.depl$key<-paste(res.data.depl$Biochemical.Feature, res.data.depl$Chain)
    res2use<-merge(res.data.depl,res.data.enrich, by="key")
    res2use$enrich.score<-res2use$Pval.y+res2use$Pval.x
    res2use<-res2use[,c(2,4,10)]
    colnames(res2use)<-c("Biochemical.Feature","Chain","Biochemical.Enrichment.Score")
    df.all[[count]]<-res2use
}
data2plot<-do.call(rbind, df.all)
bp2remove<-c("Tiny","Aliphatic.Index","Non.Polar","HydrophobicMoment.Index","Small","Net.Charge")
data2plot<-data2plot[!(data2plot$Biochemical.Feature %in% bp2remove),]
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="Hydrophobicity.Index"]<-"Hydrophobic"
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="Instability.Index"]<-"Instability"
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="PPI.Index"]<-"Protein-Protein Interaction Index"
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="Molecular.Weight"]<-"Molecular weight"
data2plot$Biochemical.Feature[data2plot$Biochemical.Feature=="Isoelectric.Point"]<-"Isoelectric point"

# plot if BEP (Biochemical Enrichment Score is > -log10(0.05))

data2plot$Biochemical.Enrichment.Score[(data2plot$Biochemical.Enrichment.Score<1.3 & data2plot$Biochemical.Enrichment.Score>-1.3)]<-NA
colnames(data2plot)[3]<-"BES"

p<-ggplot(data=data2plot, aes(x=factor(Chain), y=factor(Biochemical.Feature), size=abs(BES), fill=BES)) +
  geom_point(alpha=1, shape=21, color="black") +
  scale_size(range=c(.1, 8), name="Significance") +
  labs(y="", x="") +
  scale_fill_gradientn(colours = rev(colorspace::diverge_hcl(5))) +
  theme_bw() + theme(plot.title=element_text(hjust=0.5, face="bold", size=15), axis.text=element_text(size=11))

jpeg("Output_Data/SingleClone_BiochemicalProfile_SigChains.jpeg", res=300, height=2000, width=2000)
plot(p)
dev.off()

#   ggplot(data=data2plot, aes(x=factor(TNFeff), y=factor(Biochemical.Feature), size=abs(BES), fill=BES)) +
#   geom_point(alpha=1, shape=21, color="black") +
#   scale_size(range=c(.1, 8), name="Significance") +
#   facet_wrap(~Chain) + labs(y="", x="") +
#   scale_fill_gradientn(colours = rev(colorspace::diverge_hcl(15))) +
#   theme_bw() + theme(plot.title=element_text(hjust=0.5, face="bold", size=16), axis.text.x=element_text(vjust=0.5, size=13), strip.text=element_text(size=14, face="bold"))


