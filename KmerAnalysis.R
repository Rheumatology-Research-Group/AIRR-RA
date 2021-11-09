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


# Case-control and case-case association analysis

chains<-c("TRG","TRD","TRA","TRB","IGH","IGL","IGK")
phenos<-c("IMID","RespGM","Resp","ActB","ACPA","RF")

for (chain in chains) {
  
  # Loading Data
  
  print(paste0("Loading Data ", chain))
  
  grr<-readRDS(paste0('Output_Data/Clones_By_CDR3aa/',chain,'.rds'))
  grr$meta$RespGM<-grr$meta$Resp
  grr$meta$RespGM[grr$meta$RespGM=="MOD"]<-"GOOD"
  grr$meta$Resp[grr$meta$Resp=="MOD"]<-NA
  
  for (pheno in phenos) {
    
    # Preparing Pheno Data
    
    if (pheno=="IMID") {major.pheno<-"CTRL"; minor.pheno<-"RA"}
    if (pheno=="RespGM") {major.pheno<-"GOOD"; minor.pheno<-"NON_RESP"}
    if (pheno=="Resp") {major.pheno<-"GOOD"; minor.pheno<-"NON_RESP"}
    if (pheno=="ActB") {major.pheno<-"High"; minor.pheno<-"Low"}
    if (pheno=="ACPA") {major.pheno<-"Pos"; minor.pheno<-"Neg"}
    if (pheno=="RF") {major.pheno<-"Pos"; minor.pheno<-"Neg"}
    
    pheno.col<-which(colnames(grr$meta)==pheno)
    grr$meta$Pheno<-as.data.frame(grr$meta[,pheno.col])[,1]
    
    na.sample.exc<-as.data.frame(grr$meta[is.na(grr$meta$Pheno),"Sample"])[,1]
    meta_col<-dplyr::filter(grr$meta, Week!="RA");
    meta_col1<-meta_col[!(meta_col$Sample %in% unique(na.sample.exc)),]
    pheno.samples2keep<-as.data.frame(grr$meta[((grr$meta$Pheno==major.pheno | grr$meta$Pheno==minor.pheno) & !is.na(grr$meta$Pheno)),"Sample"])[,1]
    
    meta_col2<-meta_col1[meta_col1$Sample %in% unique(c(pheno.samples2keep)),]
    samples_col<-meta_col2$Sample;  data_col<-grr$data[match(samples_col,names(grr$data))];  grr_pheno<-list(data_col,meta_col2);  names(grr_pheno)<-c("data","meta")
    
    # Computing Kmers
    # Alternative k-mer based clustering methods involve enumerating all k-letter words in a sequence through a sliding window of length k.
    # The n×4k matrix of k-mer counts (where n is the number of sequences) can then be used in place of a multiple sequence
    kmers.pheno<-as.data.frame(getKmers(grr_pheno$data, 4))
    cols.interest<-sapply(kmers.pheno[,c(2:ncol(kmers.pheno))], as.numeric);  rownames(cols.interest)<-kmers.pheno$Kmer;
    cols.interest[is.na(cols.interest)]<-0
    kmer.count<-as.data.frame(cols.interest)
    
    #  Reads Info
    reads.data<-pubRep(grr_pheno$data, "aa", .verbose=F)
    reads.data<-as.data.frame(reads.data[,-c(1,2)])
    
    #  Kmer Info
    kmer.count$Samples<-apply(kmer.count, 1, function(x) length(which(x>0)))
    kmer.info<-data.frame("primerid"=rownames(kmer.count), "Samples"=kmer.count[,"Samples"])

    #  Kmer Reads
    kmer.count.info<-as.matrix(kmer.count[,-ncol(kmer.count)])
    kmer.count.info.cpm<-apply(kmer.count.info, 2, function(x) (x/sum(x, na.rm=T))*1000000)
    kmer.count.info.cpm.plus<-kmer.count.info.cpm+1
    kmer.count.info.cpm.plus[is.na(kmer.count.info.cpm.plus)]<-1
    kmer.count.info.cpm.log<-log2(kmer.count.info.cpm.plus)
    
    #  Sample Data
    df.reads.wk0<-readRDS("Data2analyze/ChainReadsCount.rds")
    chain.col<-which(colnames(df.reads.wk0)==chain)
    sample.info<-df.reads.wk0[(df.reads.wk0$Sample %in% grr_pheno$meta$Sample),c(1,chain.col)]
    colnames(sample.info)<-c("wellKey","Total.Reads")
    
    #  Creating SingleCellAssay object
    vbeta<-FromMatrix(kmer.count.info.cpm.log, cData=sample.info, fData=kmer.info)
    colData(vbeta)$CDR<-unname(scale(colSums(assay(vbeta)>0, na.rm=T))[,1])
    
    #  Creating vbeta by the phenotype to compare
    metadata<-readRDS(paste0("Data2analyze/Metadata_",chain,".txt"))
    pheno.values<-as.factor(grr_pheno$meta$Pheno)
    gender.values<-as.factor(as.character(metadata[match(colData(vbeta)$wellKey, metadata$Sample),"Gender"]))
    age.values<-as.numeric(metadata[match(colData(vbeta)$wellKey, metadata$Sample),"Age"])
    colData(vbeta)$Pheno<-pheno.values
    colData(vbeta)$Pheno[is.na(colData(vbeta)$Pheno)]<-NA
    colData(vbeta)$Age<-age.values
    colData(vbeta)$Gender<-gender.values
    
    #  Filtering out kmers detected in <1% of samples
    filt.value<-round(ncol(assay(vbeta))*0.01)
    n.kmers.pheno<-length(which(rowSums(assay(vbeta)>0)!=0))
    kmers2keep<-rowData(vbeta)[which(rowData(vbeta)$Samples>filt.value),"primerid"]
    vbeta<-vbeta[kmers2keep,]
    
    #  Statistical analysis
    print(paste0("\n\n",chain," - ",pheno," - Statistical analysis"))
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
    colnames(dfm)[1]<-"Kmer"
    
    saveRDS(dfm, paste0("Output_Data/Kmer_",pheno,"/",chain,"_statistics.rds"))
    
    #  Summary Table
    features1<-c("Chain", "Clinical Phenotype", "Reads (n)", "Kmers (n)", "Reads per sample (m)", "CDR per sample (m)", "Filtered Kmers (>1% samples, n)")
    features2<-c("Hurdle", "Continuous", "Discrete")
    total.chain.reads<-sum(reads.data[,colData(vbeta)$wellKey], na.rm=T)
    total.chain.kmers<-n.kmers.pheno
    reads.sample<-round(sum(reads.data[,colData(vbeta)$wellKey], na.rm=T)/nrow(colData(vbeta)),2)
    cdr.sample<-round(mean(colSums(assay(vbeta)>0)),2)
    chain.kmers.filt<-length(names(pvals))
    h.sigP<-length(which(dfm$Pval.hurdle<0.05))
    h.sigFDR<-length(which(dfm$FDR.hurdle<0.05))
    d.sigP<-length(which(dfm$Pval.disc<0.05))
    d.sigFDR<-length(which(dfm$FDR.disc<0.05))
    c.sigP<-length(which(dfm$Pval.cont<0.05))
    c.sigFDR<-length(which(dfm$FDR.cont<0.05))
    h.sig<-paste0(h.sigP," (",h.sigFDR,")")
    c.sig<-paste0(c.sigP," (",c.sigFDR,")")
    d.sig<-paste0(d.sigP," (",d.sigFDR,")")
  
    values1<-c(chain, paste0(pheno,": ",major.pheno," vs ",minor.pheno), total.chain.reads, total.chain.kmers, reads.sample, cdr.sample, chain.kmers.filt)
    values2<-c(h.sig, c.sig, d.sig)
    summ.table.1<-data.frame("Feature iRep Analysis"=features1, "Value"=values1)
    summ.table.2<-data.frame("Significant clones"=features2, "P<0.05 (FDR<0.05)"=values2)
    colnames(summ.table.1)<-c("Feature iRep Analysis", "Value")
    colnames(summ.table.2)<-c("Significant kmers", "P<0.05 (FDR<0.05)")
    st1<-tableGrob(summ.table.1, rows=NULL)
    st2<-tableGrob(summ.table.2, rows=NULL)
    st<-gtable_combine(st1, st2, along=2)
    
    jpeg(paste0("Output_Data/Kmer_",pheno,"/",chain,"_summary_table.jpeg"), res=300, width=1700, height=1175)
    grid.arrange(st)
    dev.off()
  }
}


# Kmer characterization
# This code should be modified according to the kmers of interest (e.g. kmers associated with RA)

chains.sig<-c("TRA","TRB","IGK","IGL")
chains<-chains.sig
models<-c("Hurdle","Cont","Disc")

df.models<-vector("list",length(models)); names(df.models)<-models
df.sig<-vector("list",length(chains)); names(df.sig)<-chains

for (ch in names(df.sig)) {
  df.sig[[ch]]<-df.models
  df.results<-readRDS(paste0("Output_Data/Kmer_IMID_",ch,"_statistics.rds"))
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


pdf("Output_Data/KmerCharacterization.pdf", width=10)
for (ch in names(df.sig)) {
  print(ch)
  if (ch=="TRA") {grr<-grr.tra}
  if (ch=="TRB") {grr<-grr.trb}
  if (ch=="TRD") {grr<-grr.trd}
  if (ch=="TRG") {grr<-grr.trg}
  if (ch=="IGH") {grr<-grr.igh}
  if (ch=="IGL") {grr<-grr.igl}
  if (ch=="IGK") {grr<-grr.igk}
  
  pheno<-"IMID";
  major.pheno<-"CTRL";
  minor.pheno<-"RA";
  
  pheno.col<-which(colnames(grr$meta)==pheno)
  grr$meta$Pheno<-as.data.frame(grr$meta[,pheno.col])[,1]
  
  na.sample.exc<-as.data.frame(grr$meta[is.na(grr$meta$Pheno),"Sample"])[,1]
  meta_col<-dplyr::filter(grr$meta, Week!="RA");
  meta_col1<-meta_col[!(meta_col$Sample %in% unique(na.sample.exc)),]
  pheno.samples2keep<-as.data.frame(grr$meta[((grr$meta$Pheno==major.pheno | grr$meta$Pheno==minor.pheno) & !is.na(grr$meta$Pheno)),"Sample"])[,1]
  
  meta_col2<-meta_col1[meta_col1$Sample %in% unique(c(pheno.samples2keep)),]
  samples_col<-meta_col2$Sample;  data_col<-grr$data[match(samples_col,names(grr$data))];  grr_pheno<-list(data_col,meta_col2);  names(grr_pheno)<-c("data","meta")
  
  # Computing Kmers
  kmers.pheno<-as.data.frame(getKmers(grr_pheno$data, 4))
  cols.interest<-sapply(kmers.pheno[,c(2:ncol(kmers.pheno))], as.numeric);  rownames(cols.interest)<-kmers.pheno$Kmer;
  cols.interest[is.na(cols.interest)]<-0
  kmer.count<-as.data.frame(cols.interest)
  
  #  Reads Info
  reads.data<-pubRep(grr_pheno$data, "aa", .verbose=F)
  reads.data<-as.data.frame(reads.data[,-c(1,2)])
  
  #  Kmer Info
  kmer.count$Samples<-apply(kmer.count, 1, function(x) length(which(x>0)))
  kmer.info<-data.frame("primerid"=rownames(kmer.count), "Samples"=kmer.count[,"Samples"])
  
  #  Kmer Reads
  kmer.count.info<-as.matrix(kmer.count[,-ncol(kmer.count)])
  kmer.count.info.cpm<-apply(kmer.count.info, 2, function(x) (x/sum(x, na.rm=T))*1000000)
  kmer.count.info.cpm.plus<-kmer.count.info.cpm+1
  kmer.count.info.cpm.plus[is.na(kmer.count.info.cpm.plus)]<-1
  kmer.count.info.cpm.log<-log2(kmer.count.info.cpm.plus)
  
  #  Sample Data
  df.reads.wk0<-readRDS("Data2analyze/ChainReadsCount.rds")
  chain.col<-which(colnames(df.reads.wk0)==ch)
  sample.info<-df.reads.wk0[(df.reads.wk0$Sample %in% grr_pheno$meta$Sample),c(1,chain.col)]
  colnames(sample.info)<-c("wellKey","Total.Reads")
  
  #  Creating SingleCellAssay object
  vbeta<-FromMatrix(kmer.count.info.cpm.log, cData=sample.info, fData=kmer.info)
  colData(vbeta)$CDR<-unname(scale(colSums(assay(vbeta)>0, na.rm=T))[,1])
  
  #  Creating vbeta by the phenotype to compare
  metadata<-readRDS("Data2analyze/Metadata_",ch,".rds")
  pheno.values<-as.factor(grr_pheno$meta$Pheno)
  gender.values<-as.factor(as.character(metadata[match(colData(vbeta)$wellKey, metadata$Sample),"Gender"]))
  age.values<-as.numeric(metadata[match(colData(vbeta)$wellKey, metadata$Sample),"Age"])
  colData(vbeta)$Pheno<-pheno.values
  colData(vbeta)$Pheno[is.na(colData(vbeta)$Pheno)]<-NA
  colData(vbeta)$Age<-age.values
  colData(vbeta)$Gender<-gender.values
  
  #  Filtering out kmers detected in <1% of samples
  filt.value<-round(ncol(assay(vbeta))*0.01)
  n.kmers.pheno<-length(which(rowSums(assay(vbeta)>0)!=0))
  kmers2keep<-rowData(vbeta)[which(rowData(vbeta)$Samples>filt.value),"primerid"]
  vbeta<-vbeta[kmers2keep,]
  
  for (mod in models) {
    print(mod)
    data<-df.sig[[ch]][[mod]]
    colnames(data)<-c("Kmer","Pval.Hurdle","FDR.Hurdle","Pval.Disc","FDR.Disc","Pval.Cont","FDR.Cont")
    # Summary Plot
    if (nrow(data)>1) {
      kmers<-rep(data$Kmer,3*2)
      mods<-c(rep("Hurdle",2*nrow(data)), rep("Discrete",2*nrow(data)), rep("Continuous",2*nrow(data)))
      stats<-c(rep("Pval",nrow(data)), rep("FDR",nrow(data)), rep("Pval",nrow(data)), rep("FDR",nrow(data)), rep("Pval",nrow(data)), rep("FDR",nrow(data)))
      values<-c(data[,"Pval.Hurdle"], data[,"FDR.Hurdle"], data[,"Pval.Disc"], data[,"FDR.Disc"],  data[,"Pval.Cont"],data[,"FDR.Cont"])
      df2plot<-data.frame("Kmer"=kmers, "Model"=mods, "Significance"=stats, "Value"=values)
      
      for (i in 1:nrow(data)) {
        kmer<-data[i,"Kmer"];    pheno<-"IMID";    pval<-data[i,paste0("Pval.",mod)];
        
        col.pheno<-which(colnames(metadata)==pheno)
        meta2use<-metadata[,c(1,col.pheno)]
        
        pheno.values<-meta2use[match(colData(vbeta)$wellKey, meta2use$Sample),]
        colData(vbeta)$Pheno<-pheno.values[,pheno]
        
        vbeta.subset<-vbeta[as.character(kmer),]
        
        vbeta.pheno<-subset(vbeta.subset, (Pheno=="CTRL" | Pheno=="RA"))
        flat.dat.pheno<-as(vbeta.pheno[as.character(kmer),], 'data.table')
        
        pval.pheno.text<-paste0("P=",formatC(as.numeric(pval), format="e", digits=2))
        
        if (mod=="Disc") {
          flat.dat.pheno.count<-flat.dat.pheno
          flat.dat.pheno.count$Detection<-as.character(flat.dat.pheno.count$Pheno)
          flat.dat.pheno.count$Detection[flat.dat.pheno.count$value>0]<-"Yes"
          flat.dat.pheno.count$Detection[flat.dat.pheno.count$value==0]<-"No"
          flat.dat.pheno.count$Pheno<-droplevels(flat.dat.pheno.count$Pheno)
          
          tab.ph<-as.data.frame(table(flat.dat.pheno.count$Pheno, flat.dat.pheno.count$Detection))
          
          exp.ph<-ggplot(data=tab.ph, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat="identity", position=position_dodge()) +
            ylab("\nIndividuals with Kmer detected (N)\n") + xlab("") + theme_classic() + labs (fill="Detection") +
            ggtitle("\nKmer Detection\n") + theme(plot.title=element_text(hjust=0.5, face="bold")) +
            annotate(geom="text", x=1.5, y=max(tab.ph$Freq)*1.1, size=5, label=pval.pheno.text, color="black")
          
          kdr.ph<-ggplot(flat.dat.pheno, aes(x=CDR, y=value, color=Pheno)) + theme_classic() +
            ggtitle("\nAbundance by KDR\n") + theme(plot.title=element_text(hjust=0.5, face="bold"), legend.title=element_blank()) +
            geom_jitter() + xlab('\nStandardized Kmer Detection Rate (KDR)\n') +
            ylab('\nLog2(CPM)\n') + scale_color_manual(name=pheno, values=c("#D8EBD8","#F69191"))

          grid.arrange(exp.ph, kdr.ph, ncol=2, top=textGrob(paste0("\n",kmer," from ",ch," chain significant in ",mod," model\n"),gp=gpar(fontsize=18,fontface="bold")))
          
        } else {
          
          exp.ph<-ggplot(flat.dat.pheno, aes(x=Pheno, y=value, fill=Pheno)) + 
            theme_classic() + geom_jitter(shape=20, position=position_jitter(0)) + 
            theme(plot.title=element_text(hjust=0.5, face="bold"), legend.title=element_blank()) + ylab("\nLog2(CPM)\n") + xlab("") +
            ggtitle("\nKmer Expression\n") + geom_violin(trim=F) +
            stat_summary(fun.data=mean_sdl, geom="pointrange", color="red", size=0.1) +
            scale_fill_manual(values=c("#D8EBD8","#F69191")) +
            annotate(geom="text", x=1.5, y=max(flat.dat.pheno$value)*1.1, size=5, label=pval.pheno.text, color="black")
          
          kdr.ph<-ggplot(flat.dat.pheno, aes(x=CDR, y=value, color=Pheno)) + theme_classic() +
            ggtitle("\nAbundance by KDR\n") + theme(plot.title=element_text(hjust=0.5, face="bold"), legend.title=element_blank()) +
            geom_jitter() + xlab('\nStandardized Kmer Detection Rate (KDR)\n') + 
            ylab('\nLog2(CPM)\n') + scale_color_manual(name=pheno, values=c("#D8EBD8","#F69191"))
          
          grid.arrange(exp.ph, kdr.ph, ncol=2, top=textGrob(paste0("\n",kmer," from ",ch," chain significant in ",mod," model\n"),gp=gpar(fontsize=18,fontface="bold")))
        }
      }
    }
  }
}
dev.off()
  
# Getting LogFC

for (chain in names(df.sig)) {
  
  print(paste0("Loading Data ", chain))
  
  grr<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',chain,'.rds'))
  grr$meta$RespGM<-grr$meta$Resp
  grr$meta$RespGM[grr$meta$RespGM=="MOD"]<-"GOOD"
  grr$meta$Resp[grr$meta$Resp=="MOD"]<-NA
  
  pheno="IMID";
  major.pheno<-"CTRL";
  minor.pheno<-"RA";
  
  pheno.col<-which(colnames(grr$meta)==pheno)
  grr$meta$Pheno<-as.data.frame(grr$meta[,pheno.col])[,1]
  
  na.sample.exc<-as.data.frame(grr$meta[is.na(grr$meta$Pheno),"Sample"])[,1]
  meta_col<-dplyr::filter(grr$meta, Week!="RA");
  meta_col1<-meta_col[!(meta_col$Sample %in% unique(na.sample.exc)),]
  pheno.samples2keep<-as.data.frame(grr$meta[((grr$meta$Pheno==major.pheno | grr$meta$Pheno==minor.pheno) & !is.na(grr$meta$Pheno)),"Sample"])[,1]
  
  meta_col2<-meta_col1[meta_col1$Sample %in% unique(c(pheno.samples2keep)),]
  samples_col<-meta_col2$Sample;  data_col<-grr$data[match(samples_col,names(grr$data))];  grr_pheno<-list(data_col,meta_col2);  names(grr_pheno)<-c("data","meta")
  
  # Computing Kmers
  kmers.pheno<-as.data.frame(getKmers(grr_pheno$data, 4))
  cols.interest<-sapply(kmers.pheno[,c(2:ncol(kmers.pheno))], as.numeric);  rownames(cols.interest)<-kmers.pheno$Kmer;
  cols.interest[is.na(cols.interest)]<-0
  kmer.count<-as.data.frame(cols.interest)
  
  #  Reads Info
  reads.data<-pubRep(grr_pheno$data, "aa", .verbose=F)
  reads.data<-as.data.frame(reads.data[,-c(1,2)])
  
  #  Kmer Info
  kmer.count$Samples<-apply(kmer.count, 1, function(x) length(which(x>0)))
  kmer.info<-data.frame("primerid"=rownames(kmer.count), "Samples"=kmer.count[,"Samples"])
  
  #  Kmer Reads
  kmer.count.info<-as.matrix(kmer.count[,-ncol(kmer.count)])
  kmer.count.info.cpm<-apply(kmer.count.info, 2, function(x) (x/sum(x, na.rm=T))*1000000)
  kmer.count.info.cpm.plus<-kmer.count.info.cpm+1
  kmer.count.info.cpm.plus[is.na(kmer.count.info.cpm.plus)]<-1
  kmer.count.info.cpm.log<-log2(kmer.count.info.cpm.plus)
  
  #  Sample Data
  df.reads.wk0<-readRDS("Data2analyze/ChainReadsCount.rds")
  chain.col<-which(colnames(df.reads.wk0)==chain)
  sample.info<-df.reads.wk0[(df.reads.wk0$Sample %in% grr_pheno$meta$Sample),c(1,chain.col)]
  colnames(sample.info)<-c("wellKey","Total.Reads")
  
  #  Creating SingleCellAssay object
  vbeta<-FromMatrix(kmer.count.info.cpm.log, cData=sample.info, fData=kmer.info)
  colData(vbeta)$CDR<-unname(scale(colSums(assay(vbeta)>0, na.rm=T))[,1])
  
  #  Creating vbeta by the phenotype to compare
  metadata<-readRDS("Data2analyze/Metadata_",chain,".rds")
  pheno.values<-as.factor(grr_pheno$meta$Pheno)
  gender.values<-as.factor(as.character(metadata[match(colData(vbeta)$wellKey, metadata$Sample),"Gender"]))
  age.values<-as.numeric(metadata[match(colData(vbeta)$wellKey, metadata$Sample),"Age"])
  colData(vbeta)$Pheno<-pheno.values
  colData(vbeta)$Pheno[is.na(colData(vbeta)$Pheno)]<-NA
  colData(vbeta)$Age<-age.values
  colData(vbeta)$Gender<-gender.values
  
  #  Filtering out kmers detected in <1% of samples
  filt.value<-round(ncol(assay(vbeta))*0.01)
  n.kmers.pheno<-length(which(rowSums(assay(vbeta)>0)!=0))
  kmers2keep<-rowData(vbeta)[which(rowData(vbeta)$Samples>filt.value),"primerid"]
  vbeta<-vbeta[kmers2keep,]
  
  #  Statistical analysis
  print(paste0("\n\n",chain," - ",pheno," - Statistical analysis"))
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
  colnames(dfm)[1]<-"Kmer"
  
  df1<-as.data.frame(summary(zlm.output)$datatable)
  df2<-df1[(df1$component=="logFC" & df1$contrast==comp),]
  
  sig.clones.h<-as.character(dfm[dfm$FDR.hurdle<0.05,"Kmer"])
  if (length(sig.clones.h)>0) {
    df.clones<-df2[df2$primerid %in% sig.clones.h, c(1,6,5,4)]
    df.clones<-df.clones[order(df.clones$coef),]
    colnames(df.clones)<-c("Kmer","LogFC","CI.Low","CI.High")
    saveRDS(df.clones, paste0("Output_Data/Kmer_LogFC_",chain,"_Hurdle.rds"))
  }
  
  sig.clones.c<-as.character(dfm[dfm$FDR.cont<0.05,"Kmer"])
  if (length(sig.clones.c)>0) {
    df.clones<-df2[df2$primerid %in% sig.clones.c, c(1,6,5,4)]
    df.clones<-df.clones[order(df.clones$coef),]
    colnames(df.clones)<-c("Kmer","LogFC","CI.Low","CI.High")
    saveRDS(df.clones, paste0("Output_Data/Kmer_LogFC_",chain,"_Cont.rds"))
  }
  
  sig.clones.d<-as.character(dfm[dfm$FDR.disc<0.05,"Kmer"])
  if (length(sig.clones.d)>0) {
    df.clones<-df2[df2$primerid %in% sig.clones.d, c(1,6,5,4)]
    df.clones<-df.clones[order(df.clones$coef),]
    colnames(df.clones)<-c("Kmer","LogFC","CI.Low","CI.High")
    saveRDS(df.clones, paste0("Output_Data/Kmer_LogFC_",chain,"_Disc.rds"))
  }
}

# Similarity - Clustering

get.perm.pval<-function(kmers.query, chain.kmers) {
  obs.dist.query<-mean(stringdistmatrix(kmers.query, method="lv"))
  set.dist.query.n<-length(kmers.query)
  if (set.dist.query.n>0) {
    count.low=0;  count.hi=0;  nperm<-10000;
    perm.dist.query<-list()
    for (i in 1:nperm) {
      x<-sample(chain.kmers, size=set.dist.query.n, replace=F)
      perm.d<-mean(stringdistmatrix(x, method="lv"))
      perm.dist.query[[i]]<-perm.d
    }
    perm.d<-unlist(perm.dist.query)
    df.out<-data.frame("Similarity.Metric"="Levenshtein Distance", "Obs.Value"=obs.dist.query, "Perm.Value"=mean(perm.d),
                       "PermValue.LT.ObsValue"=length(which(perm.d<obs.dist.query)),"PermValue.GT.ObsValue"=length(which(perm.d>=obs.dist.query)),
                       "Pval"=(length(which(perm.d<obs.dist.query))+1)/(nperm+1), stringsAsFactors=F)
  }
}
files<-grep("Kmer_LogFC", list.files("Output_Data/"), value=T)
for (f in files){

  print(f)

  chain<-strsplit(f,'_')[[1]][3]
  mod<-strsplit(f,'_')[[1]][4]; mod<-gsub('.rds','',mod)
  f1<-readRDS(paste0("Output_Data/",f))
  f1$RA<-f1$Kmer
  f1$RA[f1$LogFC<0]<-"Down"
  f1$RA[f1$LogFC>0]<-"Up"
  f1$RA[is.na(f1$LogFC)]<-"RA-specific" # or "Up"
  kmers.sig.up<-f1[f1$RA=="Up","Kmer"]
  kmers.sig.down<-f1[f1$RA=="Down","Kmer"]
  kmers.sig.all<-f1$Kmer
  chain.kmers<-as.character(readRDS(paste0("Output_Data/Kmer_IMID_",chain,"_statistics.rds"))$Kmer)
  
  # Permutations
  res.sig.up<-get.perm.pval(kmers.sig.up,chain.kmers)
  res.sig.down<-get.perm.pval(kmers.sig.down,chain.kmers)
  res.sig.all<-get.perm.pval(kmers.sig.all,chain.kmers)
  
  if(is.null(res.sig.up)) {res.sig.up<-data.frame("Similarity.Metric"="Levenshtein Distance", "Obs.Value"=NA, "Perm.Value"=NA, "PermValue.LT.ObsValue"=NA,"PermValue.GT.ObsValue"=NA, "Pval"=NA)}
  if(is.null(res.sig.down)) {res.sig.down<-data.frame("Similarity.Metric"="Levenshtein Distance", "Obs.Value"=NA, "Perm.Value"=NA, "PermValue.LT.ObsValue"=NA,"PermValue.GT.ObsValue"=NA, "Pval"=NA)}
  if(is.null(res.sig.all)) {res.sig.all<-data.frame("Similarity.Metric"="Levenshtein Distance", "Obs.Value"=NA, "Perm.Value"=NA, "PermValue.LT.ObsValue"=NA,"PermValue.GT.ObsValue"=NA, "Pval"=NA)}
  
  df.out<-rbind(res.sig.up, res.sig.down, res.sig.all)
  df.out$Kmer.Type<-c("Overrepresentation in RA","Underrepresentation in RA","All Kmers")
  saveRDS(df.out[,c(7,1,2,3,4,5,6)],paste0("Output_Data/Kmer_Similarity_",chain,"_",mod,".rds"))
  
  # Clustering
  m<-as.matrix(stringdistmatrix(f1$Kmer, useNames =T, method="lv"))
  table(rownames(m)==f1$Kmer)

  ha<-HeatmapAnnotation(RA=f1$RA, col=list(RA=c("Down"="lightsteelblue1","Up"="royalblue1","RA-specific"="royalblue4")))
  heatmap<-Heatmap(m, name="Dissimilarity\nIndex", top_annotation=ha, cluster_rows=T,
                   show_column_names=F, cluster_columns=T, row_names_gp=gpar(fontsize=8),
                   col=circlize::colorRamp2(c(0,1,2,3,4), rev(brewer.pal(n=5, name="OrRd"))))
  
  jpeg(paste0("Output_Data/Kmer_Similarity_Heatmap_",chain,"_",mod,".jpeg"), res=300, width=3000, height=3000)
  print(heatmap)
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
get.perm.features<-function(set.n, chain.kmers, chain, get.aa.features) {
  count.low<-0;  count.hi<-0;  nperm<-1000;
  df.perm.all<-data.frame("Chain"="TEST", "Kmer.Type"="TEST", "Perm"=0,
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
    perm.kmers<-sample(chain.kmers, size=set.n, replace=F)
    df.perm.kmers<-data.frame("Chain"="TEST", "Kmer"="TEST", "Kmer.Type"="TEST",
                                    "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                                    "Aliphatic"=0, "Aromatic"=0,
                                    "Non.Polar"=0, "Polar"=0, "Charged"=0,
                                    "Basic"=0, "Acidic"=0,
                                    "Net.Charge"=0, "Isoelectric.Point"=0,
                                    "Aliphatic.Index"=0, "Instability.Index"=0,
                                    "PPI.Index"=0, "Hydrophobicity.Index"=0,
                                    "HydrophobicMoment.Index"=0, stringsAsFactors=F)
    count.perm=0
    for (km in perm.kmers) {
        count.perm=count.perm+1
        cl.perm.prop<-get.aa.features(km)
        df.perm.kmers[count.perm,]<-c(chain,km,"Perm.Kmer",cl.perm.prop)
    }
    cols.interest.perm<-sapply(df.perm.kmers[,c(4:ncol(df.perm.kmers))], as.numeric)
    perm.prop<-colMeans(cols.interest.perm)
    df.perm.all[i,]<-c(chain,"Perm.Kmer",i,perm.prop)
  }
  perm.prop<-df.perm.all
  return(perm.prop)
}
get.final.out<-function(obs.sig,perm.sig,nperm){
  aa.propierties<-names(obs.sig)
  df.out<-data.frame("Biochemical.Feature"="TEST","Obs.Value"="TEST","Perm.Value"="TEST","PermValue.LT.ObsValue"="TEST","PermValue.GT.ObsValue"="TEST","Pval.PermValue.LT.Obs.Value"="TEST","Pval.PermValue.GT.ObsValue"="TEST", stringsAsFactors=F)
  count=0
  for (aap in aa.propierties) {
    obs.value<-obs.sig[aap]
    perm.value<-as.numeric(perm.sig[,aap])
    pval1<-(length(which(perm.value<obs.value))+1)/(nperm+1)
    pval2<-(length(which(perm.value>=obs.value))+1)/(nperm+1)
    count=count+1
    df.out[count,]<-c(aap, unname(obs.value), mean(perm.value), length(which(perm.value<obs.value)), length(which(perm.value>obs.value)), pval1, pval2)
  }
  return(df.out)
}

files<-grep("Kmer_LogFC", list.files("Output_Data/"), value=T)

for (f in files){
  
  print(f)
  
  chain<-strsplit(f,'_')[[1]][2]
  mod<-strsplit(f,'_')[[1]][3]; mod<-gsub('.rds','',mod)
  f1<-readRDS(paste0("Output_Data/",f))
  f1$RA<-f1$Kmer
  f1$RA[f1$LogFC<0]<-"Down"
  f1$RA[f1$LogFC>0]<-"Up"
  kmers.sig.up<-f1[f1$RA=="Up","Kmer"]
  kmers.sig.down<-f1[f1$RA=="Down","Kmer"]
  kmers.sig.all<-f1$Kmer
  chain.kmers<-as.character(readRDS(paste0("Output_Data/Kmer_IMID_",chain,"_statistics.rds"))$Kmer)
  
  df.sig.down<-data.frame("Chain"="TEST", "Kmer"="TEST", "Kmer.Type"="TEST",
                                 "Molecular.Weight"=0, "Tiny"=0, "Small"=0,
                                 "Aliphatic"=0, "Aromatic"=0,
                                 "Non.Polar"=0, "Polar"=0, "Charged"=0,
                                 "Basic"=0, "Acidic"=0,
                                 "Net.Charge"=0, "Isoelectric.Point"=0,
                                 "Aliphatic.Index"=0, "Instability.Index"=0,
                                 "PPI.Index"=0, "Hydrophobicity.Index"=0,
                                 "HydrophobicMoment.Index"=0, stringsAsFactors=F)
  
  df.sig.up<-df.sig.down
  df.sig.all<-df.sig.down
  
  if (length(kmers.sig.down)>0) {
    count=0
    for (km in kmers.sig.down) {
      count=count+1
      cl.prop<-get.aa.features(km)
      df.sig.down[count,]<-c(chain,km,"Sig.Kmer",cl.prop)
    }
    cols.interest<-sapply(df.sig.down[,c(4:ncol(df.sig.down))], as.numeric)
    obs.sig.down<-colMeans(cols.interest)
  } else {df.sig.down[1,]<-rep("NA",20)}
  
  if (length(kmers.sig.up)>0) {
    count=0
    for (km in kmers.sig.up) {
      count=count+1
      cl.prop<-get.aa.features(km)
      df.sig.up[count,]<-c(chain,km,"Sig.Kmer",cl.prop)
    }
    cols.interest<-sapply(df.sig.up[,c(4:ncol(df.sig.up))], as.numeric)
    obs.sig.up<-colMeans(cols.interest)  
  } else {df.sig.up[1,]<-rep("NA",20)}
  
  count=0
  for (km in kmers.sig.all) {
    count=count+1
    cl.prop<-get.aa.features(km)
    df.sig.all[count,]<-c(chain,km,"Sig.Kmer",cl.prop)
  }
  cols.interest<-sapply(df.sig.all[,c(4:ncol(df.sig.all))], as.numeric)
  obs.sig.all<-colMeans(cols.interest)
  
  if (length(kmers.sig.up)>0) {perm.sig.up<-get.perm.features(length(kmers.sig.up), chain.kmers, chain, get.aa.features)} else {perm.sig.up<-data.frame("Chain"=NA, "Kmer.Type"=NA, "Perm"=NA,
                                                                                                                                                        "Molecular.Weight"=NA, "Tiny"=NA, "Small"=NA,
                                                                                                                                                        "Aliphatic"=NA, "Aromatic"=NA,
                                                                                                                                                        "Non.Polar"=NA, "Polar"=NA, "Charged"=NA,
                                                                                                                                                        "Basic"=NA, "Acidic"=NA,
                                                                                                                                                        "Net.Charge"=NA, "Isoelectric.Point"=NA,
                                                                                                                                                        "Aliphatic.Index"=NA, "Instability.Index"=NA,
                                                                                                                                                        "PPI.Index"=NA, "Hydrophobicity.Index"=NA,
                                                                                                                                                        "HydrophobicMoment.Index"=NA, stringsAsFactors=F)}
  if (length(kmers.sig.down)>0) {perm.sig.down<-get.perm.features(length(kmers.sig.down), chain.kmers, chain, get.aa.features)} else {perm.sig.down<-data.frame("Chain"=NA, "Kmer.Type"=NA, "Perm"=NA,
                                                                                                                                                              "Molecular.Weight"=NA, "Tiny"=NA, "Small"=NA,
                                                                                                                                                              "Aliphatic"=NA, "Aromatic"=NA,
                                                                                                                                                              "Non.Polar"=NA, "Polar"=NA, "Charged"=NA,
                                                                                                                                                              "Basic"=NA, "Acidic"=NA,
                                                                                                                                                              "Net.Charge"=NA, "Isoelectric.Point"=NA,
                                                                                                                                                              "Aliphatic.Index"=NA, "Instability.Index"=NA,
                                                                                                                                                              "PPI.Index"=NA, "Hydrophobicity.Index"=NA,
                                                                                                                                                              "HydrophobicMoment.Index"=NA, stringsAsFactors=F)}
  if (length(kmers.sig.all)>0) {perm.sig.all<-get.perm.features(length(kmers.sig.all), chain.kmers, chain, get.aa.features)} else {perm.sig.all<-data.frame("Chain"=NA, "Kmer.Type"=NA, "Perm"=NA,
                                                                                                                                                           "Molecular.Weight"=NA, "Tiny"=NA, "Small"=NA,
                                                                                                                                                           "Aliphatic"=NA, "Aromatic"=NA,
                                                                                                                                                           "Non.Polar"=NA, "Polar"=NA, "Charged"=NA,
                                                                                                                                                           "Basic"=NA, "Acidic"=NA,
                                                                                                                                                           "Net.Charge"=NA, "Isoelectric.Point"=NA,
                                                                                                                                                           "Aliphatic.Index"=NA, "Instability.Index"=NA,
                                                                                                                                                           "PPI.Index"=NA, "Hydrophobicity.Index"=NA,
                                                                                                                                                           "HydrophobicMoment.Index"=NA, stringsAsFactors=F)}
  
  if (length(kmers.sig.up)>0) {final.out.sig.up<-get.final.out(obs.sig.up, perm.sig.up, nperm)} else {
    final.out.sig.up<-data.frame("Biochemical.Feature"=NA,"Obs.Value"=NA,"Perm.Value"=NA,"PermValue.LT.ObsValue"=NA,"PermValue.GT.ObsValue"=NA,"Pval.PermValue.LT.Obs.Value"=NA,"Pval.PermValue.GT.ObsValue"=NA, stringsAsFactors=F)
  }
  if (length(kmers.sig.down)>0) {final.out.sig.down<-get.final.out(obs.sig.down, perm.sig.down, nperm)} else {
    final.out.sig.down<-data.frame("Biochemical.Feature"=NA,"Obs.Value"=NA,"Perm.Value"=NA,"PermValue.LT.ObsValue"=NA,"PermValue.GT.ObsValue"=NA,"Pval.PermValue.LT.Obs.Value"=NA,"Pval.PermValue.GT.ObsValue"=NA, stringsAsFactors=F)
  }
  if (length(kmers.sig.all)>0) {final.out.sig.all<-get.final.out(obs.sig.all, perm.sig.all, nperm)} else {
    final.out.sig.all<-data.frame("Biochemical.Feature"=NA,"Obs.Value"=NA,"Perm.Value"=NA,"PermValue.LT.ObsValue"=NA,"PermValue.GT.ObsValue"=NA,"Pval.PermValue.LT.Obs.Value"=NA,"Pval.PermValue.GT.ObsValue"=NA, stringsAsFactors=F)
  }
  saveRDS(final.out.sig.up,paste0("Output_Data/Kmer_AminoAcid_Prop_",chain,"_",mod,"_Up.rds"))
  saveRDS(final.out.sig.down,paste0("Output_Data/Kmer_AminoAcid_Prop_",chain,"_",mod,"_Down.rds"))
  saveRDS(final.out.sig.all,paste0("Output_Data/Kmer_AminoAcid_Prop_",chain,"_",mod,"_All.rds"))
}

chains<-c("TRA","TRB","IGL","IGK") # Chains conatining significant k-mers
df.all<-list()
count=0
for (c in chains) {
    count=count+1
    if (c!="TRB") {mod<-"Hurdle"} else {mod<-"Cont"}
    file.name<-paste0("Output_Data/Kmer_AminoAcid_Prop_",chain,"_",mod,"_All.rds")
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

# only plot if BEP (Biochemical Enrichment Score is > -log10(0.05))

data2plot$Biochemical.Enrichment.Score[(data2plot$Biochemical.Enrichment.Score<1.3 & data2plot$Biochemical.Enrichment.Score>-1.3)]<-NA
colnames(data2plot)[3]<-"BES"

p.kmer.biochemical<-ggplot(data=data2plot, aes(x=factor(Chain), y=factor(Biochemical.Feature), size=abs(BES), fill=BES)) +
                  geom_point(alpha=1, shape=21, color="black") +
                  scale_size(range=c(.1, 8), name="Significance") +
                  labs(y="", x="") +
                  scale_fill_gradientn(colours = rev(colorspace::diverge_hcl(15))) +
                  theme_bw() + theme(plot.title=element_text(hjust=0.5, face="bold", size=16), axis.text.x=element_text(vjust=0.5, size=13), strip.text=element_text(size=14, face="bold"))

jpeg("Output_Data/BiochemicalProgile_RAassoc_kmers.jpeg", res=300, height=2000, width=1800)
grid.arrange(p.kmer.biochemical)
dev.off()
