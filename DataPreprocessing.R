#################################  0. Loading libraries  #################################

library(data.table)
library(tidyverse)
library(immunarch)
library(purrr)
library(dplyr)

#################################  1. Converting mixcr to immunarch format  #################################

dir.create("Data2analyze/Immunarch_Format/")
chains<-c("IGK","IGH","IGL","TRA","TRB","TRG","TRD")
count=0
for (i in 1:length(chains)){
  c<-chains[i]
  dir.path<-paste0("Input_Data/",c,"/")
  inds<-list.files(dir.path); inds<-unique(sort(inds[grep("tsv",inds,invert=F)]))
  dir.create(paste0("Data2analyze/Immunarch_Format/",c))
  for (j in 1:length(inds)){
    print(paste0(c,"-",j))
    sample<-inds[j]
    sample2<-gsub('.tsv','',sample)
    file.name<-paste0(dir.path,sample)
    mixcr<-immunarch::repLoad(file.name, .format='mixcr')
    command<-paste0("mixcr$data$",sample2)
    data<-eval(parse(text=command))
    df<-as.data.frame(data)
    outname<-paste0("Data2analyze/Immunarch_Format/",c,"/",gsub(paste0('__',c),'',sample2),'.csv')
    write.table(df,outname,row.names=F,col.names=T,quote=F,sep="\t")
  }
}

#################################  2. Keeping best VDJ  #################################

dir.create(paste0("Data2analyze/Immunarch_Format_BestVDJ/"))

chains<-c("TRG","TRD","TRB","TRA","IGK","IGH","IGL")
count=0
for (c in chains){
  dir.create(paste0("Data2analyze/Immunarch_Format_BestVDJ/",c))
  inds<-gsub(".csv","", grep(".csv", list.files(paste0("Data2analyze/Immunarch_Format/",c,"/")), value=T))
  for (i in inds){
    print(paste0(c," ",i))
    count=count+1
    
    f1<-fread(paste0("Data2analyze/Immunarch_Format/",c,"/",i,".csv"));
      
    f1Vname.splt<-strsplit(as.character(f1$V.name),", ");
    f1Vname.splt[lapply(f1Vname.splt,length)<1]<-NA;
    for (k in 1:length(f1Vname.splt)){
      data2query<-f1Vname.splt[[k]]
      value<-data2query[grep(c,data2query)]
      if (length(value)>0) {data2save<-value} else {data2save<-NA}
      f1Vname.splt[[k]]<-data2save
    }
    f1Vname.best<-unlist(purrr::map(f1Vname.splt, 1))
    
    f1Dname.splt<-strsplit(as.character(f1$D.name),", ");
    f1Dname.splt[lapply(f1Dname.splt,length)<1]<-NA;
    for (k in 1:length(f1Dname.splt)){
      data2query<-f1Dname.splt[[k]]
      value<-data2query[grep(c,data2query)]
      if (length(value)>0) {data2save<-value} else {data2save<-NA}
      f1Dname.splt[[k]]<-data2save
    }
    f1Dname.best<-unlist(purrr::map(f1Dname.splt, 1))
    
    f1Jname.splt<-strsplit(as.character(f1$J.name),", ");
    f1Jname.splt[lapply(f1Jname.splt,length)<1]<-NA;
    for (k in 1:length(f1Jname.splt)){
      data2query<-f1Jname.splt[[k]]
      value<-data2query[grep(c,data2query)]
      if (length(value)>0) {data2save<-value} else {data2save<-NA}
      f1Jname.splt[[k]]<-data2save
    }
    f1Jname.best<-unlist(purrr::map(f1Jname.splt, 1))
      
    f1$V.name<-f1Vname.best
    f1$D.name<-f1Dname.best
    f1$J.name<-f1Jname.best
    saveRDS(f1, paste0("Data2analyze/Immunarch_Format_BestVDJ/",c,"/",i,".rds"))
  }
}

#################################  3. Merging sample-replicate  #################################

chains<-c("IGK","IGH","IGL","TRB","TRA","TRG","TRD")
count=0
dir.create("Data2analyze/Immunarch_Format_Merged/")
for (i in 1:length(chains)){
  c<-chains[i]
  system(paste0("mkdir Data2analyze/Immunarch_Format_Merged/",c))
  dir.path<-paste0("Data2analyze/Immunarch_Format_BestVDJ/",c)
  inds<-list.files(dir.path); inds<-unique(sort(inds[grep("_R",inds,invert=T)]))
  for (j in 1:length(inds)){
    count=count+1
    ind1<-inds[j]; ind<-gsub(".rds","",ind1)
    files<-Sys.glob(paste0("Data2analyze/Immunarch_Format_BestVDJ/",c,"/",ind,"*"))
    f1<-readRDS(files[1]); f2<-readRDS(files[2]);
    f3<-rbind(f1,f2)
    f4 <- f3 %>% dplyr::group_by(Sequence) %>% dplyr::summarise(Clones=sum(Clones), CDR3.nt=paste(unlist(unique(CDR3.nt)), collapse=','), CDR3.aa=paste(unlist(unique(CDR3.aa)), collapse=','), V.name=paste(unlist(unique(V.name)), collapse=','), D.name=paste(unlist(unique(D.name)), collapse=','), J.name=paste(unlist(unique(J.name)), collapse=','), V.end=paste(unlist(unique(V.end)), collapse=','), D.start=paste(unlist(unique(D.start)), collapse=','), D.end=paste(unlist(unique(D.end)), collapse=','), J.start=paste(unlist(unique(J.start)), collapse=','), VJ.ins=paste(unlist(unique(VJ.ins)), collapse=','), VD.ins=paste(unlist(unique(VD.ins)), collapse=','), DJ.ins=paste(unlist(unique(DJ.ins)), collapse=','))
    f4<-as.data.frame(f4)
    f4$Proportion<-f4$Clones/sum(f4$Clones)
    final.file<-f4[,c(2,15,3,4,5,6,7,8,9,10,11,12,13,14,1)]
    print(paste0(c,"-",count,"-",ind,"-",nrow(final.file)))
    out.filename<-paste0("Data2analyze/Immunarch_Format_Merged/",c,"/",ind)
    saveRDS(final.file, file=out.filename)
  }
}



#################################  3. Grouping clones #################################

dir.create(paste0("Data2analyze/Clones_By_CDR3aa/"))
chains<-c("TRG","TRD","TRB","TRA","IGK","IGH","IGL")
for (c in chains){
  dir.create(paste0("Data2analyze/Clones_By_CDR3aa/",c))
  inds<-list.files(paste0("Data2analyze/Immunarch_Format_Merged/",c,"/"))
  print(paste0("Grouping by CDR3aa ",c))
  for (i in inds){
    data<-readRDS(file=paste0("Data2analyze/Immunarch_Format_Merged/",c,"/",i))
    data2<- data %>% dplyr::group_by(CDR3.aa) %>% dplyr::summarise(Clones=sum(Clones), Proportion=sum(Proportion), CDR3.nt=paste(unlist(unique(CDR3.nt)), collapse=','), V.name=paste(unlist(unique(V.name)), collapse=','), D.name=paste(unlist(unique(D.name)), collapse=','), J.name=paste(unlist(unique(J.name)), collapse=','), V.end=paste(unlist(unique(V.end)), collapse=','), D.start=paste(unlist(unique(D.start)), collapse=','), D.end=paste(unlist(unique(D.end)), collapse=','), J.start=paste(unlist(unique(J.start)), collapse=','), VJ.ins=paste(unlist(unique(VJ.ins)), collapse=','), VD.ins=paste(unlist(unique(VD.ins)), collapse=','), DJ.ins=paste(unlist(unique(DJ.ins)), collapse=','), Sequence=paste(unlist(unique(Sequence)), collapse=','))
    final.file<-data2[,c(2,3,4,1,5,6,7,8,9,10,11,12,13,14,15)]
    out.filename<-paste0("Data2analyze/Clones_By_CDR3aa/",c,"/",i,".rds")
    saveRDS(final.file, file=out.filename)
  }
}
    
dir.create(paste0("Data2analyze/Clones_By_Vname/"))
chains<-c("TRG","TRD","TRB","TRA","IGK","IGH","IGL")
count=0
for (i in 1:length(chains)){
  c<-chains[i]
  inds<-list.files(paste0("Data2analyze/Immunarch_Format_Merged/",c))
  dir.create(paste0("Data2analyze/Clones_By_Vname/",c))
  for (j in 1:length(inds)){
    ind<-inds[j];    count=count+1;   print(paste0(c,"-",count,"-",ind))
    filename=paste0("Data2analyze/Immunarch_Format_Merged/",c,"/",ind)
    data<-readRDS(file=filename)
    data1 <- tibble(Clones=integer(),CDR3.nt=character(),CDR3.aa=character(),V.name=character(),D.name=character(),J.name=character(),V.end=character(),D.start=character(),D.end=character(),J.start=character(),VJ.ins=character(),VD.ins=character(),DJ.ins=character(),Sequence=character())
    for (k in 1:dim(data)[1]){
      line<-data[k,]
      line<-line[,-2]
      lineX<-gsub(" ","",data[k,5][[1]])
      Vnames<-unlist(strsplit(lineX,","))
      for (l in 1:length(Vnames)){
        Vn<-Vnames[l]
        line$V.name<-Vn
        data1<-data1 %>% bind_rows(line)
      }
    }
    data2<- data1 %>% dplyr::group_by(V.name) %>% dplyr::summarise(Clones=sum(Clones), CDR3.aa=paste(unlist(unique(CDR3.aa)), collapse=','), CDR3.nt=paste(unlist(unique(CDR3.nt)), collapse=','), D.name=paste(unlist(unique(D.name)), collapse=','), J.name=paste(unlist(unique(J.name)), collapse=','), V.end=paste(unlist(unique(V.end)), collapse=','), D.start=paste(unlist(unique(D.start)), collapse=','), D.end=paste(unlist(unique(D.end)), collapse=','), J.start=paste(unlist(unique(J.start)), collapse=','), VJ.ins=paste(unlist(unique(VJ.ins)), collapse=','), VD.ins=paste(unlist(unique(VD.ins)), collapse=','), DJ.ins=paste(unlist(unique(DJ.ins)), collapse=','), Sequence=paste(unlist(unique(Sequence)), collapse=','))
    data2$Proportion<-data2$Clones/sum(data2$Clones)
    final.file<-data2[,c(2,15,4,3,1,5,6,7,8,9,10,11,12,13,14)]
    out.filename<-paste0("Data2analyze/Clones_By_Vname/",c,"/",ind,".rds")
    saveRDS(final.file, file=out.filename)
  }
}

dir.create(paste0("Data2analyze/Clones_By_Jname/"))
chains<-c("TRG","TRD","TRB","TRA","IGK","IGH","IGL")
count=0
for (i in 1:length(chains)){
  c<-chains[i]
  inds<-list.files(paste0("Data2analyze/Immunarch_Format_Merged/",c))
  dir.create(paste0("Data2analyze/Clones_By_Jname/",c))
  for (j in 1:length(inds)){
    ind<-inds[j];    count=count+1;   print(paste0(c,"-",count,"-",ind))
    filename=paste0("Data2analyze/Immunarch_Format_Merged/",c,"/",ind)
    data<-readRDS(file=filename)
    data1 <- tibble(Clones=integer(),CDR3.nt=character(),CDR3.aa=character(),V.name=character(),D.name=character(),J.name=character(),V.end=character(),D.start=character(),D.end=character(),J.start=character(),VJ.ins=character(),VD.ins=character(),DJ.ins=character(),Sequence=character())
    
    for (k in 1:dim(data)[1]){
      line<-data[k,]
      line<-line[,-2]
      lineX<-gsub(" ","",data[k,7][[1]])
      Jnames<-unlist(strsplit(lineX,","))
      for (l in 1:length(Jnames)){
        Jn<-Jnames[l]
        line$J.name<-Jn
        data1<-data1 %>% bind_rows(line)
      }
    }
    data2<- data1 %>% dplyr::group_by(J.name) %>% dplyr::summarise(Clones=sum(Clones), CDR3.aa=paste(unlist(unique(CDR3.aa)), collapse=','), CDR3.nt=paste(unlist(unique(CDR3.nt)), collapse=','), D.name=paste(unlist(unique(D.name)), collapse=','), V.name=paste(unlist(unique(V.name)), collapse=','), V.end=paste(unlist(unique(V.end)), collapse=','), D.start=paste(unlist(unique(D.start)), collapse=','), D.end=paste(unlist(unique(D.end)), collapse=','), J.start=paste(unlist(unique(J.start)), collapse=','), VJ.ins=paste(unlist(unique(VJ.ins)), collapse=','), VD.ins=paste(unlist(unique(VD.ins)), collapse=','), DJ.ins=paste(unlist(unique(DJ.ins)), collapse=','), Sequence=paste(unlist(unique(Sequence)), collapse=','))
    data2$Proportion<-data2$Clones/sum(data2$Clones)
    final.file<-data2[,c(2,15,4,3,1,5,6,7,8,9,10,11,12,13,14)]
    out.filename<-paste0("Data2analyze/Clones_By_Jname/",c,"/",ind,".rds")
    saveRDS(final.file, file=out.filename)
  }
}
    
dir.create(paste0("Data2analyze/Clones_By_VJname/"))
chains<-c("TRG","TRD","TRB","TRA","IGK","IGH","IGL")
count=0
for (i in 1:length(chains)){
  c<-chains[i]
  inds<-list.files(paste0("Data2analyze/Immunarch_Format_Merged/",c))
  dir.create(paste0("Data2analyze/Clones_By_VJname/",c))
  for (j in 1:length(inds)){
    ind<-inds[j];    count=count+1;   print(paste0(c,"-",count,"-",ind))
    filename=paste0("Data2analyze/Immunarch_Format_Merged/",c,"/",ind)
    data<-readRDS(file=filename)
    
    data1 <- tibble(Clones=integer(),CDR3.nt=character(),CDR3.aa=character(),V.name=character(),D.name=character(),J.name=character(),V.end=character(),D.start=character(),D.end=character(),J.start=character(),VJ.ins=character(),VD.ins=character(),DJ.ins=character(),Sequence=character())
    
    for (k in 1:nrow(data)){
      line<-data[k,-2]
      line$VJ.name<-"ToFill"
      lineJ<-gsub(" ","",data[k,7][[1]])
      lineV<-gsub(" ","",data[k,5][[1]])
      Jnames<-unlist(strsplit(lineJ,","))
      Vnames<-unlist(strsplit(lineV,","))
      VJnames.1<-expand.grid(Vnames,Jnames)
      VJnames<-paste0(VJnames.1[,"Var1"],".",VJnames.1[,"Var2"])
      for (VJn in VJnames){
        line$VJ.name<-VJn
        data1<-data1 %>% bind_rows(line)
      }
    }
    data1<-unique(data1)
    data2<- data1 %>% dplyr::group_by(VJ.name) %>% dplyr::summarise(Clones=sum(Clones), CDR3.aa=paste(unlist(unique(CDR3.aa)), collapse=','), CDR3.nt=paste(unlist(unique(CDR3.nt)), collapse=','), J.name=paste(unlist(unique(J.name)), collapse=','), D.name=paste(unlist(unique(D.name)), collapse=','), V.name=paste(unlist(unique(V.name)), collapse=','), V.end=paste(unlist(unique(V.end)), collapse=','), D.start=paste(unlist(unique(D.start)), collapse=','), D.end=paste(unlist(unique(D.end)), collapse=','), J.start=paste(unlist(unique(J.start)), collapse=','), VJ.ins=paste(unlist(unique(VJ.ins)), collapse=','), VD.ins=paste(unlist(unique(VD.ins)), collapse=','), DJ.ins=paste(unlist(unique(DJ.ins)), collapse=','), Sequence=paste(unlist(unique(Sequence)), collapse=','))
    data2$Proportion<-data2$Clones/sum(data2$Clones)
    final.file<-data2[,c(2,15,4,3,1,5,6,7,8,9,10,11,12,13,14)]
    out.filename<-paste0("Data2analyze/Clones_By_VJname/",c,"/",ind,".rds")
    saveRDS(final.file, file=out.filename)
  }
}  
    
#################################  4. Creating metadata files
#################################  5. Merging clone info and metadata in a single object  #################################

chains<-c("TRG","TRD","TRB","TRA","IGK","IGH","IGL")

# By CDR3aa
for (chain in chains){
  filenames<-list.files(paste0("Data2analyze/Clones_By_CDR3aa/",chain,"/"));
  data<-list(); names(data)<-vector()
  for (i in 1:length(filenames)){
    d1<-readRDS(paste0("Data2analyze/Clones_By_CDR3aa/",chain,"/",filenames[i]))
    data[[i]]<-d1;    names(data)[i]<-gsub(".rds", "", filenames[i])
  }
  metadata<-readRDS(paste0("Data2analyze/Metadata_",chain,".rds"))
  grr<-list(data, metadata); names(grr)<-c("data","meta")
  checking<-table(names(grr$data)==grr$meta$Sample)
  print(checking)
  names(grr$data)<-as.character(grr$meta$Sample)
  saveRDS(grr, file=paste0("Data2analyze/Clones_By_CDR3aa/",chain,".rds"))
  rm(grr,metadata)
  gc()
}

# By Vname
for (chain in chains){
  filenames<-list.files(paste0("Data2analyze/Clones_By_Vname/",chain,"/"));
  data<-list(); names(data)<-vector()
  for (i in 1:length(filenames)){
    d1<-readRDS(paste0("Data2analyze/Clones_By_Vname/",chain,"/",filenames[i]))
    data[[i]]<-d1;    names(data)[i]<-gsub(".rds", "", filenames[i])
  }
  metadata<-readRDS(paste0("Data2analyze/Metadata_",chain,".rds"))
  grr<-list(data, metadata); names(grr)<-c("data","meta")
  checking<-table(names(grr$data)==grr$meta$Sample)
  print(checking)
  names(grr$data)<-as.character(grr$meta$Sample)
  saveRDS(grr, file=paste0("Data2analyze/Clones_By_Vname/",chain,".rds"))
  rm(grr,metadata)
  gc()
}

# By Jname
for (chain in chains){
  filenames<-list.files(paste0("Data2analyze/Clones_By_Jname/",chain,"/"));
  data<-list(); names(data)<-vector()
  for (i in 1:length(filenames)){
    d1<-readRDS(paste0("Data2analyze/Clones_By_Jname/",chain,"/",filenames[i]))
    data[[i]]<-d1;    names(data)[i]<-gsub(".rds", "", filenames[i])
  }
  metadata<-readRDS(paste0("Data2analyze/Metadata_",chain,".rds"))
  grr<-list(data, metadata); names(grr)<-c("data","meta")
  checking<-table(names(grr$data)==grr$meta$Sample)
  print(checking)
  names(grr$data)<-as.character(grr$meta$Sample)
  saveRDS(grr, file=paste0("Data2analyze/Clones_By_Jname/",chain,".rds"))
  rm(grr,metadata)
  gc()
}

# By VJname
for (chain in chains){
  filenames<-list.files(paste0("Data2analyze/Clones_By_VJname/",chain,"/"));
  data<-list(); names(data)<-vector()
  for (i in 1:length(filenames)){
    d1<-readRDS(paste0("Data2analyze/Clones_By_VJname/",chain,"/",filenames[i]))
    data[[i]]<-d1;    names(data)[i]<-gsub(".rds", "", filenames[i])
  }
  metadata<-readRDS(paste0("Data2analyze/Metadata_",chain,".rds"))
  grr<-list(data, metadata); names(grr)<-c("data","meta")
  checking<-table(names(grr$data)==grr$meta$Sample)
  print(checking)
  names(grr$data)<-as.character(grr$meta$Sample)
  saveRDS(grr, file=paste0("Data2analyze/Clones_By_VJname/",chain,".rds"))
  rm(grr,metadata)
  gc()
}

