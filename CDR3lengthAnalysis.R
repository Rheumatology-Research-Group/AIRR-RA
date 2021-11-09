
# Loading libraries
    
    suppressPackageStartupMessages({
      library(immunarch)
      library(knitr)
      library(data.table)
      library(dtplyr)
      library(dplyr)
      library(tidyverse)
      library(ggpubr)
      library(forecast)
      library(MASS)
      library(vegan)
      library(gridExtra)
      library(grid)
      library(RColorBrewer)
      library(car)
      library(pheatmap)
      library(venneuler)
      library(eulerr)
      library(VennDiagram)
      library(gplots)
      library(stringdist)
      library(viridis)
      library(pals)
      library(cowplot)
      library(motifStack)
      library(protr)
      library(msa)
      library(Biostrings)
      library(ggseqlogo)
      library(ggpubr)
      library(dgof)
    })
      
    get.pheno.data<-function(pheno) {
        if (pheno=="IMID") {pheno.out<-"IMID"; pheno.name<-"IMID"; pheno1<-"RA"; pheno2<-"CTRL"; week<-"wk0"; analysis<-"bin"}
        if (pheno=="RespGM") {pheno.out<-"RespGM"; pheno.name<-"Resp"; pheno1<-"NON_RESP"; pheno2<-"GOOD"; week<-"wk0"; analysis<-"bin"}
        if (pheno=="ActB") {pheno.out<-"ActB"; pheno.name<-"ActB"; pheno1<-"High"; pheno2<-"Low"; week<-"wk0"; analysis<-"bin"}
        if (pheno=="ACPA") {pheno.out<-"ACPA"; pheno.name<-"ACPA"; pheno1<-"Pos"; pheno2<-"Neg"; week<-"wk0"; analysis<-"bin"}
        if (pheno=="RF") {pheno.out<-"RF"; pheno.name<-"RF"; pheno1<-"Pos"; pheno2<-"Neg"; week<-"wk0"; analysis<-"bin"}
        return(c(pheno.out,pheno.name,pheno1,pheno2,week,analysis))
      }
    make.chain.plots.baseline<-function(c,ph,pheno1,pheno2,week,analysis){
        grr<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'.rds'))
        grr$meta$RespGM<-grr$meta$Resp;   grr$meta$RespGM[grr$meta$RespGM=="MOD"]<-"GOOD";   grr$meta$Resp[grr$meta$Resp=="MOD"]<-NA;
        meta_col<-dplyr::filter(grr$meta, Week!="wk12");  samples_col<-meta_col$Sample;
        data_col<-grr$data[match(samples_col,names(grr$data))];  grr<-list(data_col,meta_col);  names(grr)<-c("data","meta")
        colnames(grr$meta)[which(colnames(grr$meta)==ph)]<-"Pheno"
        
        meta_col<-dplyr::filter(grr$meta, (Pheno==pheno1 & !is.na(Pheno)));   data_col<-grr$data[match(meta_col$Sample,names(grr$data))];    grr_p1<-list(data_col,meta_col);  names(grr_p1)<-c("data","meta")
        meta_col<-dplyr::filter(grr$meta, (Pheno==pheno2 & !is.na(Pheno)));   data_col<-grr$data[match(meta_col$Sample,names(grr$data))];    grr_p2<-list(data_col,meta_col);  names(grr_p2)<-c("data","meta")
        
        clones.p1<-as.data.frame(pubRep(grr_p1$data, "aa", .verbose=F))
        clones.p2<-as.data.frame(pubRep(grr_p2$data, "aa", .verbose=F))
        
        clones.p1.abundance<-rowSums(clones.p1[,-c(1,2)], na.rm=T);  clones.p1.df<-data.frame("Clone"=clones.p1$CDR3.aa, "Length"=str_length(clones.p1$CDR3.aa), "Abundance"=clones.p1.abundance)
        clones.p2.abundance<-rowSums(clones.p2[,-c(1,2)], na.rm=T);  clones.p2.df<-data.frame("Clone"=clones.p2$CDR3.aa, "Length"=str_length(clones.p2$CDR3.aa), "Abundance"=clones.p2.abundance)
        
        #  Wilcoxon
        
        # Weighted by abundance
        p1.vect<-rep(clones.p1.df$Length,clones.p1.df$Abundance); p1.ab.df<-data.frame("Length"=p1.vect, "Pheno"=rep(pheno1,length(p1.vect)))
        p2.vect<-rep(clones.p2.df$Length,clones.p2.df$Abundance); p2.ab.df<-data.frame("Length"=p2.vect, "Pheno"=rep(pheno2,length(p2.vect)))
        p1p2.ab.df<-rbind(p1.ab.df,p2.ab.df)
        p1.prop<-table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)[,pheno1];  p1.prop.df<-data.frame("Length"=as.numeric(names(p1.prop)), "Abundance"=unname(p1.prop), "Perc"=100*unname(p1.prop)/sum(unname(p1.prop)), "Pheno"=rep(pheno1,length(p1.prop)))
        p2.prop<-table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)[,pheno2];  p2.prop.df<-data.frame("Length"=as.numeric(names(p2.prop)), "Abundance"=unname(p2.prop), "Perc"=100*unname(p2.prop)/sum(unname(p2.prop)), "Pheno"=rep(pheno2,length(p2.prop)))
        p1p2.prop.df<-rbind(p1.prop.df,p2.prop.df)
        pval.w<-wilcox.test(sort(p1.vect),sort(p2.vect))$p.value
        pval.ks<-ks.test(sort(p1.vect),sort(p2.vect))$p.value
        
        df.fisher.w<-data.frame("Length"=rownames(table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)),
                                "Pheno1.Yes"=as.numeric(table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)[,pheno1]),
                                "Pheno1.No"=as.numeric(sum(table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)[,pheno1])-table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)[,pheno1]),
                                "Pheno2.Yes"=as.numeric(table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)[,pheno2]),
                                "Pheno2.No"=as.numeric(sum(table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)[,pheno2])-table(p1p2.ab.df$Length,p1p2.ab.df$Pheno)[,pheno2]), stringsAsFactors=F)
        
        df.fisher.w$Pval<-rep(0,nrow(df.fisher.w))
        for (i in 1:nrow(df.fisher.w)) {
          p1.yes<-df.fisher.w[i,2];   p1.no<-df.fisher.w[i,3];   p2.yes<-df.fisher.w[i,4];   p2.no<-df.fisher.w[i,5];
          mat2analyze<-matrix(c(p1.yes,p1.no,p2.yes,p2.no), nrow=2)
          pval<-fisher.test(mat2analyze)$p.value
          df.fisher.w$Pval[i]<-pval
        }
        
        fisher.sig<-df.fisher.w[df.fisher.w$Pval<0.05,]
        fisher.sig$lab<-rep("*",nrow(fisher.sig));   fisher.sig$lab[fisher.sig$Pval<0.005]<-"**";
        
        summ.p1<-c(summary(p1.vect), sd(p1.vect))
        names(summ.p1)[7]="Std."
        summ.p2<-c(summary(p2.vect), sd(p2.vect))
        names(summ.p2)[7]="Std."
        summ.p1.plot<-paste0(pheno1,"\n",names(summ.p1)[1]," = ", summ.p1[[1]],"\n",names(summ.p1)[2]," = ", summ.p1[[2]],"\n",
                             names(summ.p1)[3]," = ", round(summ.p1[[3]],2),"\n",names(summ.p1)[4]," = ", round(summ.p1[[4]],2),"\n",
                             names(summ.p1)[5]," = ", summ.p1[[5]],"\n",names(summ.p1)[6]," = ", summ.p1[[6]],"\n", names(summ.p1)[7]," = ", round(summ.p1[[7]],2),"\n")
        summ.p2.plot<-paste0(pheno2,"\n",names(summ.p2)[1]," = ", summ.p2[[1]],"\n",names(summ.p2)[2]," = ", summ.p2[[2]],"\n",
                             names(summ.p2)[3]," = ", round(summ.p2[[3]],2),"\n",names(summ.p2)[4]," = ", round(summ.p2[[4]],2),"\n",
                             names(summ.p2)[5]," = ", summ.p2[[5]],"\n",names(summ.p2)[6]," = ", summ.p2[[6]],"\n", names(summ.p2)[7]," = ", round(summ.p2[[7]],2),"\n")
        
        if (pval.w==0){pval.text<-"(P<1.00e-300)\n"}
        if (pval.w>0){pval.text<-paste0("(P=",formatC(as.numeric(pval.w), format="e", digits=2),")\n")}
        plot.w<-ggplot(p1p2.prop.df, aes(x=Length, y=Perc, fill=Pheno), color="black") + geom_bar(stat="identity", position=position_dodge()) +
          theme_classic() + scale_fill_manual(values=c("#F69191","#7DBB7D")) + theme_classic() +
          labs(title=paste0("\n",c,": CDR3 amino acid length distribution  ",pval.text),
               y="\nFrequency (%)\n", x="\nCDR3 amino acid length (N)\n") +
          theme(legend.position=c(0.9,0.9),legend.title=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5)) +
          scale_x_continuous(labels=as.character(p1p2.prop.df$Length),breaks=p1p2.prop.df$Length) + 
          annotate(geom="text", x=as.numeric(rownames(fisher.sig)), y=rep(max(p1p2.prop.df$Perc)*1.1,nrow(fisher.sig)), label=fisher.sig$lab, color="black", size=3) +
          annotate(geom="text", x=max(p1p2.prop.df$Length)-2, y=max(p1p2.prop.df$Perc)*0.7, label=summ.p1.plot, size=3, col="#F69191") +
          annotate(geom="text", x=max(p1p2.prop.df$Length)-2, y=max(p1p2.prop.df$Perc)*0.32, label=summ.p2.plot, size=3, col="#7DBB7D")
        
        plot.cum.w<-ggplot(p1p2.ab.df, aes(Length, color=Pheno)) + stat_ecdf(geom="step") + theme_classic() + scale_color_manual(values=c("#F69191","#7DBB7D")) +
          labs(title=paste0("\n",c,": Empirical Cumulative Density Funcion\n"), y="\nCumulative Probability\n", x="\nCDR3 amino acid length (N)\n") +
          theme(legend.position=c(0.9,0.1),axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5))
        
        plot.qqplot.w<-ggplot(mapping=aes(x=quantile(p1p2.ab.df[p1p2.ab.df$Pheno==pheno1,"Length"], seq(0,1,0.01)), y=quantile(p1p2.ab.df[p1p2.ab.df$Pheno==pheno2,"Length"], seq(0,1,0.01)))) + 
          geom_point() + geom_abline(aes(slope=1, intercept=0), linetype=2)  + theme_classic() +
          theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5)) +
          labs(x=paste("\n",pheno1, " Quantiles\n"), y=paste("\n",pheno2, "Quantiles\n"), title=paste0("\n",c,": Empirical QQ-plot\n"))
        
        return(list(plot.w, plot.cum.w, plot.qqplot.w))
      }
    make.chain.plots.longAll<-function(c,ph,pheno1,pheno2,week,analysis) {
        ctrl<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'_CTRL.rds')) # This file needs to be manually generated so that it includes only AIRR-seq data for controls 
        ra.wk0<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'_WK0.rds')) # This file needs to be manually generated so that it includes only AIRR-seq data for RA patients at baseline 
        ra.wk12<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'_WK12.rds')) # This file needs to be manually generated so that it includes only AIRR-seq data for RA patients after 12 weeks of TNFi therapy 
        
        x1<-read.table("Data2analyze/Donor2RNA_Codes.txt", head=T)
        samples.wk0<-as.character(x1$RNA1)
        samples.wk12<-x1[(as.character(x1$RNA1) %in% samples.wk0),"RNA2"]
        
        meta_col<-dplyr::filter(ra.wk0$meta, Week=="wk0");  meta_col<-meta_col[(meta_col$Sample %in% samples.wk0),];   samples_col<-meta_col$Sample;
        data_col<-ra.wk0$data[match(samples_col,names(ra.wk0$data))]; grr.ra.wk0<-list(data_col,meta_col);                       names(grr.ra.wk0)<-c("data","meta")
        
        meta_col<-dplyr::filter(ra.wk12$meta, Week=="wk12");              meta_col<-meta_col[(meta_col$Sample %in% samples.wk12),];   samples_col<-meta_col$Sample;
        data_col<-ra.wk12$data[match(samples_col,names(ra.wk12$data))];   grr.ra.wk12<-list(data_col,meta_col);                       names(grr.ra.wk12)<-c("data","meta")    
        
        clones.ctrl<-as.data.frame(pubRep(ctrl$data, "aa", .verbose=F))
        clones.wk0<-as.data.frame(pubRep(grr.ra.wk0$data, "aa", .verbose=F))
        clones.wk12<-as.data.frame(pubRep(grr.ra.wk12$data, "aa", .verbose=F))
        
        clones.ctrl.abundance<-rowSums(clones.ctrl[,-c(1,2)], na.rm=T);   clones.ctrl.df<-data.frame("Clone"=clones.ctrl$CDR3.aa, "Length"=str_length(clones.ctrl$CDR3.aa), "Abundance"=clones.ctrl.abundance)
        clones.wk0.abundance<-rowSums(clones.wk0[,-c(1,2)], na.rm=T);     clones.wk0.df<-data.frame("Clone"=clones.wk0$CDR3.aa, "Length"=str_length(clones.wk0$CDR3.aa), "Abundance"=clones.wk0.abundance)
        clones.wk12.abundance<-rowSums(clones.wk12[,-c(1,2)], na.rm=T);   clones.wk12.df<-data.frame("Clone"=clones.wk12$CDR3.aa, "Length"=str_length(clones.wk12$CDR3.aa), "Abundance"=clones.wk12.abundance)
        
        wk0.vect<-rep(clones.wk0.df$Length,clones.wk0.df$Abundance);     wk0.ab.df<-data.frame("Length"=wk0.vect,"Pheno"=rep("Week0",length(wk0.vect)))
        wk12.vect<-rep(clones.wk12.df$Length,clones.wk12.df$Abundance);  wk12.ab.df<-data.frame("Length"=wk12.vect, "Pheno"=rep("Week12",length(wk12.vect)))
        wk0wk12.ab.df<-rbind(wk0.ab.df,wk12.ab.df)
        
        wk0.prop<-table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)[,"Week0"];    wk0.prop.df<-data.frame("Length"=as.numeric(names(wk0.prop)), "Abundance"=unname(wk0.prop), "Perc"=100*unname(wk0.prop)/sum(unname(wk0.prop)), "Pheno"=rep("Week0",length(wk0.prop)))
        wk12.prop<-table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)[,"Week12"];  wk12.prop.df<-data.frame("Length"=as.numeric(names(wk12.prop)), "Abundance"=unname(wk12.prop), "Perc"=100*unname(wk12.prop)/sum(unname(wk12.prop)), "Pheno"=rep("Week12",length(wk12.prop)))
        wk0wk12.prop.df<-rbind(wk0.prop.df,wk12.prop.df)
        
        pval.w<-wilcox.test(sort(wk0.vect),sort(wk12.vect))$p.value
        pval.ks<-ks.test(sort(wk0.vect),sort(wk12.vect))$p.value
        
        df.fisher.w<-data.frame("Length"=rownames(table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)),
                                "Pheno1.Yes"=as.numeric(table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)[,"Week0"]),
                                "Pheno1.No"=as.numeric(sum(table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)[,"Week0"])-table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)[,"Week0"]),
                                "Pheno2.Yes"=as.numeric(table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)[,"Week12"]),
                                "Pheno2.No"=as.numeric(sum(table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)[,"Week12"])-table(wk0wk12.ab.df$Length,wk0wk12.ab.df$Pheno)[,"Week12"]), stringsAsFactors=F)
        
        df.fisher.w$Pval<-rep(0,nrow(df.fisher.w))
        for (i in 1:nrow(df.fisher.w)) {
          p1.yes<-df.fisher.w[i,2];   p1.no<-df.fisher.w[i,3];   p2.yes<-df.fisher.w[i,4];   p2.no<-df.fisher.w[i,5];
          mat2analyze<-matrix(c(p1.yes,p1.no,p2.yes,p2.no), nrow=2)
          pval<-fisher.test(mat2analyze)$p.value
          df.fisher.w$Pval[i]<-pval
        }
        
        fisher.sig<-df.fisher.w[df.fisher.w$Pval<0.05,]
        fisher.sig$lab<-rep("*",nrow(fisher.sig));   fisher.sig$lab[fisher.sig$Pval<0.005]<-"**";
        
        summ.p1<-c(summary(wk0.vect), sd(wk0.vect))
        names(summ.p1)[7]="Std."
        summ.p2<-c(summary(wk12.vect), sd(wk12.vect))
        names(summ.p2)[7]="Std."
        summ.p1.plot<-paste0("BASELINE\n",names(summ.p1)[1]," = ", summ.p1[[1]],"\n",names(summ.p1)[2]," = ", summ.p1[[2]],"\n",
                             names(summ.p1)[3]," = ", round(summ.p1[[3]],2),"\n",names(summ.p1)[4]," = ", round(summ.p1[[4]],2),"\n",
                             names(summ.p1)[5]," = ", summ.p1[[5]],"\n",names(summ.p1)[6]," = ", summ.p1[[6]],"\n", names(summ.p1)[7]," = ", round(summ.p1[[7]],2),"\n")
        summ.p2.plot<-paste0("WEEK.12\n",names(summ.p2)[1]," = ", summ.p2[[1]],"\n",names(summ.p2)[2]," = ", summ.p2[[2]],"\n",
                             names(summ.p2)[3]," = ", round(summ.p2[[3]],2),"\n",names(summ.p2)[4]," = ", round(summ.p2[[4]],2),"\n",
                             names(summ.p2)[5]," = ", summ.p2[[5]],"\n",names(summ.p2)[6]," = ", summ.p2[[6]],"\n", names(summ.p2)[7]," = ", round(summ.p2[[7]],2),"\n")
        
        if (pval.w==0){pval.text<-"(P<1.00e-300)\n"}
        if (pval.w>0){pval.text<-paste0("(P=",formatC(as.numeric(pval.w), format="e", digits=2),")\n")}
        
        plot.w<-ggplot(wk0wk12.prop.df, aes(x=Length, y=Perc, fill=Pheno), color="black") + geom_bar(stat="identity", position=position_dodge()) +
          theme_classic() + scale_fill_manual(values=c("red", "blue")) + theme_classic() +
          labs(title=paste0("\n",c,": CDR3 amino acid length distribution  ",pval.text),
               y="\nFrequency (%)\n", x="\nCDR3 amino acid length (N)\n") +
          theme(legend.position=c(0.9,0.9), legend.title=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5)) +
          scale_x_continuous(labels=as.character(wk0wk12.prop.df$Length),breaks=wk0wk12.prop.df$Length) + 
          annotate(geom="text", x=as.numeric(rownames(fisher.sig)), y=rep(max(wk0wk12.prop.df$Perc)*1.1,nrow(fisher.sig)), label=fisher.sig$lab, color="black", size=3) +
          annotate(geom="text", x=max(wk0wk12.prop.df$Length)-2, y=max(wk0wk12.prop.df$Perc)*0.67, label=summ.p1.plot, size=3, col="red") +
          annotate(geom="text", x=max(wk0wk12.prop.df$Length)-2, y=max(wk0wk12.prop.df$Perc)*0.32, label=summ.p2.plot, size=3, col="blue")
        
        ctrl.vect<-rep(clones.ctrl.df$Length,clones.ctrl.df$Abundance);     ctrl.ab.df<-data.frame("Length"=ctrl.vect,"Pheno"=rep("CTRL",length(ctrl.vect)))
        wk0wk12ctrl.ab.df<-rbind(wk0wk12.ab.df,ctrl.ab.df)
        
        plot.cum.w<-ggplot(wk0wk12ctrl.ab.df, aes(Length, color=Pheno)) + stat_ecdf(geom="step") + theme_classic() + scale_color_manual(values=c("red", "blue", "chartreuse4")) +
          labs(title=paste0("\n",c,": Empirical Cumulative Density Funcion\n"), y="\nCumulative Probability\n", x="\nCDR3 amino acid length (N)\n") +
          theme(legend.position=c(0.9,0.2),axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5))
        
        plot.qqplot.w<-ggplot(mapping=aes(x=quantile(wk0wk12.ab.df[wk0wk12.ab.df$Pheno=="Week0","Length"], seq(0,1,0.01)), y=quantile(wk0wk12.ab.df[wk0wk12.ab.df$Pheno=="Week12","Length"], seq(0,1,0.01)))) + 
          geom_point() + geom_abline(aes(slope=1, intercept=0), linetype=2)  + theme_classic() +
          theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5)) +
          labs(x=paste("\nWeek.0 Quantiles\n"), y=paste("\nWeek.12 Quantiles\n"), title=paste("\n",c,": Empirical QQ-plot\n"))
        
        return(list(plot.w, plot.cum.w, plot.qqplot.w))
    }
    make.chain.plots.longByResp<-function(c,ph,pheno1,pheno2,week,analysis) {
      
      ctrl<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'_CTRL.rds'))  # This file needs to be manually generated so that it includes only AIRR-seq data for controls 
      
      x1<-read.table("Data2analyze/Donor2RNA_Codes.txt", header=T)
      
      ra.wk0<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'_WK0.rds')) # This file needs to be manually generated so that it includes only AIRR-seq data for RA patients at baseline
      ra.wk0$meta$RespGM<-ra.wk0$meta$Resp;   ra.wk0$meta$RespGM[ra.wk0$meta$RespGM=="MOD"]<-"GOOD";   ra.wk0$meta$Resp[ra.wk0$meta$Resp=="MOD"]<-NA;
      meta_col<-dplyr::filter(ra.wk0$meta, IMID=="RA");
      meta_col<-meta_col[(meta_col$Sample %in% as.character(x1$RNA1)),]
      samples_col<-meta_col$Sample;  data_col<-ra.wk0$data[match(samples_col,names(ra.wk0$data))];  ra.wk0<-list(data_col,meta_col);  names(ra.wk0)<-c("data","meta")
      
      ra.wk12<-readRDS(paste0('Data2analyze/Clones_By_CDR3aa/',c,'_WK12.rds')) # This file needs to be manually generated so that it includes only AIRR-seq data for RA patients after 12 weeks of TNFi therapy
      ra.wk12$meta$RespGM<-ra.wk12$meta$Resp;   ra.wk12$meta$RespGM[ra.wk12$meta$RespGM=="MOD"]<-"GOOD";   ra.wk12$meta$Resp[ra.wk12$meta$Resp=="MOD"]<-NA;
      sample.act.wk0<-as.data.frame(ra.wk0$meta[,c("Sample","ActB")]);     colnames(sample.act.wk0)[1]<-"RNA1"
      meta_col<-dplyr::filter(ra.wk12$meta, IMID=="RA");
      meta_col<-meta_col[(meta_col$Sample %in% as.character(x1$RNA2)),]
      samples_col<-meta_col$Sample;  data_col<-ra.wk12$data[match(samples_col,names(ra.wk12$data))];  ra.wk12<-list(data_col,meta_col);  names(ra.wk12)<-c("data","meta")
      data.act<-merge(x1,sample.act.wk0, by="RNA1")
      ra.wk12$meta$ActB<-data.act$ActB[match(data.act$RNA2,names(ra.wk12$data))]
      ra.wk12$meta$RespGM<-ra.wk12$meta$Resp
      ra.wk12$meta$RespGM[ra.wk12$meta$RespGM=="MOD"]<-"GOOD"
      ra.wk12$meta$Resp[ra.wk12$meta$Resp=="MOD"]<-NA
      
      # Preparing Data for Pheno
      colnames(ra.wk0$meta)[which(colnames(ra.wk0$meta)==ph)]<-"Pheno"
      colnames(ra.wk12$meta)[which(colnames(ra.wk12$meta)==ph)]<-"Pheno"
      
      meta_col<-dplyr::filter(ra.wk0$meta, (Pheno==pheno1 & !is.na(Pheno)));    data_col<-ra.wk0$data[match(meta_col$Sample,names(ra.wk0$data))];      ra.wk0.p1<-list(data_col,meta_col);  names(ra.wk0.p1)<-c("data","meta")
      meta_col<-dplyr::filter(ra.wk0$meta, (Pheno==pheno2 & !is.na(Pheno)));    data_col<-ra.wk0$data[match(meta_col$Sample,names(ra.wk0$data))];      ra.wk0.p2<-list(data_col,meta_col);  names(ra.wk0.p2)<-c("data","meta")    
      
      meta_col<-dplyr::filter(ra.wk12$meta, (Pheno==pheno1 & !is.na(Pheno)));   data_col<-ra.wk12$data[match(meta_col$Sample,names(ra.wk12$data))];    ra.wk12.p1<-list(data_col,meta_col);  names(ra.wk12.p1)<-c("data","meta")
      meta_col<-dplyr::filter(ra.wk12$meta, (Pheno==pheno2 & !is.na(Pheno)));   data_col<-ra.wk12$data[match(meta_col$Sample,names(ra.wk12$data))];    ra.wk12.p2<-list(data_col,meta_col);  names(ra.wk12.p2)<-c("data","meta") 
      
      # Computing clones
      clones.ctrl<-as.data.frame(pubRep(ctrl$data, "aa", .verbose=F))
      clones.wk0.p1<-as.data.frame(pubRep(ra.wk0.p1$data, "aa", .verbose=F))
      clones.wk0.p2<-as.data.frame(pubRep(ra.wk0.p2$data, "aa", .verbose=F))
      clones.wk12.p1<-as.data.frame(pubRep(ra.wk12.p1$data, "aa", .verbose=F))
      clones.wk12.p2<-as.data.frame(pubRep(ra.wk12.p2$data, "aa", .verbose=F))
      
      clones.ctrl.abundance<-rowSums(clones.ctrl[,-c(1,2)], na.rm=T);  clones.ctrl.df<-data.frame("Clone"=clones.ctrl$CDR3.aa, "Length"=str_length(clones.ctrl$CDR3.aa), "Abundance"=clones.ctrl.abundance) 
      clones.wk0.p1.abundance<-rowSums(clones.wk0.p1[,-c(1,2)], na.rm=T);  clones.wk0.p1.df<-data.frame("Clone"=clones.wk0.p1$CDR3.aa, "Length"=str_length(clones.wk0.p1$CDR3.aa), "Abundance"=clones.wk0.p1.abundance)
      clones.wk0.p2.abundance<-rowSums(clones.wk0.p2[,-c(1,2)], na.rm=T);  clones.wk0.p2.df<-data.frame("Clone"=clones.wk0.p2$CDR3.aa, "Length"=str_length(clones.wk0.p2$CDR3.aa), "Abundance"=clones.wk0.p2.abundance)
      clones.wk12.p1.abundance<-rowSums(clones.wk12.p1[,-c(1,2)], na.rm=T);  clones.wk12.p1.df<-data.frame("Clone"=clones.wk12.p1$CDR3.aa, "Length"=str_length(clones.wk12.p1$CDR3.aa), "Abundance"=clones.wk12.p1.abundance)
      clones.wk12.p2.abundance<-rowSums(clones.wk12.p2[,-c(1,2)], na.rm=T);  clones.wk12.p2.df<-data.frame("Clone"=clones.wk12.p2$CDR3.aa, "Length"=str_length(clones.wk12.p2$CDR3.aa), "Abundance"=clones.wk12.p2.abundance)    
      
      # Weighted by Abundance
      
      # Data Generation
      
      wk0.p1.vect<-rep(clones.wk0.p1.df$Length,clones.wk0.p1.df$Abundance);        wk0.p1.ab.df<-data.frame("Length"=wk0.p1.vect,"Pheno"=rep("Week0",length(wk0.p1.vect)))
      wk12.p1.vect<-rep(clones.wk12.p1.df$Length,clones.wk12.p1.df$Abundance);     wk12.p1.ab.df<-data.frame("Length"=wk12.p1.vect, "Pheno"=rep("Week12",length(wk12.p1.vect)))
      wk0wk12.p1.ab.df<-rbind(wk0.p1.ab.df,wk12.p1.ab.df)
      
      wk0.p2.vect<-rep(clones.wk0.p2.df$Length,clones.wk0.p2.df$Abundance);        wk0.p2.ab.df<-data.frame("Length"=wk0.p2.vect,"Pheno"=rep("Week0",length(wk0.p2.vect)))
      wk12.p2.vect<-rep(clones.wk12.p2.df$Length,clones.wk12.p2.df$Abundance);     wk12.p2.ab.df<-data.frame("Length"=wk12.p2.vect, "Pheno"=rep("Week12",length(wk12.p2.vect)))
      wk0wk12.p2.ab.df<-rbind(wk0.p2.ab.df,wk12.p2.ab.df)  
      
      wk0.p1.prop<-table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)[,"Week0"];    wk0.p1.prop.df<-data.frame("Length"=as.numeric(names(wk0.p1.prop)), "Abundance"=unname(wk0.p1.prop), "Perc"=100*unname(wk0.p1.prop)/sum(unname(wk0.p1.prop)), "Pheno"=rep("Week0",length(wk0.p1.prop)))
      wk12.p1.prop<-table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)[,"Week12"];  wk12.p1.prop.df<-data.frame("Length"=as.numeric(names(wk12.p1.prop)), "Abundance"=unname(wk12.p1.prop), "Perc"=100*unname(wk12.p1.prop)/sum(unname(wk12.p1.prop)), "Pheno"=rep("Week12",length(wk12.p1.prop)))
      wk0wk12.p1.prop.df<-rbind(wk0.p1.prop.df,wk12.p1.prop.df)
      
      wk0.p2.prop<-table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)[,"Week0"];    wk0.p2.prop.df<-data.frame("Length"=as.numeric(names(wk0.p2.prop)), "Abundance"=unname(wk0.p2.prop), "Perc"=100*unname(wk0.p2.prop)/sum(unname(wk0.p2.prop)), "Pheno"=rep("Week0",length(wk0.p2.prop)))
      wk12.p2.prop<-table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)[,"Week12"];  wk12.p2.prop.df<-data.frame("Length"=as.numeric(names(wk12.p2.prop)), "Abundance"=unname(wk12.p2.prop), "Perc"=100*unname(wk12.p2.prop)/sum(unname(wk12.p2.prop)), "Pheno"=rep("Week12",length(wk12.p2.prop)))
      wk0wk12.p2.prop.df<-rbind(wk0.p2.prop.df,wk12.p2.prop.df)   
      
      # Statistical Tests
      
      pval.wilc.p1.w<-wilcox.test(sort(wk0.p1.vect),sort(wk12.p1.vect))$p.value
      pval.ks.p1.w<-ks.test(sort(wk0.p1.vect),sort(wk12.p1.vect))$p.value
      pval.wilc.p2.w<-wilcox.test(sort(wk0.p2.vect),sort(wk12.p2.vect))$p.value
      pval.ks.p2.w<-ks.test(sort(wk0.p2.vect),sort(wk12.p2.vect))$p.value    
      
      df.fisher.p1.w<-data.frame("Length"=rownames(table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)),
                                 "Week0.Yes"=as.numeric(table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)[,"Week0"]),
                                 "Week0.No"=as.numeric(sum(table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)[,"Week0"])-table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)[,"Week0"]),
                                 "Week12.Yes"=as.numeric(table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)[,"Week12"]),
                                 "Week12.No"=as.numeric(sum(table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)[,"Week12"])-table(wk0wk12.p1.ab.df$Length,wk0wk12.p1.ab.df$Pheno)[,"Week12"]), stringsAsFactors=F)
      
      df.fisher.p1.w$Pval<-rep(0,nrow(df.fisher.p1.w))
      for (i in 1:nrow(df.fisher.p1.w)) {
        p1.yes<-df.fisher.p1.w[i,2];   p1.no<-df.fisher.p1.w[i,3];   p2.yes<-df.fisher.p1.w[i,4];   p2.no<-df.fisher.p1.w[i,5];
        mat2analyze<-matrix(c(p1.yes,p1.no,p2.yes,p2.no), nrow=2)
        pval<-fisher.test(mat2analyze)$p.value
        df.fisher.p1.w$Pval[i]<-pval
      }
      
      fisher.sig.p1<-df.fisher.p1.w[df.fisher.p1.w$Pval<0.05,];   fisher.sig.p1$lab<-rep("*",nrow(fisher.sig.p1));   fisher.sig.p1$lab[fisher.sig.p1$Pval<0.005]<-"**";
      
      summ.wk0.p1<-c(summary(wk0.p1.vect), sd(wk0.p1.vect))
      names(summ.wk0.p1)[7]="Std."
      summ.wk12.p1<-c(summary(wk12.p1.vect), sd(wk12.p1.vect))
      names(summ.wk12.p1)[7]="Std."
      
      summ.wk0.p1.plot<-paste0("BASELINE\n",names(summ.wk0.p1)[1]," = ", summ.wk0.p1[[1]],"\n",names(summ.wk0.p1)[2]," = ", summ.wk0.p1[[2]],"\n",
                               names(summ.wk0.p1)[3]," = ", round(summ.wk0.p1[[3]],2),"\n",names(summ.wk0.p1)[4]," = ", round(summ.wk0.p1[[4]],2),"\n",
                               names(summ.wk0.p1)[5]," = ", summ.wk0.p1[[5]],"\n",names(summ.wk0.p1)[6]," = ", summ.wk0.p1[[6]],"\n", names(summ.wk0.p1)[7]," = ", round(summ.wk0.p1[[7]],2),"\n")
      summ.wk12.p1.plot<-paste0("WEEK.12\n",names(summ.wk12.p1)[1]," = ", summ.wk12.p1[[1]],"\n",names(summ.wk12.p1)[2]," = ", summ.wk12.p1[[2]],"\n",
                                names(summ.wk12.p1)[3]," = ", round(summ.wk12.p1[[3]],2),"\n",names(summ.wk12.p1)[4]," = ", round(summ.wk12.p1[[4]],2),"\n",
                                names(summ.wk12.p1)[5]," = ", summ.wk12.p1[[5]],"\n",names(summ.wk12.p1)[6]," = ", summ.wk12.p1[[6]],"\n", names(summ.wk12.p1)[7]," = ", round(summ.wk12.p1[[7]],2),"\n")
      
      
      df.fisher.p2.w<-data.frame("Length"=rownames(table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)),
                                 "Week0.Yes"=as.numeric(table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)[,"Week0"]),
                                 "Week0.No"=as.numeric(sum(table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)[,"Week0"])-table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)[,"Week0"]),
                                 "Week12.Yes"=as.numeric(table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)[,"Week12"]),
                                 "Week12.No"=as.numeric(sum(table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)[,"Week12"])-table(wk0wk12.p2.ab.df$Length,wk0wk12.p2.ab.df$Pheno)[,"Week12"]), stringsAsFactors=F)
      
      df.fisher.p2.w$Pval<-rep(0,nrow(df.fisher.p2.w))
      for (i in 1:nrow(df.fisher.p2.w)) {
        p1.yes<-df.fisher.p2.w[i,2];   p1.no<-df.fisher.p2.w[i,3];   p2.yes<-df.fisher.p2.w[i,4];   p2.no<-df.fisher.p2.w[i,5];
        mat2analyze<-matrix(c(p1.yes,p1.no,p2.yes,p2.no), nrow=2)
        pval<-fisher.test(mat2analyze)$p.value
        df.fisher.p2.w$Pval[i]<-pval
      }
      
      fisher.sig.p2<-df.fisher.p2.w[df.fisher.p2.w$Pval<0.05,];   fisher.sig.p2$lab<-rep("*",nrow(fisher.sig.p2));   fisher.sig.p2$lab[fisher.sig.p2$Pval<0.005]<-"**";
      
      summ.wk0.p2<-c(summary(wk0.p2.vect), sd(wk0.p2.vect))
      names(summ.wk0.p2)[7]="Std."
      summ.wk12.p2<-c(summary(wk12.p2.vect), sd(wk12.p2.vect))
      names(summ.wk12.p2)[7]="Std."
      
      summ.wk0.p2.plot<-paste0("BASELINE\n",names(summ.wk0.p2)[1]," = ", summ.wk0.p2[[1]],"\n",names(summ.wk0.p2)[2]," = ", summ.wk0.p2[[2]],"\n",
                               names(summ.wk0.p2)[3]," = ", round(summ.wk0.p2[[3]],2),"\n",names(summ.wk0.p2)[4]," = ", round(summ.wk0.p2[[4]],2),"\n",
                               names(summ.wk0.p2)[5]," = ", summ.wk0.p2[[5]],"\n",names(summ.wk0.p2)[6]," = ", summ.wk0.p2[[6]],"\n", names(summ.wk0.p2)[7]," = ", round(summ.wk0.p2[[7]],2),"\n")
      summ.wk12.p2.plot<-paste0("WEEK.12\n",names(summ.wk12.p2)[1]," = ", summ.wk12.p2[[1]],"\n",names(summ.wk12.p2)[2]," = ", summ.wk12.p2[[2]],"\n",
                                names(summ.wk12.p2)[3]," = ", round(summ.wk12.p2[[3]],2),"\n",names(summ.wk12.p2)[4]," = ", round(summ.wk12.p2[[4]],2),"\n",
                                names(summ.wk12.p2)[5]," = ", summ.wk12.p2[[5]],"\n",names(summ.wk12.p2)[6]," = ", summ.wk12.p2[[6]],"\n", names(summ.wk12.p2)[7]," = ", round(summ.wk12.p2[[7]],2),"\n")
      
      # Plots
      if (pval.wilc.p1.w>0) {pval.text.p1<-paste0("(P=",formatC(as.numeric(pval.wilc.p1.w), format="e", digits=2),")")}
      if (pval.wilc.p1.w==0) {pval.text.p1<-"(P<1e-300)"}
      plot.p1.w<-ggplot(wk0wk12.p1.prop.df, aes(x=Length, y=Perc, fill=Pheno), color="black") + geom_bar(stat="identity", position=position_dodge()) +
        theme_classic() + scale_fill_manual(values=c("red", "blue")) + theme_classic() +
        labs(title=paste0("\n",c," - ",pheno1,":\nCDR3aa length distribution ",pval.text.p1,"\n"),
             y="\nFrequency (%)\n", x="\nCDR3 amino acid length (N)\n") +
        theme(legend.title=element_blank(), legend.position=c(0.9,0.9),axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5)) +
        scale_x_continuous(labels=as.character(wk0wk12.p1.prop.df$Length),breaks=wk0wk12.p1.prop.df$Length) + 
        annotate(geom="text", x=as.numeric(rownames(fisher.sig.p1)), y=rep(max(wk0wk12.p1.prop.df$Perc)*1.1,nrow(fisher.sig.p1)), label=fisher.sig.p1$lab, color="black", size=3) +
        annotate(geom="text", x=max(wk0wk12.p1.prop.df$Length)-2, y=max(wk0wk12.p1.prop.df$Perc)*0.67, label=summ.wk0.p1.plot, size=3, col="red") +
        annotate(geom="text", x=max(wk0wk12.p1.prop.df$Length)-2, y=max(wk0wk12.p1.prop.df$Perc)*0.32, label=summ.wk12.p1.plot, size=3, col="blue")
      
      if (pval.wilc.p2.w>0) {pval.text.p2<-paste0("(P=",formatC(as.numeric(pval.wilc.p2.w), format="e", digits=2),")")}
      if (pval.wilc.p2.w==0) {pval.text.p2<-"(P<1e-300)"}
      plot.p2.w<-ggplot(wk0wk12.p2.prop.df, aes(x=Length, y=Perc, fill=Pheno), color="black") + geom_bar(stat="identity", position=position_dodge()) +
        theme_classic() + scale_fill_manual(values=c("red", "blue")) + theme_classic() +
        labs(title=paste0("\n",c," - ",pheno2,":\nCDR3aa length distribution ",pval.text.p2,"\n"),
             y="\nFrequency (%)\n", x="\nCDR3 amino acid length (N)\n") +
        theme(legend.title=element_blank(), legend.position=c(0.9,0.9), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5)) +
        scale_x_continuous(labels=as.character(wk0wk12.p2.prop.df$Length),breaks=wk0wk12.p2.prop.df$Length) + 
        annotate(geom="text", x=as.numeric(rownames(fisher.sig.p2)), y=rep(max(wk0wk12.p2.prop.df$Perc)*1.1,nrow(fisher.sig.p2)), label=fisher.sig.p2$lab, color="black", size=3) +
        annotate(geom="text", x=max(wk0wk12.p2.prop.df$Length)-2, y=max(wk0wk12.p2.prop.df$Perc)*0.67, label=summ.wk0.p2.plot, size=3, col="red") +
        annotate(geom="text", x=max(wk0wk12.p2.prop.df$Length)-2, y=max(wk0wk12.p2.prop.df$Perc)*0.32, label=summ.wk12.p2.plot, size=3, col="blue")    
      
      ctrl.vect<-rep(clones.ctrl.df$Length,clones.ctrl.df$Abundance);     ctrl.ab.df<-data.frame("Length"=ctrl.vect,"Pheno"=rep("CTRL",length(ctrl.vect)))
      wk0wk12ctrl.p1.ab.df<-rbind(wk0wk12.p1.ab.df,ctrl.ab.df)
      wk0wk12ctrl.p2.ab.df<-rbind(wk0wk12.p2.ab.df,ctrl.ab.df)
      
      plot.cum.p1.w<-ggplot(wk0wk12ctrl.p1.ab.df, aes(Length, color=Pheno)) + stat_ecdf(geom="step") + theme_classic() + scale_color_manual(values=c("red", "blue", "chartreuse4")) +
        labs(title=paste0("\n",c," - ",pheno1,":\nEmpirical Cumulative Density Function\n"), y="\nCumulative Probability\n", x="\nCDR3 amino acid length (N)\n") +
        theme(legend.title=element_blank(), legend.position=c(0.9,0.2),axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5))
      
      plot.cum.p2.w<-ggplot(wk0wk12ctrl.p2.ab.df, aes(Length, color=Pheno)) + stat_ecdf(geom="step") + theme_classic() + scale_color_manual(values=c("red", "blue", "chartreuse4")) +
        labs(title=paste0("\n",c," - ",pheno2,":\nEmpirical Cumulative Density Function\n"), y="\nCumulative Probability\n", x="\nCDR3 amino acid length (N)\n") +
        theme(legend.title=element_blank(), legend.position=c(0.9,0.2),axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5))    
      
      plot.qqplot.p1.w<-ggplot(mapping=aes(x=quantile(wk0wk12.p1.ab.df[wk0wk12.p1.ab.df$Pheno=="Week0","Length"], seq(0,1,0.01)), y=quantile(wk0wk12.p1.ab.df[wk0wk12.p1.ab.df$Pheno=="Week12","Length"], seq(0,1,0.01)))) + 
        geom_point() + geom_abline(aes(slope=1, intercept=0), linetype=2)  + theme_classic() +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5)) +
        labs(x=paste("\nWeek.0 Quantiles\n"), y=paste("\nWeek.12 Quantiles\n"), title=paste0("\n",c," - ",pheno1,":\nEmpirical QQ-plot\n"))
      
      plot.qqplot.p2.w<-ggplot(mapping=aes(x=quantile(wk0wk12.p2.ab.df[wk0wk12.p2.ab.df$Pheno=="Week0","Length"], seq(0,1,0.01)), y=quantile(wk0wk12.p2.ab.df[wk0wk12.p2.ab.df$Pheno=="Week12","Length"], seq(0,1,0.01)))) + 
        geom_point() + geom_abline(aes(slope=1, intercept=0), linetype=2) + theme_classic() +
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), plot.title=element_text(size=14, face="bold", hjust=0.5)) +
        labs(x=paste("\nWeek.0 Quantiles\n"), y=paste("\nWeek.12 Quantiles\n"), title=paste0("\n",c," - ",pheno2,":\nEmpirical QQ-plot\n"))
      
      return(list(plot.p1.w, plot.p2.w, plot.cum.p1.w, plot.cum.p2.w, plot.qqplot.p1.w, plot.qqplot.p2.w))
    }
      
    pdf("Output_Data/CDR3length_Comparisons.pdf", width=30, height=52)
    phenos<-c("IMID", "RespGM", "ActB", "ACPA", "RF")
    chains<-c("TRA", "TRB", "TRD", "TRG", "IGH", "IGL", "IGK")
    for (ph in phenos) {
        pheno.data<-get.pheno.data(ph)
        pheno.out<-pheno.data[1]; pheno.name<-pheno.data[2]; pheno1<-pheno.data[3]; pheno2<-pheno.data[4]; week<-pheno.data[5]; analysis<-pheno.data[6]
        print(paste0(ph," - TRA"))
        tra.plot<-make.chain.plots.baseline("TRA",ph,pheno1,pheno2,week,analysis)
        trb.plot<-make.chain.plots.baseline("TRB",ph,pheno1,pheno2,week,analysis)
        trd.plot<-make.chain.plots.baseline("TRD",ph,pheno1,pheno2,week,analysis)
        trg.plot<-make.chain.plots.baseline("TRG",ph,pheno1,pheno2,week,analysis)
        igh.plot<-make.chain.plots.baseline("IGH",ph,pheno1,pheno2,week,analysis)
        igl.plot<-make.chain.plots.baseline("IGL",ph,pheno1,pheno2,week,analysis)
        igk.plot<-make.chain.plots.baseline("IGK",ph,pheno1,pheno2,week,analysis)
        if (ph=="IMID") {title="\nCASE-CONTROL ANALYSIS\n"}
        if (ph=="RespGM") {title="\nCASE-CASE ANALYSIS: RESPONSE TO TNFi THERAPY\n"}  
        if (ph=="ActB") {title="\nCASE-CASE ANALYSIS: DISEASE ACTIVITY\n"}  
        if (ph=="ACPA") {title="\nCASE-CASE ANALYSIS: ACPA PHENOTYPE\n"} 
        if (ph=="RF") {title="\nCASE-CASE ANALYSIS: RF PHENOTYPE\n"} 
        
        grid.arrange(tra.plot[[1]], tra.plot[[2]], tra.plot[[3]],
                     trb.plot[[1]], trb.plot[[2]], trb.plot[[3]],
                     trd.plot[[1]], trd.plot[[2]], trd.plot[[3]],
                     trg.plot[[1]], trg.plot[[2]], trg.plot[[3]],
                     igh.plot[[1]], igh.plot[[2]], igh.plot[[3]],
                     igl.plot[[1]], igl.plot[[2]], igl.plot[[3]],
                     igk.plot[[1]], igk.plot[[2]], igk.plot[[3]], ncol=3,
                     top=textGrob(title,gp=gpar(fontsize=22,fontface="bold")))
    }
    
    phenos<-c("IMID")
    chains<-c("TRA", "TRB", "TRD", "TRG", "IGH", "IGL", "IGK")
    for (ph in phenos) {
      pheno.data<-get.pheno.data(ph)
      pheno.out<-pheno.data[1]; pheno.name<-pheno.data[2]; pheno1<-pheno.data[3]; pheno2<-pheno.data[4]; week<-pheno.data[5]; analysis<-pheno.data[6]
    
      print(paste0(ph," - TRA"))
      tra.plot<-make.chain.plots.longAll("TRA",ph,pheno1,pheno2,week,analysis)
      trb.plot<-make.chain.plots.longAll("TRB",ph,pheno1,pheno2,week,analysis)
      trd.plot<-make.chain.plots.longAll("TRD",ph,pheno1,pheno2,week,analysis)
      trg.plot<-make.chain.plots.longAll("TRG",ph,pheno1,pheno2,week,analysis)
      igh.plot<-make.chain.plots.longAll("IGH",ph,pheno1,pheno2,week,analysis)
      igl.plot<-make.chain.plots.longAll("IGL",ph,pheno1,pheno2,week,analysis)
      igk.plot<-make.chain.plots.longAll("IGK",ph,pheno1,pheno2,week,analysis)
      title="\nLONGITUDINAL ANALYSIS: BASELINE vs. TNFi THERAPY\n"
        
      grid.arrange(tra.plot[[1]], tra.plot[[2]], tra.plot[[3]],
                     trb.plot[[1]], trb.plot[[2]], trb.plot[[3]],
                     trd.plot[[1]], trd.plot[[2]], trd.plot[[3]],
                     trg.plot[[1]], trg.plot[[2]], trg.plot[[3]],
                     igh.plot[[1]], igh.plot[[2]], igh.plot[[3]],
                     igl.plot[[1]], igl.plot[[2]], igl.plot[[3]],
                     igk.plot[[1]], igk.plot[[2]], igk.plot[[3]], ncol=3,
                     top=textGrob(title,gp=gpar(fontsize=22,fontface="bold")))
    }
    
    phenos<-c("RespGM")
    chains<-c("TRA", "TRB", "TRD", "TRG", "IGH", "IGL", "IGK")
    for (ph in phenos) {
      pheno.data<-get.pheno.data(ph)
      pheno.out<-pheno.data[1]; pheno.name<-pheno.data[2]; pheno1<-pheno.data[3]; pheno2<-pheno.data[4]; week<-pheno.data[5]; analysis<-pheno.data[6]
    
      print(paste0(ph," - TRA"))
      tra.plot<-make.chain.plots.longByResp("TRA",ph,pheno1,pheno2,week,analysis)
      trb.plot<-make.chain.plots.longByResp("TRB",ph,pheno1,pheno2,week,analysis)
      trd.plot<-make.chain.plots.longByResp("TRD",ph,pheno1,pheno2,week,analysis)
      trg.plot<-make.chain.plots.longByResp("TRG",ph,pheno1,pheno2,week,analysis)
      igh.plot<-make.chain.plots.longByResp("IGH",ph,pheno1,pheno2,week,analysis)
      igl.plot<-make.chain.plots.longByResp("IGL",ph,pheno1,pheno2,week,analysis)
      igk.plot<-make.chain.plots.longByResp("IGK",ph,pheno1,pheno2,week,analysis)
      if (ph=="RespGM") {title="\nLONGITUDINAL ANALYSIS STRATIFIED BY CLINICAL RESPONSE TO TNFi THERAPY\n"}
    
      grid.arrange(tra.plot[[1]], tra.plot[[2]], tra.plot[[3]], tra.plot[[4]], tra.plot[[5]], tra.plot[[6]],
                   trb.plot[[1]], trb.plot[[2]], trb.plot[[3]], trb.plot[[4]], trb.plot[[5]], trb.plot[[6]],
                   trd.plot[[1]], trd.plot[[2]], trd.plot[[3]], trd.plot[[4]], trd.plot[[5]], trd.plot[[6]],
                   trg.plot[[1]], trg.plot[[2]], trg.plot[[3]], trg.plot[[4]], trg.plot[[5]], trg.plot[[6]],
                   igh.plot[[1]], igh.plot[[2]], igh.plot[[3]], igh.plot[[4]], igh.plot[[5]], igh.plot[[6]],
                   igl.plot[[1]], igl.plot[[2]], igl.plot[[3]], igl.plot[[4]], igl.plot[[5]], igl.plot[[6]],
                   igk.plot[[1]], igk.plot[[2]], igk.plot[[3]], igk.plot[[4]], igk.plot[[5]], igk.plot[[6]], ncol=6,
                   top=textGrob(title,gp=gpar(fontsize=22,fontface="bold")))
    }
    dev.off()
    

