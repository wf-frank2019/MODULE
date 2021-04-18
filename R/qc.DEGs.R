#' Data pretreatment and differentially expressed genes identification
#'
#' @title qc.DEGs function
#' @param expr data.frame, Expression profile,the row name is the gene name,the column name is the sample name.
#' @param group data.frame, two columns:sample name, sample label,for example:GSM2718949,control.
#' @param lowEXP default 1,Low quality value cutoff,lowEXP must > 0.The expression value less than 'lowEXP' will be seen as missing values.
#' @param filter default 0.3,The proportion of the number of missing values in all samples in one row,the interval is (0,1).A row will be removed if the proportion of missing values in this row is greater than 'filter'.
#' @param fc default 2,Fold Change cutoff,fc must be greater than 1,generally fc is 1.5 or 2 or greater.
#' @param p default 0.05,Adjust p value cutoff,the interval is (0,1),generally p is 0.05 or 0.01 or 0.001.
#' @param func default FALSE,if 'func=TRUE',GO and KEGG pathways Enrichment analysis for DEGs(human genes) will be performed.
#' @import limma Hmisc ggplot2 clusterProfiler org.Hs.eg.db
#' @export qc.DEGs
#' @author Fan Wang

qc.DEGs<- function(expr,
                   group,
                   lowEXP=1,
                   filter=0.3,
                   fc=2,
                   p=0.05,
                   func=FALSE){
  #check the format of parameters
  if(!class(expr)=="data.frame" || !class(expr[1,1]) == "numeric")
  {
    stop("First param 'Expression profile' input error!
         Please input dataframe with gene/protein.. expression value!")
  }
  else if(!class(group)=="data.frame" || ncol(group) <2)
  {
    stop("Second param 'group' input error!
         Please input dataframe with two columns(sample name;sample label)!")
  }
  else if(lowEXP <=0)
  {
    stop("Third param 'lowEXP' input error!
         Please input a number > 0 !")
  }
  else if(filter <=0 || filter>=1)
  {
    stop("Fourth param 'filter' input error!
         Please input a number in (0,1) !")
  }
  else if(fc <=1)
  {
    stop("Fifth param 'fc' input error!
         Please input a number > 1 !")
  }
  else if(p <=0 || p>=1)
  {
    stop("Sixth param 'p' input error!
         Please input a number in (0,1) !")
  }
  time1<-Sys.time()
  expr.q <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  contrName=(deparse(substitute(group)))
  if(lowEXP>=expr.q[2])
  {
    #Control low quality value cutoff to be less than the first quartile
    stop(paste(c("the third param 'lowEXP' cannot > ",expr.q[2])))
  }
  #Check whether the expression profile has been logarithmic
  Log.test <- (expr.q[5] > 100) || (expr.q[6]-expr.q[1] > 50 && expr.q[2] > 0) ||
    (expr.q[2] > 0 && expr.q[2] < 1 && expr.q[4] > 1 && expr.q[4] < 2)
  if (Log.test){
    #Screening low-abundance data
    expr[expr<lowEXP]<-NA
    expr <- log2(expr)
  }
  else{
    expr[expr<log2(lowEXP)]<-NA
  }
  #Imputation for missing values
  filter.count<-length(which(apply(expr, 1, function(temp) sum(is.na(temp)))
                             > ncol(expr)*filter))
  if(filter.count>nrow(expr)*(1-filter)){
    message("Tips: too much objects show low expression abundance", appendLF = T)
  }
  if(filter.count > 0){
    expr <- expr[-which(apply(expr, 1, function(temp) sum(is.na(temp)))
                        > ncol(expr)*filter),]
  }
  expr=as.data.frame(t(apply(expr, 1, function(x){impute(x,mean)})))
  message(paste(c("row counts after quality control:",
                  as.character(nrow(expr)), collapse = "")))
  #Export data after pretreatment
  write.table(expr,"Pre_expr.txt",sep="\t",quote=F,col.names = NA)
  message("Data pretreatment is complete !",collapse = "")
  #Difference analysis in limma
  sample_label <- model.matrix(~0+factor(group[,2]))
  colnames(sample_label)=levels(factor(group[,2]))
  rownames(sample_label)=colnames(expr)
  if(!contrName=="group")
  {
    message(" Error! Please rename second input parameter as 'group'",collapse = "")
  }
  de.matrix<-makeContrasts(paste0(unique(group[,2]),collapse = "-"),levels = sample_label)
  de.model <- lmFit(expr,sample_label)
  de.model2 <- contrasts.fit(de.model,de.matrix)
  de.model2 <- eBayes(de.model2)
  de.Output = topTable(de.model2, coef=1, n=Inf,adjust="BH")
  de.all = na.omit(de.Output)
  #Export genes with different significance level
  write.table(de.all, "DEGs_all.txt",sep="\t",quote=F,col.names=NA)
  de.temp = de.all[de.all$adj.P.Val < p & abs(de.all$logFC)>log2(fc),]
  de=de.temp[,1:2]
  de[,2]=de.temp$adj.P.Val
  colnames(de)=c("logFC","adj.P.Val")
  write.table(de, paste('DEGs_p',p,'_FC',fc,'.txt',sep='')
              ,sep="\t",quote=F,col.names=NA)
  #volcano plot of DEGs
  de.all$label[(de.all$adj.P.Val > p|de.all$adj.P.Val=="NA")|
                 (de.all$logFC <= log2(fc))& de.all$logFC >= -log2(fc)] <- "no significant"
  de.all$label[de.all$adj.P.Val <= p & de.all$logFC > log2(fc)] <- "up regulated"
  de.all$label[de.all$adj.P.Val <= p & de.all$logFC < -log2(fc)] <- "down regulated"
  x_lim <- max(de.all$logFC,-de.all$logFC)
  de.all$label<-factor(de.all$label,levels = c("up regulated","down regulated","no significant"))
  p_DEGs = ggplot(de.all,aes(logFC,-1*log10(adj.P.Val),color=label)) +
    geom_point(size=1.2) + xlim(-x_lim,x_lim) +
    labs(x="log2(Fold Change)",y="-log10(adjust.p)",title="DEGs in experimental / control")+
    scale_color_manual(values =c("red","dark green","grey"))+
    geom_hline(aes(yintercept=-1*log10(p)),colour="black", linetype="dashed") +
    geom_vline(xintercept=c(-log2(fc),log2(fc)),colour="black", linetype="dashed")+
    theme(legend.text=element_text(size=14),
          legend.title=element_text(size=15),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+theme(legend.key = element_blank())+
    theme(legend.background = element_blank())+theme(plot.title =
                                                       element_text(hjust = 0.5,size = 16, face = "bold"),axis.text=
                                                       element_text(size=12,face = "bold"),axis.title.x=
                                                       element_text(size=14),axis.title.y=element_text(size=14))
  #GO and KEGG enrichment analysis for human differentially expressed genes
  if(func)
  {
    symbol=as.data.frame(rownames(de))
    colnames(symbol)="SYMBOL"
    symbol_func <- bitr(symbol[,1], fromType = "SYMBOL",
                        toType = "ENTREZID",OrgDb = org.Hs.eg.db)
    #GO enrichment analysis with p value <0.05
    MF <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = symbol_func$ENTREZID,
                   pvalueCutoff = 0.05,
                   ont = "MF",
                   readable=TRUE)
    CC <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = symbol_func$ENTREZID,
                   pvalueCutoff = 0.05,
                   ont = "CC",
                   readable=TRUE)
    BP <- enrichGO(OrgDb="org.Hs.eg.db",
                   gene = symbol_func$ENTREZID,
                   pvalueCutoff = 0.05,
                   ont = "BP",
                   readable=TRUE)
    def<-rbind(as.data.frame(MF)[1:30, ],as.data.frame(CC)[1:30, ],
               as.data.frame(BP)[1:30, ])
    #Export top 30 GO terms
    write.table(def,"GO_top30.txt",sep="\t",quote=F,row.names =F)
    #KEGG enrichment analysis with p value <0.05
    ekk<- enrichKEGG(gene = symbol_func$ENTREZID,
                     organism   = 'hsa',
                     pvalueCutoff = 0.05)
    ekk<-setReadable(ekk,OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    #Export all KEGG pathways
    write.table(ekk@result,"KEGG_pathways.txt",row.names =F,quote=F,sep="\t")
  }
  #Time-consuming feedback
  time2<-Sys.time()
  message("", appendLF = T)
  message(paste(c("You start at ", as.character(time1)), collapse = ""), appendLF = T)
  message(paste(c("You finish at ", as.character(time2)), collapse = ""), appendLF = T)
  message("", appendLF = T)
  return(p_DEGs)
  }
