#' Node and edge Prioritization-based Community Analysis framework
#'
#' @title ne.PCA function
#' @param nodes data.frame, two columns:object,label(0:candidates/1:seeds).
#' @param edges data.frame, two columns:gene symbol A,gene symbol B.
#' @param r default 0.8,Restart probability of random walk,interval is (0,1).
#' @param core default 10,Top percentile of ranking for object,interval is (0,100].
#' @import igraph dplyr dnet GOSemSim reshape2 ggplot2
#' @export ne.PCA
#' @author Fan Wang <justin15751821901@gmail.com>

ne.PCA <- function(nodes,
                   edges,
                   r=0.8,
                   core=10) {
  #Get the core network
  RWRs.CN(nodes,edges,r,core)
  coreNet<-read.table(paste('top',core,'%_CoreNet.txt',sep=''),header=F)
  colnames(coreNet)=c("name1","name2")
  #Check input format
  if(!ncol(edges) == 2)
  {
    stop("Param 'edges' input error!
         Only two columns needed!
         for example: TP53,EGFR
                      TP53,KRAS
                      EGFR,KRAS")
  }
  #Get the edge weight TFC
  TFCs(edges)
  ecount=nrow(edges)
  weight<-read.table(paste('TFC_',ecount,'interactions.txt',sep=''),header=T)
  WCN=merge(coreNet,weight,by=c("name1","name2"))
  #Export weight core network(WCN)
  write.table(WCN,"WCN.txt",sep="\t",quote=F,
              row.names=F)
  #Identify robust modules in WCN
  MODULEs(coreNet)
}
