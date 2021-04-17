#' Measure the correlation of candidate nodes and known nodes
#'
#' @title RWRs.CN function
#' @param nodes data.frame, two columns:object,object label(0/1),0 represents candidate nodes and 1 represents known nodes(seeds).
#' @param edges data.frame, two columns:interactorA,interactorB,A and B have a certain interaction relationship.
#' @param r default 0.8,Restart probability of random walk,interval is (0,1),r is preferably above 0.5.
#' @param core default 10,Top percentile of ranking for candidates,interval is (0,100].Top core percent nodes and their edges will be extracted as core network.
#' @import igraph dplyr dnet
#' @export RWRs.CN
#' @author Fan Wang

RWRs.CN <- function(nodes,
                 edges,
                 r=0.8,
                 core=10) {
  #check the format of parameters
  if(length(which(nodes[,2]==1))*length(which(nodes[,2]==0)) <1)
  {
    stop("First param 'nodes' input error!
         Please input dataframe with two columns!
         for example: TP53,1
                      EGFR,0
                      KRAS,0")
  }
  else if(!class(edges) == "data.frame" ||
          ncol(edges) <2)
  {
    stop("Second param 'edges' input error!
         Please input dataframe with two columns!
         for example: TP53,EGFR
                      TP53,KRAS
                      EGFR,KRAS")
  }
  else if(r >=1 || r<=0)
  {
    stop("Third param 'r' input error!
         'r' must be in (0,1)!")
  }
  else if(core >100 || core<=0)
  {
    stop("Fourth param 'core' input error!
         'core' must be in (0,100]!")
  }
  options(scipen=200)
  #Set seeds list
  nodes <-nodes[order(nodes[,2], decreasing= T), ]
  seedcount=as.numeric(length(which(nodes[,2]==1)))
  topper=ceiling((nrow(nodes)-seedcount)*core*0.01)
  g <- graph_from_data_frame(d = edges,vertices = nodes,directed = F)
  Seeds=as.data.frame(matrix(nc=1,nr=vcount(g)))
  Seeds$V1=0
  rownames(Seeds) <- nodes[,1]
  Seeds[1:seedcount,1]=1
  result=matrix(nrow = vcount(g), ncol = seedcount)
  result<-as.data.frame(result)
  #Random walk with restart from each seed
  for(i in 1:seedcount)
  {
    Seeds[1:seedcount,1]=0
    Seeds[i,1]=1
    PTmatrix<- dRWR(g=g, normalise="laplacian",
                    setSeeds=Seeds, restart=r, parallel=FALSE)
    result[,i]=PTmatrix[,1]
  }
  colnames(result)=nodes[1:seedcount,1]
  rownames(result)=nodes[,1]
  result$RWRs=apply(result[,c(rep(1:seedcount,1))],1,sum,na.rm=T)
  result <- mutate(result,rn=row.names(result))
  temp<- result
  result2=as.data.frame(temp$rn)
  result2$RWRs=temp$RWRs
  #Export random walk score of all nodes
  write.table(result2,paste('RWR_',r,'_score.txt',sep=''),quote=F,
              sep="\t",row.names=F,col.names=c("name","RWRscore"))
  result$label=0
  result$label[1:seedcount]=1
  #Random walk result sorting
  result<- arrange(result,desc(RWRs))
  result_noseed<- result[result$label==0,]
  result_noseed<-as.data.frame(result_noseed$rn[1:topper])
  result_seed<- as.data.frame(result[result$label==1,seedcount+2])
  result_topten<-as.data.frame(union(result_seed[,1],
                                     result_noseed[,1]))
  colnames(result_topten)="V1"
  colnames(edges)=c("V1","V2")
  #Extract core network
  temp1<-merge(result_topten,edges,all.x=TRUE,sort=TRUE)
  temp2=edges
  temp2[,1]<-edges[,2]
  temp2[,2]<-edges[,1]
  temp3<-merge(result_topten,temp2,all.x=TRUE,sort=TRUE)
  temp4=temp3
  temp3[,1]<-temp4[,2]
  temp3[,2]<-temp4[,1]
  temp5<-rbind(temp1,temp3)
  temp6<-temp5[complete.cases(temp5),]
  coreNet<-temp6[!duplicated(temp6, fromLast=TRUE), ]
  #Overview of core network
  g_core<-graph_from_data_frame(coreNet, directed=FALSE)
  plot(g_core)
  #Export core network
  write.table(coreNet,paste('top',core,'%_CoreNet.txt',sep=''),
              quote=F,sep="\t",row.names=F,col.names=F)
}
