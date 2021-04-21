#' Measure the role of interactions as an information medium
#'
#' @title TFCs function
#' @param edges data.frame, two columns,for exmaple:gene symbol A, gene symbol B.
#' @import igraph GOSemSim
#' @export TFCs
#' @author Fan Wang

TFCs <- function(edges) {
  #check the format of parameters
  if(!class(edges) == "data.frame")
  {
    stop("Param 'edges' input error!
         Please input the dataframe!")
  }
  else if(ncol(edges) <2)
  {
    stop("Param 'edges' input error!
         Please input two columns!
         For example:protein-protein interactions")
  }
  options(scipen=200)
  time1<-Sys.time()
  ecount<-nrow(edges)
  #Create the network
  g = graph_from_data_frame(edges, directed=FALSE)
  nodes=as.data.frame(rbind(as.matrix(edges[,1]),as.matrix(edges[,2])))
  nodes=unique(nodes)
  #Calculate the edge betweenness
  bet<-edge_betweenness(g, e = E(g), directed = F)
  bet<-cbind(edges,bet)
  colnames(bet)=c("name1","name2","bet")
  #Calculate functional similarity
  hsGO <- godata('org.Hs.eg.db', keytype = "SYMBOL",
                  ont="BP", computeIC=FALSE)
  funcsim<-mgeneSim(nodes[,1], semData=hsGO, measure="Wang",
                  combine="BMA", verbose=FALSE)
  upsim <- upper.tri(funcsim)
  mergesim<-data.frame(row = rownames(funcsim)[row(funcsim)[upsim]],
                    column = rownames(funcsim)[col(funcsim)[upsim]],
                    cor =(funcsim)[upsim] )
  temp=mergesim
  temp[,1]=mergesim[,2]
  temp[,2]=mergesim[,1]
  finalsim<-rbind(mergesim,temp)
  colnames(finalsim)=c("name1","name2","cor")
  colnames(edges)=c("name1","name2")
  GO<-merge(finalsim,edges,by=c("name1","name2"),sort=F)
  #TFC of edges without functional similarity to be zero
  nasim=dplyr::anti_join(edges[,1:2], GO[,1:2],by = c("name1", "name2"))
  nasim$cor=0
  GO=rbind(GO,nasim)
  #Obtain topological-functional connection score
  NET=merge(GO,bet,by=c("name1","name2"))
  scaleA<-NET$bet-min(NET$bet)
  scaleB<-max(NET$bet)-min(NET$bet)
  #Edge betweenness normalization
  NET$scaleTN <- scaleA/scaleB
  TF<-NET$scaleTN+NET$cor
  NET$TFC  <- 100*TF/(2-TF)
  TFC=NET[,1:2]
  TFC$TFCs=NET$TFC
  #Top TFC score interactions distribution
  #Zoom in to see details
  TFCplot<-TFC[order(TFC$TFCs,decreasing=T),]
  TFCplot$name=paste(TFCplot$name1, TFCplot$name2, sep="-")
  TFCplot=TFCplot[1:50,]
  maxcontr<-ceiling(max(TFCplot$TFCs))
  barplot(rev(TFCplot$TFCs),horiz=T,xlim=c(-35,maxcontr+20),
          axes=F,col = "orange")
  text(seq(from=0.7,length.out=nrow(TFCplot),by=1.2),x=-16,
       label=rev(TFCplot$name))
  axis(3,c(0,ceiling(0.25*max(TFCplot$TFCs)),ceiling(0.5*max(TFCplot$TFCs)),
           ceiling(0.75*max(TFCplot$TFCs)),maxcontr))
  title(xlab = "TFC score",ylab = "Interactions Name",line=0,col.lab = 9,
        cex.lab = 1.5,font.lab = 4)
  #Export TFC score
  write.table(TFC,paste('TFC_',ecount,'interactions.txt',sep=''),
              sep="\t",quote=F,row.names=F)
  #Time-consuming feedback
  time2<-Sys.time()
  message("", appendLF = T)
  message(paste(c("TFCs start at ", as.character(time1)), collapse = ""), appendLF = T)
  message(paste(c("TFCs finish at ", as.character(time2)), collapse = ""), appendLF = T)
  message("", appendLF = T)
}
