###############################
setwd("*:/***/")
library(igraph)
library(dnet)
library(dplyr)
library(GOSemSim)
library(igraph)
library(ggplot2)
library(reshape2)
options(scipen=200)
i=j=1
w=f=0
links <-read.csv("links.csv",header=T) #PPIs
nodes <- read.csv("nodes.csv",header=T) #name,label(0/1)
seedcount=as.numeric(length(which(nodes[,2]==1)))
topten=ceiling((nrow(nodes)-seedcount)/10)
g = graph_from_data_frame(links, directed=FALSE)
Seeds=as.data.frame(matrix(nc=1,nr=vcount(g)))
Seeds$V1=0
rownames(Seeds) <- nodes[,1]
Seeds[1:seedcount,1]=1
#rwr
result=matrix(nrow = vcount(g), ncol = seedcount)
result<-as.data.frame(result)
for(i in 1:seedcount)
{
  Seeds[1:seedcount,1]=0
  Seeds[i,1]=1
  PTmatrix<- dRWR(g=g, normalise="laplacian", setSeeds=Seeds, restart=0.8, parallel=FALSE) 
  result[,i]=PTmatrix[,1]
}
colnames(result)=nodes[1:seedcount,1]
rownames(result)=nodes[,1]
result$sum=apply(result[,c(rep(1:seedcount,1))],1,sum,na.rm=T)
result <- mutate(result,rn=row.names(result))
result$label=0;result$label[1:seedcount]=1
result<- arrange(result,desc(sum))  #write.table(result,"RWR.txt",quote=F,sep="\t",col.names=NA)
result_noseed<- result[result$label==0,]
result_noseed<-as.data.frame(result_noseed$rn[1:topten])
result_seed<- as.data.frame(result[result$label==1,seedcount+2])
result_topten<-as.data.frame(union(result_seed[,1],result_noseed[,1])) #phenotype-related genes
colnames(result_topten)="V1";colnames(links)=c("V1","V2")
temp1<-merge(result_topten,links,all.x=TRUE,sort=TRUE)
temp2=links
temp2[,1]<-links[,2]
temp2[,2]<-links[,1]
temp3<-merge(result_topten,temp2,all.x=TRUE,sort=TRUE)
temp4=temp3
temp3[,1]<-temp4[,2]
temp3[,2]<-temp4[,1]
temp5<-rbind(temp1,temp3)
temp6<-temp5[complete.cases(temp5),]
coreNet<-temp6[!duplicated(temp6, fromLast=TRUE), ] #write.table(coreNet,"coreNet.txt",quote=F,sep="\t",row.names=F,col.names=F)
#TFC
bet<-edge_betweenness(g, e = E(g), directed = F)
bet<-cbind(links,bet)
colnames(coreNet)=c("name1","name2");colnames(bet)=c("name1","name2","bet")
subbet<-merge(bet,coreNet,by=c("name1","name2"),sort=F)  #core NET edge bet
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)
temp7<-mgeneSim(nodes[,1], semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)
temp8 <- upper.tri(temp7)
temp9<-data.frame(row = rownames(temp7)[row(temp7)[temp8]],column = rownames(temp7)[col(temp7)[temp8]], cor =(temp7)[temp8] ) 
temp10=temp9
temp10[,1]=temp9[,2]
temp10[,2]=temp9[,1]
temp11<-rbind(temp9,temp10);colnames(temp11)=c("name1","name2","cor")
subGO<-merge(temp11,coreNet,by=c("name1","name2"),sort=F)
temp12=dplyr::anti_join(subbet[,1:2], subGO[,1:2],by = c("name1", "name2"))
temp12$cor=0
subGO=rbind(subGO,temp12)
subNET=merge(subGO,subbet,by=c("name1","name2"))
subNET$scaleTN <- (subNET$bet-min(subNET$bet))/(max(subNET$bet)-min(subNET$bet))
subNET$TFC  <- 100*(subNET$scaleTN+subNET$cor)/(2-subNET$scaleTN-subNET$cor) #write.table(subNET,"WCN.txt",sep="\t",quote=F,row.names=F)
#modules please remove outliers in advance
fan=subNET[,1:2];fan$TFC=subNET$TFC
g_fan = graph_from_data_frame(fan, directed=FALSE) #E(g_fan)$weight = fan$TFC
fan_GN<-cluster_edge_betweenness(g_fan,weights=NULL) #sizes(fan_GN) view quantity
fan_LP<-cluster_louvain(g_fan,weights=NULL)  
fan_GN_label=matrix(nrow = length(get.vertex.attribute(g_fan)[[1]]), ncol = 2)
fan_GN_label<-as.data.frame(fan_GN_label)
fan_GN_label[,1]=get.vertex.attribute(g_fan)[[1]]
fan_GN_label[,2]=fan_GN$membership
fan_LP_label=matrix(nrow = length(get.vertex.attribute(g_fan)[[1]]), ncol = 2)
fan_LP_label<-as.data.frame(fan_LP_label)
fan_LP_label[,1]=get.vertex.attribute(g_fan)[[1]]
fan_LP_label[,2]=fan_LP$membership
#write.table(fan_GN_label,"./module/WCN_GN.txt",sep="\t",quote=F,row.names=F)
#write.table(fan_LP_label,"./module/WCN_LP.txt",sep="\t",quote=F,row.names=F)
fan_GN_label_list <- list()
for(i in 1:length(fan_GN)){
  fan_GN_label_list[[i]]<-fan_GN_label[fan_GN_label$V2==i,1]
}
fan_LP_label_list <- list()
for(i in 1:length(fan_LP)){
  fan_LP_label_list[[i]]<-fan_LP_label[fan_LP_label$V2==i,1]
}#each nodes in each modules
fan_GN_label_list2<- list()
fan_GN_label_list2<-fan_GN_label_list[(lengths(fan_GN_label_list) >1)]
fan_LP_label_list2<- list()
fan_LP_label_list2<-fan_LP_label_list[(lengths(fan_LP_label_list) >1)]
fan_GN_LP2<-as.data.frame(matrix(nrow=length(fan_LP_label_list2),ncol=length(fan_GN_label_list2)));fan_GN_LP_union2=fan_GN_LP2;fan_GN_LP_phyper2=fan_GN_LP2
for(i in 1:length(fan_GN_label_list2)){
  for(j in 1:length(fan_LP_label_list2)){
    fan_GN_LP2[j,i]=length(intersect(fan_GN_label_list2[[i]],fan_LP_label_list2[[j]]))
    fan_GN_LP_union2[j,i]=length(union(fan_GN_label_list2[[i]],fan_LP_label_list2[[j]]))
  }
}
for(i in 1:ncol(fan_GN_LP_phyper2)){
  for(j in 1:nrow(fan_GN_LP_phyper2)){
    fan_GN_LP_phyper2[j,i]=1-phyper(fan_GN_LP2[j,i],length(fan_GN_label_list2[[i]]),fan_GN_LP_union2[j,i],length(fan_LP_label_list2[[j]]))
  }
}
colnames(fan_GN_LP_phyper2)<-paste("GN",seq(from=1,to=length(fan_GN_label_list2),by=1),sep="_")
rownames(fan_GN_LP_phyper2)<-paste("LP",seq(from=1,to=length(fan_LP_label_list2),by=1),sep="_")
cormat2=as.data.frame(t(fan_GN_LP_phyper2))
cormat2$ID <-  row.names(cormat2)
fanplot2 <- melt(cormat2, id.vars=c("ID"))
p_fan2 <- ggplot(fanplot2, aes(y=variable,x=ID))+geom_tile(aes(fill=value))+ scale_fill_gradient(low = "red", high = "skyblue") +labs(x='GN',y='LPA',title='Module Differences Significance')+theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45)) 
sig_fan2=as.data.frame(which(fan_GN_LP_phyper2<0.01, arr.ind = TRUE))
sig_fan_list2<- list()
for(i in 1:nrow(sig_fan2)){
  sig_fan_list2[[i]]<-intersect(fan_GN_label_list2[[sig_fan2[i,2]]],fan_LP_label_list2[[sig_fan2[i,1]]])
}
for(i in 1:length(sig_fan_list2))
{
  w=w+length(sig_fan_list2[[i]])
}
modules=as.data.frame(matrix(nc=2,nr=w))
for(i in 1:length(sig_fan_list2))
{
  w=f+1
  f=f+length(sig_fan_list2[[i]])
  modules[w:f,1]<- unlist(sig_fan_list2[[i]])
  modules[w:f,2]<- i
}
write.table(modules,"modules.txt",sep="\t",quote=F,row.names=F,col.names=c("name","label"))
###############################
