# COMMUNITY DETECTION
The node and edge Prioritization based Community Analysis (ne.PCA) for Protein-Protein Interaction Network<br>
ref: https://www.mdpi.com/2073-4409/10/2/402

##  Usage
	# install.packages("devtools")
	library(devtools)
	devtools::install_github("wf-frank2019/ne.PCA")
	library(ne.PCA) 
	# ? or help(" ") for description of five functions：qc.DEGs();RWRs.CN();TFCs();MODULEs();ne.PCA() 
	
	# If you wanna identify key biological network components and robust modules, just type：
	ne.PCA(nodes,edges)

##  Overall Presentation
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/overview.PNG "Outline")
   
### 1 Preliminary omics research
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/GSE.PNG "Microarray")
        DEGs were identified respectively:588 common DEGs

### 2 NSCLC seed genes collection
	NSCLC seed genes were collected from seven database with different focus
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/seeds.PNG "Seeds")
#### details:
#### DisGeNet(http://www.disgenet.org/downloads) Pinero J, Ramirez-Anguita JM, Sauch-Pitarch J et al. The DisGeNET knowledge platform for disease genomics: 2019 update,NU-CLEIC ACIDS RESEARCH 2020.
#### KEGG(https://www.kegg.jp/kegg/)   Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes,NUCLEIC ACIDS RESEARCH 2000.
#### Cancer Gene Census(https://cancer.sanger.ac.uk/census) Sondka,Z.,Bamford, S.,Cole,C.G.et al.The COSMIC Cancer Gene Census: describing genetic dysfunction across all human cancers,Nat Rev Cancer 2018.
#### OMIM(https://www.omim.org/downloads)  Amberger JS, Bocchini CA, Scott AF, Hamosh A. OMIM.org: leveraging knowledge across phenotype-gene relationships,Nucleic Acids Res 2019.
#### Diseases(https://diseases.jensenlab.org/Downloads) Pletscher-Frankild S, Pallejà, Albert, Tsafou K, et al. DISEASES: Text mining and data integration of disease–gene associations,Methods 2015.
#### IntAct(https://www.ebi.ac.uk/intact/downloads) Orchard S1, Ammari M, Aranda B, et al. The MIntAct project--IntAct as a common curation platform for 11 molecular interaction databases,Nucleic Acids Research 2014.
#### OncoPPI(http://oncoppi.emory.edu) Li,Z.,Ivanov,A.,Su,R.et al.The OncoPPi network of cancer-focused protein–protein interactions to inform biological insights and therapeutic strategies,Nat Commun 2017.
	
### 3 PPIN Construction
	Protein-protein interaction network based on STRING(https://string-db.org/)
   	DEGs include 12 known genes from KEGG(https://www.kegg.jp/kegg/),COSMIC(https://cancer.sanger.ac.uk/cosmic/),DisGenet(https://www.disgenet.org/)
	KeyNetwork：nodes and edges are weighted by RWRscore and TFscore respectively
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/git2.PNG "PPN_WCN")

### 4 Communities
	Communities are detected by GN_LPA Model
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/community.PNG "Module")


