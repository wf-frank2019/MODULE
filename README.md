# COMMUNITY DETECTION
The node and edge Prioritization based Community Analysis (ne-PCA) for NSCLC Protein-Protein Interaction Network 

##  Overall Presentation
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/Outline1.png "Outline")
   
### 1 NSCLC seeds genes collection
#### 1 DisGnet(http://www.disgenet.org/downloads)  158 in total   
	Screening："Non-Small Cell Lung Carcinoma":111  "Adenocarcinoma of lung":60
#### 2 KEGG(https://www.kegg.jp/kegg/)  81 in total
	Network variation(https://www.kegg.jp/kegg-bin/show_network?id=nt06266) nt06266:44
	PATHWAY(https://www.kegg.jp/kegg-bin/show_pathway?hsa05223+H00014) map05223:64
  	DISEASE(https://www.kegg.jp/dbget-bin/www_bget?ds:H00014) H00014:9 
#### 3 Cancer Gene Census(https://cancer.sanger.ac.uk/census)  42 in total
     	Screening："NSCLC“:37  “lung adenocarcinoma":5
#### 4 OMIM(https://www.omim.org/downloads)  24 in total
	Screening："Nonsmall cell lung cancer":14  "Adenocarcinoma of lung":12 "lung adenocarcinoma":2
#### 5 Diseases(https://diseases.jensenlab.org/Downloads)  2 in total
	Screening："non-small cell lung carcinoma":2
#### 6 IntAct(https://www.ebi.ac.uk/intact/downloads)  18 in total
	Screening："non-small-cell lung cancer":9 "non-small cell lung cancer":4 "NSCLC":6
#### 7 OncoPPI(http://oncoppi.emory.edu/)  78 in total

### 2 Preliminary omics research
	GSE19804 & GSE101929 (94 NSCLCs vs 92 normal samples)
   	DEGs were identified respectively(588 common DEGs)

### 3 PPIN Construction
	Protein-protein interaction network based on STRING(https://string-db.org/)
   	DEGs include 12 known genes from KEGG(https://www.kegg.jp/kegg/),COSMIC(https://cancer.sanger.ac.uk/cosmic/),DisGenet(https://www.disgenet.org/)
	KeyNetwork：nodes and edges are weighted by RWRscore and TFscore respectively
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/git2.PNG "PPN_WCN")

### 4 Communities
	Communities are detected by GN_LPA Model
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/community.PNG "Module")


