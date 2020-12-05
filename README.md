# COMMUNITY DETECTION
The node and edge Prioritization based Community Analysis (ne-PCA) for NSCLC Protein-Protein Interaction Network 

## 一 Overall Presentation
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/Outline1.png "Outline")
### 1 Preliminary omics research
	GSE19804 & GSE101929 (94 NSCLCs vs 92 normal samples)
   	DEGs were identified respectively(588 common DEGs)

### 2 PPIN Construction
	Protein-protein interaction network based on STRING(https://string-db.org/)
   	DEGs include 12 known genes from KEGG(https://www.kegg.jp/kegg/),COSMIC(https://cancer.sanger.ac.uk/cosmic/),DisGenet(https://www.disgenet.org/)
	KeyNetwork：nodes and edges are weighted by RWRscore and TFscore respectively
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/git2.PNG "PPN_WCN")

### 3 Communities
	Communities are detected by GN_LPA Model
   ![frank](https://github.com/wf-frank2019/-storehouse/blob/master/res/community.PNG "Module")


