# hallmarks
R package hallmarks allows to perform over-representative analysis on a vector of genes with the hallmarks gene sets







## Script used in the package

### import data

> data(custom)

### or transform limma results in gene vector

> custom<-geneqvector(limma_res,updown="up",thres=0.5,q=0.05)


### compute of over-representative analysis on custum vector of genes

> res<-computehallmarks(custom,q=0.05)


### plot barplot of the analysis

> plothallmarks(res, font=14, label=5, title="Enriched hallmarks")

![hallmarks](https://github.com/cdesterke/hallmarks/blob/main/hallmarks.png)


### version 1.1.1 plot heatmap of jaccard index between significant gene set

> heatjaccard(custom,res,xmarg=22,ymarg=22)

![jaccard](https://github.com/cdesterke/hallmarks/blob/main/jaccard.png)


