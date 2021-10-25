ident<-read.table("id.txt",h=T,sep="\t")

genes<-read.table("hallmarks.csv",h=F)
colnames(genes)<-c("genes")
genesplit<-strsplit(as.character(genes$genes),split=';')


genesplit<-lapply(genesplit, function(z){ z[!is.na(z) & z != ""]})


custom<-c("PCNA","PSMD8","PSMD7","SET","SNRPA1","RAN","SRSF2","G3BP1","STARD7","NPM1","BUB3","EIF3D",
"XPO1","FBL","EIF4A1","CANX","NAP1L1","ANGPTL4","ITGA2","SPRY2","AURKA","BRCA2","CCP110","CENPE","CKS2","DCLRE1B","DNMT1",
"PTPRC","IL12B","TGFB1","IL12A","CD3E","CD3D","CD28","LYN","HCLS1","IL18","CRTAM","IFNG","CD3G","CD86","IL10","FOS",
"TAP1","E2F5","RAB27A","CCND3","SLC25A4","SPR","CDKN2B","SHOX2","OLFM1","IGFBP2","LYN","FOXO3","BTG1","TXNIP")

hallmarks <- ident$pathway
hallmarks <- as.data.frame(hallmarks)

for(i in 1:length(genesplit)){
hallmarks$intersect[i]<-length(intersect(custom,genesplit[[i]]))}
hallmarks$input<-length(custom)
for(i in 1:length(genesplit)){
hallmarks$geneset[i]<-length(genesplit[[i]])}
hallmarks$totaldb<-length(unique(unlist(genesplit)))

hallmarks


hallmarks<-hallmarks[(hallmarks$intersect != "0"),]
df<-hallmarks[with(hallmarks,order(-intersect)),]

row.names(df)<-df$hallmarks
df$hallmarks<-NULL

res1 <- NULL
for (i in 1:nrow(df)){
  table <- matrix(c(df[i,1], df[i,2], df[i,3], df[i,4]), ncol = 2, byrow = TRUE)
  o <- fisher.test(table, alternative="two.sided")$estimate
  # save all odds in a vector
  res1 <- c(res1,o)
}
df$ES <- res1
df


res2 <- NULL
for (i in 1:nrow(df)){
  table <- matrix(c(df[i,1], df[i,2], df[i,3], df[i,4]), ncol = 2, byrow = TRUE)
  p <- fisher.test(table, alternative="two.sided")$p.value
  # save all p values in a vector
  res2 <- c(res2,p)
}
df$pvalues <- res2
df

df$qvalues<-p.adjust(df$pvalues,method="fdr")

df<-df[with(df,order(qvalues)),]
df<-df[(df$qvalues <= 0.05),]

library(dplyr)

ident%>%select(-link)->ident

df$pathway<-row.names(df)

ident%>%right_join(df,by="pathway")->df

library(ggplot2)
library(pals)

p=ggplot(data=df,aes(x=reorder(pathway,ES),y=ES,fill=class))+geom_bar(stat="identity")+
coord_flip()+theme_minimal()
p +  xlab("Genesets") + ylab("Normalized enrichment score")
p + scale_fill_manual(values=cols25())+
geom_text(aes(label=round(ES,2)),hjust=0, vjust=0.5,color="black",position= position_dodge(0),size=6,angle=0)+
xlab("Gene sets") + ylab("Enrichment score")+
ggtitle("Enriched Hallmarks") +theme(text = element_text(size = 16))

save(ident,file="ident.rda")
save(genesplit,file="genesplit.rda")

