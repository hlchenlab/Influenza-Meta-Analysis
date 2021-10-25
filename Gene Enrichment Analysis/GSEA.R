#############################

##Author Name: Conor Cremin, Sara Ballouz
##GO Analysis
##Dataset:Influenza

#############################

#Packages:
if( !require("lemon")){
  install.packages("lemon")
}
if( !require("gplots")){
  install.packages("gplots")
}
if( !require("dplyr")){
  install.packages("dplyr")
}
if( !require("DESeq2")){
  install.packages("DESeq2")
}
if( !require("gridExtra")){
  install.packages("gridExtra")
}
if( !require("ggplotify")){
  install.packages("ggplotify")
}
if( !require("cowplot")){
  install.packages("cowplot")
}
if( !require("GEOquery")){
  BiocManager::install("GEOquery")
}
if( !require("pheatmap")){
  install.packages("pheatmap")
}
if( !require("circlize")){
  install.packages("circlize")
}
if( !require("rlist")){
  install.packages("rlist")
}
if( !require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}

#Gene enrichment function:

gene_set_enrichment <- function(genes, genes.labels, voc){
  
  genes.names = rownames(genes.labels)
  labels.names = colnames(genes.labels)
  genes.counts = rowSums(genes.labels)
  labels.counts = colSums(genes.labels)              			# p
  
  m = match ( genes, genes.names )
  filt.genes  = !is.na(m)
  filt.labels = m[filt.genes]
  
  
  labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g
  
  m = match (labels.names, voc[,1])
  v.f = !is.na(m)
  v.g = m[v.f]
  
  universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
  if(  length(filt.labels) == 1 ) { genes.counts.set = genes.labels[filt.labels,] }
  else { genes.counts.set = colSums(genes.labels[filt.labels,]) }             ## does weird things with 0 sets
  
  test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
  pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
  pvals.adj = p.adjust( pvals, method="BH")
  
  results = cbind(voc[v.g,1:2], test[v.f,c(1)]+1, test[v.f,c(2)] , pvals[v.f], pvals.adj[v.f] )
  colnames(results) = c("term", "descrp","p", "q", "pvals", "padj" )
  return (results)
  
}

# Step 1: Load the GO matrices and the gene annotation files
#For Human
load("GO.human.Rdata")
load("gene_annotations_v29.Rdata")
#Or Mouse:
load("GO.mouse.Rdata")
attr = read.csv("Mouse_Annotation_Ensembl_v99.csv", header = T, stringsAsFactors = F) # Can get this in Differential Expression directory

#Step 2: Select GO group by sizes

gosums = colSums(GO.human.nonIEA )
filtGO = gosums <=200 & gosums>=100 
annot =  GO.human.nonIEA[,filtGO]

#Or for Mouse:
gosums = colSums(GO.mouse.nonIEA )
filtGO = gosums <=150 & gosums>=20
annot =  GO.mouse.nonIEA[,filtGO]

#Change label sizes:

m = match(rownames(annot), attr$entrezID)  #Convert EntrezIDs to Ensembl IDs  (Colnames may be different "attr" between files, so change accordingly)
f.a = !is.na(m)          
f.t = m[f.a]
rownames(annot)[f.a] = as.character(attr[f.t,]$ensemblID)

#Choose GO Catagory:

filt = voc$V3 == "biological_process"
#Or
filt = voc$V3 == "cellular_component"
#Or
filt = voc$V3 == "molecular_function" 

# Option 1: Extract significant enriched GO groups from enriched conditions

res = list()
g.up = geneset.up   # This dataframe was made from the last section of the Differential Expression.R code
n = ncol(g.up)
i = 1

for (i in 1:n){
  genes = rownames(g.up[g.up[,i] == 1,])
  m = match (genes, rownames(annot))
  f.a = !is.na(m)
  genes = genes[f.a]
  resMG = gene_set_enrichment(genes, annot, voc[filt,])
  filt2 = resMG[,"padj"] <= 0.05
  res[[i]] = resMG
}
n = length(res)
TFterms = matrix(0,nrow = nrow(voc), ncol = n)
rownames(TFterms) = voc[,1]
colnames(TFterms) = colnames(g.up)

# Option 2 for Clusters: (The Clusters dataframe corresponds to the clust dataframe generated from the get_cluster_ids function in the Clustering Analysis folder)

res = list()
n = length(unique(clusters$unmergedLabels.mod))
m = match(rownames(annot), clusters$V1)
f.a = !is.na(m)
genes = rownames(annot)[f.a]
resMG = gene_set_enrichment(genes, annot, voc1[filt,])
filt2 = resMG[,"padj"] < 0.05
res[[1]] = resMG[filt2,,drop = F]

for (i in 1:n){
  genes = as.character(clusters[clusters$unmergedLabels.mod ==i,]$V1)
  m = match (genes, rownames(annot))
  f.a = !is.na(m)
  genes = genes[f.a]
  resMG = gene_set_enrichment(genes, annot, voc1[filt,])
  filt2 = resMG[,"padj"] < 0.05
  res[[i+1]] = resMG[filt2,,drop = F]
}
n = length(res)
TFterms = matrix(0,nrow = nrow(voc), ncol = n)
rownames(TFterms) = voc[,1]
colnames(TFterms) = c("All",unique(clusters$unmergedLabels.mod))

i =1
for (i in 1:n) {
  m = match(rownames(TFterms), as.matrix(res[[i]][,1]))
  f.a = !is.na(m)
  f.t = m[f.a]
  TFterms[f.a, i] = -log10(as.numeric(res[[i]][f.t,5]))
}
filt4 = rowSums(TFterms) > 0
g = TFterms[filt4,]

#Select the top 5 GO groups in each condition:

topVarGenes <- sapply(1:ncol(g), function(x) head(order(g[,x], decreasing = T),5)) %>% as.numeric() %>% unique(.)

#Convert IDs:
m = match(rownames(g), voc[,1])
f.a = !is.na(m)
f.t = m[f.a]
rownames(g)[f.a] = as.character(voc[f.t,2])

#Colours:
#Greens:
col_fun = colorRamp2(c(seq(from = 0, to = 12, length.out = 15)), colorpanel(15,low="white", mid = "Green3",high="#003300"))

Heatmap(as.matrix(g[topVarGenes,]),
        name = "Gene Enrichment \n -log10(P-value)",
        col = col_fun,
        rect_gp = gpar(col = "black", lwd = 1),
        column_title_gp = gpar(fontsize = 7, fontface = "bold"),
        column_title = "Conditions",
        row_title_gp = gpar(fontsize = 7, fontface = "bold"),
        row_names_rot = 0,
        column_gap = unit(2,"mm"),
        heatmap_legend_param = list(
          at = c(0,4,8,12),
          color_bar="continuous",
          legend_height = unit(2, "cm"),
          legend_width = unit(2, "cm"),
          title_position = "topleft",
          border = "black",
          labels_gp = gpar(fontsize = 6),
          title_gp=gpar(fontsize=6, fontface="bold"),
          legend_direction = "horizontal"),
        cluster_columns = F,
        column_dend_side = "bottom",
        show_column_dend = F,
        cluster_rows = T,
        show_column_names = F,
        row_names_gp = gpar(fontsize = 7),
        row_names_side = "right", row_dend_side = "left", 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 4),
        column_names_rot = 90,
        height = unit(12,"cm"),
        width = unit(3,"cm")
)
