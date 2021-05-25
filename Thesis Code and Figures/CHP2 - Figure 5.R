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

# Step 1: Load the GO matrices and the gene annotation files
setwd('E:/')
load('Coexpression-influenza-human-network and datasets.rdata')
load("E:/gene_annotations_v29.Rdata")
load("GO Annotation Jan2020/GO.human.Rdata")

#Filter1: GO group sizes
gosums = colSums(GO.human.nonIEA )
filtGO = gosums <=500 & gosums>=100
annot =  GO.human.nonIEA[,filtGO]

#Change label sizes:

m = match(rownames(annot), attr$entrezID)
f.a = !is.na(m)          
f.t = m[f.a]
rownames(annot)[f.a] = as.character(attr[f.t,]$ensemblID)

#Voc filter:

filt = voc$V3 == "biological_process"

res = list()
n = length(unique(clusters$unmergedLabels.mod))
m = match(rownames(annot), clusters$V1)
f.a = !is.na(m)
genes = rownames(annot)[f.a]
resMG = gene_set_enrichment(genes, annot, voc[filt,])
filt2 = resMG[,"padj"] < 0.05
res[[1]] = resMG[filt2,,drop = F]

for (i in 1:n){
  genes = as.character(clusters[clusters$unmergedLabels.mod ==i,]$V1)
  m = match (genes, rownames(annot))
  f.a = !is.na(m)
  genes = genes[f.a]
  resMG = gene_set_enrichment(genes, annot, voc[filt,])
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

m = match(rownames(g), voc[,1])
f.a = !is.na(m)
f.t = m[f.a]
rownames(g)[f.a] = as.character(voc[f.t,2])
gh =g

prots = c("white", unique(clusters$V3))
names(prots) = colnames(gh)
cc = data.frame(colnames(gh))
colnames(cc) = "Clusters"
hb = columnAnnotation(df = cc,
                      col = list(Clusters = prots),
                      show_legend = F,
                      show_annotation_name = F,
                      annotation_label = NULL,
                      border = T)
col_fun = colorRamp2(c(seq(from = 0, to = 12, length.out = 30)), colorpanel(30,low="white", mid = "#FF3333",high="#660000"))

an=hb
topVarGenes <- sapply(1:ncol(gh), function(x) head(order(gh[,x], decreasing = T),10)) %>% as.numeric() %>% unique(.)
ht =Heatmap(head(gh[topVarGenes,],25),
        name = "Upregulated GO-BP\n-log10(P-value)",
        col = col_fun,
        use_raster = F,
        rect_gp = gpar(col = "black", lwd = 2),
        row_names_rot = 0,
        top_annotation = an,
        heatmap_legend_param = list(
          at = c(0,4,8,12),
          color_bar="continuous",
          legend_direction = "vertical",
          legend_height = unit(2, "cm"),
          title_position = "leftcenter-rot",
          title_gp=gpar(fontsize=8, fontface="bold"),
          border = T),
        cluster_columns = F,
        show_column_dend = F,
        cluster_rows = T,
        show_row_dend = F,
        row_names_gp = gpar(fontsize = 10),
        row_names_side = "left", 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 0,
        column_names_centered = T,
        width = unit(3,"cm"),
        height = unit(10,"cm"))
ht = grid::grid.grabExpr(draw(ht,heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE))

#Mouse - Influenza

load('Coexpression-influenza-mouse-network and datasets.rdata')

#Filter1: GO group sizes
load("GO Annotation Jan2020/GO.mouse.Rdata")
load('mouse_gene_annotations_v100.Rdata')
gosums = colSums(GO.mouse.nonIEA )
filtGO = gosums <=500 & gosums>=100
annot =  GO.mouse.nonIEA[,filtGO]

#Change label sizes:
attrm = read.csv("E:/mouse_ext_IDs.csv", header = T, stringsAsFactors = F)
m = match(rownames(annot), attrm$NCBI.gene.ID)
f.a = !is.na(m)          
f.t = m[f.a]
rownames(annot)[f.a] = as.character(attrm[f.t,]$Gene.stable.ID)

#Voc filter:

filt = voc$V3 == "biological_process"

res = list()
n = length(unique(clusters$unmergedLabels.mod))
m = match(rownames(annot), clusters$V1)
f.a = !is.na(m)
genes = rownames(annot)[f.a]
resMG = gene_set_enrichment(genes, annot, voc[filt,])
filt2 = resMG[,"padj"] < 0.05
res[[1]] = resMG[filt2,,drop = F]

for (i in 1:n){
  genes = as.character(clusters[clusters$unmergedLabels.mod ==i,]$V1)
  m = match (genes, rownames(annot))
  f.a = !is.na(m)
  genes = genes[f.a]
  resMG = gene_set_enrichment(genes, annot, voc[filt,])
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

m = match(rownames(g), voc[,1])
f.a = !is.na(m)
f.t = m[f.a]
rownames(g)[f.a] = as.character(voc[f.t,2])


prots = c("white", unique(clusters$V3))
names(prots) = colnames(TFterms)
cc = data.frame(colnames(g))
colnames(cc) = "Clusters"
hb = columnAnnotation(df = cc,
                      col = list(Clusters = prots),
                      show_legend = F,
                      show_annotation_name = F,
                      annotation_label = NULL,
                      border = T)
gm =g
anm =hb
topVarGenes <- sapply(1:ncol(gm), function(x) head(order(gm[,x], decreasing = T),10)) %>% as.numeric() %>% unique(.)
col_fun = colorRamp2(c(seq(from = 0, to = 12, length.out = 30)), colorpanel(30,low="white", mid = "#FF3333",high="#660000"))
htm = Heatmap(head(gm[topVarGenes,],25),
        name = "Upregulated GO-BP\n -log10(P-value)",
        col = col_fun,
        use_raster = F,
        rect_gp = gpar(col = "black", lwd = 2),
        row_names_rot = 0,
        top_annotation = anm,
        heatmap_legend_param = list(
          at = c(0,4,8,12),
          color_bar="continuous",
          legend_direction = "vertical",
          legend_height = unit(2, "cm"),
          title_position = "leftcenter-rot",
          title_gp=gpar(fontsize=8, fontface="bold"),
          border = T),
        cluster_columns = F,
        show_column_dend = F,
        cluster_rows = T,
        show_row_dend = F,
        row_names_gp = gpar(fontsize = 10),
        row_names_side = "left", 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 0,
        column_names_centered = T,
        width = unit(2.8,"cm"),
        height = unit(10,"cm"))

htm = grid::grid.grabExpr(draw(htm,heatmap_legend_side = "right",merge_legend = TRUE))

#Filter1: GO group sizes
gosums = colSums(GO.human.nonIEA )
filtGO = gosums <=500 & gosums>=100
annot =  GO.human.nonIEA[,filtGO]

#Change label sizes:

m = match(rownames(annot), attr$entrezID)
f.a = !is.na(m)          
f.t = m[f.a]
rownames(annot)[f.a] = as.character(attr[f.t,]$ensemblID)

#Voc filter:

filt = voc$V3 == "biological_process"
load('Coexpression-Sars-CoV-2-human-network and datasets.rdata')
clusters = clustersh
res = list()
n = length(unique(clusters$unmergedLabels.mod))
m = match(rownames(annot), clusters$V1)
f.a = !is.na(m)
genes = rownames(annot)[f.a]
resMG = gene_set_enrichment(genes, annot, voc[filt,])
filt2 = resMG[,"padj"] < 0.05
res[[1]] = resMG[filt2,,drop = F]

for (i in 1:n){
  genes = as.character(clusters[clusters$unmergedLabels.mod ==i,]$V1)
  m = match (genes, rownames(annot))
  f.a = !is.na(m)
  genes = genes[f.a]
  resMG = gene_set_enrichment(genes, annot, voc[filt,])
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

m = match(rownames(g), voc[,1])
f.a = !is.na(m)
f.t = m[f.a]
rownames(g)[f.a] = as.character(voc[f.t,2])
gcov =g

prots = c("white", unique(clusters$V3))
names(prots) = colnames(gcov)
cc = data.frame(colnames(gcov))
colnames(cc) = "Clusters"
hb = columnAnnotation(df = cc,
                      col = list(Clusters = prots),
                      show_legend = F,
                      show_annotation_name = F,
                      annotation_label = NULL,
                      border = T)
col_fun = colorRamp2(c(seq(from = 0, to = 12, length.out = 30)), colorpanel(30,low="white", mid = "#FF3333",high="#660000"))

acov=hb
topVarGenes <- sapply(1:ncol(gcov), function(x) head(order(gcov[,x], decreasing = T),10)) %>% as.numeric() %>% unique(.)
hcov =Heatmap(head(gcov[topVarGenes,],25),
            name = "Upregulated GO-BP\n-log10(P-value)",
            col = col_fun,
            use_raster = F,
            rect_gp = gpar(col = "black", lwd = 2),
            row_names_rot = 0,
            top_annotation = acov,
            heatmap_legend_param = list(
              at = c(0,4,8,12),
              color_bar="continuous",
              legend_direction = "vertical",
              legend_height = unit(2, "cm"),
              title_position = "leftcenter-rot",
              title_gp=gpar(fontsize=8, fontface="bold"),
              border = T),
            cluster_columns = F,
            show_column_dend = F,
            cluster_rows = T,
            show_row_dend = F,
            row_names_gp = gpar(fontsize = 10),
            row_names_side = "left", 
            column_names_side = "top",
            column_names_gp = gpar(fontsize = 10),
            column_names_rot = 0,
            column_names_centered = T,
            width = unit(1.5,"cm"),
            height = unit(10,"cm"))


pdf("Thesis figures/Chp2-Figure 5.pdf", height = 8,width = 7, useDingbats =  F)
ht
htm
hcov
dev.off()


