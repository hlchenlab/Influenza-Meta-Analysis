
#Packages:
if( !require("lemon")){
  install.packages("lemon")
}
if( !require("circlize")){
  install.packages("circlize")
}
if( !require("gplots")){
  install.packages("gplots")
}
if( !require("dplyr")){
  install.packages("dplyr")
}
if( !require("EGAD")){
  install.packages("EGAD")
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
if( !require("ggplot2")){
  install.packages("ggplot2")
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
if( !require("viridis")){
  BiocManager::install("viridis")
}

##Set-up

setwd("D:/")

################For Human#####################

load("GO Annotation Jan2020/GO.human.Rdata")
load("gene_annotations_v29.Rdata")
load("genesets_for_human_recurrence.rdata")
#Change genetype####
attr$type = as.character(attr$type)
attr[grep("pseudogene", attr$type),]$type = "Pseudogene"
attr[grep("3prime_overlapping", attr$type),]$type = "Other ncRNA"
attr[grep("bidirectional_promoter", attr$type),]$type = "lncRNA"
attr[grep("sense", attr$type),]$type = "lncRNA"
attr[grep("lincRNA", attr$type),]$type = "lncRNA"
attr[grep("IG_", attr$type),]$type = "Other ncRNA"
attr[grep("protein_coding", attr$type),]$type = "Protein coding"
attr[grep("TEC", attr$type),]$type = "Other ncRNA"
attr[grep("snoRNA", attr$type),]$type = "snRNA"
attr[grep("scaRNA", attr$type),]$type = "snRNA"
attr[grep("TR_", attr$type),]$type = "Other ncRNA"
attr[grep("misc_RNA", attr$type),]$type = "Other ncRNA"
attr[grep("Mt_", attr$type),]$type = "Other ncRNA"
attr[grep("macro", attr$type),]$type = "Other ncRNA"
attr[grep("processed", attr$type),]$type = "Other ncRNA"
attr[grep("ribozyme", attr$type),]$type = "Other ncRNA"
######################################################

#For Mouse
load("recurrence_matrix_MOUSE.Rdata")
load('mouse_gene_annotations_v100.Rdata')  ##Entrez ids to ensebl id conversion dataframe: can be acquired from BioMart

#Filter: GO group sizes
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
#filt = voc$V3 == "cellular_component"
#filt = voc$V3 == "molecular_function" 

# For Studies
res = list()
g.up = genesets.down  #Change
n = ncol(g.up)
i = 1
for (i in 1:n){
  genes = rownames(g.up[g.up[,i] == 1,])
  m = match (genes, rownames(annot))
  f.a = !is.na(m)
  genes = genes[f.a]
  resMG = gene_set_enrichment(genes, annot, voc[filt,])
  filt2 = resMG[,"padj"] <= 0.05
  res[[i]] = resMG[filt2,]
}
n = length(res)
TFterms = matrix(0,nrow = nrow(voc), ncol = n)
rownames(TFterms) = voc[,1]
colnames(TFterms) = colnames(g.up)
i =1
for (i in 1:n) {
  m = match(rownames(TFterms), as.matrix(res[[i]][,1]))
  f.a = !is.na(m)
  f.t = m[f.a]
  TFterms[f.a, i] = -log10(as.numeric(res[[i]][f.t,5]))
}
filt4 = rowSums(TFterms) > 0
g = TFterms[filt4,]

#Filter to get top 10 most significant terms
topVarGenes <- sapply(1:ncol(g), function(x) head(order(g[,x], decreasing = T),10)) %>% as.numeric() %>% unique(.)
m = match(rownames(g), voc[,1])
f.a = !is.na(m)
f.t = m[f.a]
rownames(g)[f.a] = as.character(voc[f.t,2])

#Greens:
col_fun = colorRamp2(c(seq(from = 10, to = 60, length.out = 30)), colorpanel(30,low="white", mid = "Green3",high="#003300"))
#Red:
col_fun = colorRamp2(c(seq(from = 10, to = 60, length.out = 30)), colorpanel(30,low="white", mid = "#FF3333",high="#660000"))
#Blue
col_fun = colorRamp2(c(seq(from = 0, to = 15, length.out = 20)), colorpanel(20,low="white", mid = "#33CCFF",high="#0000FF"))

colnames(g) = c("GSE84204","GSE89008", "GSE97672", "GSE103477", "GSE103604", "GSE104168", "GSE156060", "GSE156152")
save(h.up,h.down,m.up,m.down, file = "study gene enrichment.rdata")

ht =list()

g = h.up
topVarGenes <- sapply(1:ncol(g), function(x) head(order(g[,x], decreasing = T),10)) %>% as.numeric() %>% unique(.)
col_fun = colorRamp2(c(seq(from = 0, to = 12, length.out = 30)), colorpanel(30,low="white", mid = "#FF3333",high="#660000"))
ht[[1]] = Heatmap(head(g[topVarGenes,],25),
                  name = "Upregulated GO-BP\n -log10(P-value)",
                  col = col_fun,
                  use_raster = F,
                  rect_gp = gpar(col = "black", lwd = 2),
                  row_names_rot = 0,
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
                  column_names_rot = 90,
                  width = unit(4,"cm"),
                  height = unit(12,"cm"))
g = h.down
topVarGenes <- sapply(1:ncol(g), function(x) head(order(g[,x], decreasing = T),10)) %>% as.numeric() %>% unique(.)
col_fun = colorRamp2(c(seq(from = 0, to = 12, length.out = 30)), colorpanel(30,low="white", mid = "#33CCFF",high="#0000FF"))
ht[[2]] = Heatmap(head(g[topVarGenes,],25),
                  name = "Downregulated GO-BP\n -log10(P-value)",
                  col = col_fun,
                  use_raster = F,
                  rect_gp = gpar(col = "black", lwd = 2),
                  row_names_rot = 0,
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
                  column_names_rot = 90,
                  width = unit(4,"cm"),
                  height = unit(12,"cm"))
g = m.up
topVarGenes <- sapply(1:ncol(g), function(x) head(order(g[,x], decreasing = T),10)) %>% as.numeric() %>% unique(.)
col_fun = colorRamp2(c(seq(from = 10, to = 60, length.out = 30)), colorpanel(30,low="white", mid = "#FF3333",high="#660000"))
ht[[3]] = Heatmap(head(g[topVarGenes,],25),
                  name = "Upregulated GO-BP\n -log10(P-value)",
                  col = col_fun,
                  use_raster = F,
                  rect_gp = gpar(col = "black", lwd = 2),
                  row_names_rot = 0,
                  heatmap_legend_param = list(
                    at = c(0,20,40,60),
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
                  column_names_rot = 90,
                  width = unit(3.6,"cm"),
                  height = unit(12,"cm"))
g = m.down
topVarGenes <- sapply(1:ncol(g), function(x) head(order(g[,x], decreasing = T),10)) %>% as.numeric() %>% unique(.)
col_fun = colorRamp2(c(seq(from = 0, to = 15, length.out = 20)), colorpanel(20,low="white", mid = "#33CCFF",high="#0000FF"))
ht[[4]] = Heatmap(head(g[topVarGenes,],25),
                  name = "Downregulated GO-BP\n -log10(P-value)",
                  col = col_fun,
                  use_raster = F,
                  rect_gp = gpar(col = "black", lwd = 2),
                  row_names_rot = 0,
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
                  column_names_rot = 90,
                  width = unit(3.6,"cm"),
                  height = unit(12,"cm"))

pdf("E:/Thesis figures/Figure 3.pdf", height = 15,width = 10, useDingbats = F)
draw(ht[[1]] %v% ht[[2]], ht_gap = unit(0.5, "cm"))
draw(ht[[3]] %v% ht[[4]], ht_gap = unit(0.5, "cm"))
dev.off()

##Plot 2 - Genetypes per study

gmat = genesets.down
i=1
df = list()
for (i in 1:ncol(gmat)) {
  genes = data.frame(rownames(gmat[gmat[,i] == 1,]))
  colnames(genes) = "Genes"
  m = match(attr$ensemblID, genes$Genes)
  f.a = !is.na(m)
  f.t = m[f.a]
  genes$Genetypes = rep("Unknown", nrow(genes))
 genes[f.t,]$Genetypes = attr[f.a,]$type
  df[[i]] = (table(genes$Genetypes)/colSums(gmat)[i]) %>% as.data.frame()
colnames(df[[i]])[2] = colnames(gmat)[i]
}

tmat = matrix(0, nrow = length(unique(attr$type))+1, ncol = length(df))
rownames(tmat) = c(names(table(attr$type)[order(table(attr$type), decreasing = T)]),"Unknown")
colnames(tmat) = colnames(gmat)
for(i in 1:length(df)){
  m = match(df[[i]]$Var1, rownames(tmat))
  f.a =!is.na(m)
  f.t =m[f.a]
  tmat[f.t,i] = df[[i]][f.a,2]
}
tmat = tmat[rowSums(tmat)>0,]
df = reshape2::melt(tmat, value.name = c("Percentage"))
colnames(df)[1:2] = c("Genetype","Study")
df$Study  = as.character(df$Study)

save(df.down,df.up,m.up,m.down , file = 'human_mouse_genetype proportions.rdata')

#For Human
df[grep("HKUStudy_293", df$Study),]$Study = "GSE156152"
df[grep("HKUStudy", df$Study),]$Study = "GSE156060"
ids1 = c("GSE84204","GSE89008", "GSE97672", "GSE103477", "GSE103604", "GSE104168", "GSE156060", "GSE156152")

#For Mouse
ids2 = c("GSE49933", "GSE52405", "GSE100522", "GSE107488", "GSE117029", "ERP020504", "SRP061303") # list real studies here 

m.up$`Differential Expression` = 'Upregulated'
m.down$`Differential Expression` = 'Downregulated'
df.up$`Differential Expression` = 'Upregulated'
df.down$`Differential Expression` = 'Downregulated'
m.up$Species = 'Mouse'
m.down$Species = 'Mouse'
df.up$Species = 'Human'
df.down$Species = 'Human'

load('human_mouse_genetype proportions.rdata')

res = rbind.data.frame(df.up, df.down,m.up,m.down)
res$`Differential Expression` = factor(res$`Differential Expression`, levels = c('Upregulated', 'Downregulated'))
p = ggplot(data = res, aes(x = factor(Study, levels = c(ids1,ids2)), y = Percentage, fill =Genetype)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = 'PuBu',
                    direction = -1) +
  theme_classic() +
  labs(x = "Study", y = "Gene-type proportion") +
  facet_grid(`Differential Expression`~ Species, scales="free_x") +
  theme(axis.text.x.bottom = element_text(angle = 90,hjust=0.95,vjust=0.2),
        legend.position = "bottom",
        legend.key = element_rect(colour = 'black'))

g2 <- ggplot_gtable(ggplot_build(p))

stripr <- which(grepl('strip', g2$layout$name))

fills <- c("Red1","Green3", "#FF3333", "#33CCFF")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("Thesis figures/Figure 3d.pdf" ,height = 5, width = 7, useDingbats = F)
grid::grid.draw(g2)
dev.off()

##Plot 3 - Barplots of differential expression gene totals

load("genesets_for_human_recurrence.rdata")
colnames(genesets.down) = ids1
colnames(genesets.up) = ids1
ids = colnames(genesets.up)
df1 = cbind.data.frame(colSums(genesets.up, na.rm = TRUE),rep("Upregulated", length(colSums(genesets.up, na.rm = TRUE))), names(colSums(genesets.up, na.rm = TRUE)))
df2 = cbind.data.frame(colSums(genesets.down, na.rm = TRUE),rep("Downregulated", length(colSums(genesets.down, na.rm = TRUE))), names(colSums(genesets.down, na.rm = TRUE)))
colnames(df1) = c('Totals', 'Class','Study')
colnames(df2) = c('Totals', 'Class','Study')
df = rbind.data.frame(df1,df2)
df$Totals = as.numeric(df$Totals)
for (i in 1:nrow(df)) {
  if(df[i,]$Class == "Downregulated"){
    df[i,1] = df[i,1]*-1
  }
}
df$Species = 'Human'
m.df$Species = 'Mouse'
res1 = rbind.data.frame(m.df,df)
res1$Class = factor(res1$Class, levels = c("Upregulated", "Downregulated"))
save(res1, file = 'Differential_gene_totals.rdata')

load('Differential_gene_totals.rdata')

p = ggplot(res1,aes(x=factor(Study,levels = c(ids1,ids2) ),y=Totals,fill=Class)) +
  geom_col(color = "black")+
  labs(y = "Differentially expressed genes", x = "Studies") +
  theme_classic() + 
  scale_fill_manual("Class", values = c("Upregulated" = "#FF3333","Downregulated" = "#33CCFF")) +
  geom_text(
    aes(y = Totals+(Totals/abs(Totals)*200),label = abs(Totals),group = Study), size = 3) +
  theme_classic() +
  facet_grid(~Species, scales = 'free_x') +
  theme(axis.text.x.bottom = element_text(angle = 90,hjust=0.95,vjust=0.2),
        legend.position = "bottom",
        legend.key = element_rect(colour = 'black'))

g <- ggplot_gtable(ggplot_build(p))

stripr <- which(grepl('strip-t', g$layout$name))

fills <- c("Red1","Green3")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

b = plot.new()
g1 = plot_grid(g, b, ncol = 2, rel_widths = c(1.93, 0.07))
pdf("Thesis figures/Figure 3c+d.pdf" ,height = 9, width = 7, useDingbats = F)
plot_grid(g1, g2, nrow = 2, 
          labels = LETTERS[1:2], 
          label_size = 18,
          label_fontfamily = 'serif',rel_heights = c(1, 1))
dev.off()

##Plot 4 - Multifunctionality plots with DE prior ranks

load("Mouse_DE_prior.rdata")  ##See Crow et al, for permission
load("recurrence_matrix_MOUSE.Rdata")

#Human
load("DE_prior.rdata")   ##See Crow et al, for permission
load("gene_annotations_v29.Rdata")
load("genesets_for_human_recurrence.rdata")
colnames(genesets.down) = ids1
colnames(genesets.up) = ids1

genesets = genesets.up
genesets = genesets[rownames(genesets) %in% rownames(DE_list),, drop = F] # if using a specific annotation set e.g. clusters etc
m = match(rownames(DE_list), rownames(genesets)) 
f.go = !is.na(m)
DE_list_reviewed = DE_list[f.go,]
genesets <- genesets[rownames(DE_list_reviewed),]
DE_list_reviewed[,4] = rank(DE_list_reviewed[,4])

roc = list()
aucs.all = auc_multifunc( genesets, DE_list_reviewed[,4])
names(aucs.all) = colnames(genesets)

i=1
for (i in 1:dim(genesets)[[2]]) {
  roc[[i]] <- get_roc(DE_list_reviewed[,4], genesets[,i])
  roc[[i]] <- data.frame(roc[[i]])
  roc[[i]]$Studies = colnames(genesets)[i]
}

rocm= do.call(rbind.data.frame, roc)
labels = sapply(1:length(ids2), function(x) paste0(ids2[x],'\n(AUC =',round(aucs.all[x], digits = 2),')'))
rocm$Studies = factor(rocm$Studies, levels = ids2, labels = labels)

#DE prior multifunctionality scores per study saved in:
save(roc1, rocm, file = 'Study specific multifunc scores-de-prior.rdata')

p =ggplot(data = rocm ) +
  geom_line(aes(x = fpr, y = tpr)) +
  facet_wrap(~ Studies,nrow = 2) +
  geom_abline(color = 'Red3') +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = 'grey', linetype = 2, size = 0.2),
        panel.grid.major.y = element_line(color = 'grey', linetype = 2, size = 0.2))
g2 <- ggplot_gtable(ggplot_build(p))


stripr <- which(grepl('strip-t', g$layout$name))
#For Mouse:
stripr = stripr[-4]

fills <- rep("Red1", length(stripr))
fills <- rep("Green3", length(stripr))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#Print output:

pdf("Thesis figures/Figure 3b.pdf", height = 9, width = 8, useDingbats = F)
plot_grid(g, g2, nrow = 2,
          labels = LETTERS[1:2], 
          label_size = 18,
          label_fontfamily = 'serif')
dev.off()

