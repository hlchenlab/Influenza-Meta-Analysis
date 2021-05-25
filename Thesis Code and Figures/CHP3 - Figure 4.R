#Packages:

if( !require("dplyr")){
  BiocManager::install("dplyr")
}
if( !require("ggrepel")){
  BiocManager::install("ggrepel")
}
if( !require("cowplot")){
  BiocManager::install("cowplot")
}
if( !require("Seurat")){
  BiocManager::install("Seurat")
}
if( !require("patchwork")){
  BiocManager::install("patchwork")
}
if( !require("gtable")){
  BiocManager::install("gtable")
}
if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}
if( !require("gplots")){
  BiocManager::install("gplots")
}
if( !require("grid")){
  BiocManager::install("grid")
}
if( !require("gridExtra")){
  BiocManager::install("gridExtra")
}
if( !require("future")){
  BiocManager::install("future")
}
if( !require("future.apply")){
  BiocManager::install("future.apply")
}
if( !require("sctransform")){
  BiocManager::install("sctransform")
}
if( !require("ggpubr")){
  BiocManager::install("ggpubr")
}


#Set-up

setwd("E:/")
load("single-cell data/SCTransform/wsn.integrated.rdata")
lvls = c('B Cell','CD4+ T Cell','Th2 CD4+ Cell','CD8+ T Cell','Cytotoxic CD8+ T Cell','NK Cell','Neutrophil', 'Macrophage',
         'Epithelial Cell', 'Type I Pneumocyte', 'Type II Pneumocyte','Lipofibroblast', 'Myofibroblast', 'Fibroblast', 'Matrix Fibroblast',
         'Endothelial Cell')

#Getting UMI Totals per celltype

df = matrix(0, nrow = length(lvls), ncol = length(unique(wsn.integrated$Condition)))
colnames(df) = unique(wsn.integrated$Condition)
rownames(df) = lvls
geneUMI = list()
k =1
for (i in 1:nrow(df)) {
  for (j in 1:ncol(df)) {
    df[i,j] = sum(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[j]])
  }
}

#Calculate normalized gene expression and relative
k=1
i=1

for (i in 1:nrow(df)) {
  cells = names(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[1]])
  geneUMI[[k]] = (rowSums(wsn.integrated[["RNA"]]@counts[,cells])/df[i,1])*100000
  print(paste(paste0(colnames(df)[1],"-",rownames(df)[i]),":",median(geneUMI[[k]][geneUMI[[k]] >0 ])))
  geneUMI[[k]] = data.frame(geneUMI[[k]][geneUMI[[k]] >0]/median(geneUMI[[k]][geneUMI[[k]] >0 ]))
  colnames(geneUMI[[k]]) = paste0(colnames(df)[1],"-",rownames(df)[i])
  
  k = k+1
  cells = names(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[2]])
  geneUMI[[k]] = (rowSums(wsn.integrated[["RNA"]]@counts[,cells])/df[i,2])*100000
  print(paste(paste0(colnames(df)[2],"-",rownames(df)[i]),":",median(geneUMI[[k]][geneUMI[[k]] >0 ])))
  geneUMI[[k]] = data.frame(geneUMI[[k]][geneUMI[[k]] >0]/median(geneUMI[[k]][geneUMI[[k]] >0 ]))
  colnames(geneUMI[[k]]) = paste0(colnames(df)[2],"-",rownames(df)[i])
  
  k = k+1
  cells = names(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[3]])
  geneUMI[[k]] = (rowSums(wsn.integrated[["RNA"]]@counts[,cells])/df[i,3])*100000
  print(paste(paste0(colnames(df)[3],"-",rownames(df)[i]),":",median(geneUMI[[k]][geneUMI[[k]] >0 ])))
  geneUMI[[k]] = data.frame(geneUMI[[k]][geneUMI[[k]] >0]/median(geneUMI[[k]][geneUMI[[k]] >0 ]))
  colnames(geneUMI[[k]]) = paste0(colnames(df)[3],"-",rownames(df)[i])
  
  k = k+1
  cells = names(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[4]])
  geneUMI[[k]] = (rowSums(wsn.integrated[["RNA"]]@counts[,cells])/df[i,4])*100000
  print(paste(paste0(colnames(df)[4],"-",rownames(df)[i]),":",median(geneUMI[[k]][geneUMI[[k]] >0 ])))
  geneUMI[[k]] = data.frame(geneUMI[[k]][geneUMI[[k]] >0]/median(geneUMI[[k]][geneUMI[[k]] >0 ]))
  colnames(geneUMI[[k]]) = paste0(colnames(df)[4],"-",rownames(df)[i])
  
  k = k+1
}

pcc = matrix(0, nrow = length(geneUMI), ncol = length(geneUMI))
colnames(pcc) = sapply(geneUMI, colnames)
rownames(pcc) = sapply(geneUMI, colnames)

i=1
j=1
for (i in 1:length(geneUMI)) {
  for (j in 1:length(geneUMI)) {
    m = match(rownames(geneUMI[[i]]), rownames(geneUMI[[j]]))
    f.a = !is.na(m)
    f.t = m[f.a]
    x = geneUMI[[i]][f.a,, drop = F]
    y = geneUMI[[j]][f.t,,drop = F]
    pcc[colnames(geneUMI[[i]]),colnames(geneUMI[[j]])] = as.numeric(cor.test(x[,1], y[,1], method = "pearson")["estimate"])
  }
}

df = colnames(pcc) %>% data.frame()
df$Condition = matrix(unlist(strsplit(as.character(df[,1]), "-")),ncol=2, byro=T)[,1]
df$Celltype = matrix(unlist(strsplit(as.character(df[,1]), "-")),ncol=2, byro=T)[,2]
rownames(df) = df[,1]
df$Condition = factor(df$Condition, levels = c("Mock", "Wildtype","3M","dNS1"))
df$Celltype = factor(df$Celltype, levels = lvls)
df = df[,-1]
names(col_fun) = lvls
col_fun = col_fun[1:length(unique(df$Celltype))]
hb = columnAnnotation(df = df,
                      col =list(Condition = colour_sample,
                                Celltype = col_fun),
                      annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                      show_legend = T,
                      border = T,
                      gap = unit(1, "mm"),
                      annotation_name_side = "left",
                      annotation_legend_param = list(title_gp=gpar(fontsize=10, fontface="bold"),
                                                     labels_gp = gpar(fontsize = 8),
                                                     legend_height = unit(5, "cm")))
mat = pcc

ht = Heatmap(mat,
             name = "PCC",
             col = plasma(100),
             column_title_gp = gpar(fontsize = 18, fontface = "bold"),
             row_title_gp = gpar(fontsize = 18, fontface = "bold"),
             heatmap_legend_param = list(
               at = c(0,0.5,1),
               labels_gp = gpar(fontsize = 8),
               legend_height = unit(2, "cm"),
               title_position = "topleft",
               title_gp=gpar(fontsize=10, fontface="bold"),
               border = T),
             top_annotation = hb,
             cluster_columns = T,
             column_dend_side = "top",
             show_row_dend = F,
             show_column_dend = F,
             cluster_rows = T,
             row_names_gp = gpar(fontsize = 5),
             show_row_names = T,
             row_names_side = "left",
             use_raster = T,
             show_column_names = F,
             column_names_rot = 90,
             height = unit(12,"cm"),
             width = unit(12,"cm"))
ht = grid::grid.grabExpr(draw(ht,heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE))

#Plot Barplot of DEGs

load("single-cell data/SCTransform/deg_totals_df.rdata")
df = df[df$Celltype %in% lvls,]
colnames(df)[4] = 'Differential expression'
df$Celltype = factor(df$Celltype, levels = lvls)
df$Condition = factor(df$Condition, levels = c("Wildtype","3M","dNS1"))
df$`Differential expression` = factor(df$`Differential expression`, levels = c("Upregulated", "Downregulated"))
df$Total = as.numeric(df$Total)
for (i in 1:nrow(df)) {
  if(df[i,]$`Differential expression` == "Downregulated"){
    df[i,5] = df[i,5]*-1
  }
}

p = ggplot(df,aes(x= factor(Celltype,levels = rev(lvls)),y=Total,fill=`Differential expression`)) +
  geom_col(color = "black")+
  labs(y = "Differentially expressed genes", x = NULL) +
  theme_classic() + 
  facet_grid(~Condition) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_fill_manual("Differential expression", values = c("Upregulated" = "#FF3333","Downregulated" = "#33CCFF")) +
  geom_text(aes(y = Total+(Total/abs(Total)*200),label = abs(Total),group = Celltype), size = 2.5) +
  theme_classic() +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = 2,size = 0.5, fill =NA), 
        panel.spacing.x = unit(0.5,"cm"),
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))

g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- colour_sample[2:4]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("Thesis figures/Chp3-Figure 4.pdf", width = 9, height = 11, paper = 'a4', useDingbats = F)
plot_grid(ht,g, nrow = 2,label_y = c(1,1.06), rel_heights = c(1.3,1),labels = c(LETTERS[1:2]), label_fontfamily = 'serif', label_size = 18)
dev.off()

####################Expression plots####################

load("single-cell data/SCTransform/allgenes_mouse.rdata")
orthos = read.csv("mouse_human.csv") %>% .[.$Mouse.homology.type == "ortholog_one2one",]
load('cross-species.clusters.rdata')
col_fun = rainbow(n = length(lvls), start = 0, end = 0.7, s = 0.7, v= 0.85) 
names(col_fun) = lvls
colour_sample = c("Mock" = "tan1","Wildtype" = "firebrick1", "3M" = "limegreen", "dNS1"="blue2")

m = match(h.clusters$genes, orthos$Gene.stable.ID)
f.a =!is.na(m)
f.t =m[f.a]
h.clusters$ids = rep(0, nrow(h.clusters))
h.clusters[f.a,]$ids = paste0("mmc9-",orthos[f.t,]$Mouse.gene.name)
h.clusters = h.clusters[f.a,]

m = match(m.clusters$genes, orthos$Gene.stable.ID)
f.a =!is.na(m)
f.t =m[f.a]
m.clusters[f.a,]$genes = orthos[f.t,]$Mouse.gene.stable.ID
m = match(m.clusters$genes, allgenes$Gene.stable.ID)
f.a =!is.na(m)
f.t =m[f.a]
m.clusters$ids = rep(0, nrow(m.clusters))
m.clusters[f.a,]$ids = paste0("mmc9-",allgenes[f.t,]$Gene.name)
m.clusters = m.clusters[f.a,]

#Scale and normalize by "LogNormalize"

DefaultAssay(wsn.integrated) = "RNA"
wsn.integrated = NormalizeData(wsn.integrated)
wsn.integrated = FindVariableFeatures(wsn.integrated, nfeatures = 5000)
wsn.integrated <- ScaleData(wsn.integrated, vars.to.regress = c("percent.mt","nCount_RNA", "nFeature_RNA"),verbose = FALSE)

res = list()
cells = wsn.integrated@assays$RNA@data[row.names(wsn.integrated@assays$RNA@data) %in% m.clusters$ids,] %>% as.data.frame()
i=1

#Average cluster expression per cell
for (i in 1:length(unique(m.clusters$labels))){
  m = match(rownames(cells), m.clusters[m.clusters$labels == i,]$ids)
  f.a = !is.na(m)
  res[[i]] = data.frame(sapply(1:ncol(cells[f.a,]), function(x) mean(cells[f.a,x])))
  colnames(res[[i]]) = "Average Expression"
  res[[i]]$Condition = wsn.integrated$Condition
  res[[i]]$Cluster = i
  res[[i]]$Celltype = wsn.integrated$Cell.type
}

#Or 
cells = as.matrix(cells)

i =1
res =list()
for (i in 1:length(unique(m.clusters$labels))){
  m = match(rownames(cells), m.clusters[m.clusters$labels == i,]$ids)
  f.a = !is.na(m)
  res[[i]] = reshape2::melt(cells[f.a,],value.name = 'Normalized Expression')
  colnames(res[[i]]) = c("gene","cell","Normalized Expression")
  m = match(res[[i]]$cell, names(wsn.integrated$orig.ident))
  f.a = !is.na(m)
  f.t = m[f.a]
  res[[i]]$Condition = rep(0, nrow(res[[i]]))
  res[[i]]$Celltype = rep(0, nrow(res[[i]]))
  res[[i]][f.a,]$Condition = wsn.integrated$Condition[f.t] %>% as.character()
  res[[i]]$Cluster = i
  res[[i]][f.a,]$Celltype = wsn.integrated$Cell.type[f.t] %>% as.character()
}


head(res[[i]])
res1 = do.call(rbind.data.frame, res)
res1$Condition = factor(res1$Condition, levels = c('Mock','Wildtype','3M','dNS1'))
res1 = res1[res1$Celltype %in% lvls,]
#res1 = res1[res1$`Normalized Expression` > 0,]
ht = ggplot(data =res1,aes(x=Condition, y =`Average Expression`,fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = colour_sample, 
                    name="Condition") +
  facet_grid(Celltype ~ Cluster) +
  stat_summary(fun=mean, colour="black", geom="text", size = 1.8,show.legend = FALSE, 
               vjust = -4,aes(label=paste(round(..y.., digits=2)))) +
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05),
              textsize=2, 
              comparisons = list(c("Mock", "Wildtype"),
                                 c("Wildtype", "3M"),
                                 c("Wildtype", "dNS1")),
              show.legend = T, 
              test = "t.test",
              step_increase = 0.1,
              margin_top = -0.3) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = 2,size = 0.3, fill = NA), 
        panel.spacing.x = unit(0.2,"cm"),
        strip.text.y = element_text(size = 5),
        axis.text.x.bottom = element_text(angle = 45,hjust = 0.95))

g <- ggplot_gtable(ggplot_build(ht))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(m.clusters$colors)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
stripr <- which(grepl('strip-r', g$layout$name))
fills <- col_fun
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
pdf("E:/Thesis figures/Chp3-Figure 5b.pdf", width = 13, height = 20, paper = 'a4',useDingbats = F)
grid::grid.draw(g)
dev.off()

####################Getting UMI Totals per celltype####################

df = matrix(0, nrow = length(lvls), ncol = length(unique(wsn.integrated$Condition)))
colnames(df) = unique(wsn.integrated$Condition)
rownames(df) = lvls
geneUMI = list()
k =1
for (i in 1:nrow(df)) {
  for (j in 1:ncol(df)) {
    #For total UMI count per sample per cell type/cluster
    df[i,j] = sum(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[j]])
    #For cell number per sample per cell type/cluster
    #df[i,j] = length(wsn.integrated$orig.ident[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$orig.ident == colnames(df)[j]]) %>% as.numeric()
  }
}

#Calculate normalized gene expression
k=1
i=1

for (i in 1:nrow(df)) {
  cells = names(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[1]])
  if(length(cells) > 0 ){
    geneUMI[[k]] =  data.frame((rowSums(as.matrix(wsn.integrated[["RNA"]]@counts[,cells]))/df[i,1])*100000)
  } else {
    geneUMI[[k]] = data.frame(rep(0,nrow(wsn.integrated[["RNA"]]@counts)))
    rownames(geneUMI[[k]]) = rownames(wsn.integrated[["RNA"]]@counts)
  }
  colnames(geneUMI[[k]]) = paste0(colnames(df)[1],"-",rownames(df)[i])
  
  k = k+1
  cells = names(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[2]])
  if(length(cells) > 0 ){
    geneUMI[[k]] =  data.frame((rowSums(as.matrix(wsn.integrated[["RNA"]]@counts[,cells]))/df[i,2])*100000)
  } else {
    geneUMI[[k]] = data.frame(rep(0,nrow(wsn.integrated[["RNA"]]@counts)))
    rownames(geneUMI[[k]]) = rownames(wsn.integrated[["RNA"]]@counts)
  }
  colnames(geneUMI[[k]]) = paste0(colnames(df)[2],"-",rownames(df)[i])
  
  k = k+1
  cells = names(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[3]])
  if(length(cells) > 0 ){
    geneUMI[[k]] =  data.frame((rowSums(as.matrix(wsn.integrated[["RNA"]]@counts[,cells]))/df[i,3])*100000)
  } else {
    geneUMI[[k]] = data.frame(rep(0,nrow(wsn.integrated[["RNA"]]@counts)))
    rownames(geneUMI[[k]]) = rownames(wsn.integrated[["RNA"]]@counts)
  }
  colnames(geneUMI[[k]]) = paste0(colnames(df)[3],"-",rownames(df)[i])
  
  k = k+1
  cells = names(wsn.integrated$nCount_RNA[wsn.integrated$Cell.type == rownames(df)[i] & wsn.integrated$Condition == colnames(df)[4]])
  if(length(cells) > 0 ){
    geneUMI[[k]] =  data.frame((rowSums(as.matrix(wsn.integrated[["RNA"]]@counts[,cells]))/df[i,4])*100000)
  } else {
    geneUMI[[k]] = data.frame(rep(0,nrow(wsn.integrated[["RNA"]]@counts)))
    rownames(geneUMI[[k]]) = rownames(wsn.integrated[["RNA"]]@counts)
  }
  colnames(geneUMI[[k]]) = paste0(colnames(df)[4],"-",rownames(df)[i])
  
  k = k+1
  
}

genexp = do.call(cbind.data.frame, geneUMI)
save(genexp, file = "norm_exp.rdata")

####################################################################################################

load("norm_exp.rdata")
tcounts = genexp

m = match(rownames(tcounts), m.clusters$ids)
f.a = !is.na(m)
tcounts = tcounts[f.a,]

##Plot - Heatmap

clustmat = matrix(0, nrow = length(m.clusters$ids), ncol = dim(tcounts)[[2]]) %>% as.data.frame()
rownames(clustmat) = as.character(m.clusters$ids)
colnames(clustmat) = colnames(tcounts)
m = match(rownames(clustmat), rownames(tcounts))
f.a = !is.na(m)
f.t = m[f.a]
clustmat[f.a,] = tcounts[f.t,]
m = match(m.clusters$ids, rownames(clustmat))
clustmat = clustmat[m,]

COL = c("Mock" = "tan1","Wildtype" = "firebrick1", "dNS1"="blue2", "3M" = "limegreen")

samples = as.data.frame(colnames(clustmat))
samples$Conditions = matrix(unlist(strsplit(as.character(samples$`colnames(clustmat)`),"-")), byrow = T, ncol = 2)[,1]
samples$Celltype = matrix(unlist(strsplit(as.character(samples$`colnames(clustmat)`),"-")), byrow = T, ncol = 2)[,2]
rownames(samples) = samples$`colnames(clustmat)`
samples = samples[,-1]
samples$Conditions = factor(samples$Conditions, levels = c("Mock", "Wildtype","3M", "dNS1"))
samples$Celltype = factor(samples$Celltype, levels = lvls)
cols = as.character(unique(m.clusters$colors))
names(cols) = as.character(unique(m.clusters$labels))
col_fun1 = cividis(100)
prot = m.clusters$labels %>% as.data.frame()
rownames(prot) = m.clusters$ids
colnames(prot) ="Clusters"
prot = prot[rownames(clustmat),,drop = F]

mat_scaled =t(scale(t(clustmat)))
mat_scaled = na.omit(mat_scaled)

m = match(rownames(prot), rownames(mat_scaled))
f.a = !is.na(m)
prot = prot[f.a,,drop = F]
prot$Clusters = factor(prot$Clusters, levels = 1:5)

ha = rowAnnotation(df = prot,
                   col = list(Clusters = cols),
                   show_annotation_name = F,
                   border = T,
                   annotation_legend_param = list(grid_height = unit(0.5, "cm"),
                                                  title_gp=gpar(fontsize=10, fontface="bold"),
                                                  labels_gp = gpar(fontsize = 8),
                                                  legend_direction = "horizontal",
                                                  nrow = 1)
)
hb = columnAnnotation(df = samples,
                      col = list(Conditions = COL,
                                 Celltype = col_fun),
                      show_annotation_name = F,
                      show_legend = T,
                      border = T,
                      annotation_legend_param = list(grid_height = unit(0.5, "cm"),
                                                     title_gp=gpar(fontsize=10, fontface="bold"),
                                                     labels_gp = gpar(fontsize = 8),
                                                     legend_direction = "horizontal",
                                                     nrow = 4)
)

#Scaling:

ht = Heatmap(mat_scaled,
             name = "Z-Score",
             col = viridis(200),
             column_title_gp = gpar(fontsize = 10, fontface = "bold"),
             row_title_gp = gpar(fontsize = 10, fontface = "bold" ),
             heatmap_legend_param = list(
               at = c(-2,0,2, 4,6),
               color_bar="continuous",
               labels_gp = gpar(fontsize = 8),
               legend_width = unit(2, "cm"),
               title_position = "topleft",
               legend_direction = 'horizontal',
               title_gp=gpar(fontsize=10, fontface="bold"),
               border = T),
             top_annotation = hb,
             left_annotation = ha,
             cluster_columns = F,
             column_dend_side = "top",
             row_split = prot,
             show_column_dend = F,
             cluster_rows = T,
             row_names_gp = gpar(fontsize = 8),
             show_row_names = F,
             show_column_names = F,
             row_names_rot = 90,
             row_dend_side = "right", 
             row_names_side = "left",
             column_names_side = "bottom",
             column_names_gp = gpar(fontsize = 8),
             column_names_rot = 45,
             height = unit(15,"cm"),
             width = unit(15,"cm"))

ht = grid::grid.grabExpr(draw(ht,heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE))

pdf("Thesis figures/Chp3-Figure 5e.pdf", height = 9, width = 7, paper = 'a4', useDingbats = F)
grid.draw(ht)
dev.off()


tcounts = as.matrix(tcounts)
tcounts = reshape2::melt(tcounts, value.name ="Relative Expression")
colnames(tcounts)[1:2] = c('gene','sample')
tcounts$Condition = matrix(unlist(strsplit(as.character(tcounts$sample), "-")),ncol=2, byro=T)[,1]
tcounts$Celltype = matrix(unlist(strsplit(as.character(tcounts$sample), "-")),ncol=2, byro=T)[,2]
tcounts$Cluster = rep(0, nrow(tcounts))
m =match(tcounts$gene,m.clusters$ids)
f.a =!is.na(m)
f.t =m[f.a]
tcounts[f.a,]$Cluster = m.clusters[f.t,]$labels
tcounts$Condition = factor(tcounts$Condition, levels = c('Mock','Wildtype','3M','dNS1'))
tcounts = tcounts[tcounts$Celltype %in% lvls,]

ht =ggplot(data =tcounts[tcounts$Celltype %in% lvls[1:8],],aes(x=Condition, y =log(`Relative Expression`+1),fill = Condition)) +
  geom_violin() +
  scale_fill_manual(values = colour_sample, 
                    name="Condition") +
  facet_grid(Celltype ~ Cluster) +
  labs(x =NULL) +
  stat_summary(fun=mean, colour="black", geom="text", size = 2.5,show.legend = FALSE, 
               vjust = -11,aes(label=paste(round(..y.., digits=2)))) +
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05) , 
              textsize=2, 
              comparisons = list(c("Mock", "Wildtype"),
                                 c("Wildtype", "3M"),
                                 c("Wildtype", "dNS1")),
              show.legend = T,
              vjust = 0.2,
              test = "wilcox.test",
              step_increase = 0.1,
              margin_top = 0.2) +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))
g <- ggplot_gtable(ggplot_build(ht))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(m.clusters$colors)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
stripr <- which(grepl('strip-r', g$layout$name))
fills <- col_fun[1:8]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
g1 = g

ht =ggplot(data =tcounts[tcounts$Celltype %in% lvls[9:16],],aes(x=Condition, y =log(`Relative Expression`+1),fill = Condition)) +
  geom_violin() +
  scale_fill_manual(values = colour_sample, 
                    name="Condition") +
  facet_grid(Celltype ~ Cluster) +
  labs(x =NULL) +
  stat_summary(fun=mean, colour="black", geom="text", size = 2.5,show.legend = FALSE, 
               vjust = -11,aes(label=paste(round(..y.., digits=2)))) +
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05) , 
              textsize=2, 
              comparisons = list(c("Mock", "Wildtype"),
                                 c("Wildtype", "3M"),
                                 c("Wildtype", "dNS1")),
              show.legend = T,
              vjust = 0.2,
              test = "wilcox.test",
              step_increase = 0.1,
              margin_top = 0.2) +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8))
g <- ggplot_gtable(ggplot_build(ht))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(m.clusters$colors)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
stripr <- which(grepl('strip-r', g$layout$name))
fills <- col_fun[9:16]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
g2 = g

pdf("E:/Thesis figures/Chp3-Figure 5a.pdf", width = 7, height = 14, paper = 'a4',useDingbats = F)
plot_grid(g1,labels = c(LETTERS[1]), label_fontfamily = 'serif', label_size = 18)
dev.off()
pdf("E:/Thesis figures/Chp3-Figure 5b.pdf", width = 7, height = 14, paper = 'a4',useDingbats = F)
plot_grid(g2,labels = c(LETTERS[2]), label_fontfamily = 'serif', label_size = 18)
dev.off()

