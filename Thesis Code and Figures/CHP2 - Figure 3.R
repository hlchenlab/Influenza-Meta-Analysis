#Packages:

if( !require("dynamicTreeCut")){
  BiocManager::install("dynamicTreeCut")
}
if( !require("viridisLite")){
  BiocManager::install("viridisLite")
}
if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}
if( !require("cowplot")){
  BiocManager::install("cowplot")
}
if( !require("gplots")){
  BiocManager::install("gplots")
}
if( !require("dplyr")){
  BiocManager::install("dplyr")
}
if( !require("stringr")){
  BiocManager::install("stringr")
}
if( !require("ggplotify")){
  BiocManager::install("ggplotify")
}
if( !require("EGAD")){
  BiocManager::install("EGAD")
}
if( !require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}
if( !require("tibble")){
  BiocManager::install("tibble")
}
if( !require("ggdendro")){
  install.packages("ggdendro")
}


# Set up:
setwd('E:/')

load("recurrence_matrix_MOUSE.Rdata")
load("mouse_gene_annotations_v100.Rdata")
load("coexpression network/Mouse Co-expression network SB/agg.rank.ENSEMBL_IDs.ranked.Rdata")

m = match(rownames(reranked),rownames(genesets.up[rowSums(genesets.up) >=3,]))
f.a = !is.na(m)
local = reranked[f.a,f.a]
local = local[!duplicated(rownames(local)), !duplicated(colnames(local))]
dim(local)

##Plot 1: Co-expression heatmap
temp = local

ids = rownames(temp)
consTree = hclust( dist(temp, method = "euclidean"), method = "average")

ggdendrogram(consTree,rotate = TRUE,labels = F, theme_dendro = FALSE) +
  geom_hline(aes(yintercept = max(dist(temp))*0.29)) +
  theme_classic()

consDend = as.dendrogram(consTree)

unmergedLabels3 = cutree(consTree,h = max(dist(temp))*0.29) ######Change when needed

unmergedColors3 = sample(cividis( max(unmergedLabels3) + 1 ))[ as.numeric(unmergedLabels3)+1]

# Label and count modules
i.prev = ""
ki = 1
ji = 0
unmergedLabels.mod = as.numeric(unmergedLabels3[consTree$order]) * 0
for( ii in as.numeric(unmergedLabels3[consTree$order]) ){
  if( ii == i.prev){
    unmergedLabels.mod[ki] = ji
  } else {
    i.prev = ii
    ji = ji + 1
    unmergedLabels.mod[ki] = ji
    
  }
  ki = ki + 1
}
clust =  cbind(ids[consTree$order], unmergedLabels.mod, unmergedColors3[consTree$order] )
clusters = as.data.frame(clust)

cols = as.data.frame(unique(clust[,2:3]),row.names = 1) %>% tibble::deframe()
col_fun = cividis(100)
prot = data.frame(clust[,1:2])
rownames(prot) = prot$V1
prot = prot[,-1, drop = F]
colnames(prot) ="Clusters"
prot[,1] = factor(prot[,1], levels = unique(prot[,1])[order(as.numeric(unique(prot[,1])), decreasing = F)])
prot = prot[rownames(temp),,drop = F]
ha = rowAnnotation(df = prot,
                   col = list(Clusters = cols),
                   annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                   show_legend = F,
                   border = T
)
hb = columnAnnotation(df = prot,
                      col = list(Clusters = cols),
                      annotation_name_side = "right",
                      show_legend = T,
                      annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
                      annotation_legend_param = list(title_gp=gpar(fontsize=10, fontface="bold"),
                                                     legend_height = unit(2, "cm"),
                                                     border = "black"),
                      annotation_label = NULL,
                      border = T)
ht =Heatmap(temp, 
            name = "Ranked\nco-expression",  
            col=viridis(100), 
            border = T,
            right_annotation = ha,
            heatmap_legend_param = list(at = c(0,0.5,1),
                                        #legend_height = unit(2, "cm"),
                                        legend_width = unit(2, "cm"),
                                        title_position = "topleft",
                                        border = "black",
                                        labels_gp = gpar(fontsize = 8),
                                        title_gp=gpar(fontsize=10, fontface="bold"),
                                        legend_direction = "vertical"),
            top_annotation = hb,
            column_title_gp = gpar(fontsize = 12),
            column_title = "Genes",
            column_title_side = "bottom",
            column_dend_side = "top",
            column_dend_gp = gpar(col = "black"),
            column_dend_height = unit(1, "cm"),
            row_title = "Genes",
            use_raster = T,
            row_title_side = "left",
            row_title_gp = gpar(fontsize = 12),
            row_dend_side = "left",
            row_dend_gp = gpar(col = "black"),
            show_column_dend = F,
            show_row_dend = F,
            show_row_names = F,
            show_column_names = F,
            cluster_columns = consDend,
            cluster_rows = rev(consDend),
            height = unit(8,"cm"),
            width = unit(8,"cm")
)

ht = grid::grid.grabExpr(draw(ht,heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE))


#Plot 2: Coexpression boxplots

res = list()
for (i in 1:length(unique(clusters$unmergedLabels.mod))) {
  res[[i]] = as.numeric(local[clusters[clusters$unmergedLabels.mod == i,]$V1,clusters[clusters$unmergedLabels.mod == i,]$V1])
  res[[i]] = as.data.frame(res[[i]])
  res[[i]]$Cluster = i
}
res1 = do.call(rbind.data.frame, res)
colnames(res1)[1] = 'Co-expression'
res1$Cluster = factor(res1$Cluster, levels = 1:max(res1$Cluster))
cols = unique(clusters[,2:3]) %>% data.frame(.,row.names = 1) %>% t

c1 = ggplot()+
  geom_boxplot(data = res1, aes(Cluster,`Co-expression`, fill = Cluster), width = 0.5,
               show.legend = F,outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = cols, 
                    name="Cluster") 


#Plot 3 - DE prior ranks v clusters

load("Mouse_DE_prior.rdata")
gmat = matrix(0, ncol = 2, nrow = nrow(DE_list))
rownames(gmat) = rownames(DE_list)
m = match(rownames(gmat), clusters$V1)
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,1] = clusters[f.t,2]
gmat[,2] = DE_list[,4] 
gmat = data.frame(gmat[gmat[,1] >= 1,], row.names = NULL)
gmat[,2] = as.numeric(gmat[,2])
gmat[,2] = gmat[,2]/max(gmat[,2])

colnames(gmat) = c('Cluster', 'DE prior')
gmat$Cluster = factor(gmat$Cluster, levels = 1:max(gmat$Cluster))
gmat$`DE prior` = as.numeric(gmat$`DE prior`)

c2 =ggplot() + 
  geom_boxplot(data = gmat, aes( Cluster,`DE prior`, fill = Cluster), 
               width = 0.5,show.legend = F,outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = cols, 
                    name="Cluster") 


#Plot 4 - global v cluster node degrees

bes = list()
for (k in 1:length(unique(clusters$unmergedLabels.mod))) {
  #local node degree
  network = agg.rank[rownames(agg.rank) %in% clusters[clusters$unmergedLabels.mod == k,]$V1, colnames(agg.rank) %in% clusters[clusters$unmergedLabels.mod == k,]$V1]
  n = length(unique(clusters[clusters$unmergedLabels.mod == k,]$V1))
  local = (colSums(network))/(nrow(network))
  
  #global node degree
  m = match(rownames(agg.rank), clusters[clusters$unmergedLabels.mod == k,]$V1)
  f.a = is.na(m)
  network = agg.rank[f.a, clusters[clusters$unmergedLabels.mod == k,]$V1]
  global = (colSums(network))/(nrow(network))
  bes[[k]] = cbind.data.frame(local, global)
  bes[[k]]$Cluster = k
}

bes1 = do.call(rbind.data.frame, bes)
bes1$Cluster = factor(bes1$Cluster , levels = 1:max(bes1$Cluster))
ggplot(data = bes1) +
  geom_point(aes(global,local)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  theme_classic() +
  facet_grid(~Cluster)

#Plot 5 - Genetypes per cluster

i=1
df = list()
for (i in 1:length(unique(clusters$unmergedLabels.mod))) {
  genes = data.frame(clusters[clusters$unmergedLabels.mod ==i,]$V1)
  colnames(genes) = "Genes"
  m = match(attrm$ensemblID, genes$Genes)
  f.a = !is.na(m)
  f.t = m[f.a]
  genes$Genetypes = rep("Unknown", nrow(genes))
  genes[f.t,]$Genetypes = attrm[f.a,]$type
  df[[i]] = (table(genes$Genetypes)/nrow(genes)) %>% as.data.frame()
  colnames(df[[i]])[2] = i
}

tmat = matrix(0, nrow = length(unique(attrm$type))+1, ncol = length(df))
rownames(tmat) = c(names(table(attrm$type)[order(table(attrm$type), decreasing = T)]),"Unknown")
colnames(tmat) = 1:length(unique(clusters$unmergedLabels.mod))
for(i in 1:length(df)){
  m = match(df[[i]]$Var1, rownames(tmat))
  f.a =!is.na(m)
  f.t =m[f.a]
  tmat[f.t,i] = df[[i]][f.a,2]
}
tmat = tmat[rowSums(tmat)>0,]
df = reshape2::melt(tmat, value.name = c("Percentage"))
colnames(df)[1:2] = c("Genetype","Cluster")
df$Cluster  = factor(df$Cluster , levels = 1:max(df$Cluster))

c3 =ggplot(data = df, aes(x = Cluster, y = Percentage, fill =Genetype)) +
  geom_bar(position="stack", stat="identity", color ='black') +
  scale_fill_brewer(palette = 'PuBu',
                    direction = -1) +
  theme_classic() +
  labs(x = "Cluster", y = "Percentage") +
  theme(legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.title = element_text(size=10, face = "bold"),
        legend.text = element_text(size=10))

gtps = plot_grid(c1,c2,c3, ncol = 3,axis = 'bottom',rel_widths = c(1,1,1.7), labels = LETTERS[2:4], label_fontfamily = 'serif', label_size = 18)

pdf("Thesis figures/Chp2-Figure 3.pdf", height = 7,width = 8, useDingbats =  F)
plot_grid(ht,gtps, nrow = 2, rel_heights = c(2,1), labels = LETTERS[1], label_fontfamily = 'serif', label_size = 18)
dev.off()

save(df, gmat, res1, local, clusters,file = 'Coexpression-influenza-mouse-network and datasets.rdata')
