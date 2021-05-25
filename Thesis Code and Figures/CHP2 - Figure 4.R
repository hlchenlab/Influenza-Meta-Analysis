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


orthos = read.csv("orthos.csv", stringsAsFactors = F)
orthos = orthos[orthos$Ferret.homology.type == "ortholog_one2one",]
orthos = orthos[orthos$Golden.Hamster.homology.type == "ortholog_one2one",]
orthos = orthos[orthos$Mouse.homology.type == "ortholog_one2one",]

load("COVID/Clinical studies/cross.species.genesets.Rdata")

load("mouse_gene_annotations_v100.Rdata")
load("coexpression network/Mouse Co-expression network SB/agg.rank.ENSEMBL_IDs.ranked.Rdata")
load("gene_annotations_v29.Rdata")
load("E:/coexpression network/agg.reranked.Rdata")

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

#Get gene orthologue ids
genes = rownames(genesets.up[rowSums(genesets.up) >=6,])
m = match(orthos$Gene.stable.ID, genes)
f.a = !is.na(m)
genes = orthos[f.a,]
m = match(genes$Gene.stable.ID, rownames(agg.rank))
f.a = !is.na(m)
genes = genes[f.a,]
m = match(genes$Mouse.gene.stable.ID, rownames(reranked))
f.a = !is.na(m)
genes = genes[f.a,]

#Get both networks
local = reranked[genes$Mouse.gene.stable.ID, genes$Mouse.gene.stable.ID]
localm = local[!duplicated(rownames(local)), !duplicated(colnames(local))]
dim(localm)
local = agg.rank[genes$Gene.stable.ID, genes$Gene.stable.ID]
localh = local[!duplicated(rownames(local)), !duplicated(colnames(local))]
dim(localh)

##Plot 1: Co-expression heatmap
temp = localh

ids = rownames(temp)
consTree = hclust( dist(temp, method = "euclidean"), method = "average")
consDend = as.dendrogram(consTree)

ggdendrogram(consTree,rotate = TRUE,labels = F, theme_dendro = FALSE) +
  geom_hline(aes(yintercept = max(dist(temp))*0.5)) +
  theme_classic()


unmergedLabels3 = cutree(consTree,h = max(dist(temp))*0.5) ######Change when needed

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
clustersh = as.data.frame(clust)

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
                   show_annotation_name = F,
                   border = T
)
hb = columnAnnotation(df = prot,
                      col = list(Clusters = cols),
                      annotation_name_side = "right",
                      show_legend = T,
                      show_annotation_name = F,
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
            height = unit(5,"cm"),
            width = unit(5,"cm")
)

ht = grid::grid.grabExpr(draw(ht,heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE))


#Plot 2: Coexpression boxplots

clusters = clustersh
local = localh
res = list()
for (i in 1:length(unique(clusters$unmergedLabels.mod))) {
  res[[i]] = as.numeric(local[clusters[clusters$unmergedLabels.mod == i,]$V1,clusters[clusters$unmergedLabels.mod == i,]$V1])
  res[[i]] = as.data.frame(res[[i]])
  res[[i]]$Cluster = i
}
res1 = do.call(rbind.data.frame, res)
colnames(res1)[1] = 'Co-expression'
res1$Cluster = factor(res1$Cluster, levels = c(paste0(1:2)))
res1$Network = 'Human'
cols = unique(clusters[,2:3]) %>% data.frame(.,row.names = 1) 
cols$Cluster = paste0(rownames(cols),'H')
cols$Network = "Human"
resh = res1
colh = cols
clusters = clustersm
local = localm
res = list()
for (i in 1:length(unique(clusters$unmergedLabels.mod))) {
  res[[i]] = as.numeric(local[clusters[clusters$unmergedLabels.mod == i,]$V1,clusters[clusters$unmergedLabels.mod == i,]$V1])
  res[[i]] = as.data.frame(res[[i]])
  res[[i]]$Cluster = i
}
res1 = do.call(rbind.data.frame, res)
colnames(res1)[1] = 'Co-expression'
res1$Cluster = factor(res1$Cluster, levels = 1:max(res1$Cluster))
res1$Network = 'Mouse'
cols = unique(clusters[,2:3]) %>% data.frame(.,row.names = 1)
cols$Cluster = rownames(cols)
cols$Network = "Mouse"
c1 = ggplot(data = resh,aes(Cluster,`Co-expression`, fill = Cluster))+
  geom_boxplot(width = 0.5,show.legend = F,outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = colh$V3)
c1m = ggplot(data = res1,aes(Cluster,`Co-expression`, fill = Cluster))+
  geom_boxplot(width = 0.5,show.legend = F,outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = cols$V3)

#Plot 3 - DE prior ranks v clusters

clusters = clustersh
load("DE_prior.rdata")

gmat = matrix(0, ncol = 2, nrow = nrow(DE_list))
rownames(gmat) = rownames(DE_list)
m = match(rownames(gmat), clusters$V1)
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,1] = clusters[f.t,2]
gmat[,2] = DE_list[,4]
gmat = data.frame(gmat[gmat[,1] >= 1,], row.names = NULL)

colnames(gmat) = c('Cluster', 'DE prior')
gmat$Cluster = factor(gmat$Cluster, levels =  c(1:2))
gmat$`DE prior` = as.numeric(gmat$`DE prior`)
gmat$Network = "Human"
gmah = gmat

clusters = clustersm
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
gmat$Network = "Mouse"

c2 =ggplot() + 
  geom_boxplot(data = gmah, aes( Cluster,`DE prior`, fill = Cluster), 
               width = 0.5,show.legend = F,outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = colh$V3, 
                    name="Cluster")
c2m =ggplot() + 
  geom_boxplot(data = gmat, aes( Cluster,`DE prior`, fill = Cluster), 
               width = 0.5,show.legend = F,outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = cols$V3, 
                    name="Cluster")

##Plot Cluster overlaps and Jaccard plots:

clusters = clustersh
hnet = localh
consTree = hclust( dist(hnet, method = "euclidean"), method = "average")
consDend = as.dendrogram(consTree)
clusterids <- data.frame( genes = clusters$V1, labels = as.numeric(clusters$unmergedLabels.mod),
                          labels_unmerged =  clusters$unmergedLabels.mod,
                          colors = clusters$V3)

m = match(clusters$V1, rownames(hnet))

clust_hnet <- list(as.matrix(hnet), consTree, consDend, m, clusterids)
names(clust_hnet) <- c( "distance_matrix", "tree", "dendrogram", "order", "clusters")

clusters = clustersm
hnet = localm
m = match(clusters$V1, orthos$Mouse.gene.stable.ID)
f.a = !is.na(m)
f.t = m[f.a]
clusters[f.a,]$V1 = orthos[f.t,]$Gene.stable.ID
m = match(rownames(hnet), orthos$Mouse.gene.stable.ID)
f.a = !is.na(m)
f.t = m[f.a]
rownames(hnet)[f.a] = orthos[f.t,]$Gene.stable.ID
colnames(hnet)[f.a] = orthos[f.t,]$Gene.stable.ID

consTree = hclust( dist(hnet, method = "euclidean"), method = "average")
consDend = as.dendrogram(consTree)
clusterids <- data.frame( genes = clusters$V1, labels = as.numeric(clusters$unmergedLabels.mod),
                          labels_unmerged =  clusters$unmergedLabels.mod,
                          colors = clusters$V3)

m = match(clusters$V1, rownames(hnet))

clust_mnet <- list(as.matrix(hnet), consTree, consDend, m, clusterids)
names(clust_mnet) <- c( "distance_matrix", "tree", "dendrogram", "order", "clusters")

#Use Jaccard.R for heatmap generation.

human = plot_grid(ht,c1,c2, ncol = 3,nrow = 1,axis = 'top',rel_widths = c(1.5,1,1), labels = LETTERS[1:3], label_fontfamily = 'serif', label_size = 18)
mouse = plot_grid(htm,c1m,c2m, ncol = 3,nrow = 1,axis = 'top',rel_widths = c(1.5,1,1), labels = LETTERS[4:6], label_fontfamily = 'serif', label_size = 18)
blk = plot.new()
ott = plot_grid(blk,htj,blk, ncol = 3,nrow = 1,axis = 'bottom',rel_widths = c(1.3,1,1.3), labels = c("",LETTERS[7],""), label_fontfamily = 'serif', label_size = 18)
pdf("Thesis figures/Chp2-Figure 4.pdf", height = 8,width = 9, useDingbats =  F)
plot_grid(human, mouse, ott, nrow = 3, rel_heights = c(1,1,0.9))
dev.off()

save(clust_hnet,clust_mnet, localh, localm, clustersh, clustersm,file = 'Coexpression-Sars-CoV-2-human-network and datasets.rdata')
