#Packages:

library(dynamicTreeCut)
library(viridisLite)
library(ggplot2)
library(gplots)
library(dplyr)
library(stringr)
library(ggplotify)
library(EGAD)
library(ComplexHeatmap)
if( !require("cluster")){
  BiocManager::install("cluster")
}
if( !require("pvclust")){
  BiocManager::install("pvclust")
}

# Set up:

source("http://peterhaschke.com/Code/multiplot.R") #load multiplot function
load("E:/genesets_for_human_recurrence.rdata")
load("E:/gene_annotations_v29.Rdata")
allgenes = read.csv("E:/allgenes.csv", header = T,stringsAsFactors = F)
load("agg.reranked.Rdata")

network = agg.rank
m = match(rownames(network),rownames(genesets.up[rowSums(genesets.up) >=4,]))
f.a = !is.na(m)
local = network[f.a,f.a]
dim(local)

# The Clustering Function: heatmap and cluster table included


get_cluster_ids <- function( temp, filtMin=2){
  library(ComplexHeatmap)
  if( !require("tibble")){
    BiocManager::install("tibble")
  }
  if( !require("dynamicTreeCut")){
    BiocManager::install("dynamicTreeCut")
  }
  ids = rownames(temp)
  consTree = hclust( dist(temp, method = "euclidean"), method = "complete")
  consDend = as.dendrogram(consTree)
  #unmergedLabels3 = cutreeDynamic(dendro = consTree, distM = temp,
  #                                deepSplit = 0, minAbsSplitHeight = max(dist(temp))*0.60 , minClusterSize = 400,                                  # Change if needed.
  #                               pamRespectsDendro = FALSE, respectSmallClusters = TRUE);
  
  unmergedLabels3 = cutree(consTree, 
                           h = max(dist(temp))*0.7)
  
  unmergedColors3 = sample(plasma( max(unmergedLabels3) + 1 ))[ as.numeric(unmergedLabels3)+1]
  
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
  clust =  cbind(ids[consTree$order], unmergedLabels.mod, unmergedColors3[consTree$order], consTree$order )
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
                     annotation_name_gp = gpar(fontsize = 6, fontface = "bold"),
                     show_legend = F
  )
  hb = columnAnnotation(df = prot,
                        col = list(Clusters = cols),
                        annotation_name_side = "right",
                        show_legend = T,
                        annotation_name_gp = gpar(fontsize = 6, fontface = "bold"),
                        annotation_legend_param = list(title_gp=gpar(fontsize=6, fontface="bold"),
                                                       legend_direction = "horizontal",
                                                       legend_height = unit(2, "cm"),
                                                       labels_gp = gpar(fontsize = 8),
                                                       nrow = 1,
                                                       legend_width = unit(2, "cm")),
                        annotation_label = NULL,
                        border = T)
  ht = Heatmap(temp, 
               name = "Co-expression",  
               col=viridis(100), 
               border = T,
               right_annotation = ha,
               heatmap_legend_param = list(at = c(0,0.5,1),
                                           legend_height = unit(2, "cm"),
                                           legend_width = unit(2, "cm"),
                                           title_position = "topleft",
                                           border = "black",
                                           labels_gp = gpar(fontsize = 8),
                                           title_gp=gpar(fontsize=6, fontface="bold"),
                                           legend_direction = "horizontal"),
               top_annotation = hb,
               column_title_gp = gpar(fontsize = 12, fontface = "bold"),
               column_title = "Genes",
               column_title_side = "bottom",
               column_dend_side = "top",
               column_dend_gp = gpar(col = "black"),
               column_dend_height = unit(1, "cm"),
               row_title = "Genes",
               use_raster = T,
               row_title_side = "right",
               row_title_gp = gpar(fontsize = 12, fontface = "bold"),
               row_dend_side = "left",
               row_dend_gp = gpar(col = "black"),
               row_dend_width = unit(1, "cm"),
               show_row_names = F,
               show_column_names = F,
               cluster_columns = consDend,
               cluster_rows = rev(consDend),
               height = unit(10,"cm"),
               width = unit(10,"cm")
  )
  pdf("Coexpression network.pdf", height = 7,width = 7, useDingbats =  F)
  ht
  dev.off()
  print(consDend)
  return(clust)
} 
clusters = as.data.frame(get_cluster_ids(local))
