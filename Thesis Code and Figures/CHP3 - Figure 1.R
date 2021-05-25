#Packages
if( !require("gplots")){
  BiocManager::install("gplots")
}
if( !require("dplyr")){
  BiocManager::install("dplyr")
}
if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}
if( !require("DESeq2")){
  BiocManager::install("DESeq2")
}
if( !require("RColorBrewer")){
  BiocManager::install("RColorBrewer")
}
if( !require("limma")){
  BiocManager::install("limma")
}
if( !require("pheatmap")){
  BiocManager::install("pheatmap")
}
if( !require("gridExtra")){
  BiocManager::install("gridExtra")
}
if( !require("grid")){
  BiocManager::install("grid")
}
if( !require("fields")){
  BiocManager::install("fields")
}
if( !require("ggplotify")){
  BiocManager::install("ggplotify")
}
if( !require("EnhancedVolcano")){
  BiocManager::install("EnhancedVolcano")
}
if( !require("cowplot")){
  BiocManager::install("cowplot")
}
if( !require("heatmap3")){
  BiocManager::install("heatmap3")
}
if( !require("circlize")){
  BiocManager::install("circlize")
}
if( !require("viridis")){
  BiocManager::install("viridis")
}
if( !require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}
if( !require("tidyr")){
  BiocManager::install("tidyr")
}
if( !require("data.table")){
  BiocManager::install("data.table")
}


#Set directory:
setwd("D:/")
load("HKUStudy.DE.Rdata")
load('cross-species.clusters.rdata')

dds_obj = dds_all
COL = c("Mock" = "#ffaf00","Wildtype" = "#d70000", "dNS1"="#5f5fff", "3M" = "#00af00")

#Plot 1 - Volplots for conditions

ids = as.data.frame(colData(dds_all)) %>% .[,'infection'] %>% unique() %>% as.character() %>% .[-1]
mylist = list()
for (i in 1:length(ids)) {
  id = ids[i]
  res <- results(dds_all, contrast = c("infection", id, "Mock"), pAdjustMethod = "BH")
  res <- na.omit(res, cols = c("log2FoldChange", "padj"))
  mylist[[i]] = as.data.frame(res)
  mylist[[i]]$Study = id
  mylist[[i]]$genes = rownames(mylist[[i]])
}  
mylist = do.call(rbind.data.frame, mylist)
mylist$Study = factor(mylist$Study, levels = ids)

#Custom colours to highlight cluster 6
keyvals.colour = rep('black', nrow(mylist))
keyvals.colour[(mylist$log2FoldChange >= 1.5 | mylist$log2FoldChange <= -1.5 ) & mylist$pvalue <= 0.05] = 'red3'
keyvals.colour[(mylist$log2FoldChange >= 1.5 | mylist$log2FoldChange <= -1.5 ) & mylist$pvalue > 0.05] = 'green3'
keyvals.colour[mylist$log2FoldChange > -1.5 & mylist$log2FoldChange < 1.5 & mylist$pvalue <= 0.05] = 'blue'
keyvals.colour[mylist$genes %in% h.clusters[h.clusters$labels_unmerged== 6,]$genes] = 'gold'

names(keyvals.colour)[keyvals.colour == 'gold'] <- 'Cluster 6'
names(keyvals.colour)[keyvals.colour == 'black'] <- 'Not significant'
names(keyvals.colour)[keyvals.colour == 'green3'] <- '-1.5 > Log2FC > 1.5'
names(keyvals.colour)[keyvals.colour == 'red3'] <- 'Significant'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'p <= 0.05'
names(keyvals.colour) = factor(names(keyvals.colour) , levels = c("Not significant",
                                                                  "-1.5 > Log2FC > 1.5",
                                                                  "p <= 0.05",
                                                                  "Significant",
                                                                  'Cluster 6') )
myplots = EnhancedVolcano(mylist,
                          lab = NA,
                          selectLab = NA,
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlab = expression(paste('Log'[2],'FC')),
                          ylab = expression(paste('-Log'[10],'P')),
                          pointSize = 0.5,
                          pCutoff = 0.05,
                          FCcutoff = 1.5,
                          title = NULL,
                          cutoffLineType = "dashed",
                          raster = TRUE,
                          subtitle = NULL,
                          colCustom = keyvals.colour,
                          legendLabels = factor(names(keyvals.colour) , levels = c("Not significant",
                                                                                   "-1.5 > Log2FC > 1.5",
                                                                                   "p <= 0.05",
                                                                                   "Significant",
                                                                                   'Cluster 6') ),
                          legendIconSize = 6,
                          caption = NULL,
                          colAlpha = 1,
                          border = 'full',
                          legendPosition = 'none') +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.direction = 'vertical',
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Study, ncol = 3) +
  guides(colour=guide_legend(ncol=2,byrow=F))

g1 <- ggplot_gtable(ggplot_build(myplots))

stripr <- which(grepl('strip-t', g1$layout$name))

fills <- COL[c(2,4,3)]
k =1
for (i in stripr) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#Plot2 - Expression boxplots

genes = as.character(h.clusters$genes)
m = match(genes, rownames(dds_obj))
f.a =!is.na(m)
genes = genes[f.a]

tcounts <- t(log2((counts(dds_obj[genes, ], normalized=TRUE, replaced=FALSE)+1))) %>%
  merge(colData(dds_obj), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(genes)+1):ncol(.))
tcounts = tcounts %>% 
  dplyr::select(Row.names, infection, gene, expression)
test = list()
##Averaging genes in clusters per sample:
for (i in 1:length(unique(h.clusters$labels_unmerged))) {
  samps = tcounts[tcounts$gene %in% h.clusters[h.clusters$labels_unmerged == i,]$genes,]
  trs = list()
  for (j in 1:length(unique(samps$Row.names))) {
    bd = samps[samps$Row.names == unique(samps$Row.names)[j],]
    md = mean(bd$expression)
    trs[[j]] = cbind(as.character(unique(bd$Row.names)), as.character(unique(bd$infection)),i,md)
  }
  test[[i]] = do.call(rbind.data.frame,trs)
}
test = do.call(rbind.data.frame, test)
colnames(test) = c("SampleID",'Condition', 'Cluster',"expression")
test$expression = as.numeric(test$expression)
test$Condition = factor(test$Condition, levels = c("Mock","Wildtype", "3M","dNS1"))
test$Cluster = factor(test$Cluster, levels = c(1:6))
m = match(test$Condition, names(COL)) 
f.a =!is.na(m)
f.t = m[f.a]
test$colours = rep(0, nrow(test))
test[f.a,]$colours = COL[f.t]
COL = c("Mock" = "#ffaf00","Wildtype" = "#d70000", "dNS1"="#5f5fff", "3M" = "#00af00")
p = ggplot(test, aes(Condition, expression, fill = Condition)) + 
  geom_boxplot()  + 
  labs(x=NULL, 
       y="Average Gene Expression") +
  facet_grid(~ Cluster) +
  theme_classic()+
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05) , textsize=3, 
              comparisons = list(c("Mock", "Wildtype"),
                                 c("Wildtype", "3M"),
                                 c("Wildtype", "dNS1")),
              show.legend = T, 
              test = "t.test",
              vjust = 0.5,
              hjust = 0.4,
              step_increase = 0.1) +
  scale_fill_manual(values = COL, 
                    name="Condition") +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white"), 
        legend.position = "bottom",
        legend.direction = "horizontal")

g2 <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g2$layout$name))
fills <- unique(h.clusters$colors)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#Plot2 - heatmap
vsd <- vst(dds_all, blind = TRUE) 

clustmat = matrix(0, nrow = length(h.clusters$genes), ncol = dim(assay(vsd))[[2]])
rownames(clustmat) = as.character(h.clusters$genes)
colnames(clustmat) = colnames(assay(vsd))
m = match(rownames(clustmat), rownames(assay(vsd)))
f.a = !is.na(m)
f.t = m[f.a]
clustmat[f.a,] = assay(vsd)[f.t,]
m = match(clusters$V1, rownames(clustmat))
clustmat = clustmat[m,]

COL = c("Mock" = "tan1","Wildtype" = "firebrick1", "dNS1"="blue2", "3M" = "limegreen")

samples = as.data.frame(colData(dds_all))
m = match(samples$infection, names(COL))
f.a = !is.na(m)
f.t = m[f.a]
conds = as.data.frame(samples[,"infection"])
rownames(conds) = samples$sample
colnames(conds) ="Conditions"
conds$Conditions = factor(conds$Conditions, levels = c("Mock", "Wildtype","3M", "dNS1"))

cols = as.character(unique(h.clusters$colors))
names(cols) = as.character(unique(h.clusters$labels))
col_fun = cividis(100)
prot = h.clusters$labels %>% as.data.frame()
rownames(prot) = h.clusters$genes
colnames(prot) ="Clusters"
prot = prot[rownames(clustmat),,drop = F]
ha = rowAnnotation(df = prot,
                   col = list(Clusters = cols),
                   show_annotation_name = F,
                   border = T,
                   annotation_legend_param = list(grid_height = unit(0.5, "cm"),
                                                  title_gp=gpar(fontsize=10, fontface="bold"),
                                                  labels_gp = gpar(fontsize = 8))
)
hb = columnAnnotation(df = conds,
                      col = list(Conditions = COL),
                      show_annotation_name = F,
                      show_legend = T,
                      border = T,
                      annotation_legend_param = list(grid_height = unit(0.5, "cm"),
                                                     title_gp=gpar(fontsize=10, fontface="bold"),
                                                     labels_gp = gpar(fontsize = 8))
)

#Scaling:
mat_scaled =t(scale(t(clustmat)))
colnames(mat_scaled) =c('Mock-1','Mock-2','Mock-3','Wildtype-1','Wildtype-2','Wildtype-3','3M-1','3M-2','3M-3','dNS1-1','dNS1-2','dNS1-3')
ht = Heatmap(mat_scaled,
        name = "Z-Score",
        col = viridis(200),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        row_title_gp = gpar(fontsize = 10, fontface = "bold" ),
        heatmap_legend_param = list(
          at = c(-2,-1,0,1,2),
          color_bar="continuous",
          labels_gp = gpar(fontsize = 8),
          legend_height = unit(2, "cm"),
          title_position = "topleft",
          title_gp=gpar(fontsize=10, fontface="bold"),
          border = T),
        top_annotation = hb,
        left_annotation = ha,
        cluster_columns = F,
        column_dend_side = "top",
        row_split = prot,
        show_column_dend = T,
        cluster_rows = F,
        row_names_gp = gpar(fontsize = 8),
        show_row_names = F,
        show_column_names = T,
        row_names_rot = 90,
        row_dend_side = "right", 
        row_names_side = "left",
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 45,
        height = unit(9,"cm"),
        width = unit(6,"cm"))

ht = grid::grid.grabExpr(draw(ht,heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE))

r1 = plot_grid(g1, ht, nrow = 1, ncol = 2, rel_heights = c(2,1), labels = LETTERS[1:2], label_fontfamily = 'serif', label_size = 18)

pdf("Thesis figures/Chp3-Figure 1.pdf", height = 8,width = 9, paper = "a4", useDingbats =  F)
plot_grid(r1,g2, nrow = 2,  rel_heights = c(1.5,1),labels = c('',LETTERS[3]), label_fontfamily = 'serif', label_size = 18)
dev.off()



