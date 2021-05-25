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

setwd("E:/")
load("single-cell data/SCTransform/wsn.integrated.rdata")
lvls = c('B Cell','CD4+ T Cell','Th2 CD4+ Cell','CD8+ T Cell','Cytotoxic CD8+ T Cell','NK Cell','Neutrophil', 'Macrophage',
         'Epithelial Cell', 'Type I Pneumocyte', 'Type II Pneumocyte','Lipofibroblast', 'Myofibroblast', 'Fibroblast', 'Matrix Fibroblast',
         'Endothelial Cell')
colour_sample = c("Mock" = "tan1","Wildtype" = "firebrick1", "3M" = "limegreen", "dNS1"="blue2")
lvls2 =data.frame(wsn.integrated@meta.data[,c("Clusters","Cell.type")], row.names = NULL) %>% .[order(.[,1]),] %>% unique() %>% .[,2] %>% as.character()

#Plot 1 - UMAP-plots

p1 = DimPlot(wsn.integrated, 
             reduction = "umap", 
             pt.size = 0.6,
             ncol = 1,
             raster = TRUE,
             group.by = "Cell.type", 
             label = F) +
  theme_classic() +
  theme(legend.position = 'none',
        plot.title = element_blank())

p1 = LabelClusters(p1,id ='Cell.type', repel = T,size = 3,
                   arrow = arrow(length = unit(0.02, "npc"),
                                 type = 'open',
                                 ends = 'last'),
                   box.padding = 1)
pdf("Thesis figures/Chp3-Figure 3a.pdf", width = 6, height = 5, useDingbats = F)
plot_grid(p1)
dev.off()

p2 = DimPlot(wsn.integrated, 
             reduction = "umap", 
             pt.size = 0.6,
             ncol = 1,
             raster = TRUE,
             group.by = "orig.ident", 
             label = F) +
  theme_classic() +
  theme(legend.position = 'bottom',
        plot.title = element_blank(),
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.title = element_text(size=10, face = "bold"),
        legend.text = element_text(size=8)) +
  guides(color = guide_legend(ncol = 2))


# Plot 2 - Plotting markers for clusters

i = c("Cd4,Cd8a,Cd14,Serping1,H2-DMb2,Il2ra,Bst2,Hhip,Lpl,Ms4a7,Pecam1,Col4a3,Ncr1,Hmgb2,Slc34a2,Cxcr2") %>% strsplit(.,",")
i = sapply(i, function(x) paste0("mmc9-",x))
p4 =FeaturePlot(wsn.integrated, 
                reduction = "umap",
                features = c(i),
                ncol = 4,
                raster = T,
                pt.size = 0.5,
                cols = c("lightgrey","#0000ff","#00005f"),
                combine = F,
                min.cutoff = 0) 
for (i in 1:length(p4)) {
  p4[[i]]$labels$title = gsub("mmc9-","", p4[[i]]$labels$title)
  p4[[i]] = p4[[i]] +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank())
}
p3 =plot_grid(plotlist = p4, ncol = 4)


#Plot 3 - Genetype proportions per sample/condition

res = list()
j=1
for (i in unique(wsn.integrated$Condition)) {
  total = length(wsn.integrated$Cell.type[wsn.integrated$Condition == i])
  res[[j]] = (table(wsn.integrated$Cell.type[wsn.integrated$Condition == i])/total)*100
  res[[j]] = as.data.frame(res[[j]])
  res[[j]]$Condition = i  
  j = j+1
}
res = do.call(rbind.data.frame, res)
lvls2 =data.frame(wsn.integrated@meta.data[,c("Clusters","Cell.type")], row.names = NULL) %>% .[order(.[,1]),] %>% unique() %>% .[,2] %>% as.character()
colnames(res) = c("Celltype", "Percentage", "Condition")
res = res[res$Celltype %in% lvls2,]
col_fun2 = rainbow(n = length(unique(res$Celltype)), start = 0, end = 0.7, s = 0.7, v= 0.85) 

res$Condition = factor(res$Condition, levels = c("Mock", "Wildtype", "3M", "dNS1"))
res$Celltype = factor(as.character(res$Celltype), levels = lvls2)
p4 = ggplot(data = res, aes(x = Condition, y = Percentage, fill =Celltype)) +
  geom_bar(position="stack", stat="identity", colour = "black", width = 0.5) +
  scale_fill_manual(values = col_fun2, drop = T) +
  theme_classic() +
  labs(x = "Condition", y = "Gene-type proportion (%)") +
  theme(axis.text.x.bottom = element_text(angle = 45,hjust = 0.95),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.title = element_text(size=10, face = "bold"),
        legend.text = element_text(size=8))
#Or
res = list()
j=1
for (i in unique(wsn.integrated$orig.ident)) {
  total = length(wsn.integrated$Cell.type[wsn.integrated$orig.ident == i])
  res[[j]] = (table(wsn.integrated$Cell.type[wsn.integrated$orig.ident == i])/total)*100
  res[[j]] = as.data.frame(res[[j]])
  res[[j]]$Condition = i  
  j = j+1
}
res = do.call(rbind.data.frame, res)
colnames(res) = c("Celltype", "Percentage", "Condition")
res$Condition = as.character(res$Condition)
res$Celltype = as.character(res$Celltype)
res[res$Condition == 'mock-1',]$Condition = 'Mock-1'
res[res$Condition == 'mock-2',]$Condition = 'Mock-2'
res[res$Condition == 'WSN-1',]$Condition = 'Wildtype-1'
res[res$Condition == 'WSN-2',]$Condition = 'Wildtype-2'
res[res$Condition == 'WSN3M-1',]$Condition = '3M-1'
res[res$Condition == 'WSN3M-2',]$Condition = '3M-2'
res[res$Condition == 'WSNdelNS1-1',]$Condition = 'dNS1-1'
res[res$Condition == 'WSNdelNS1-2',]$Condition = 'dNS1-2'
res = res[res$Celltype %in% lvls2,]

res$Condition = factor(res$Condition, levels = c("Mock-1","Mock-2", 
                                                 "Wildtype-1","Wildtype-2", 
                                                 "3M-1","3M-2", 
                                                 "dNS1-1","dNS1-2"))
res$Celltype = factor(as.character(res$Celltype), levels = lvls2)
p4 =ggplot(data = res, aes(x = Condition, y = Percentage, fill =Celltype)) +
  geom_bar(position="stack", stat="identity", colour = "black", width = 0.5) +
  scale_fill_manual(values = col_fun2) +
  theme_classic() +
  labs(x = NULL, y = "Celltype (%)") +
  theme(axis.text.x.bottom = element_text(angle = 45,hjust = 0.95),
        legend.position = "bottom",
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,0,-10),
        legend.title = element_text(size=10),
        legend.text = element_text(size=8)) +
        guides(fill = guide_legend(ncol = 3, title.position = 'top'))

r1 =plot_grid(p2,p4, ncol = 2, rel_widths = c(1,1.7),align = 'b',labels = c(LETTERS[2:3]), label_fontfamily = 'serif', label_size = 18)

pdf("Thesis figures/Chp3-Figure 3b.pdf", width = 8, height = 11, paper = 'a4',useDingbats = F)
plot_grid(p3, r1, nrow = 2,  rel_heights = c(1.8,1),labels = c(LETTERS[1],''), label_fontfamily = 'serif', label_size = 18)
dev.off()

#Plot 4 - Gene type enrichment

res = list()
j=1
for (i in unique(wsn.integrated$orig.ident)) {
  total = length(wsn.integrated$Cell.type[wsn.integrated$orig.ident == i])
  res[[j]] = (table(wsn.integrated$Cell.type[wsn.integrated$orig.ident == i])/total)
  res[[j]] = as.data.frame(res[[j]])
  res[[j]]$Sample = i
  res[[j]]$Condition = wsn.integrated$Condition[wsn.integrated$orig.ident == i] %>% unique()
  j = j+1
}
res = do.call(rbind.data.frame, res)
colnames(res) = c("Celltype", "Percentage", "Sample", "Condition")
res = res[res$Celltype %in% lvls,]
res$Sample = factor(res$Sample, levels = c("mock-1","mock-2", 
                                                 "WSN-1","WSN-2", 
                                                 "WSN3M-1","WSN3M-2", 
                                                 "WSNdelNS1-1","WSNdelNS1-2"))
res$Condition = factor(res$Condition, levels = c("Mock", "Wildtype", "3M", "dNS1"))
res$Celltype = factor(as.character(res$Celltype), levels = lvls)
res$`Fold Enrichment` = rep(0,nrow(res))
res$Cond = factor(res$Cond, levels = c("Mock", "Wildtype", "3M", "dNS1"))
res$Pval = rep(0, nrow(res))
res$Average = rep(0, nrow(res))
for (i in levels(res$Celltype)) {
  
  g = mean(res[res$Celltype == i & grepl("Mock", res$Condition),2])     #Get average of celltype % per condition
  
  #Generate Fold enrichment:
  for (j in as.character(res[res$Celltype == i & !grepl("Mock", res$Condition),3])) {
    res[res$Celltype == i & grepl(j, res$Sample),]$`Fold Enrichment` = res[res$Celltype == i & grepl(j, res$Sample),2]/g
  }
  
  #Generate P-values:
  for(k in as.character(unique(res[res$Celltype == i & !grepl("Mock", res$Condition),6]))){
    x = res[res$Celltype == i & grepl(k, res$Cond),2]
    y = res[res$Celltype == i & grepl("Mock", res$Cond),2]
    res[res$Celltype == i & grepl(k, res$Cond),]$Average = mean(log2(res[res$Celltype == i & grepl(k, res$Cond),]$`Fold Enrichment`))
    res[res$Celltype == i & grepl(k, res$Cond),]$Pval = t.test(x, y, alternative = "two.sided")[["p.value"]]
  }
  print(i)
}

res = res[!grepl("Mock", res$Condition),]
res$P.signif = rep("", nrow(res))
res[res$Pval <= 0.001,]$P.signif =c("***")
res[res$Pval <= 0.01 & res$Pval > 0.001,]$P.signif = c("**")
res[res$Pval <= 0.05 & res$Pval > 0.01,]$P.signif = c("*")

res2 = unique(res[,-c(2,3,5,6,7)])
res2 = res2[res2$P.signif != "",]
res = res[,-c(8:10)]
col_fun = rainbow(n = length(unique(res$Celltype))+2, start = 0, end = 0.7, s = 0.7, v= 0.85) 

p =ggplot(data = res, aes(x = factor(Celltype, levels = rev(lvls)), y = log2(`Fold Enrichment`), fill =Celltype)) +
  geom_boxplot() +
  facet_grid(~ Cond) +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  theme_classic() +
  scale_fill_manual(values = col_fun) +
  #geom_text(data = res2, aes(x = Celltype, y = Average+0.5, label = P.signif, group = Condition)) +
  coord_flip() +
  labs(y = "Log2(Fold-change)\n(Infected/Mock)", x = "Celltype") +
  scale_x_discrete(position = "top") +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = 2,size = 0.5, fill =NA), 
        panel.spacing.x = unit(0.5,"cm")) 


g3 <- ggplot_gtable(ggplot_build(p))
fills = c('red3','green3','blue3')
stripr <- which(grepl('strip-t', g3$layout$name))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g3$grobs[[i]]$grobs[[1]]$childrenOrder))
  g3$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("Thesis figures/Chp3-Figure 3c.pdf", width = 7, height = 6, paper = 'a4',useDingbats = F)
grid::grid.draw(g3)
dev.off()
