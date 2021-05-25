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
if( !require("circlize")){
  BiocManager::install("circlize")
}
if( !require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}
if( !require("gridBase")){
  BiocManager::install("gridBase")
}

orthos = read.csv("orthos.csv", stringsAsFactors = F)
orthos = orthos[orthos$Mouse.homology.type == "ortholog_one2one",]
orthos = orthos[orthos$Ferret.homology.type == "ortholog_one2one",]
orthos = orthos[orthos$Golden.Hamster.homology.type == "ortholog_one2one",]

##Plot 1 - Cluster overlaps:

load('Coexpression-Sars-CoV-2-human-network and datasets.rdata')

clusters = clustersh
hnet = localh

consTree = hclust( dist(hnet, method = "euclidean"), method = "average")
consDend = as.dendrogram(consTree)
clusterids <- data.frame( genes = clusters$V1, labels = as.numeric(clusters$unmergedLabels.mod),
                          labels_unmerged =  clusters$unmergedLabels.mod,
                          colors = clusters$V3)

m = match(clusters$V1, rownames(hnet))

clust_cnet <- list(as.matrix(hnet), consTree, consDend, m, clusterids)
names(clust_cnet) <- c( "distance_matrix", "tree", "dendrogram", "order", "clusters")


load('Coexpression-influenza-human-network and datasets.rdata')

clusters = clusters
hnet = local
consTree = hclust( dist(hnet, method = "euclidean"), method = "average")
consDend = as.dendrogram(consTree)
clusterids <- data.frame( genes = clusters$V1, labels = as.numeric(clusters$unmergedLabels.mod),
                          labels_unmerged =  clusters$unmergedLabels.mod,
                          colors = clusters$V3)
m = match(clusters$V1, rownames(hnet))
clust_hnet <- list(as.matrix(hnet), consTree, consDend, m, clusterids)
names(clust_hnet) <- c( "distance_matrix", "tree", "dendrogram", "order", "clusters")

load('Coexpression-influenza-mouse-network and datasets.rdata')
clusters = clusters
hnet = local

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

##Practice!!!!!!

m.clusters = clust_mnet$clusters
h.clusters = clust_hnet$clusters
c.clusters = clust_cnet$clusters
save(m.clusters,h.clusters,c.clusters, file = 'cross-species.clusters.rdata')

load('cross-species.clusters.rdata')

res = list()
res2 = list()
res3 = list()

for (i in 1:length(unique(h.clusters$labels))) {
  m = match(m.clusters$genes, c.clusters$genes)
  f.b =is.na(m)
  m = match(m.clusters[f.b,]$genes,h.clusters[h.clusters$labels == i,]$genes)
  f.a =!is.na(m)
  g= table(m.clusters[f.b,][f.a,]$labels)
  human_cluster = rep(paste0("H",i),length(g))
  mouse_cluster = as.list(names(g)) %>% sapply(., function(x) paste0("M",x))
  res[[i]] = cbind.data.frame(human_cluster, mouse_cluster, g) %>% .[,-3]
}
for (i in 1:length(unique(h.clusters$labels))) {
  m = match(c.clusters$genes, m.clusters$genes)
  f.b =is.na(m)
  m = match(c.clusters[f.b,]$genes,h.clusters[h.clusters$labels == i,]$genes)
  f.a =!is.na(m)
  g= table(c.clusters[f.b,][f.a,]$labels)
  human_cluster = rep(paste0("H",i),length(g))
  mouse_cluster = as.list(names(g)) %>% sapply(., function(x) paste0("C",x))
  res2[[i]] = cbind.data.frame(human_cluster, mouse_cluster, g) %>% .[,-3]
}
for (i in 1:length(unique(m.clusters$labels))) {
  m = match(c.clusters$genes, h.clusters$genes)
  f.b =is.na(m)
  m = match(c.clusters[f.b,]$genes,m.clusters[m.clusters$labels == i,]$genes)
  f.a =!is.na(m)
  g= table(c.clusters[f.b,][f.a,]$labels)
  human_cluster = rep(paste0("M",i),length(g))
  mouse_cluster = as.list(names(g)) %>% sapply(., function(x) paste0("C",x))
  res3[[i]] = cbind.data.frame(human_cluster, mouse_cluster, g) %>% .[,-3]
}

edges = do.call(rbind.data.frame,res)
edges2 = do.call(rbind.data.frame,res2)
edges3 = do.call(rbind.data.frame,res3)
nf = rbind.data.frame(edges,edges2,edges3)
colnames(nf) =c("from","to","value")

res = list()
res2 = list()
res3 = list()
for (i in 1:length(unique(h.clusters$labels))) {
  m = match(m.clusters$genes, c.clusters$genes)
  f.b =!is.na(m)
  m = match(m.clusters[f.b,]$genes,h.clusters[h.clusters$labels == i,]$genes)
  f.a =!is.na(m)
  g= table(m.clusters[f.b,][f.a,]$labels)
  human_cluster = rep(paste0("H",i),length(g))
  mouse_cluster = as.list(names(g)) %>% sapply(., function(x) paste0("M",x))
  res[[i]] = cbind.data.frame(human_cluster, mouse_cluster, g) %>% .[,-3]
}
for (i in 1:length(unique(h.clusters$labels))) {
  m = match(c.clusters$genes, m.clusters$genes)
  f.b =!is.na(m)
  m = match(c.clusters[f.b,]$genes,h.clusters[h.clusters$labels == i,]$genes)
  f.a =!is.na(m)
  g= table(c.clusters[f.b,][f.a,]$labels)
  human_cluster = rep(paste0("H",i),length(g))
  mouse_cluster = as.list(names(g)) %>% sapply(., function(x) paste0("C",x))
  res2[[i]] = cbind.data.frame(human_cluster, mouse_cluster, g) %>% .[,-3]
}
for (i in 1:length(unique(m.clusters$labels))) {
  m = match(c.clusters$genes, h.clusters$genes)
  f.b =!is.na(m)
  m = match(c.clusters[f.b,]$genes,m.clusters[m.clusters$labels == i,]$genes)
  f.a =!is.na(m)
  g= table(c.clusters[f.b,][f.a,]$labels)
  human_cluster = rep(paste0("M",i),length(g))
  mouse_cluster = as.list(names(g)) %>% sapply(., function(x) paste0("C",x))
  res3[[i]] = cbind.data.frame(human_cluster, mouse_cluster, g) %>% .[,-3]
}

edges = do.call(rbind.data.frame,res)
edges2 = do.call(rbind.data.frame,res2)
edges3 = do.call(rbind.data.frame,res3)
df = rbind.data.frame(edges,edges2,edges3)
colnames(df) =c("from","to","value")

set.seed(999)
res = list(c(table(m.clusters$labels)),
           c(table(h.clusters$labels)),
           c(table(c.clusters$labels)))
names(res[[1]]) = paste0('M',names(res[[1]]))           
names(res[[2]]) = paste0('H',names(res[[2]]))           
names(res[[3]]) = paste0('C',names(res[[3]]))           

res1 = data.frame(to = rep(0,3), from = c(sum(res[[1]]), sum(res[[2]]),sum(res[[3]])), row.names = c('Mouse','Human','SARS-CoV-2'))
sectors = rownames(res1)

b1 = c(0,44,115,281,363,944)
b2 = c(0,200, 236,460, 514,780,984)
b3 = c(0,87,101)

#Plot
circlize_plot = function() {
  circos.par("track.height" = 0.1, gap.after = c(10,10,10),cell.padding = c(0.02, 0, 0.02, 0))
  circos.initialize(xlim = res1)
  circos.track(sectors, ylim = c(0,1),
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, 
                             CELL_META$cell.ylim[2] +mm_y(5),
                             CELL_META$sector.index)
                 circos.axis(labels.cex = 0.6)
                 if(CELL_META$sector.index == "Mouse"){
                   circos.rect(0, 0,
                               CELL_META$xlim, 1,
                               col = 'green3', border = 'black')
                 } else if(CELL_META$sector.index == "Human") {
                   circos.rect(0, 0,
                               CELL_META$xlim, 1,
                               col = 'red3', border = 'black')
                 } else {
                   circos.rect(0, 0,
                               CELL_META$xlim, 1,
                               col = 'purple', border = 'black')
                   
                 }
               })
  
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    if(CELL_META$sector.index == "Mouse"){
      breaks = b1
      n_breaks = length(breaks)
      col1 = as.character(unique(m.clusters$colors))
      circos.rect(breaks[-n_breaks], 0,
                  breaks[-1], 1,
                  col = col1, border = 'black')
    } else if(CELL_META$sector.index == "Human") {
      breaks = b2
      n_breaks = length(breaks)
      col1 = c(unique(h.clusters$colors))
      circos.rect(breaks[-n_breaks], 0,
                  breaks[-1], 1,
                  col = col1, border = 'black')
    } else {
      breaks = b3
      n_breaks = length(breaks)
      col1 = c(unique(c.clusters$colors))
      circos.rect(breaks[-n_breaks], 0,
                  breaks[-1], 1,
                  col = col1, border = 'black')
      
    }
  })
  
  #save(df,nf,file = 'cluster_overlaps_in_recurrence.rdata')
  
  #Genes overlaps across three recurrence sets
  
  #Human            #Mouse             #SARS-CoV-2
  #H1 = 0, 200      M1 = 0, 44         C1 = 0, 87
  #H2 = 200, 236    M2 = 44, 115       C2 = 87, 101
  #H3 = 236, 460    M3 = 115, 281
  #H4 = 460, 514    M4 = 281, 363
  #H5 = 514, 780    M5 = 363, 944
  #H6 = 780, 984
  
  #Non-specific links:
  #H1 (MAX - 180)
  circos.link('Human',c(170,171),'SARS-CoV-2', c(28,29),col = '#af875f')
  circos.link('Human',c(160,162),'Mouse', c(100,102),col = '#af875f')
  circos.link('Human',c(155,157),'Mouse', c(270,272),col = '#af875f')
  circos.link('Human',c(130,135),'Mouse', c(330,335),col = '#af875f')
  circos.link('Human',c(40,47),'Mouse', c(800, 807),col = '#af875f')
  
  #H2
  circos.link('Human',c(214,217),'Mouse', c(20,21),col = '#af875f')
  circos.link('Human',c(212,213),'Mouse', c(95,96),col = '#af875f')
  circos.link('Human',c(205,206),'Mouse', c(790,791),col = '#af875f')
  
  #H3
  circos.link('Human',c(310,311),'Mouse', c(89,90),col = '#af875f')
  circos.link('Human',c(295,296),'Mouse', c(321,322),col = '#af875f')
  circos.link('Human',c(290,291),'Mouse', c(720,721),col = '#af875f')
  
  #H4
  circos.link('Human', c(470,480), 'Mouse', c(700,710),col = '#af875f')
  
  #H5
  circos.link('Human', c(710,711), 'Mouse', c(250,251),col = '#af875f')
  circos.link('Human', c(705,707), 'Mouse', c(308,309),col = '#af875f')
  circos.link('Human', c(550,552), 'Mouse', c(670,672),col = '#af875f')
  
  #H6
  circos.link('Human', c(950,954), 'SARS-CoV-2',c(4,8),col = '#af875f')
  circos.link('Human', c(945,946), 'Mouse', c(130,131),col = '#af875f')
  circos.link('Human', c(920,922), 'Mouse', c(285,287),col = '#af875f')
  circos.link('Human', c(785,807), 'Mouse', c(645,667),col = '#af875f')
  
  #Mouse clusters
  
  circos.link('Mouse',c(160,162), 'SARS-CoV-2',c(80, 82),col = '#af875f')
  circos.link('Mouse', c(420,450),'SARS-CoV-2',c(52, 82),col = '#af875f')
  circos.link('Mouse',c(282,283), 'SARS-CoV-2',c(100,101),col = '#af875f')
  circos.link('Mouse', c(364,370),'SARS-CoV-2',c(95,101),col = '#af875f')
  
  
  #Paired links
  #M4
  #H2-M4
  circos.link('Human',c(206,208),'Mouse', c(325,327),col = "#5f0087")
  #H2-C2
  circos.link('Human',c(206,208),'SARS-CoV-2', c(90,92),col = "#5f0087")
  #M4-C2
  circos.link('Mouse',c(325,327),'SARS-CoV-2', c(90,92),col = "#5f0087")
  #H3-M4
  circos.link('Human',c(298,299),'Mouse', c(318,319),col = "#5f0087")
  #H3-C2
  circos.link('Human',c(298,299),'SARS-CoV-2', c(93,94),col = "#5f0087")
  #H6-M4
  circos.link('Human',c(780,780+1),'Mouse', c(281+2+1,281+2+1+1),col = "#5f0087")
  #H6-C1 (sPECIFIC to M4)
  circos.link('Human',c(780,780+1),'SARS-CoV-2', c(0,1),col = "#5f0087")    #(1 gene)
  #M4-C1
  circos.link('Mouse',c(281,281+1),'SARS-CoV-2', c(0,1),col = "#5f0087")
  
  #M2
  #H4-M2
  circos.link('Human',c(490, 491),'Mouse', c(45,46),col = "#5f0087")
  #H4-C1 (SPECIFIC TO M2)
  circos.link('Human',c(490, 491),'SARS-CoV-2', c(86,87),col = "#5f0087")
  #M2-C1
  circos.link('Mouse',c(45, 46),'SARS-CoV-2', c(85,87),col = "#5f0087") #Two genes (For H4,H6) in this connection (mAY NEED TO INCREASE)
  #H6-M2
  circos.link('Human',c(830,831),'Mouse', c(44,45),col = "#5f0087")
  #H6-C1 (SPECIFIC TO M2)
  circos.link('Human',c(830,831),'SARS-CoV-2', c(85,86),col = "#5f0087") #(1 gene)
  
  #M3
  #H6-M3
  circos.link('Human',c(870, 872),'Mouse', c(190,192),col = "#5f0087")
  #M3-C1
  circos.link('Mouse',c(190,192),'SARS-CoV-2', c(40,42),col = "#5f0087")  # 2 genes outside core network
  #H6-C1 (SPECIFIC TO M3)
  circos.link('Human',c(870, 872),'SARS-CoV-2', c(n,n+2),col = "#5f0087") 
  
  #M5
  #H2-M5
  circos.link('Human',c(210,212),'Mouse', c(785,787),col = "#5f0087")
  #M5-C1 (SPECIFIC TO H2)
  circos.link('Mouse',c(785,787),'SARS-CoV-2', c(47,49),col = "#5f0087")
  #H2-C1  (SPECIFIC TO M5)
  circos.link('Human',c(210,212),'SARS-CoV-2', c(47,49),col = "#5f0087") ###?
  #H4-M5
  circos.link('Human',c(482,486),'Mouse', c(680,684),col = "#5f0087")
  #M5-C1 (SPECIFIC TO H4)
  circos.link('Mouse',c(680,684),'SARS-CoV-2', c(40,44),col = "#5f0087")
  #H4-C1 (SPECIFIC TO M5)
  circos.link('Human',c(482,486),'SARS-CoV-2', c(40,44),col = "#5f0087")
  #H6-M5
  circos.link('Human',c(820,850),'Mouse', c(540,570),col = "#5f0087")
  #M5-C1
  circos.link('Mouse', c(540,570),'SARS-CoV-2', c(16,46),col = "#5f0087")
  #H6-C1
  circos.link('Human',c(820,850),'SARS-CoV-2', c(16,46),col = "#5f0087")
  
  circos.clear()
  
}

lg1 = Legend(labels = c('Influenza (Human)', 'Influenza (Mouse)','SARS-CoV-2'),
             type = "grid",
             legend_gp = gpar(fill = c('red3', 'green3', 'purple')),
             title_position = "topleft",
             border = 'black',
             title = "Recurrence Set")
lg2 = Legend(labels = c(1:6),
             type = "grid",
             legend_gp = gpar(fill = c(unique(h.clusters$colors))),
             title_position = "topleft",
             border = 'black',
             title = "Human Clusters")
lg3 = Legend(labels = c(1:5),
             type = "grid",
             legend_gp = gpar(fill = c(unique(m.clusters$colors))),
             title_position = "topleft",
             border = 'black',
             title = "Mouse Clusters")
lg4 = Legend(labels = c(1:2),
             type = "grid",
             legend_gp = gpar(fill = c(unique(c.clusters$colors))),
             title_position = "topleft",
             border = 'black',
             title = "SARS-CoV-2 Clusters")
lg5 = Legend(labels = c('2 Sets','3 Sets'),
             type = "grid",
             legend_gp = gpar(fill = c('#af875f', '#5f0087')),
             title_position = "topleft",
             border = 'black',
             title = "Overlaps")


circle_size = unit(1.3, "snpc") # snpc unit gives you a square region
lgd_list = packLegend(lg1, lg2, lg3, lg4, lg5)

pdf("Thesis figures/Chp2-Figure 6.pdf", height = 5,width = 8, useDingbats =  F)
circlize_plot()
draw(lgd_list, x = circle_size, just = "left")
dev.off()
