##SARS-CoV-2 figures


##Packages:

if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}
if( !require("gplots")){
  BiocManager::install("gplots")
}
if( !require("EnhancedVolcano")){
  BiocManager::install("EnhancedVolcano")
}
if( !require("DESeq2")){
  BiocManager::install("DESeq2")
}
if( !require("EGAD")){
  BiocManager::install("EGAD")
}
if( !require("dplyr")){
  BiocManager::install("dplyr")
}

####################Volplots for human####################

##Set-up:

setwd('D:/')
i = 1 
orthos = read.csv("orthos.csv", stringsAsFactors = F)
orthos = orthos[orthos$Ferret.homology.type == "ortholog_one2one",]
orthos = orthos[orthos$Golden.Hamster.homology.type == "ortholog_one2one",]
orthos = orthos[orthos$Mouse.homology.type == "ortholog_one2one",]


setwd("COVID/Clinical studies")
#Get DESeq2 output and formatting:
ids = c("GSE147507.A549","GSE147507.NHBE", "GSE147507.Calu-3","GSE147507.biopsy","GSE150316", "GSE150819",
        "COVID1","COVID2","GSE147507.ferret","GSE154104") # list real studies here 
species = c(rep('Human',6), rep('Hamster',2), 'Ferret','Mouse')
mylist = list()

for (i in 1:length(ids)) {
  id = ids[i]
  filename = paste0(id,".DE.Rdata") 
  load(filename)
  mylist[[i]] = as.data.frame(res)
  mylist[[i]]$GeneID = rownames(mylist[[i]])
  mylist[[i]]$Study = id
  mylist[[i]]$Species = species[i]
  rownames(mylist[[i]]) = NULL
}  
mylist = do.call(rbind.data.frame, mylist)
mylist$Study = factor(mylist$Study, levels = c("GSE150316","GSE150819","GSE147507.biopsy","GSE147507.NHBE","GSE147507.A549", "GSE147507.Calu-3","GSE147507.ferret","GSE154104","COVID2","COVID1"), 
                      labels = c("GSE150316","GSE150819","Biopsies","NHBE","A549","Calu-3","Trachea","GSE154104","GSE156005", 'Unpublished'))
mylist$Species = factor(mylist$Species, levels = c('Human', 'Ferret', 'Mouse', 'Hamster'))

#Facetting Volplots:
covnames = mylist[grep('ENS', mylist$GeneID, invert = T),]$GeneID %>% unique()
cols = c(rep('red3',6),rep('LightSkyBlue1',2), 'LightGoldenrod2','green3')
i=1

p = EnhancedVolcano(mylist,
                    lab = 'GeneID',
                    selectLab = covnames,
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    xlab = expression(paste('Log'[2],'FC')),
                    ylab = expression(paste('-Log'[10],'P')),
                    pointSize = 0.5,
                    pCutoff = 0.05,
                    FCcutoff = 0.5,
                    title = NULL,
                    cutoffLineType = "dashed",
                    raster = TRUE,
                    subtitle = NULL,
                    col=c('black', 'green3', 'blue', 'red3'),
                    legendLabels = c("Not significant",
                                     "-0.5 > Log2FC > 0.5",
                                     "p <= 0.05",
                                     "Significant"),
                    legendIconSize = 3,
                    caption = NULL,
                    colAlpha = 1,
                    border = 'full',
                    legendPosition = 'none') +
  theme_classic() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ Study,nrow = 2, drop = T)
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c('red3','LightGoldenrod2','green3','LightSkyBlue1', 'LightSkyBlue1', rep('red3', 5))
k <- 1
for (f in stripr) {
  j <- which(grepl('rect', g$grobs[[f]]$grobs[[1]]$childrenOrder))
  g$grobs[[f]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("D:/Thesis figures/Figure 5A.pdf", height = 6, width = 7, useDingbats = F)
grid::grid.draw(g)
dev.off()
save(mylist, file = 'COVID/Clinical studies/COVID-DESeq2-output-alldatasets.Rdata')

####Plot 2 - Barplots of DEGs and genetypes ###

load("COVID/Clinical studies/cross.species.genesets.Rdata")

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
df$Class = factor(df$Class, levels = c("Upregulated", "Downregulated"))
df$Species = rep(c(rep('Human', 6),rep('Hamster',2), 'Ferret', 'Mouse'),2)
df$Species = factor(df$Species, levels = c('Human', 'Ferret','Mouse','Hamster'))
df$Study = factor(df$Study, levels = c('GSE150316','GSE150819','Biopsies','NHBE','A549',
                                       'Calu-3','Trachea','GSE154104','GSE156005','Unpublished'))
p = ggplot(df,aes(x=Study,y=Totals,fill=factor(Class, levels = c("Upregulated","Downregulated")))) +
  geom_point(size=0,shape =NA) +
  geom_bar(stat="identity",color = "black", width = 0.7) +
  labs(y = "Totals", x = NULL) +
  theme_classic() + 
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual("Class", values = c("Upregulated" = "#FF3333","Downregulated" = "#33CCFF")) +
  facet_grid(~Species, scale = 'free', space = 'free_x')+
  geom_text(aes(y = Totals+(Totals/abs(Totals)*190),
                label = abs(Totals),
                group = Study), 
            vjust = 0.4) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 0.95, size = 10),
        plot.margin = unit(c(2,1,2,1), "cm"),
        legend.position = 'bottom',
        legend.direction = 'horizontal')   

g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- c('red3','LightGoldenrod2','green3','LightSkyBlue1')
k <- 1
for (f in stripr) {
  j <- which(grepl('rect', g$grobs[[f]]$grobs[[1]]$childrenOrder))
  g$grobs[[f]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
 
#Human
load("DE_prior.rdata") 
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
roch= do.call(rbind.data.frame, roc)

#Mouse
load("Mouse_DE_prior.rdata")
genesets = genesets.up[,c(7:10)]
m = match(rownames(genesets), orthos$Gene.stable.ID)
f.a = !is.na(m)
f.t = m[f.a]
rownames(genesets)[f.a] = orthos[f.t,]$Mouse.gene.stable.ID
genesets = genesets[rownames(genesets) %in% rownames(DE_list),, drop = F] # if using a specific annotation set e.g. clusters etc
m = match(rownames(DE_list), rownames(genesets)) 
f.go = !is.na(m)
DE_list_reviewed = DE_list[f.go,]
genesets <- genesets[rownames(DE_list_reviewed),]
DE_list_reviewed[,4] = rank(DE_list_reviewed[,4])

roc = list()
aucs.allm = auc_multifunc( genesets, DE_list_reviewed[,4])
names(aucs.allm) = colnames(genesets)
aucs = c(aucs.all, aucs.allm)
aucs = aucs[ids]
i=1
for (i in 1:dim(genesets)[[2]]) {
  roc[[i]] <- get_roc(DE_list_reviewed[,4], genesets[,i])
  roc[[i]] <- data.frame(roc[[i]])
  roc[[i]]$Studies = colnames(genesets)[i]
}
rocm= do.call(rbind.data.frame, roc)

roc = rbind.data.frame(rocm, roch)
ids = c("GSE150316","GSE150819","Biopsies","NHBE","A549","Calu-3","Trachea","GSE154104","GSE156005", 'Unpublished')
labels = sapply(1:length(ids), function(x) paste0(ids[x],'\n(AUC =',round(aucs[x], digits = 2),')'))
roc$Studies = factor(roc$Studies, levels = ids, labels = labels)

#DE prior multifunctionality scores per study saved in:
save(roc1, rocm, file = 'Study specific multifunc scores-de-prior.rdata')

p = ggplot(data = roc ) +
  geom_line(aes(x = fpr, y = tpr)) +
  facet_wrap(~ Studies,nrow = 2) +
  geom_abline(color = 'Red3') +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = 'grey', linetype = 2, size = 0.2),
        panel.grid.major.y = element_line(color = 'grey', linetype = 2, size = 0.2))
g2 <- ggplot_gtable(ggplot_build(p))

stripr <- which(grepl('strip-t', g2$layout$name))
fills <- c('red3','LightGoldenrod2','green3','LightSkyBlue1', 'LightSkyBlue1', rep('red3', 5))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

b = plot.new()
h = plot_grid(b,g,b,ncol = 3, rel_widths = c(0.9,3.2,0.9))
pdf("D:/Thesis figures/Figure 5C.pdf", height = 11, width = 9,useDingbats = F)
plot_grid(h, g2, nrow = 2,
          labels = LETTERS[1:2], 
          label_size = 18,
          label_fontfamily = 'serif',
          rel_heights = c(1.2,0.8))

dev.off()

####Plot 3 - GO Enrichment #####

df.res = mylist

m = match(df.res$GeneID, orthos$Ferret.gene.stable.ID)
f.a = !is.na(m)
f.t = m[f.a]
df.res[f.a,]$GeneID = orthos[f.t,]$Gene.stable.ID
m = match(df.res$GeneID, orthos$Mouse.gene.stable.ID)
f.a = !is.na(m)
f.t = m[f.a]
df.res[f.a,]$GeneID = orthos[f.t,]$Gene.stable.ID
m = match(df.res$GeneID, orthos$Golden.Hamster.gene.stable.ID)
f.a = !is.na(m)
f.t = m[f.a]
df.res[f.a,]$GeneID = orthos[f.t,]$Gene.stable.ID


load("GO Annotation Jan2020/GO.human.Rdata")
load("gene_annotations_v29.Rdata")
load("DE_prior.rdata") 
allgenes = read.csv("allgenes_Ensembl_v97.csv", header = T, stringsAsFactors = F)

load("GO Annotation Jan2020/GO.mouse.Rdata")
load('mouse_gene_annotations_v100.Rdata')
load("recurrence_matrix_MOUSE.Rdata")
load("Mouse_DE_prior.rdata")

setwd("E:/")
i =1
go.up = list()
go.down = list()
count.up = list()
count.down = list()

for (i in 1:length(unique(df.res$Species))) {

    #Loading human GO for human studies only.
    load("GO Annotation Jan2020/GO.human.Rdata")
    load("gene_annotations_v29.Rdata")
    gosums = colSums(GO.human.nonIEA )
    filtGO = gosums <=500 & gosums>=100
    annot =  GO.human.nonIEA[,filtGO]
    filt = voc$V3 == "biological_process"
    
    #Change label sizes:
    m = match(rownames(annot), attr$entrezID)
    f.a = !is.na(m)          
    f.t = m[f.a]
    rownames(annot)[f.a] = as.character(attr[f.t,]$ensemblID)
    
    #Get DEGs
    res.up = df.res[df.res$Species == unique(df.res$Species)[i] & df.res$log2FoldChange >= 0.5 & df.res$pvalue <= 0.05,]
    res.down = df.res[df.res$Species == unique(df.res$Species)[i] & df.res$log2FoldChange <= -0.5 & df.res$pvalue <= 0.05,]
    
    m = match(res.up$GeneID, orthos$Gene.stable.ID)
    f.a = !is.na(m)
    res.up = res.up[f.a,]
    m = match(res.down$GeneID, orthos$Gene.stable.ID)
    f.a = !is.na(m)
    res.down = res.down[f.a,]
    
    #Binary matrices for GSEA
    g.up = matrix(0, ncol = length(unique(res.up$Study)), nrow = length(orthos$Gene.stable.ID))
    g.down = matrix(0, ncol = length(unique(res.down$Study)), nrow = length(orthos$Gene.stable.ID))
    rownames(g.up) = orthos$Gene.stable.ID
    rownames(g.down) = orthos$Gene.stable.ID
    colnames(g.up) = unique(res.up$Study)
    colnames(g.down) = unique(res.down$Study)
    
    #Matching DEGs and scoring
    for (j in 1:length(unique(res.up$Study))) {
      genes = res.up[res.up$Study == unique(res.up$Study)[j],]$GeneID
      g.up[rownames(g.up) %in% genes,as.character(unique(res.up$Study)[j])] = 1
      genes = res.down[res.down$Study == unique(res.down$Study)[j],]$GeneID
      g.down[rownames(g.down) %in% genes,as.character(unique(res.down$Study)[j])] = 1
    }
    count.up[[i]] = g.up
    count.down[[i]] = g.down
    ##########Upregulated GO Enrichment##########
    res = list()
    n = ncol(g.up)
    k = 1
    for (k in 1:n){
      genes = rownames(g.up[g.up[,k] == 1,,drop = F])
      m = match (genes, rownames(annot))
      f.a = !is.na(m)
      genes = genes[f.a]
      resMG = gene_set_enrichment(genes, annot, voc[filt,])
      filt2 = resMG[,"padj"] <= 0.05
      res[[k]] = resMG[filt2,]
    }
    n = length(res)
    TFterms = matrix(0,nrow = nrow(voc), ncol = n)
    rownames(TFterms) = voc[,1]
    colnames(TFterms) = colnames(g.up)
    l =1
    for (l in 1:n) {
      m = match(rownames(TFterms), as.matrix(res[[l]][,1]))
      f.a = !is.na(m)
      f.t = m[f.a]
      TFterms[f.a, l] = -log10(as.numeric(res[[l]][f.t,5]))
    }
    filt4 = rowSums(TFterms) > 0
    g = TFterms[filt4,,drop = F]
    m = match(rownames(g), voc[,1])
    f.a = !is.na(m)
    f.t = m[f.a]
    rownames(g)[f.a] = as.character(voc[f.t,2])
    g1 = g
    
    res = list()
    n = ncol(g.down)
    k = 1
    for (k in 1:n){
      genes = rownames(g.down[g.down[,k] == 1,,drop = F])
      m = match (genes, rownames(annot))
      f.a = !is.na(m)
      genes = genes[f.a]
      resMG = gene_set_enrichment(genes, annot, voc[filt,])
      filt2 = resMG[,"padj"] <= 0.05
      res[[k]] = resMG[filt2,]
    }
    n = length(res)
    TFterms = matrix(0,nrow = nrow(voc), ncol = n)
    rownames(TFterms) = voc[,1]
    colnames(TFterms) = colnames(g.up)
    l =1
    for (l in 1:n) {
      m = match(rownames(TFterms), as.matrix(res[[l]][,1]))
      f.a = !is.na(m)
      f.t = m[f.a]
      TFterms[f.a, l] = -log10(as.numeric(res[[l]][f.t,5]))
    }
    filt4 = rowSums(TFterms) > 0
    g = TFterms[filt4,,drop = F]
    m = match(rownames(g), voc[,1])
    f.a = !is.na(m)
    f.t = m[f.a]
    rownames(g)[f.a] = as.character(voc[f.t,2])
    g2 = g
 
  go.up[[i]] = g1
  go.down[[i]] = g2
}

gonames = unique(c(rownames(go.down[[1]]),rownames(go.down[[2]]),rownames(go.down[[3]]),rownames(go.down[[4]])))
idnames = unique(c(colnames(go.down[[1]]),colnames(go.down[[2]]),colnames(go.down[[3]]),colnames(go.down[[4]])))
g = matrix(0, nrow = length(gonames), ncol = length(idnames))
rownames(g) = gonames
colnames(g) = idnames

m = match(rownames(g), rownames(go.down[[1]]))
f.a = !is.na(m)
f.t = m[f.a]
n = match(colnames(g), colnames(go.down[[1]]))
f.b = !is.na(n)
f.u = n[f.b]
g[f.a,f.b] = go.down[[1]][f.t,f.u]
m = match(rownames(g), rownames(go.down[[2]]))
f.a = !is.na(m)
f.t = m[f.a]
n = match(colnames(g), colnames(go.down[[2]]))
f.b = !is.na(n)
f.u = n[f.b]
g[f.a,f.b] = go.down[[2]][f.t,f.u]
m = match(rownames(g), rownames(go.down[[3]]))
f.a = !is.na(m)
f.t = m[f.a]
n = match(colnames(g), colnames(go.down[[3]]))
f.b = !is.na(n)
f.u = n[f.b]
g[f.a,f.b] = go.down[[3]][f.t,f.u]
m = match(rownames(g), rownames(go.down[[4]]))
f.a = !is.na(m)
f.t = m[f.a]
n = match(colnames(g), colnames(go.down[[4]]))
f.b = !is.na(n)
f.u = n[f.b]
g[f.a,f.b] = go.down[[4]][f.t,f.u]

topVarGenes <- sapply(1:ncol(g.up), function(x) head(order(g.up[,x], decreasing = T),5)) %>% as.numeric() %>% unique(.)
prot = c("Ferret" = "LightGoldenrod2", "Hamster" = "LightSkyBlue1","Human" = "LightSalmon1", "Mouse" = "DarkSeaGreen1" )
df = data.frame(c(rep("Human", 6),rep("Hamster", 2),"Ferret", "Mouse" ))
colnames(df) = "Species"
df$Species = factor(df$Species, levels = c('Human','Ferret','Mouse', 'Hamster'))
rownames(df) = colnames(g.up)

hb = columnAnnotation(df = df,
                      col = list(Species = prot),
                      show_legend = TRUE,
                      annotation_name_gp = gpar(fontsize=9, fontface="bold"),
                      annotation_name_side = "left",
                      show_annotation_name = FALSE,
                      annotation_legend_param = list(title_gp=gpar(fontsize=9, fontface="bold"),
                                                     legend_direction = "horizontal",
                                                     legend_height = unit(2, "cm"),
                                                     labels_gp = gpar(fontsize = 8),
                                                     legend_width = unit(2, "cm")),
                      annotation_label = NULL,
                      border = T)
col_fun = colorRamp2(c(seq(from = 0, to = 12, length.out = 30)), colorpanel(30,low="white", mid = "#FF3333",high="#660000"))
ht = list()
ht[[1]] = Heatmap(head(g.up[topVarGenes,][order(rowSums(g.up[topVarGenes,]),decreasing = T),],25),
        name = "Upregulated GO-BP\n -log10(P-value)",
        col = col_fun,
        use_raster = F,
        rect_gp = gpar(col = "black", lwd = 2),
        row_names_rot = 0,
        top_annotation = hb,
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
        column_split = df,
        column_gap = unit(2, "mm"),
        show_row_dend = F,
        row_names_gp = gpar(fontsize = 10),
        row_names_side = "left", 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 90,
        width = unit(4,"cm"),
        height = unit(12,"cm"))

topVarGenes <- sapply(1:ncol(g.down), function(x) head(order(g.down[,x], decreasing = T),5)) %>% as.numeric() %>% unique(.)
col_fun = colorRamp2(c(seq(from = 0, to = 12, length.out = 30)), colorpanel(30,low="white", mid = "#33CCFF",high="#0000FF"))
ht[[2]] =Heatmap(head(g.down[topVarGenes,][order(rowSums(g.down[topVarGenes,]),decreasing = T),],25),
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
        column_split = df,
        column_gap = unit(2, "mm"),
        show_row_dend = F,
        row_names_gp = gpar(fontsize = 10),
        row_names_side = "left", 
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 90,
        width = unit(4,"cm"),
        height = unit(12,"cm"))

#save(g.up, g.down, file = 'COVID_Studies GO Human Enrichment.Rdata')

pdf("E:/Thesis figures/Figure 6.pdf", height = 15,width = 10, useDingbats = F)
draw(ht[[1]] %v% ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

