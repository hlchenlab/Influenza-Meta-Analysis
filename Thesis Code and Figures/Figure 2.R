

#Packages:
if( !require("EGAD")){
  BiocManager::install("EGAD")
}
if( !require("grid")){
  BiocManager::install("grid")
}
if( !require("dplyr")){
  BiocManager::install("dplyr")
}
if( !require("gridExtra")){
  BiocManager::install("gridExtra")
}
if( !require("ggplotify")){
  BiocManager::install("ggplotify")
}
if( !require("cowplot")){
  BiocManager::install("cowplot")
}
if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}

#Make Human/Mouse Ranked List:

setwd("")
load("mouse_human_high_confidence_orthos.rdata") # File of Human and Mouse IDs.
mouse_ortho = mouse_ortho[!duplicated(mouse_ortho$Human.gene.stable.ID),]
mouse_ortho = mouse_ortho[!duplicated(mouse_ortho$Gene.stable.ID),]

#Human
ids = c("GSE84204","GSE89008", "GSE97672", "GSE103477", "GSE103604", "GSE104168", "HKUStudy", "HKUStudy_293")  ##DESeq2 results output for each study
ngenes = mouse_ortho$Human.gene.stable.ID # For Human Ranking
n = length(ids)
N = length(ngenes)
geneset.rank = matrix(0, ncol=n, nrow= N ) 

for( i in 1:n ){
  id = ids[i]
  filename = paste0(id,".DE.Rdata") 
  load(filename)
  m = match(ngenes, rownames(res) )  
  f.a = !is.na(m)
  f.r = m[f.a] 
  genes.both = res[f.r,]$log2FoldChange # Change to pvalue or log2FoldChange to get desired stats.
  geneset.rank[f.a,i] = genes.both
}
colnames(geneset.rank) = ids
rownames(geneset.rank) = ngenes

##Averaging log2FC
human = sapply(1:nrow(geneset.rank), function(x) mean(geneset.rank[x,]))
names(human) = rownames(geneset.rank)
HUMAN_gene_rank_FC = as.data.frame(rank(human), stringsAsFactors = F)
##And/Or Averaging p-val
human = sapply(1:nrow(geneset.rank), function(x) mean(geneset.rank[x,]))
names(human) = rownames(geneset.rank)
HUMAN_gene_rank_pval = as.data.frame(rank(human*-1), stringsAsFactors = F)

##Average and Reranking:
g = cbind(HUMAN_gene_rank_FC, HUMAN_gene_rank_pval) %>% as.data.frame(.,stringAsFactors = F)
g = sapply(1:nrow(g), function(x) mean(as.numeric(g[x,])))
names(g) = rownames(HUMAN_gene_rank_FC)
HUMAN_gene_rank = as.data.frame(rank(g))

#For mouse rankings

ids = c("GSE49933","GSE52405", "GSE100522", "GSE107488", "GSE117029", "ERP020504", "SRP061303") 
ngenes = mouse_ortho$Gene.stable.ID # For Mouse Ranking
n = length(ids)
N = length(ngenes)
geneset.rank = matrix(0, ncol=n, nrow= N ) 
for( i in 1:n ){
  id = ids[i]
  filename = paste0(id,".DE.Rdata") 
  load(filename)
  m = match(ngenes, rownames(res) )  
  f.a = !is.na(m)
  f.r = m[f.a] 
  genes.both = res[f.r,]$log2FoldChange # Change to pvalue or log2FoldChange to get desired stats.
  geneset.rank[f.a,i] = genes.both
}
colnames(geneset.rank) = ids
rownames(geneset.rank) = ngenes

##Averaging log2FC
mouse = sapply(1:nrow(geneset.rank), function(x) mean(geneset.rank[x,]))
names(mouse) = rownames(geneset.rank)
MOUSE_gene_rank_FC = as.data.frame(rank(human), stringsAsFactors = F)
##And/Or Averaging p-val
mouse = sapply(1:nrow(geneset.rank), function(x) mean(geneset.rank[x,]))
names(mouse) = rownames(geneset.rank)
MOUSE_gene_rank_pval = as.data.frame(rank(mouse*-1), stringsAsFactors = F)

##Average and Reranking:
g = cbind(MOUSE_gene_rank_FC, MOUSE_gene_rank_pval) %>% as.data.frame(.,stringAsFactors = F)
g = sapply(1:nrow(g), function(x) mean(as.numeric(g[x,])))
names(g) = rownames(MOUSE_gene_rank_FC)
MOUSE_gene_rank = as.data.frame(rank(g))

save(HUMAN_gene_rank, HUMAN_gene_rank_FC, HUMAN_gene_rank_pval, 
     MOUSE_gene_rank, MOUSE_gene_rank_FC, MOUSE_gene_rank_pval, file = 'D:/DE Probability Ranking.Rdata')

#Convert Mouse IDs to human (Only for Mouse Rank list)

m = match(rownames(MOUSE_gene_rank), mouse_ortho$Gene.stable.ID) ##Change IDs when switching between mouse/human
f.a = !is.na(m)
f.t = m[f.a]
rownames(MOUSE_gene_rank)[f.a] = as.character(mouse_ortho[f.t,]$Human.gene.stable.ID) ##Change IDs when switching between mouse/human

load('Mouse_DE_prior.rdata')
DE_listM = DE_list
m = match(rownames(DE_listM), mouse_ortho$Gene.stable.ID) ##Change IDs when switching between mouse/human
f.a = !is.na(m)
f.t = m[f.a]
rownames(DE_listM)[f.a] = as.character(mouse_ortho[f.t,]$Human.gene.stable.ID) ##Change IDs when switching between mouse/human

load("DE_prior.rdata")

#Rerank priors:
m = match(rownames(DE_list), rownames(DE_listM))
f.a = !is.na(m)
f.t =m[f.a]
DE_list = DE_list[f.a,]
DE_listM= DE_listM[f.t,]
DE_list[,4] = rank(DE_listM[,4])
DE_listM[,4] = rank(DE_listM[,4])

####For Mouse Predictibility:
ids = c("GSE84204","GSE89008", "GSE97672", "GSE103477", "GSE103604", "GSE104168", "HKUStudy", "HKUStudy_293")  
ngenes = mouse_ortho$Human.gene.stable.ID #To match mouse DE Genes
ranking_list = MOUSE_gene_rank
####For Human Predictibility:
ids = c("GSE49933","GSE52405", "GSE100522", "GSE107488", "GSE117029", "ERP020504", "SRP061303") 
ngenes = mouse_ortho$Gene.stable.ID #To match human DE Genes
ranking_list = HUMAN_gene_rank
####For Mouse Prior Ranks:
ids = c("GSE84204","GSE89008", "GSE97672", "GSE103477", "GSE103604", "GSE104168", "HKUStudy", "HKUStudy_293")  
ngenes = mouse_ortho$Human.gene.stable.ID #To match mouse DE Genes
ranking_list = DE_listM
####For Human Prior Rankings:
ids = c("GSE49933","GSE52405", "GSE100522", "GSE107488", "GSE117029", "ERP020504", "SRP061303") 
ngenes = mouse_ortho$Gene.stable.ID #To match human DE Genes
ranking_list = DE_list

###########################Make DE gene matrix for DE thresholds to be tested##############################

logFC = seq(0,4,by = 0.1)

#Get desired Gene List:
CLUST = matrix(0, nrow = length(ngenes), ncol =length(logFC))  
rownames(CLUST) = ngenes                
colnames(CLUST) = sapply(1:length(logFC), function(x) paste0('Log2FC_',logFC[x]))

for (h in 1:length(logFC)) {
  ##Get Genes
  n = length(ids)
  N = length(ngenes)
  genesets.up = matrix(0, ncol=n, nrow= N ) 
  for( i in 1:n ){
    id = ids[i]
    filename = paste0(id,".DE.Rdata") 
    load(filename)
    m = match(ngenes, rownames(res) )  
    f.a = !is.na(m)
    f.r = m[f.a] 
    genes.up = (res$log2FoldChange >= logFC[h] & res$padj < 0.05 ) *1 # Originally included the filter res$baseMean > 20, therefore plots are different then to those in paper
    genesets.up[f.a,i] = genes.up[f.r]
  }
  colnames(genesets.up) = ids
  rownames(genesets.up) = ngenes
  genes100 = rownames(genesets.up[rowSums(genesets.up) > 2,])
  m = match(rownames(CLUST), genes100)
  f.a = !is.na(m)
  CLUST[f.a,h] = 1
}

########For Mouse-Human id Conversion##########

m = match(rownames(CLUST),mouse_ortho$Gene.stable.ID,)
f.a = !is.na(m)
f.t =m[f.a]
rownames(CLUST)[f.a] = as.character(mouse_ortho[f.t,]$Human.gene.stable.ID)

########################################################################

#Make plotting dataframe

roc_CLUST = matrix(0,nrow = nrow(ranking_list),ncol = ncol(CLUST))
rownames(roc_CLUST) = rownames(ranking_list)
colnames(roc_CLUST) =colnames(CLUST)
m =match(row.names(roc_CLUST), rownames(CLUST))
f.a =!is.na(m)
f.t =m[f.a]
roc_CLUST[f.a,] = CLUST[f.t,]

#restable = list()

j =4
restable[[j]] = matrix(0,nrow = length(logFC), ncol = 3)
rownames(restable[[j]]) = logFC
colnames(restable[[j]]) =c("AUROC","Gene Totals","log2FC Thresholds")
restable[[j]][,2] = colSums(roc_CLUST)
for (i in 1:nrow(restable[[j]])) {
  roc = data.frame(roc_CLUST[,i]) %>% as.matrix()
  restable[[j]][i,1] = as.numeric(auroc_analytic(ranking_list[,4], roc))
  restable[[j]][i,3] = logFC[i]
}

#Use Mouse "mf" with human "CLUST", vice versa.

restable[[1]] = as.data.frame(restable[[1]])
restable[[2]] = as.data.frame(restable[[2]])
restable[[3]] = as.data.frame(restable[[3]])
restable[[4]] = as.data.frame(restable[[4]])

save(restable, file = "Human Mouse DE predicitibilities.rdata")

#Plots
load("Human Mouse DE predicitibilities.rdata")
color <- c("Mouse" = "green3", "Human" = "red3")

p = list()
p[[1]] = ggplot() +
  geom_point(data = restable[[1]],          
             aes(x = `log2FC Thresholds`, y = AUROC, color ="Mouse"),
             size =2) +
  geom_line(data = restable[[1]],                                       
            aes(x = `log2FC Thresholds`, y = AUROC,color ="Mouse")) +                       
  geom_point(data = restable[[2]],                                  
             aes(x = `log2FC Thresholds`, y = AUROC, colour="Human"),
             size = 2) +
  geom_line(data = restable[[2]],  
            aes(x = `log2FC Thresholds`, y = AUROC, colour="Human")) +                       
  labs(color = "DE Predictibility") + 
  ylab("AUROC") +
  geom_vline(aes(xintercept = restable[[1]][restable[[1]]$`log2FC Thresholds` == "1", 3]), colour = "green3", linetype = 2) +
  geom_vline(aes(xintercept = restable[[2]][restable[[2]]$`log2FC Thresholds` == "1.5", 3]), colour = "red3", linetype = 2) +
  scale_color_manual(values = color) +
  theme_classic() +
  theme(legend.position = "none",
        legend.direction = "horizontal")

legend = cowplot::get_legend(p[[1]])

p[[2]] = ggplot() +
  geom_point(data = restable[[1]],          
             aes(x = `Gene Totals`, y = AUROC, color ="Mouse"),
             size =2) +
  geom_line(data = restable[[1]],                                       
            aes(x = `Gene Totals`, y = AUROC,color ="Mouse")) +  
  geom_point(data = restable[[2]],                                  
             aes(x = `Gene Totals`, y = AUROC, colour="Human"),
             size = 2) +
  geom_line(data = restable[[2]],  
            aes(x = `Gene Totals`, y = AUROC, colour="Human")) +                       
  labs(color = "DE Predictibility") +  
  ylab("AUROC") +
  geom_vline(aes(xintercept = restable[[1]][restable[[1]]$`log2FC Thresholds` == "1", 2]), colour = "green3", linetype = 2) +
  geom_vline(aes(xintercept = restable[[2]][restable[[2]]$`log2FC Thresholds` == "1.5", 2]), colour = "red3", linetype = 2) +
  scale_color_manual(values = color) +
  theme_classic() +
  theme(legend.position = "none",
        legend.direction = "horizontal")

p[[3]] = ggplot() +
  geom_point(data = restable[[3]],          
             aes(x = `log2FC Thresholds`, y = AUROC, color ="Mouse"),
             size =2) +
  geom_line(data = restable[[3]],                                       
            aes(x = `log2FC Thresholds`, y = AUROC,color ="Mouse")) +                       
  geom_point(data = restable[[4]],                                  
             aes(x = `log2FC Thresholds`, y = AUROC, colour="Human"),
             size = 2) +
  geom_line(data = restable[[4]],  
            aes(x = `log2FC Thresholds`, y = AUROC, colour="Human")) +                       
  labs(color = "DE prior") + 
  ylab("AUROC") +
  geom_vline(aes(xintercept = restable[[3]][restable[[3]]$`log2FC Thresholds` == "1", 3]), colour = "green3", linetype = 2) +
  geom_vline(aes(xintercept = restable[[4]][restable[[4]]$`log2FC Thresholds` == "1.5", 3]), colour = "red3", linetype = 2) +
  scale_color_manual(values = color) +
  theme_classic() +
  theme(legend.position = "none",
        legend.direction = "horizontal")

p[[4]] = ggplot() +
  geom_point(data = restable[[3]],          
             aes(x = `Gene Totals`, y = AUROC, color ="Mouse"),
             size =2) +
  geom_line(data = restable[[3]],                                       
            aes(x = `Gene Totals`, y = AUROC,color ="Mouse")) +                       
  geom_point(data = restable[[4]],                                  
             aes(x = `Gene Totals`, y = AUROC, colour="Human"),
             size = 2) +
  geom_line(data = restable[[4]],  
            aes(x = `Gene Totals`, y = AUROC, colour="Human")) +                       
  labs(color = "DE prior") +  
  ylab("AUROC") +
  geom_vline(aes(xintercept = restable[[1]][restable[[1]]$`log2FC Thresholds` == "1", 2]), colour = "green3", linetype = 2) +
  geom_vline(aes(xintercept = restable[[2]][restable[[2]]$`log2FC Thresholds` == "1.5", 2]), colour = "red3", linetype = 2) +
  scale_color_manual(values = color) +
  theme_classic() +
  theme(legend.position = "none",
        legend.direction = "horizontal")

legend2 = cowplot::get_legend(p[[4]])
g1 = plot_grid(p[[1]], p[[2]],ncol = 2, labels = LETTERS[1:2])
g1 = plot_grid(g1,legend, nrow = 2,rel_heights = c(9,1))
g2 = plot_grid(p[[3]], p[[4]],ncol = 2, labels = LETTERS[3:4])
g2 = plot_grid(g2,legend2, nrow = 2,rel_heights = c(9,1))
g3 = as.ggplot(~ scatterhist(x = human,y = mouse, 
                             xlab = "Average Human Log2FC", 
                             ylab = "Average Mouse Log2FC", 
                             main = NULL))
b =plot.new()
g3 = plot_grid(b, g3, b, ncol = 3, rel_widths = c(1,2,1), labels = c("","E", ""))
g = plot_grid(g1,g2,g3, nrow = 3)
pdf("D:/Thesis figures/Figure 2.pdf", height = 10,width = 7, useDingbats = F)
g
dev.off()
