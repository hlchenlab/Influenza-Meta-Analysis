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
load("DE_prior.rdata") 
allgenes = read.csv("allgenes_Ensembl_v97.csv", header = T, stringsAsFactors = F)

################For Mouse#####################

load("GO Annotation Jan2020/GO.mouse.Rdata")
load('mouse_gene_annotations_v100.Rdata')
load("recurrence_matrix_MOUSE.Rdata")
load("Mouse_DE_prior.rdata")
allgenes = read.delim('mouse_allgenes.txt', header = T, stringsAsFactors = F)

##Plot 1 - Recurrence
df.fdr = as.data.frame(rep(0, 4))
colnames(df.fdr) = 'Species'
df.fdr$`Total` = rep(0,4)
df.fdr$`Differential Expression` = rep(0,4)
df.fdr$counts = rep(0,4)

res = list()

load("genesets_for_human_recurrence.rdata")
#or
load("recurrence_matrix_MOUSE.Rdata")

genesets = genesets.up
recur = rowSums(genesets, na.rm = TRUE)
fdrs = calc_fdrs_recur(genesets)
df.fdr[3,1] = 'Mouse'
df.fdr$`Differential Expression`[3] = 'Upregulated'
df.fdr$`Total`[3] = fdrs$sig
df.fdr$counts[3] = fdrs$Pt
df = data.frame(recur[recur>0])
res[[3]] = data.frame(table(df), row.names = NULL)
colnames(res[[3]]) = c("Recurrence","counts")
res[[3]]$Species = 'Mouse'
res[[3]]$`Differential Expression` = 'Upregulated'


genesets = genesets.down
recur = rowSums(genesets, na.rm = TRUE)
fdrs = calc_fdrs_recur(genesets)
df.fdr[4,1] = 'Mouse'
df.fdr$`Differential Expression`[4] = 'Downregulated'
df.fdr$`Total`[4] = fdrs$sig
df.fdr$counts[4] = fdrs$Pt
df = data.frame(recur[recur>0])
res[[4]] = data.frame(table(df), row.names = NULL)
colnames(res[[4]]) = c("Recurrence","counts")
res[[4]]$Species = 'Mouse'
res[[4]]$`Differential Expression` = 'Downregulated'
res[[4]]$counts = res[[4]]$counts*-1

res1 = do.call(rbind.data.frame, res)
res1$`Differential Expression` = factor(res1$`Differential Expression`, levels = c('Upregulated','Downregulated'))
colnames(df.fdr)[4] = 'Recurrence'
df.fdr$counts = rep(0, 4)
res1$Recurrence = as.numeric(res1$Recurrence)
p =ggplot(data = res1, aes(x = Recurrence, y = counts,fill = `Differential Expression`)) +  
  geom_rect(data = df.fdr, aes(xmin = Recurrence-.5,
                               xmax = Inf,
                               ymin = -Inf, 
                               ymax = Inf),fill="DarkSeaGreen2") +
  geom_vline(data = df.fdr, aes(xintercept=Recurrence-.5, color = 'FDR <= 0.05'), 
             linetype="dashed",
             size = 1) +
  geom_col(color = "black") +
  labs(x = "Recurrence", y = "Number of Genes" ) +
  theme_classic() +
  scale_fill_manual('Differential expression', values = c("Upregulated" = "#FF3333","Downregulated" = "#33CCFF")) +
  geom_text(aes(y = counts+(counts/abs(counts)*250),label = abs(counts),group = Recurrence), 
            size = 3) +
  facet_grid(~ Species, scales = 'free') +
  scale_color_manual('Recurrence threshold', values = c("FDR <= 0.05" = "red3")) +
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom') +
  scale_x_continuous(breaks = 1:8)

g <- ggplot_gtable(ggplot_build(p))

stripr <- which(grepl('strip-t', g$layout$name))

fills <- c("Red1","Green3")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

##Plot 2 - DEprior v recurrence

load("DE_prior.rdata")
load("genesets_for_human_recurrence.rdata")
load("gene_annotations_v29.Rdata")
#or
load("recurrence_matrix_MOUSE.Rdata")
load("Mouse_DE_prior.rdata")
load('mouse_gene_annotations_v100.Rdata')
DE_list$MF.rank = DE_list$MF.rank/max(DE_list$MF.rank)
genesets = genesets.up
recur = rowSums(genesets, na.rm = TRUE)
gmat = matrix(0, ncol = 4, nrow = nrow(DE_list))
rownames(gmat) = rownames(DE_list)
m = match(rownames(gmat), names(recur))
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,1] = recur[f.t]
gmat[,2] = (DE_list[,4])
gmat = as.data.frame(gmat[gmat[,1] >= 1,])
m = match(rownames(gmat), attrm$ensemblID)
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,3] = as.character(attrm[f.t,]$name)
gmat = as.data.frame(gmat)
colnames(gmat) = c('Recurrence', 'DE prior', 'Name', 'Species')
gmat$Species = 'Mouse'

df = rbind.data.frame(gmat, gmath)
gmath = gmat
m = df[df[,1] == 7 & df$Species == 'Mouse',]
h = df[df$Recurrence == 8,]
ggplot() + 
  geom_boxplot(data = df, aes(x= factor(Recurrence, levels = 0:max(Recurrence)), 
                                y = `DE prior`), 
               width = 0.5, outlier.shape = NA,fill = '#FF3333') +
  geom_jitter(data =df, aes(x= factor(Recurrence, levels = 0:max(Recurrence)), 
                              y = `DE prior`), 
              width = 0.2, alpha = .02) +
  geom_text(data = h,
            aes(x= factor(Recurrence, levels = 0:max(Recurrence)), 
                y = `DE prior`,label = Name), 
            check_overlap = T,
            size = 2) +
  geom_text(data = m,
            aes(x= factor(Recurrence, levels = 0:max(Recurrence)), 
                y = `DE prior`,label = Name), 
            check_overlap = T,
            size = 2) + 
  facet_grid(~ Species, scales = 'free') +
  geom_vline(data = df, aes(xintercept=6.5, color = 'FDR <= 0.05'), 
             linetype="dashed",
             size = 1) +
  scale_color_manual('Recurrence threshold', values = c("FDR <= 0.05" = "red3")) +  labs(x="Recurrence", 
       y="DE prior") +
  theme_classic() +
  theme(legend.direction = 'horizontal',
        legend.position = 'bottom')
g2 <- ggplot_gtable(ggplot_build(p))

stripr <- which(grepl('strip-t', g2$layout$name))

fills <- c("Red1","Green3")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#Plot 3 - Gene types per recurrence

genesets = genesets.up
recur = rowSums(genesets, na.rm = TRUE)
gmat = matrix(0, ncol = 2, nrow = nrow(genesets)) %>% as.data.frame()
rownames(gmat) = rownames(genesets)
colnames(gmat) = c('Recurrence','Genetype')
m = match(rownames(gmat), names(recur))
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,]$Recurrence = recur[f.t]
gmat = gmat[gmat$Recurrence != 0,]
m = match(rownames(gmat), attr$ensemblID)
f.a = !is.na(m)
f.t = m[f.a]
gmat[f.a,]$Genetype = as.character(attr[f.t,]$type)
gmat[gmat$Genetype == 0 ,]$Genetype = 'Unknown'
tf = lapply(1:max(gmat$Recurrence), function(x) data.frame(table(gmat[gmat$Recurrence == x,]$Genetype)/sum(table(gmat[gmat$Recurrence == x,]$Genetype)), row.names = NULL))

for (i in 1:length(tf)) {
  k = tf[[i]]
  k$Recurrence = i
  k$Sp = 'Human'
  colnames(k) = c('Genetype', 'Proportion','Recurrence', 'Species')
  tf[[i]] = k
}
gmat = do.call(rbind.data.frame, tf)
gmatm = gmat
lvls = c('Protein coding', 'Pseudogene', 'lncRNA', 'Other ncRNA','snRNA', 'miRNA', 'rRNA','scRNA', 'Unknown')
df1 = rbind.data.frame(gmat, gmatm)
df1$Genetype = factor(df$Genetype, levels = lvls)
q =ggplot(data = df1, aes(x = factor(Recurrence, levels = 1:8), y = Proportion, fill =Genetype)) +
  geom_bar(position="stack", stat="identity", color = 'black') +
  scale_fill_brewer(palette = 'PuBu',
                    direction = -1) +
  theme_classic() +
  labs(x = "Recurrence", y = "Gene-type proportion") +
  facet_grid(~ Species, scales="free_x") +
  geom_vline(data = df1, aes(xintercept=2.5, color = 'FDR <= 0.05'), 
             linetype="dashed",
             size = 1) +
  scale_color_manual('Recurrence threshold', values = c("FDR <= 0.05" = "red3")) +  
  labs(x="Recurrence",y="DE prior") +
  theme(legend.position = "bottom")
g3<- ggplot_gtable(ggplot_build(q))

stripr <- which(grepl('strip-t', g3$layout$name))

fills <- c("Red1","Green3")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g3$grobs[[i]]$grobs[[1]]$childrenOrder))
  g3$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf("Thesis figures/Figure 4.pdf" ,height = 11, width = 9, useDingbats = F)
plot_grid(g, g2,g3,  nrow = 3, 
          labels = LETTERS[1:3], 
          label_size = 18,
          label_fontfamily = 'serif')
dev.off()

save(res1, df, df1, file = 'meta-analysis-dataframes-fig 4.rdata')
