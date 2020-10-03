##############################

##Author Name: Conor Cremin
##Title:Meta-analysis

#Requirements:

##############################

#For Human:
allgenes = read.csv("Human_Annotation_Ensembl_V97.csv", header = T, stringsAsFactors = F)
#For mouse:
allgenes = read.csv('Mouse_Annotation_Ensembl_v99.csv', header = T, stringsAsFactors = F)

# Packages:

library(gplots)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(limma)
library(pheatmap)
library(gridExtra)
library(grid)
library(ggplotify)
library(EnhancedVolcano)
library(extrafont)
library(tidyr)
library(cowplot)
library(fields)
library(heatmap3)


# Step 2: Designate Studies for analysis and create empty matrices for results.

i = 1 
ids  = c() #List study ids here. This should be an identifier "STUDY_ID_IS_HERE.DE.Rdata" of the deseq2 output from Differential Expression.R, to extract the ids for DEGs per study.
#For human
#ids = c("GSE84204","GSE89008", "GSE97672", "GSE103477", "GSE103604", "GSE104168", "GSE156060", "GSE156152")
#Or for mouse
#ids = c("GSE49933","GSE52405", "GSE100522", "GSE107488", "GSE117029", "ERP020504", "SRP061303")  
n = length(ids)
ngenes = unique(allgenes$Gene.stable.ID)
N = length(ngenes)
genesets.up = matrix(0, ncol=n, nrow= N ) 
genesets.down = matrix(0, ncol=n, nrow= N ) 

# Step 3: Loop function for recurrence analysis across all data sets.

for( i in 1:n ){
  
  id = ids[i]
  filename = paste0(id,".DE.Rdata") #Ensure path of deseq2 output is correct!!!!!
  load(filename)
  
  m = match(ngenes, rownames(res) )  
  f.a = !is.na(m)
  f.r = m[f.a] 
  
  genes.up = (res$baseMean > 20 & res$log2FoldChange >= 1.5 & res$padj < 0.05 ) *1 #Upregulated DEGS
  genes.down = (res$log2FoldChange < -1.5 & res$padj < 0.05)*1  #Downregulated DEGS
  
  genesets.up[f.a,i] = genes.up[f.r]
  genesets.down[f.a,i] = genes.down[f.r]
}

colnames(genesets.up) = ids
rownames(genesets.up) = ngenes
colnames(genesets.down) = ids
rownames(genesets.down) = ngenes
genesets.up = genesets.up %>% as.data.frame()
genesets.down = genesets.down %>% as.data.frame()

# Step 4: Calculating Recurrence Statistics (Ensure helper_functions(Recurrence_Analysis).R is LOADED)

#Number of Significantly up-regulated genes per study:

colSums(genesets.up, na.rm = TRUE)
recur = rowSums(genesets.up, na.rm = TRUE)
fdrs = calc_fdrs_recur( genesets.up )     #Get function from https://github.com/sarbal/OutDeCo/tree/master/R

#Histogram of recurrence of Upregulated genes across studies. 

pdf("Cross-species recurrence.pdf", height=7, width=8)
df = data.frame(recur[recur>0]+0.01)
colnames(df) = "counts"
p = 0.35
ggplot(df, aes(x = counts )) +  
  labs(x = "Recurrence", y = "Number of Genes" ) +
  theme_classic() +
  geom_rect(aes(xmin = fdrs$Pt-p,xmax = Inf,ymin = -Inf, ymax = Inf),fill="DarkSeaGreen2") +
  geom_vline(aes(xintercept=fdrs$Pt-p), 
             color="black",
             linetype="dashed",
             show.legend = T) +
  annotate("text", x = 7.5, y = max(plyr::count(df[,1]))/2, 
           label = paste0("N (FDR < 0.05) = ", fdrs$sig), 
           size= 3, 
           color = 'firebrick', 
           angle=0 ) +
  annotate("text",x=fdrs$Pt-p, label=paste0("FDR < 0.05"), y=max(plyr::count(df[,1]))/2, 
           colour="firebrick",
           size=3,
           angle=90, 
           vjust = 1.2) +
  geom_bar(color = "black",width = 0.7,
           fill = "blue2") +
  scale_x_continuous(breaks=1:max(recur)) +
  theme(plot.margin = unit(c(2,1,2,1), "cm"))

dev.off()

