
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

####################Volplots for human####################

##Set-up:

setwd('E:/')
ids = c("GSE84204","GSE89008", "GSE97672", "GSE103477", "GSE103604", "GSE104168", "HKUStudy", "HKUStudy_293") # list real studies here 
ids = c("GSE49933", "GSE52405", "GSE100522", "GSE107488", "GSE117029", "ERP020504", "SRP061303") # list real studies here 
mylist = list()

#Get DESeq2 output and formatting:

for (i in 1:length(ids)) {
  id = ids[i]
  filename = paste0(id,".DE.Rdata") 
  load(filename)
  mylist[[i]] = as.data.frame(res)
  mylist[[i]]$GeneID = rownames(mylist[[i]])
  mylist[[i]]$Study = id
  rownames(mylist[[i]]) = NULL
}  
mylist = do.call(rbind.data.frame, mylist)

#For human
mylist[mylist$Study == 'HKUStudy_293',]$Study = 'GSE156152'
mylist[mylist$Study == 'HKUStudy',]$Study = 'GSE156060'
mylist$Study = factor(mylist$Study, levels = c("GSE84204","GSE89008", "GSE97672", "GSE103477", "GSE103604", "GSE104168", "GSE156060", "GSE156152"))

#Or for mouse
mylist$Study = factor(mylist$Study, levels = ids)

#Facetting Volplots:

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
                          col=c('black', 'green3', 'blue', 'red3'),
                          legendLabels = c("Not significant",
                                           "-1.5 > Log2FC > 1.5",
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
  facet_wrap(~Study, ncol = 4)

#Color Panel headers:

g <- ggplot_gtable(ggplot_build(myplots))

#For Human
stripr <- which(grepl('strip-t', g$layout$name))
#For Mouse:
stripr = stripr[-4]

fills <- rep("Red1", length(stripr))
fills <- rep("Green3", length(stripr))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

#Print output:

pdf("E:/Thesis figures/Figure 1C.pdf", height = 5, width = 7, useDingbats = F)
grid::grid.draw(g)
dev.off()
######################################################################
####################Barplot (Human)####################

load("E:/genesets_for_human_recurrence.rdata")
colnames(genesets.up)[c(7,8)] =c("GSE156060", "GSE156152")
colnames(genesets.down)[c(7,8)] =c("GSE156060", "GSE156152")

#Or Mouse
load("E:/recurrence_matrix_MOUSE.Rdata")

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

##Barplot 1: Mouse and Human influenza datasets:
g2 =ggplot(df,aes(x=factor(Study,levels = colnames(genesets.up) ),y=Totals,fill=Class)) +
  geom_col(color = "black")+
  labs(y = "Differentially expressed genes", x = "Studies") +
  theme_classic() + 
  scale_fill_manual("Class", values = c("Downregulated" = "DarkBlue","Upregulated" = "DodgerBlue1")) +
  geom_text(
    aes(y = Totals+(Totals/abs(Totals)*200),label = abs(Totals),group = Study), size = 4) +
  theme_classic() +
  theme(axis.text.x.bottom = element_text(angle = 90))


pdf("E:/Thesis figures/Figure 1A.pdf", height = 5, width = 7, useDingbats = F)
grid.draw(g)
dev.off()
pdf("E:/Thesis figures/Figure 1D.pdf", height = 5, width = 6, useDingbats = F)
grid.draw(g2)
dev.off()
