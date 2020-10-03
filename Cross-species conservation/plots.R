##########################

##Author name: Conor Cremin
##Title: Riverplots and Venndiagram
##Dataset: Influenza

##########################

#Packages:
library(limma)
library(grid)
library(dplyr)
library(gridExtra)
library(ggplotify)
library(cowplot)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
l

if( !require("tidyverse")){
  install.packages("tidyverse")
}
if( !require("riverplot")){
  install.packages("riverplot")
}

load("E:/mouse_human_high_confidence_orthos.rdata") # File of Human and Mouse IDs.

##For this example, we plotted the recurrence overlap between mouse and human recurrence.
##Get mouse recurrence and human recurrence dataframes genesets.up. 
##These are generated using the Meta-analysis.R script in the Differential Expression directory

##human.rec = rec = rownames(genesets.up[rowSums(genesets.up) >= 3,]) ##Get genes at or above recurrence threshold. Calculated in Meta-analysis.R
##mouse.rec = rec = rownames(genesets.up[rowSums(genesets.up) >= 3,]) ##Get genes at or above recurrence threshold. Calculated in Meta-analysis.R

m = match(mouse.rec = mouse_ortho$Gene.stable.ID)
f.a = !is.na(m)
f.t = m[f.a]
mouse.rec[f.a] = mouse_ortho[f.t,]$Human.gene.stable.ID  #Convert orthos to  human id form.

pdf("Figures/Venn.pdf", height = 7, width = 8)
venn.plot <- venn.diagram(
  x = list(
    Human = human.rec %>% unlist() , 
    Mouse = mouse.rec %>% unlist()
  ),
  category.names = c("Human Recurrence" , "Mouse Recurrence"),
  filename = NULL,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("red3",0.3), alpha('green3',0.3)),
  cex = 1,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 2,
  cat.default.pos = "outer",
  cat.fontface ="bold",
  cat.fontfamily = "sans",
  cat.col = c("red3", 'green3'),
)

grid.draw(venn.plot)
dev.off()


##Riverplot

m.clusters = mouse.clusters  #Get from clustering analysis
clusters = human_clusters  #Get from clustering analysis
m = match(m.clusters$V1, mouse_ortho$Gene.stable.ID)
f.a =!is.na(m)
f.t = m[f.a]
m.clusters$V1 = as.character(m.clusters$V1)
m.clusters[f.a,]$V1 = as.character(mouse_ortho[f.t,]$Human.gene.stable.ID)
res = list()
for (i in 1:length(unique(clusters$unmergedLabels.mod))) {
  m = match(m.clusters$V1,clusters[clusters$unmergedLabels.mod == i,]$V1)
  f.a =!is.na(m)
  g= table(m.clusters[f.a,]$unmergedLabels.mod)
  human_cluster = rep(paste0("H",i),length(g))
  mouse_cluster = as.list(names(g)) %>% sapply(., function(x) paste0("M",x))
  res[[i]] = cbind.data.frame(human_cluster, mouse_cluster, g) %>% .[,-3]
}

edges = do.call(rbind.data.frame,res)
colnames(edges) =c("N1","N2","Value")
nodes = data.frame(ID = unique(c(as.character(edges$N1), as.character(edges$N2))), stringsAsFactors = FALSE)
nodes$x = c(rep(1,length(unique(c(as.character(edges$N1))))), c(2,2,2))
nodes$y = c(0,0.9,0,0.9,1)
nodes$col = c(rep("red3", length(unique(c(as.character(edges$N1))))), rep("green3", length(unique(c(as.character(edges$N2))))))
rownames(nodes) = nodes$ID
palette = paste0(brewer.pal(2, "RdYlGn"), "70")
edges$col = palette[factor(edges$N1)]
#edges$edgecol = "col"
style <- list(nodestyle= "regular", 
              srt = 0,
              lty = 2,
              textcol = "black",
              textcex = 2,
              edgecol = "gradient",
              edgestyle= "sin")
river <- makeRiver( nodes, edges)
pdf("human and mouse riverplot.pdf", width = 5, height = 5, useDingbats = F)
riverplot( river,mar = c(0, 0, 0, 0), plot_area = 0.8, yscale = 0.017,default_style= style )
dev.off()

#Add Jaccard Index to plots:

jac = edges 
jac$N1 = gsub("H","", jac$N1)
jac$N2 = gsub("M","", jac$N2)
human =sapply(jac$N1, function(x) nrow(clusters[clusters$unmergedLabels.mod == x,]))
mouse =sapply(jac$N2, function(x) nrow(m.clusters[m.clusters$unmergedLabels.mod == x,]))
jac$union = human+mouse-jac$Value
jac$jaccard = jac$Value/jac$union
jac$labels = c("H1->M1","H1->M2","H1->M3","H2->M1","H2->M3")
jac$`Recurrent Proportion (%)` = (jac$Value/fdrs$sig)
pdf("Jaccard Index.pdf", height=5, width=5, useDingbats = F)
ggplot(data = jac)+
  labs(y = "Jaccard Index") +
  geom_point(aes(x =`Recurrent Proportion (%)`, y = jaccard),colour="black",pch=19, size=3) +
  geom_text(aes(x =`Recurrent Proportion (%)`, y = jaccard,label = labels), nudge_y = -0.03,nudge_x = 0.02,size = 2) +
  theme_classic() +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()
