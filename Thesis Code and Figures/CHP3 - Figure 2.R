##Date:7/5/21
##Written By: Conor Cremin
##Title: CHiPSeeker

#Set-up:

if( !require("DiffBind")){
  BiocManager::install("DiffBind")
}
if( !require("cowplot")){
  BiocManager::install("cowplot")
}
if( !require("ggplot2")){
  BiocManager::install("ggplot2")
}
if (!require("genomation")){
  BiocManager::install("genomation")
}
if (!require("ChIPseeker")){
  BiocManager::install("ChIPseeker")
}
if (!require("ggimage")){
  BiocManager::install("ggimage")
}
if (!require("dplyr")){
  BiocManager::install("dplyr")
}
if (!require("GenomicRanges")){
  BiocManager::install("GenomicRanges")
}
if (!require("clusterProfiler")){
  BiocManager::install("clusterProfiler", ask = F)
}
if (!require("ReactomePA")){
  BiocManager::install("ReactomePA", ask = F)
}
if (!require("org.Hs.eg.db")){
  BiocManager::install("org.Hs.eg.db", ask = F)
}
if (!require("TxDb.Hsapiens.UCSC.hg38.knownGene")){
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", ask = F)
}
if (!require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap", ask = F)
}  
if (!require("VennDetail")){
  BiocManager::install("VennDetail", ask = F)
}  

setwd('E:/')
##After Differential Analysis and BEDfile generation:
load("E:/gene_annotations_v29.Rdata")
allgenes = read.csv("E:/allgenes_Ensembl_v97.csv", stringsAsFactors = F)
load('cross-species.clusters.rdata')

#Plots - ChipSeeker
study = "HKUStudy_CHIP"
setwd(paste0(study))
files = list.files(pattern = paste0("narrowPeak"))
pol2 = files[c(1,12,10,11)]
names(files) <- matrix(unlist(strsplit(pol2,"_peaks")), ncol = 2, byrow = T)[,1]
names(pol2) = c('Mock','Wildtype','3M','dNS1')

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnnoList <- lapply(pol2, annotatePeak, TxDb=txdb,level = "gene", tssRegion=c(-500, 500), verbose=FALSE)

p =plotAnnoBar(peakAnnoList,xlab = "", ylab = "Percentage(%)", title = "") +
  theme_classic() +
  scale_fill_viridis_d(option = 'A', guide = guide_legend(reverse = TRUE) ) +
  theme(legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.title = element_text(size=10, face = "bold"),
        legend.text = element_text(size=8))

q = plotDistToTSS(peakAnnoList, title = "",) +
  theme_classic() +
  scale_fill_viridis_d(option = 'A', guide = guide_legend(reverse = TRUE) ) +
  theme(legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.key.height = unit(4, "mm"),
        legend.key.width = unit(4, "mm"), 
        legend.title = element_text(size=10, face = "bold"),
        legend.text = element_text(size=8))



##Plot 3

setwd("/home/hlchen/Conor/HKUStudy_CHIP")

for (i in 1:6) {
  bedinter = '/software/BEDTools/2.27.1/bin/intersectBed'
  bigSum = '/software/deeptools/3.1.3/bin/multiBigwigSummary'
  bed = paste0('/home/hlchen/Conor/BED/9K/','Clusters_',i,'.9k.bed')
  cmd1 = paste('sort -k1,1 -k2,2n',bed,">", paste0('/home/hlchen/Conor/BED/9K/','Clusters_',i,'.9k.sorted,bed'))
  system(cmd1)
  bigs = 'BigWigs2/Pol2_Mock.bw BigWigs2/Pol2_WSNWT.bw BigWigs2/Pol2_WSN3M.bw BigWigs2/Pol2_WSNdelNS1.bw'
  cmd2 = paste(bigSum,"BED-file --bwfiles", bigs,
               "--BED", paste0('/home/hlchen/Conor/BED/9K/','Clusters_',i,'.9k.sorted,bed'),
               "--numberOfProcessors 10", paste0("--outRawCounts BED_Coverages/cluster_",i,".txt"),
               "--labels Mock WSNWT WSN3M WSNdelNS1 -o BED_Coverages/allgenes-TSS.npz")
  system(cmd2)
  cmd3 = paste(bedinter, "-wa -wb -f 1.0 -r -a",paste0("BED_Coverages/cluster_",i,".txt"),
               '-b', paste0('/home/hlchen/Conor/BED/9K/','Clusters_',i,'.9k.sorted,bed'),
               ">", paste0("BED_Coverages/cluster",i,".bed"))
  system(cmd3)
}

files = list.files("E:/HKUStudy_CHIP/Bed_Coverages/Pol Densities", full.names = T)
names(files) = files %>% gsub(".bed", "", .) %>% gsub("E:/HKUStudy_CHIP/Bed_Coverages/Pol Densities/", "", .)
reads = files %>% lapply(., function(x) read.delim(x, header = F, check.names = F,stringsAsFactors = F))
#files = list.files(pattern = "TSS.bed")
#names(files) = files %>% gsub(".TSS.bed", "", .)
#tss = files %>% lapply(., function(x) read.delim(x, header = F, check.names = F,stringsAsFactors = F))

for (i in 1:length(reads)) {
  #colnames(tss[[i]]) = c('chr','start','end','Mock','WSNWT','WSN3M','WSNdelNS1','chr','start','end','id','coff','strand')
  #tss[[i]] = tss[[i]][!duplicated(tss[[i]]$id),]
  #print(nrow(tss[[i]]))
  #tss[[i]] = na.omit(tss[[i]])
  #print(nrow(tss[[i]]))
  colnames(reads[[i]]) = c('chr','start','end','Mock','WSNWT','WSN3M','WSNdelNS1','chr','start','end','id','coff','strand')
  reads[[i]] = reads[[i]][!duplicated(reads[[i]]$id),]
  print(nrow(reads[[i]]))
  reads[[i]] = na.omit(reads[[i]])
  print(nrow(reads[[i]]))
}

i=1
for (i in 1:length(reads)) {
  reads[[i]][,c("Mock", "WSNWT","WSN3M", "WSNdelNS1")][reads[[i]][,c("Mock", "WSNWT","WSN3M", "WSNdelNS1")] <=0] = 0
  #Set negative coverage values to 0 and add pseudocount of 1:
  reads[[i]] = reads[[i]][,c("id", "Mock", "WSNWT","WSN3M", "WSNdelNS1")]
  reads[[i]][,c(2:5)] = reads[[i]][,c(2:5)]
  print(nrow(reads[[i]]))
}

COL = c("Mock" = "tan1","Wildtype" = "firebrick1", "3M"="blue2", "dNS1" = "limegreen")

df = list()
for (i in 1:length(unique(clusters$unmergedLabels.mod))) {
  df[[i]] = suppressMessages(reshape2::melt(reads[[i]][,2:5]) )%>% as.data.frame() 
  colnames(df[[i]]) = c("Condition", "log2CPM")
  df[[i]]$Cluster = i
}
df = do.call(rbind.data.frame, df) 
df$Condition = gsub("WSNWT", "Wildtype", df$Condition)
df$Condition = gsub("WSN3M","3M", df$Condition)
df$Condition = gsub("WSNdelNS1","dNS1", df$Condition)
df$Condition = factor(df$Condition, levels = c("Mock","Wildtype", "3M", "dNS1"))
df$Cluster = factor(df$Cluster, levels = c(1:6))


r = ggplot(df, aes(Condition, log2CPM, fill = Condition)) + 
  geom_boxplot(outlier.shape = NA)  + 
  labs(x="Conditions", 
       y="Log2(CPM)") +
  facet_grid(~ Cluster) +
  theme_classic()+
  geom_signif(map_signif_level = c("***"=0.001, "**"=0.01, "*"=0.05) , textsize=3, 
              comparisons = list(c("Mock", "Wildtype"),
                                 c("Wildtype", "3M"),
                                 c("Wildtype", "dNS1")),
              show.legend = T, 
              test = "t.test",
              vjust = 0.1,
              step_increase = 0.1,
              margin_top = -0.7) +
  scale_fill_manual(values = COL, 
                    name="Condition") +
  theme(axis.text.x.bottom = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom',
        legend.direction = 'horizontal')

ylim1 = boxplot.stats(df$log2CPM)$stats[c(1, 5)]

# scale y limits based on ylim1
r = r + coord_cartesian(ylim = ylim1*7)


g <- ggplot_gtable(ggplot_build(r))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- unique(h.clusters$colors)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

r1 = plot_grid(p,q, nrow = 1, labels = c(LETTERS[1:2]), label_fontfamily = 'serif', label_size = 18)

pdf('Thesis figures/Chp3-Figure 2a.pdf', useDingbats = F, width = 8, height = 6, paper = 'a4')
plot_grid(r1, g, rel_heights = c(1,1.5), nrow = 2)
dev.off()

#Motif plotting.R used for Motif plots (Homer output)


