##Packages

if( !require("EGAD")){
  BiocManager::install("EGAD")
}
if( !require("ggplot2")){
  install.packages("ggplot2")
}
if( !require("ggplotify")){
  BiocManager::install("ggplotify")
}
if( !require("cowplot")){
  BiocManager::install("cowplot")
}
if( !require("dplyr")){
  BiocManager::install("dplyr")
}

##Set-up
setwd("E:/")
load("coexpression network/agg.reranked.Rdata")
load("coexpression network/Mouse Co-expression network SB/agg.rank.ENSEMBL_IDs.ranked.Rdata")

orthos = read.csv("orthos.csv", stringsAsFactors = F)
#orthos = orthos[orthos$Ferret.homology.type == "ortholog_one2one",]
#orthos = orthos[orthos$Golden.Hamster.homology.type == "ortholog_one2one",]
orthos = orthos[orthos$Mouse.homology.type == "ortholog_one2one",]

orthos = orthos[,c(1,9)]
colnames(orthos) = c('human','mouse')


load("agg.reranked.Rdata")
m = match(orthos$human, rownames(agg.rank))
f.a = !is.na(m)
orthos = orthos[f.a,]
load('agg.rank.ENSEMBL_IDs.ranked.Rdata')
m = match(orthos$mouse, rownames(reranked))
f.a = !is.na(m)
orthos = orthos[f.a,]

network = agg.rank[orthos$human, orthos$human]
networkh = matrix(rank(network ,na.last = "keep", ties.method = "average"), nrow=dim(network)[1], ncol=dim(network)[2])
rownames(networkh) = rownames(network)
colnames(networkh) = rownames(network)

network = reranked[orthos$mouse, orthos$mouse]
networkm = matrix(rank(network ,na.last = "keep", ties.method = "average"),nrow=dim(network)[1], ncol=dim(network)[2])
rownames(networkm) = rownames(network)
colnames(networkm) = rownames(network)

networkh = networkh/max(networkh, na.rm = T)
networkm = networkm/max(networkm, na.rm = T)

m = match(rownames(networkm), orthos$mouse)
f.a = !is.na(m)
f.t = m[f.a]
rownames(networkm)[f.a] = orthos[f.t,]$human
colnames(networkm)[f.a] = orthos[f.t,]$human

rm(agg.rank, reranked)

load("GO Annotation Jan2020/GO.human.Rdata")
load("gene_annotations_v29.Rdata")

gosums = colSums(GO.human.nonIEA )
filtGO = gosums <=1000 & gosums>=20
annot =  GO.human.nonIEA[,filtGO]
m = match(rownames(annot), attr$entrezID)
f.a = !is.na(m)          
f.t = m[f.a]
rownames(annot)[f.a] = as.character(attr[f.t,]$ensemblID)
voc = voc[voc$V3 =="biological_process",]
m = match(colnames(annot), voc$V1)
f.a = !is.na(m)
annot = annot[,f.a]

#Make Annotation for GO:

m = match(rownames(annot), rownames(networkh))
f.a = !is.na(m)
CLUSTH = annot[f.a,]

roc_aggm = run_GBA(networkm, CLUSTH, min = 0, max = 30000)
roc_aggh = run_GBA(networkh, CLUSTH, min = 0, max = 30000)

load('Human_GO.ortholog_networks.Rdata')

filt <- !is.na(roc_aggh[[1]][,1])
aucA <- roc_aggh[[1]][filt,1] 
filt <- !is.na(roc_aggm[[1]][,1])
aucB <- roc_aggm[[1]][filt,1]
aucA_dens <- density(aucA, adjust = 2)
aucA_dens$y <- aucA_dens$y/(max(aucA_dens$y))
aucB_dens <- density(aucB, adjust = 2)
aucB_dens$y <- aucB_dens$y/(max(aucB_dens$y))
xlim <- range(aucA_dens$x, aucB_dens$x)
ylim <- range(aucA_dens$y, aucB_dens$y)
a <- aucA_dens$x[which.max(aucA_dens$y)]
b <- aucB_dens$x[which.max(aucB_dens$y)]
box <- cbind(c(a, a, b, b), c(0, 1, 1, 0))

p1 = ggplot() +
  geom_rect(aes(xmin = b, xmax = a, ymin = 0, ymax= 1), fill = "grey84") +
  geom_line(aes(x = aucA_dens$x, y = aucA_dens$y), color = 'red3',size = 1) +
  geom_line(aes(x = aucB_dens$x, y = aucB_dens$y), color = 'green3', size = 1) + 
  geom_vline(aes(xintercept = a), linetype = 2, color = 'red3', size = 1) +
  geom_vline(aes(xintercept = b), linetype = 2, color = 'green3',size = 1) +
  labs(x = 'AUROCs', y = 'Density') +
  theme_classic()
  
df = cbind.data.frame(aucA,aucB)
df$Name = rownames(df)

p2 = ggplot(df,aes(x = df[,1], y = df[,2])) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept =0 ),size =1 ,color = "red3",show.legend = FALSE) +
  geom_hline(aes(yintercept = mean(df[,2])),linetype="dashed", size = 1,color ="green3") +
  geom_vline(aes(xintercept = mean(df[,1])),linetype="dashed", size = 1,color ="red3") +
  labs(x =  "Human AUROCs (NV)", y = "Mouse AUROCs (NV)") +
  scale_x_continuous(limits =c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_classic() 


#Or
df2 = cbind.data.frame(roc_aggh[[1]][,3], roc_aggm[[1]][,3])

p2 = ggplot(df2,aes(x = df2[,1], y = df2[,2])) +
  geom_point() +
  geom_abline(aes(slope = 1, intercept =0 ),size =1 ,color = "red3",show.legend = FALSE) +
  geom_hline(aes(yintercept = mean(df2[,2])),linetype="dashed", size = 1,color ="gray80") +
  geom_vline(aes(xintercept = mean(df2[,1])),linetype="dashed", size = 1,color ="gray80") +
  labs(x = "Human AUROCs (ND)", y = "Mouse AUROCs (ND)") +
  scale_x_continuous(limits =c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_classic() 


ndh = node_degree(agg.rank) %>% as.data.frame()
ndh$Network = 'Human'
ndm = node_degree(reranked) %>% as.data.frame()
ndm$Network = 'Mouse'
nd = rbind.data.frame(ndh,ndm)

dat = data.frame(v1 =c(10950.75, 8534.25), Network = c('Human','Mouse'))
l =ggplot(data = nd, aes(x = `.`)) +
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity") + 
  theme_classic() +
  labs(x = 'Node degree', y = 'Density') +
  geom_vline(data = dat,aes(xintercept=v1),
             color="red3", linetype = 3, size=1) +
  facet_grid(~Network) 
g2 <- ggplot_gtable(ggplot_build(l))

stripr <- which(grepl('strip', g2$layout$name))

fills <- c("Red1","Green3")
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g2$grobs[[i]]$grobs[[1]]$childrenOrder))
  g2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

r1 = plot_grid(p1, p2, labels = LETTERS[2:3], label_size = 18, label_fontfamily = 'serif')
pdf("Thesis figures/Chp2-Figure 1.pdf" ,height = 6, width = 7, useDingbats = F)
plot_grid(g2, r1, ncol = 1, labels = LETTERS[1], label_size = 18, label_fontfamily = 'serif' )
dev.off()



load("GO Annotation Jan2020/GO.mouse.Rdata")
load('agg.rank.ENSEMBL_IDs.ranked.Rdata')
attrm = read.csv("mouse_ext_IDs.csv", header = T, stringsAsFactors = F)

gosums = colSums(GO.mouse.nonIEA )
filtGO = gosums <=1000 & gosums>=20
annot =  GO.mouse.nonIEA[,filtGO]
m = match(rownames(annot), attrm$NCBI.gene.ID)
f.a = !is.na(m)          
f.t = m[f.a]
rownames(annot)[f.a] = as.character(attrm[f.t,]$Gene.stable.ID)
voc1 = voc[voc$V3 =="biological_process",]
m = match(colnames(annot), voc1$V1)
f.a = !is.na(m)
annot = annot[,f.a]

network = reranked
rm(reranked)
#Make Annotation for GO:

m = match(rownames(annot), rownames(networkm))
f.a = !is.na(m)
CLUSTM = annot[f.a,]

save(networkh, networkm, file = 'reranked.ortholog_networks.Rdata')
roc_aggh = run_GBA(networkh, CLUSTH, min = 0, max = 30000)
roc_aggm = run_GBA(networkm, CLUSTM, min = 0, max = 30000)
save(roc_aggh, roc_aggm = file = 'GO_ROC_orthologue_network.rdata')