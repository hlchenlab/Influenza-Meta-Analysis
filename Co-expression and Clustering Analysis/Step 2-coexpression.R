##############################

##Author Name: Conor Cremin, Sara Ballouz
##Title: Coexpression Network Construction-step 2

#Requirements:

#Rdata objects from Step 1, countdata and filttmat

##############################

projectsids = as.character(as.matrix(unique(filttmat[,1])))
genefilts = list() 

# low gene expression count filters

i=1
for( i in 1:length(projectsids)) {
  filt = filttmat[,1] == projectsids[i]
  subcount = countdata[,filt]
  genefilts[[i]] =  rowSums(subcount >= 10 )  >= 10 
}

genefilts.mat = do.call(cbind, genefilts)
finalgenefilt = rowSums(genefilts.mat) >= 5 
table(finalgenefilt)

# Or sapply version 

#finalgenefilt = rowSums( sapply( 1:length(projectsids), function(i) rowSums(countdata[,filttmat[,1]==projectsids[i]] >= 10 )  >= 10 ) ) >= 10 

finalgenefilt = as.matrix(finalgenefilt)
colnames(finalgenefilt) = "outcome"
genes = names(finalgenefilt[grep("TRUE", finalgenefilt),])
m = match(allgenes$Gene.stable.ID, genes)  #Can use Human_Annotation_v97.csv in Differential Expression folder
f.a = !is.na(m)
res = allgenes[f.a,]
table(droplevels(res$Gene.type))

m = match(rownames(countdata), as.character(res$Gene.stable.ID))
f.a = !is.na(m)
countdata = countdata[f.a,]
save(countdata, res, projectsids,filttmat, file="genematrix.rdata") #For step 3
