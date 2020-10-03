##############################

##Author Name: Conor Cremin, Sara Ballouz
##Title: Coexpression Network Construction-step 3

#Requirements:

#Rdata objects from Step 2, countdata and projectsids
#This step may take a long time, best to run overnight with plenty of memory!!!!
##############################

#Packages:

library(EGAD)

# coexpression
load("genematrix.rdata") #From step 2

for( i in 1:length(projectsids)) {
  filt = filttmat[,1]  == projectsids[i]
  m  = match ( rownames(countdata), res[,1] ) 
  f.c = !is.na(m)
  f.r = m[f.c]
  
  subcount = countdata[f.c,filt]
  
  net = build_coexp_network(subcount,res[f.r,1] , flag = "rank") # change flag arguement to rank to get only (+) correlations in matrix
  
  med = median(net, na.rm=T)
  net[is.na(net)] = med
  #save(net, file=paste0(projectsids[i],".coexp.rank.Rdata") ) #This will save the coexpression network of each study (Optional)
  
  if(i==1 ) {
    agg = net
  } else {
    agg = net + agg
  }
}
save(agg, file = "coexp_agg.rdata") #IMPORTANT TO SAVE THIS!!!!!!!!!!!!!!!! This is the unranked raw network!!!!!
coexp.rank =  matrix(rank(agg, na.last = "keep", ties.method = "average"), nrow=dim(agg)[1], ncol=dim(agg)[2] )
rownames(coexp.rank) = rownames(agg)
colnames(coexp.rank) = rownames(agg)
coexp.rank = coexp.rank/max(coexp.rank, na.rm=T)
save(coexp.rank, file="coexp.agg.rank.Rdata")

# OR
# absolute ranked:

agg.absrank =  matrix(rank(agg, na.last = "keep", ties.method = "average"), nrow=dim(agg)[1], ncol=dim(agg)[2] )
rownames(agg.absrank) = rownames(agg)
colnames(agg.absrank) = rownames(agg)
agg.absrank = agg.absrank/max(agg.absrank, na.rm=T)
save(agg.absrank, file="coexp.agg.absolute_rank.rdata")
