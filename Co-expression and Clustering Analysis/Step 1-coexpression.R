##############################

##Author Name: Conor Cremin, Sara Ballouz
##Title: Coexpression Network Construction-step 1

#Requirements:

#Since network construction will use a large number of studies, it is necessary that study directories have the same layout.
#Directory name should be the study id e.g. GSE id etc.,within this directory alignment output from STAR should be in a directory called, "star_BAM".
#This will make it eazier to scan read counts if study directories are organized similarly. 

##############################

#Packages:

library(dplyr)

# Step 1: Set the parameters:

logname = "Log.final.out"
dir = getwd()
files = data.frame(list.dirs(full.names = T))

#Step 1: Log-output Matrix construction (to get alignnment stats - get a feel for the data)

filedir = paste(files, "star_BAM", sep="/") # Specify file path
f1 = matrix(list.files(filedir, full.names = T, pattern = logname), ncol = 1)

# for Loop formula: Getting file paths to get SRR and SRP numbers

Ns = list()
i = 1

for( n in f1[,1] ){
  N = list()
  logfile = n
  
  if( file.exists(logfile) ) {
    logstats =read.delim(logfile, header = F)
    N=logstats
  } 
  
  Ns[[i]] = N
  
  i = i + 1
  
}
data = f1
mat = matrix( unlist(strsplit( as.character(data[,1]) , "/" ) ),ncol = 8, by=T )
mat[,8] = gsub("Log.final.out", "", mat[,8]  )
sraids = mat[,c(6,8)]

filenames = paste( sraids[,1], "/star_BAM/", sraids[,2], "Log.final.out", sep="")
projectid = sraids[,1]
sampleid = sraids[,2]

#  Extracting the appropriate STATs from alignment log files.

filt = c(8, 9, 24, 29)
f2 <- sapply(1:length(f1[,1]), function(i) Ns[[i]][filt,2] )
f2 = t(f2)
f2 <- gsub("%", "", as.matrix(f2))
mat = as.data.frame(cbind(projectid, sampleid, filenames, f2))
colnames(mat) = c("project_id", "sample_id", "filenames", "Uniquely mapped reads number", "Unique_Mapping", "Multi-Mapping", "Unmapped")
mat[,4:7] <- lapply(mat[,4:7], function(x) as.numeric(as.character(x)))
filt = rowSums(mat[,5:7]) > 0
mat = mat[filt,] 

# Plot 1: Histogram of Unique Mapping Percentages (Distribution)

hist(mat[mat$project_id == "GSE152418",]$Unique_Mapping, main = "Unique Mapping Distribution between all samples ", breaks = 20)
abline(v=mean(mat$Unique_Mapping))

# Filters to retain samples. 

filt = mat[,5] <= 50 
mat = mat[filt,]

## Step 2: Building the Expression Matrix:

m = match(as.factor(sraids[,2]),mat$sample_id)
f.a = !is.na(m)
ids = sraids[f.a,]                                      # f.t output numeric positions matching to the second arguement
filenames = paste( ids[,1], "/", ids[,2], "ReadsPerGene.out.tab", sep="")
dir = getwd()
filedir = paste(dir, "results", "", sep="/") # Specify file path
m = match(mat$sample_id, ids[,2])
f.a = !is.na(m)
table(f.a)                                              # Should be all true!!

# Counting Loop formula:

files = filenames
Ns = list()
i = 1

for( n in files ){
  N = list()
  countfile = paste0(filedir, n)
  
  if( file.exists(countfile) ) {
    print(countfile)
    counts =  read.delim(countfile, header = F)
    
    N$unmapped =  counts[1,]
    N$multimapping = counts[2,]
    N$noFeature =   counts[3,]
    N$ambiguous = counts[4,]
    N$length = dim(counts)[1]-4
    N$genes = counts[ (1:N$length)+4,1]
    N$counts1 = counts[ (1:N$length)+4,2]
    N$counts2 = counts[ (1:N$length)+4,3]
    N$counts3 = counts[ (1:N$length)+4,4]
  } else {
    N$counts1 = rep(0, length(attr$ensemblID ) )
  } 
  if( i == 1  ){
    counts_exp = rep(0, length(N$counts1)) 
  } 
  if( sum(N$counts3) > sum(N$counts2)+10){
    counts_exp = cbind(counts_exp, N$counts3)  
  } 
  else if (sum(N$counts2) > sum(N$counts3)+10){
    counts_exp = cbind(counts_exp, N$counts2)  
  } 
  else if (sum(N$counts2) == sum(N$counts3)+10 | sum(N$counts2) == sum(N$counts3)-10){
    counts_exp = cbind(counts_exp, N$counts1)  
  } else {
    counts_exp = cbind(counts_exp, N$counts1)  
  }
  Ns[[i]] = N
  print(i)
  i = i + 1
}

# Step 3: Format Row Names to remove the decimal from ENSG names and Save:

counts_exp = counts_exp[, -1] ## Get rid of the all zero column (Do only Once!!!!)

colnames(counts_exp) = as.character(ids[,2]) 
genes = matrix(unlist(strsplit( as.character(N$genes), "\\.") ) , ncol=2, byro=T)[,1] # GETS RID OF ENSG Decimal
rownames(counts_exp) = genes 

# Manual Assessment of count file matrix

head(order(rank(colSums(counts_exp))),20)
sum(counts_exp[,1931])    ## Change to select col number

filt2 = colSums(counts_exp) >= 1000000
m = match(mat$sample_id, colnames(counts_exp[,filt2]))
f.a = !is.na(m)
filttmat = as.data.frame(mat[f.a,])

# remove studies with less than ten samples

filttmat = filttmat %>%
  group_by(project_id) %>%
  filter(n() >9)     #Studies with less than 10 samples are removed

countdata = counts_exp[,as.character(filttmat$sample_id)]
