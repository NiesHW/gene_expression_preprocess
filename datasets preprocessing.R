# preprocess for raw data to read CEL files
# to get entrez gene ID to match with pathways and directGraph
# http://www.genomicscode.org/2013/02/affymetrix-microarray-analysis-using.html

# Affymetrix HGU133_PLUS2.0
source("http://bioconductor.org/biocLite.R")
# Install packages
biocLite("affy")
biocLite("annotate")
biocLite("hgu133plus2.db")
# load packages
library("affy")
library("annotate")
library("hgu133plus2.db")
# Read CEL files
affydata = ReadAffy()
# normalizing data
# rma() - 'Robust Multi-Chip' average
# use rma to background correct and normalize probe levels
rma = rma(affydata)
# expression set
ed = exprs(rma)
write.csv(ed,"probe.csv")
# hgu133plus2.db
Entrez_IDs = unlist(mget(rownames(ed), hgu133plus2ENTREZID, ifnotfound=NA))
write.csv(Entrez_IDs,"probe_ID.csv")
# combine matrix with gene Entrez ID
mRNA_matrixTrain = cbind(Entrez_IDs,ed)
row.names(mRNA_matrixTrain) = NULL
row.names(mRNA_matrixTrain) = mRNA_matrixTrain[,1]
# remove probe ID
mRNA_matrixTrain = mRNA_matrixTrain[,2:ncol(mRNA_matrixTrain)]
mode(mRNA_matrixTrain) <- "numeric"
write.csv(mRNA_matrixTrain,"RAW.csv")

# Affymetrix HGU133_A
source("http://bioconductor.org/biocLite.R")
# Install packages
biocLite("affy")
biocLite("annotate")
biocLite("hgu133a.db")
# load packages
library("affy")
library("annotate")
library("hgu133a.db")
# Read CEL files
affydata = ReadAffy()
# normalizing data
# rma() - 'Robust Multi-Chip' average
# use rma to background correct and normalize probe levels
rma = rma(affydata)
# expression set
ed = exprs(rma)
write.csv(ed,"probe.csv")
# hgu133a.db
Entrez_IDs = unlist(mget(rownames(ed), hgu133aENTREZID, ifnotfound=NA))
write.csv(Entrez_IDs,"probe_ID.csv")
# combine matrix with gene Entrez ID
mRNA_matrixTrain = cbind(Entrez_IDs,ed)
row.names(mRNA_matrixTrain) = NULL
row.names(mRNA_matrixTrain) = mRNA_matrixTrain[,1]
# remove probe ID
mRNA_matrixTrain = mRNA_matrixTrain[,2:ncol(mRNA_matrixTrain)]
mode(mRNA_matrixTrain) <- "numeric"
write.csv(mRNA_matrixTrain,"RAW.csv")

# to solve repeated gene IDs
name=unique(rownames(mRNA_matrixTrain))            # left non-repeated gene IDs
d=matrix(NA,length(name),ncol(mRNA_matrixTrain))   # create new matrix
colnames(d)=colnames(mRNA_matrixTrain)
rownames(d)=name
for(i in 1:length(rownames(d))){
    tmp=rownames(d)[i]
    id=which(rownames(mRNA_matrixTrain)==tmp)
    if(length(id)>1){
        tmp_matrix=matrix(NA,length(id),ncol(mRNA_matrixTrain))
        for(j in 1:length(id)){
            idx=id[j]
            tmp_matrix[j,]=mRNA_matrixTrain[idx,]
        }
        d[i,]=colMeans(tmp_matrix)   # average gene expression values
        
    }
    if(length(id)==1){
        d[i,]=mRNA_matrixTrain[id,]
    }
}
# print(d)  # new matrix formed
write.csv(d,"no_repeated.csv")
# remove whole row with missing gene ID
r=which(rownames(d)%in%NA)
newd=(d[-r,])
write.csv(newd,"no_missing.csv")

###########################
# sample
###########################
# c
#    m n  p
# a  1 4 17 # repeat
# a  2 5 28 # repeat
# b  3 6 59
#    4 5  6 # missing
###########################
# d
#      m   n    p
# a  1.5 4.5 22.5
# b  3.0 6.0 59.0
#    4.0 5.0  6.0 # missing
###########################
# newd
#     m   n    p
# a 1.5 4.5 22.5
# b 3.0 6.0 59.0
###########################