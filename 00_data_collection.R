######################################
### 00_data_collection.R
### Author: Adolfo López-Cerdán
### Contact: adlpecer@gmail.com
######################################

# Loading required libraries
library(ArrayExpress)
library(GEOquery)


######### Downloading GEO studies #########


# Array containing the IDs of the included studies from GEO
gse_list <- unique(c("GSE28894","GSE8397", "GSE20295",
                     "GSE7621", "GSE20146"))

GSE <- lapply(c(gse_list), function(x) getGEO(x)[[1]])
names(GSE) <- gse_list


######### Downloading ArrayExpress studies #########

#Some issues were detected with ArrayExpress library.
#Therefore, data from AE was manually downloaded.
#Here we load raw cel files and process them.

# E-MTAB-1194

# oligo package is preferred to process data from 
# Affymetrix GeneChip Human Gene 1.1 ST Array
library(oligo) 
library(hugene11sttranscriptcluster.db)
path = "./Data/E-MTAB-1194.raw.1"
df <- read.delim("./Data/E-MTAB-1194.raw.1/E-MTAB-1194.sdrf.txt")
rownames(df) <- df$Array.Data.File
fnames <- list.files(path, full=T,pattern = ".CEL")
ABatch <- oligo::read.celfiles(filenames=fnames)
array <- rma(ABatch)
ids <- rownames(array)
anot <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
                              keys=ids,
                              columns=c("ENTREZID","ENSEMBL","SYMBOL"),
                              keytype="PROBEID",
                              multiVals="first")
anot <- anot[!duplicated(anot$PROBEID),]
fData(array) <- anot
pData(array) <- df

dims = dim(Biobase::exprs(array))
sprintf("Probes: %s, Subjects: %s",dims[1],dims[2])
sprintf("Platform: %s",Biobase::annotation(array))
arraysP[["EMTAB1194"]] <- array

# E-MEXP-1416

# AFFY package is preferred to process data from 
# Affymetrix GeneChip Human X3P Array
library(affy)
library(u133x3p.db)
path = "./Data/E-MEXP-1416.raw.1"
df <- read.delim("./Data/E-MEXP-1416.raw.1/E-MEXP-1416.sdrf.txt",stringsAsFactors = F)
rownames(df) <- df$Array.Data.File
fnames <- list.celfiles(path=path)
ABatch <- ReadAffy(filenames=fnames,celfile.path =path,phenoData =df)
array <- rma(ABatch)
ids <- rownames(array)
anot <- AnnotationDbi::select(u133x3p.db,
                              keys=ids,
                              columns=c("ENTREZID","ENSEMBL","SYMBOL"),
                              keytype="PROBEID",
                              multiVals="first")
anot <- anot[!duplicated(anot$PROBEID),]
fData(array) <- anot

dims = dim(Biobase::exprs(array))
sprintf("Probes: %s, Subjects: %s",dims[1],dims[2])
sprintf("Platform: %s",Biobase::annotation(array))
arraysP[["EMEXP1416"]] <- array


######### Save lists of studies #########
dir.create("./Data", recursive = TRUE)
save(arraysP, file = "./Data/arrays.RData")
