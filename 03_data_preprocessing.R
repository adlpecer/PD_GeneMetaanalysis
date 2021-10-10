######################################
### 03_data_preproccesing.R
### Author: Adolfo López-Cerdán
### Contact: adlpecer@gmail.com
######################################

# Here we normalize expression data when necessary,
# annotate probes, filter outliers, adapt brain region
# codes (phenotipic data) and select only samples from
# the selected brain regions (FC, ST and SN).

# load libraries
library(Biobase)
library(dplyr)
library(oligo)

# Function to pick median values when duplicated genes
medianReps <- function(matriz){
  ID <- as.character(rownames(matriz))
  ID <- factor(ID, levels = unique(ID))
  df <- by(matriz, ID, function(x) apply(x, 2, stats::median))
  mat <- do.call("rbind", df)
  return(mat)
}

# Load the list of ExpressionSets
load("./Output/arrays_st.RData")

######### GSE28894 #########
# Preproccess expression matrix with log2 transform
# minimum value is set to 1 before transformation to get minimum equal to 0
array <- arraysP[["GSE28894"]]
exprs(array) <- exprs(array) + (abs(min(exprs(array))) + 1)
exprs(array) <- log2(exprs(array))

# Formatting of brain region column and filter selected regions
df <- pData(array)
df$brain_region <- as.character(df$brain_region)
df$brain_region[grepl("cortex",x = df$brain_region)] <- "FC"
df$brain_region[grepl("striatum",x = df$brain_region)] <- "ST"
pData(array) <- df
array <- array[,df$brain_region %in% c("FC","ST","SN")]

# Annotation of probes with gene symbol and Entrez ID
# Select median values when duplicated genes
genedata <- fData(array)
genedata <- genedata[,c("Symbol","Entrez_Gene_ID")]
genedata$Symbol <- sub(" ///.*","",genedata$Symbol)
genedata$Entrez_Gene_ID <- sub(" ///.*","",genedata$Entrez_Gene_ID)
exprs0 <- exprs(array)
rownames(exprs0) <- genedata$Symbol
exprs0 <- medianReps(exprs0)
exprs0 <- exprs0[rownames(exprs0) != "",]
genedata <- genedata[match(rownames(exprs0),genedata$Symbol),c("Symbol","Entrez_Gene_ID")]
rownames(genedata) <- genedata$Symbol

# Build the new preprocessed ExpressionSet and save it in list of arrays
proc_array <- new("ExpressionSet",
                     exprs=exprs0,
                     annotation=annotation(array))
fData(proc_array) <- genedata
pData(proc_array) <- pData(array)

arraysP[["GSE28894"]] <- proc_array

######### GSE8397 #########

array <- arraysP[["GSE8397"]]

# Data already log2 transformed

# Formatting of brain region column and filter selected regions
df <- pData(array)
df$brain_region <- as.character(df$brain_region)
df$brain_region[grepl("substantia nigra",x = df$brain_region)] <- "SN"
pData(array) <- df
array <- array[,df$brain_region %in% c("FC","ST","SN")]

# Annotation of probes with gene symbol and Entrez ID
# Select median values when duplicated genes
genedata <- fData(array)
genedata <- genedata[,c("Gene Symbol","ENTREZ_GENE_ID")]
colnames(genedata) <- c("Symbol","Entrez_Gene_ID")
genedata$Symbol <- sub(" ///.*","",genedata$Symbol)
genedata$Entrez_Gene_ID <- sub(" ///.*","",genedata$Entrez_Gene_ID)
exprs0 <- exprs(array)
rownames(exprs0) <- genedata$Symbol
exprs0 <- medianReps(exprs0)
exprs0 <- exprs0[rownames(exprs0) != "",]
genedata <- genedata[match(rownames(exprs0),genedata$Symbol),c("Symbol","Entrez_Gene_ID")]
rownames(genedata) <- genedata$Symbol

# Build the new preprocessed ExpressionSet and save it in list of arrays
proc_array <- new("ExpressionSet",
                  exprs=exprs0,
                  annotation=annotation(array))
fData(proc_array) <- genedata
pData(proc_array) <- pData(array)

arraysP[["GSE8397"]] <- proc_array

######### GSE20295 #########

array <- arraysP[["GSE20295"]]
# Some samples in this study don't seem to be normalized
# Here we filter and normalize by quantile method those samples who need it
# Then we log2 transform them
norm <-apply(exprs(array), 2, mean) > 10
expr0 <- exprs(array)
expr0[,norm] <- expr0[,norm] + (1 - min(expr0[,norm]))
expr0[,norm] <- log2(expr0[,norm])
expr0 <- normalize(expr0,method="quantile")
exprs(array) <- expr0

# Formatting of brain region column and filter selected regions
df <- pData(array)
df$brain_region <- as.character(df$brain_region)
df$brain_region[grepl("cortex",x = df$brain_region)] <- "FC"
df$brain_region[grepl("substantia nigra",x = df$brain_region)] <- "SN"
df$brain_region[grepl("putamen",x = df$brain_region)] <- "ST"
pData(array) <- df
array <- array[,df$brain_region %in% c("FC","ST","SN")]

# Annotation of probes with gene symbol and Entrez ID
# Select median values when duplicated genes
genedata <- fData(array)
genedata <- genedata[,c("Gene Symbol","ENTREZ_GENE_ID")]
colnames(genedata) <- c("Symbol","Entrez_Gene_ID")
genedata$Symbol <- sub(" ///.*","",genedata$Symbol)
genedata$Entrez_Gene_ID <- sub(" ///.*","",genedata$Entrez_Gene_ID)
exprs0 <- exprs(array)
rownames(exprs0) <- genedata$Symbol
exprs0 <- medianReps(exprs0)
exprs0 <- exprs0[rownames(exprs0) != "",]
genedata <- genedata[match(rownames(exprs0),genedata$Symbol),c("Symbol","Entrez_Gene_ID")]
rownames(genedata) <- genedata$Symbol

# Build the new preprocessed ExpressionSet and save it in list of arrays
proc_array <- new("ExpressionSet",
                  exprs=exprs0,
                  annotation=annotation(array))
fData(proc_array) <- genedata
pData(proc_array) <- pData(array)

arraysP[["GSE20295"]] <- proc_array

######### GSE7621 #########

array <- arraysP[["GSE7621"]]

# Preproccess expression matrix with log2 transform
# minimum value is set to 1 before transformation to get minimum equal to 0
exprs(array) <- exprs(array) + (1 - min(exprs(array)))
exprs(array) <- log2(exprs(array))

# Formatting of brain region column and filter selected regions
df <- pData(array)
df$brain_region <- as.character(df$brain_region)
df$brain_region[grepl("substantia nigra",x = df$brain_region)] <- "SN"
pData(array) <- df
array <- array[,df$brain_region %in% c("FC","ST","SN")]

# Annotation of probes with gene symbol and Entrez ID
# Select median values when duplicated genes
genedata <- fData(array)
genedata <- genedata[,c("Gene Symbol","ENTREZ_GENE_ID")]
colnames(genedata) <- c("Symbol","Entrez_Gene_ID")
genedata$Symbol <- sub(" ///.*","",genedata$Symbol)
genedata$Entrez_Gene_ID <- sub(" ///.*","",genedata$Entrez_Gene_ID)
exprs0 <- exprs(array)
rownames(exprs0) <- genedata$Symbol
exprs0 <- medianReps(exprs0)
exprs0 <- exprs0[rownames(exprs0) != "",]
genedata <- genedata[match(rownames(exprs0),genedata$Symbol),c("Symbol","Entrez_Gene_ID")]
rownames(genedata) <- genedata$Symbol

# Build the new preprocessed ExpressionSet and save it in list of arrays
proc_array <- new("ExpressionSet",
                  exprs=exprs0,
                  annotation=annotation(array))
fData(proc_array) <- genedata
pData(proc_array) <- pData(array)

arraysP[["GSE7621"]] <- proc_array

######### GSE20146 #########

array <- arraysP[["GSE20146"]]

# Preproccess expression matrix with log2 transform
# minimum value is set to 1 before transformation to get minimum equal to 0
exprs(array) <- exprs(array) + (1 - min(exprs(array)))
exprs(array) <- log2(exprs(array))

# Formatting of brain region column and filter selected regions
df <- pData(array)
df$brain_region <- as.character(df$brain_region)
df$brain_region[grepl("globus",x = df$brain_region)] <- "ST"
pData(array) <- df
array <- array[,df$brain_region %in% c("FC","ST","SN")]

# Annotation of probes with gene symbol and Entrez ID
# Select median values when duplicated genes
genedata <- fData(array)
genedata <- genedata[,c("Gene Symbol","ENTREZ_GENE_ID")]
colnames(genedata) <- c("Symbol","Entrez_Gene_ID")
genedata$Symbol <- sub(" ///.*","",genedata$Symbol)
genedata$Entrez_Gene_ID <- sub(" ///.*","",genedata$Entrez_Gene_ID)
exprs0 <- exprs(array)
rownames(exprs0) <- genedata$Symbol
exprs0 <- medianReps(exprs0)
exprs0 <- exprs0[rownames(exprs0) != "",]
genedata <- genedata[match(rownames(exprs0),genedata$Symbol),c("Symbol","Entrez_Gene_ID")]
rownames(genedata) <- genedata$Symbol

# Build the new preprocessed ExpressionSet and save it in list of arrays
proc_array <- new("ExpressionSet",
                  exprs=exprs0,
                  annotation=annotation(array))
fData(proc_array) <- genedata
pData(proc_array) <- pData(array)

arraysP[["GSE20146"]] <- proc_array



######### EMTAB1194 #########

array <- arraysP[["EMTAB1194"]]

# Data already log2 transformed

# Formatting of brain region column and filter selected regions
df <- pData(array)
df$brain_region <- as.character(df$brain_region)
df$brain_region[grepl("frontal_lobe",x = df$brain_region)] <- "FC"
pData(array) <- df
array <- array[,df$brain_region %in% c("FC","ST","SN")]

# Filtering outlier 21_G05_P2.CEL
array <- array[,rownames(df) != "21_G05_P2.CEL"]

# Annotation of probes with gene symbol and Entrez ID
# Select median values when duplicated genes
genedata <- fData(array)
colnames(genedata) <- c("PROBEID","Entrez_Gene_ID","ENSEMBL","Symbol")
exprs0 <- exprs(array)
rownames(exprs0) <- genedata$Symbol
exprs0 <- medianReps(exprs0)
genedata <- genedata[match(rownames(exprs0),genedata$Symbol),c("Symbol","Entrez_Gene_ID")]
rownames(genedata) <- genedata$Symbol

# Build the new preprocessed ExpressionSet and save it in list of arrays
proc_array <- new("ExpressionSet",
                  exprs=exprs0,
                  annotation=annotation(array))
fData(proc_array) <- genedata
pData(proc_array) <- pData(array)

arraysP[["EMTAB1194"]] <- proc_array


######### EMEXP1416 #########

array <- arraysP[["EMEXP1416"]]

# Data already log2 transformed

# Formatting of brain region column and filter selected regions
df <- pData(array)
df$brain_region <- as.character(df$brain_region)
df$brain_region[grepl("substantia nigra",x = df$brain_region)] <- "SN"
pData(array) <- df
array <- array[,df$brain_region %in% c("FC","ST","SN")]

# Annotation of probes with gene symbol and Entrez ID
# Select median values when duplicated genes
genedata <- fData(array)
colnames(genedata) <- c("PROBEID","Entrez_Gene_ID","ENSEMBL","Symbol")
exprs0 <- exprs(array)
rownames(exprs0) <- genedata$Symbol
exprs0 <- medianReps(exprs0)
genedata <- genedata[match(rownames(exprs0),genedata$Symbol),c("Symbol","Entrez_Gene_ID")]
rownames(genedata) <- genedata$Symbol

# Build the new preprocessed ExpressionSet and save it in list of arrays
proc_array <- new("ExpressionSet",
                  exprs=exprs0,
                  annotation=annotation(array))
fData(proc_array) <- genedata
pData(proc_array) <- pData(array)

arraysP[["EMEXP1416"]] <- proc_array
######### Save arrays #########
save(arraysP, file = "./Output/arrays_proc.RData")
