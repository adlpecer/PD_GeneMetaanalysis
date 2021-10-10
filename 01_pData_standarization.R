######################################
### 01_pData_standarization.R
### Author: Adolfo López-Cerdán
### Contact: adlpecer@gmail.com
######################################

# Standarization of all ExpressionSets to make easier further analyses.
# We save three columns from phenotipic data:
## sex("male" or "female")
## brain_region(Options will be set to "FC","ST" and "SN" in 03 script)
# Current categories are kept for a better understanding of data with
# data exploratory analysis.
## disease("case" or "control")

# load libraries 
library(Biobase)

# Load retrieved arrays
load("./Data/arrays.RData")



######### GSE28894 #########

array <- arraysP[["GSE28894"]]
df <- pData(array)
dplot = data.frame(sex=tolower(df$`gender:ch1`),
                   brain_region=tolower(df$`brain region:ch1`),
                   group=tolower(df$`disease state:ch1`))
rownames(dplot) <- rownames(df)
dplot$group <- as.character(dplot$group)
dplot$group[grepl("parkinson",x = dplot$group)] <- "case"
colnames(dplot)[3] <- "disease"
dplot <- dplot[,c("sex","brain_region","disease")]
pData(array) <- dplot
arraysP[["GSE28894"]] <- array


######### GSE8397 #########

array <- arraysP[["GSE8397"]]
df <- pData(array)
dplot = data.frame(sex=tolower(df$`age:ch1`),
                   brain_region=tolower(df$source_name_ch1),
                   group=tolower(df$title))
rownames(dplot) <- rownames(df)
dplot$brain_region <- gsub(dplot$brain_region,pattern ="post mortem ",replacement ="")

dplot$sex <- as.character(dplot$sex)
dplot$sex[grepl(dplot$sex,pattern ="gender: m")] <- "male"
dplot$sex[grepl(dplot$sex,pattern ="gender: f")] <- "female"

dplot$group <- as.character(dplot$group)
dplot$group[grepl(dplot$group,pattern ="control")] <- "control"
dplot$group[grepl(dplot$group,pattern ="parkinson")] <- "case"

colnames(dplot)[3] <- "disease"

dplot <- dplot[,c("sex","brain_region","disease")]
pData(array) <- dplot
arraysP[["GSE8397"]] <- array

######### GSE20295 #########

array <- arraysP[["GSE20295"]]
df <- pData(array)
dplot = data.frame(sex=tolower(df$`gender:ch1`),
                   brain_region=tolower(df$`brain region:ch1`),
                   group=tolower(df$`disease state:ch1`))
rownames(dplot) <- rownames(df)

dplot$group <- as.character(dplot$group)
dplot$group[grepl(dplot$group,pattern ="parkinson")] <- "case"

dplot$brain_region <- gsub(dplot$brain_region,pattern =" from postmortem brain",replacement ="")
colnames(dplot)[3] <- "disease"

dplot <- dplot[,c("sex","brain_region","disease")]
pData(array) <- dplot
arraysP[["GSE20295"]] <- array

######### GSE7621 #########

array <- arraysP[["GSE7621"]]
df <- pData(array)
dplot = data.frame(sex=tolower(df$characteristics_ch1.1),
                   brain_region="substantia nigra",
                   group=tolower(df$characteristics_ch1))
rownames(dplot) <- rownames(df)

dplot$group <- as.character(dplot$group)
dplot$group[grepl(pattern = "parkinson",x=dplot$group)] <- "case"
dplot$group[grepl(pattern = "control",x=dplot$group)] <- "control"

colnames(dplot)[3] <- "disease"
dplot <- dplot[,c("sex","brain_region","disease")]
pData(array) <- dplot
array <- array[complete.cases(exprs(array)),]
arraysP[["GSE7621"]] <- array

######### GSE20146 #########

array <- arraysP[["GSE20146"]]
df <- pData(array)
dplot = data.frame(sex=tolower(df$`gender:ch1`),
                   brain_region=tolower(df$`tissue:ch1`),
                   group=tolower(df$`disease state:ch1`))
rownames(dplot) <- rownames(df)

dplot$group <- as.character(dplot$group)
dplot$group[grepl(pattern = "parkinson",x=dplot$group)] <- "case"

colnames(dplot)[3] <- "disease"
dplot <- dplot[,c("sex","brain_region","disease")]
pData(array) <- dplot
arraysP[["GSE20146"]] <- array

######### EMTAB1194 #########

array <- arraysP[["EMTAB1194"]]
df <- pData(array)
dplot = data.frame(sex=tolower(df$Characteristics.Sex.),
                   brain_region=tolower(df$Characteristics.OrganismPart.),
                   group=tolower(df$Characteristics.DiseaseState.))
rownames(dplot) <- rownames(df)

dplot$group <- as.character(dplot$group)
dplot$group[grepl("parkinson",x = dplot$group)] <- "case"
dplot$group[grepl("normal",x = dplot$group)] <- "control"
dplot$brain_region <- "frontal_lobe"
colnames(dplot)[3] <- "disease"
pData(array) <- dplot
arraysP[["EMTAB1194"]] <- array

######### EMEXP1416 #########

array <- arraysP[["EMEXP1416"]]
df <- pData(array)
df <- Biobase::pData(array)
dplot = data.frame(sex=tolower(df$Characteristics..sex.),
                   brain_region=tolower(df$Characteristics.cell.type.),
                   group=tolower(df$Characteristics..disease.))
rownames(dplot) <- rownames(df)

dplot$group <- as.character(dplot$group)
dplot$group[grepl("parkinson",x = dplot$group)] <- "case"
dplot$group[grepl("normal",x = dplot$group)] <- "control"
dplot$brain_region <- "neurons from substantia nigra"
colnames(dplot)[3] <- "disease"
pData(array) <- dplot
arraysP[["EMEXP1416"]] <- array

######### Save arrays #########
dir.create("./Output", recursive = TRUE)
save(arraysP, file = "./Output/arrays_st.RData")
