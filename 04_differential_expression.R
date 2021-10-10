######################################
### 04_differential_expression.R
### Author: Adolfo López-Cerdán
### Contact: adlpecer@gmail.com
######################################

# Differential expression analysis (DEA)
# with the following contrasts:
# contrast 1 -> femalecase - femalecontrol,
# contrast 2 -> malecase - malecontrol,
# contrast 3 -> (femalecase - femalecontrol)-(malecase - malecontrol)
# DEA applied to ALL tissues, FC, ST and SN.

# load libraries
library(Biobase)
library(limma)
library(dplyr)

# Differential expression function
diffexprs <- function(X,batch=F){
  sexdisease=vector("list",ncol(X))
  for(i in seq_along(sexdisease)) {
    sexdisease[[i]]=paste(pData(X)[,"sex"][i],pData(X)[,"disease"][i],sep="")
  }
  
  sexdisease=factor(unlist(sexdisease))
  
  # Consider batch effect or not(default)
  if (batch==T){
    batch=factor(pData(X)[,"batch"])
    design <- model.matrix(~ 0 + sexdisease + batch)
  } else {
    design <- model.matrix(~ 0 + sexdisease)
  }
  
  colnames(design)[1:4] <- levels(sexdisease)
  fit <- lmFit(X,design)
  cont.matrix <- makeContrasts(dif1=femalecase - femalecontrol,
                               dif2=malecase - malecontrol,
                               dif12=(femalecase - femalecontrol)-(malecase - malecontrol),
                               levels = design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  return(fit2)
}

load("./Output/arrays_proc.RData")
wd <- "./Results/DEA/"
dir.create(wd, recursive = TRUE)

# Empty lists of results
ALLfit <- list()
FCfit<- list()
STfit <- list()
SNfit <- list()
for (study in names(arraysP)){
  array <- arraysP[[study]]
  # All tissues DEA
  fit <- diffexprs(array)
  ALLfit[[study]] <- fit
  # Write significant
  for (i in c(1:3)){
    cont <- topTable(fit,coef = i,adjust="BH",sort="none",n=Inf)
    cont <- cont[cont$adj.P.Val<0.05,]
    write.csv(cont,paste0(wd,study,"_ALL_cont",i,".csv"))
  }
  # FC differential expression
  if("FC" %in% pData(array)$brain_region) {
    filtered <- array[,pData(array)$brain_region == "FC"]
    fit <- diffexprs(filtered)
    FCfit[[study]] <- fit
    for (i in c(1:3)){
      cont <- topTable(fit,coef = i,adjust="BH",sort="none",n=Inf)
      cont <- cont[cont$adj.P.Val<0.05,]
      write.csv(cont,paste0(wd,study,"_FC_cont",i,".csv"))
    }
  }
  # ST differential expression
  if("ST" %in% pData(array)$brain_region) {
    filtered <- array[,pData(array)$brain_region == "ST"]
    fit <- diffexprs(filtered)
    STfit[[study]] <- fit
    for (i in c(1:3)){
      cont <- topTable(fit,coef = i,adjust="BH",sort="none",n=Inf)
      cont <- cont[cont$adj.P.Val<0.05,]
      write.csv(cont,paste0(wd,study,"_ST_cont",i,".csv"))
    }
  }
  # SN differential expression
  if("SN" %in% pData(array)$brain_region) {
    filtered <- array[,pData(array)$brain_region == "SN"]
    fit <- diffexprs(filtered)
    SNfit[[study]] <- fit
    for (i in c(1:3)){
      cont <- topTable(fit,coef = i,adjust="BH",sort="none",n=Inf)
      cont <- cont[cont$adj.P.Val<0.05,]
      write.csv(cont,paste0(wd,study,"_SN_cont",i,".csv"))
    }
  }
}
# Save lists of "fit" objects
save(ALLfit, file = "./Output/ALLfit.RData")
save(FCfit, file = "./Output/FCfit.RData")
save(STfit, file = "./Output/STfit.RData")
save(SNfit, file = "./Output/SNfit.RData")