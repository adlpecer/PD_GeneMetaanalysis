######################################
### 05_meta-analysis.R
### Author: Francisco García-García
### Modified by: Adolfo López-Cerdán
### Contact: fgarcia@cipf.es
###          adlpecer@gmail.com
######################################


## load libraries 
library(Biobase)
library(metafor)
library(stringr)
library(limma)

## load data
load("./Output/ALLfit.RData")

## directories
wd <- "./Results/MA/" #Output dir
metaanalysis_name <- "MA_20211008_ALL"
wd <- paste0(wd,metaanalysis_name,"/")
dir.create(wd, recursive = TRUE)
#setwd(wd)

# STEP 0. Pre-processing previous data
# ===============================================================


## Calculate SE
SE_array <- function(fit) {
  #OPTION1: https://support.bioconductor.org/p/70175/ by Gordon Smyth/January Weiner
  #The effect sizes are contained in fit$coefficients
  summary(fit$coefficients)
  head(fit$coefficients)
  #The standard errors can be obtained from 2 sub-options:
  # SE.coef <- sqrt(fit$s2.post) #JANUARY  (Here I have a problem when having several contrasts
  #                                        #because I have the same information for all contrasts)
  # head(SE.coef)
  # summary(SE.coef)
  SE.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled[,3] #GORDON
  head(SE.coef)
  summary(SE.coef)
  
  #OPTION2: https://support.bioconductor.org/p/70175/ by  Steve Lianoglou  (SE PARECE A GORDON)
  allgenes <- topTable(fit, number = "all", confint=TRUE, adjust.method = "fdr",coef =3)
  dim(allgenes)
  allgenes[, "SE"] <- (allgenes[, "CI.R"] - allgenes[, "CI.L"])/ 3.92
  head(allgenes)
  
  #final results
  table(rownames(SE.coef) == rownames(fit$coefficients))
  mat <- cbind(fit$coefficients[,3], SE.coef)
  colnames(mat) <- c("coef", "se.coef")
  head(mat)
  
  int <- intersect(rownames(allgenes), rownames(mat))
  length(int)
  res <- cbind(allgenes, mat[rownames(allgenes),])
  head(res)
  dim(res)
  return(res)
}

# Use the loaded list as function input
EDs_sel <- lapply(ALLfit, SE_array)


# STEP 1. Preparing input for meta-analysis: LOR and SE matrix
# ===============================================================

# we search a list including all unique ID genes for all studies
genes <- NULL
for (fi in EDs_sel){
  genes <- c(genes, rownames(fi))
}

length (genes)
genes <- unique (genes)
length (genes)
genes[grepl("[0-9]+-[A-Z][a-z]{2}",genes)]
genes <- genes[!grepl("[0-9]+-[A-Z][a-z]{2}",genes)]
genes[grepl("[0-9]+-[A-Z][a-z]{2}",genes)]
genes <- base::sort(genes)


### generating matrix with all logFC for all studies
mat.logFC <- matrix (NA, nrow = length (genes), ncol = length(EDs_sel))
rownames (mat.logFC) <- genes
colnames (mat.logFC) <- names(EDs_sel)
head (mat.logFC)

for (i in 1:length(EDs_sel)){
  co <- names(EDs_sel[i])
  res <- EDs_sel[[i]]
  logFC <- res$logFC
  names (logFC) <- (rownames (res))
  mat.logFC[, co] <- logFC[rownames(mat.logFC)] 
}
##### prueba
#mat.logFC <- data.frame(Symbol=genes,row.names =genes)
#for(co in names(EDs_sel1)){
#  res <- EDs_sel1[[co]]
#  res <- res[,c("Symbol","logFC")]
#  colnames(res) <- c("Symbol",co)
#  mat.logFC <- merge(mat.logFC,res,by = "Symbol",all.x = TRUE)
#}
#rownames(mat.logFC) <- genes
#mat.logFC <- mat.logFC[,-1]
head (mat.logFC)
tail(mat.logFC)
table (is.na(mat.logFC))
dim (mat.logFC)

# select genes included at least in 2 or more studies
mat.logFC.NA <- is.na(mat.logFC)
head(mat.logFC.NA)
sum.NA <-  apply(mat.logFC.NA, 1, sum)
table(sum.NA)
min.sum.NA <- sum.NA < ncol(mat.logFC) - 1
table(min.sum.NA)

# filter by min.sum.NA
mat.logFC <- mat.logFC[min.sum.NA == T, ]
dim(mat.logFC)


### generating matrix with all SE for all studies  
mat.SE <- matrix (NA, nrow = length (genes), ncol = length(EDs_sel))
rownames (mat.SE) <- genes
colnames (mat.SE) <- names(EDs_sel)
head (mat.SE)


# (SE FROM GORDON: se.coef)
for (i in 1:length(EDs_sel)){
  co <- gsub("_ED", "", names(EDs_sel[i]))
  res <- EDs_sel[[i]]
  SE <- res$se.coef
  names (SE) <- (rownames (res))
  mat.SE[, co] <- SE[rownames(mat.SE)] 
}


head (mat.SE)
tail(mat.SE)
table (is.na(mat.SE))
dim (mat.SE)

# filter by min.sum.NA
mat.SE <- mat.SE[min.sum.NA == T, ]
dim(mat.SE)




# STEP 2. Meta-analysis for genes
# ===============================================================

# suppose between-study variance is non-zero.
# there are different methods to estimate this variance:
# DL (Dersimonian-Laird), REML (Restricted maximum-likelihood, default)....
# Now we have logFC and SE  (not VARIANCE), so:
# yi -> logFC   sei -> SE
# result.lor <- rma(yi = mat.logFC[1, ], 
#                   sei = mat.SE[1, ],   #pay attention, not vi (varianze)   
#                   method = "DL") # DerSimonian-Laird.


# explore the function to do the meta-analysis
#?rma

MA <- lapply(1:length(rownames(mat.logFC)),
             function(x){rma(yi = mat.logFC[x, ],
                             sei = mat.SE[x, ],
                             method = "DL")})

# MA <- lapply(1:length(rownames(mat.logFC)),
#              function(x){rma(yi = mat.logFC[x, ],
#                              sei = mat.SE[x, ],
#                              method = "FE")})

names (MA) <- rownames(mat.logFC)
class (MA)
length(MA)
head (MA)
MA[[1]]

#result.logFC$pval      #p-value about logFC = 0
#result.logFC$ci.lb     #IC down
#result.logFC$ci.ub     #IC up
#re sult.logFC$b         #estimation of combined logFC

#data.frame including all detailed results:
result_meta <- as.data.frame(do.call("rbind",
                                     lapply(MA,
                                            function(x){
                                              c(x$ci.lb, x$b, x$ci.ub, 
                                                x$pval, x$QE, x$QEp, x$se,
                                                x$tau2, x$I2, x$H2)
                                            })))

colnames(result_meta) <- c("lower_bound", "logFC", "upper_bound",
                           "pvalue", "QE", "QEp", "SE", "tau2", "I2", "H2")

p.adjust.fdr <- stats::p.adjust(result_meta[,4], method = "fdr")
p.adjust.BY  <- stats::p.adjust(result_meta[,4], method = "BY")
result_meta <- round(cbind(result_meta, p.adjust.fdr, p.adjust.BY), 3)
head(result_meta)


# significant genes
corte = 0.05
table(result_meta[, "pvalue"] < corte)
table(result_meta[, "p.adjust.fdr"] < corte)
table(result_meta[, "p.adjust.BY"] < corte)


# add number of studies where the gene is evaluated
sum.NA <- sum.NA[sum.NA<ncol(mat.logFC)-1]
n.studies <-  ncol(mat.logFC) - sum.NA 
table(n.studies)
n.studies <- n.studies [rownames(mat.logFC)]
length(n.studies)
result_meta[, "n.studies"]  <- n.studies
head(result_meta)
summary(result_meta$p.adjust.fdr)

# no tenemos genes significativos en estos estudios y he subido el FDR para seleccioar
# al menos un grupo con los que ver el funcionamiento del script

#corte = 0.2
sig.genes.df = result_meta[result_meta$p.adjust.fdr < corte,] 
dim(sig.genes.df)

write.table(x = sig.genes.df[order(sig.genes.df$p.adjust.fdr),], file = paste0(wd,"sig.genes.tsv"), sep = "\t", quote = FALSE)
write.table(x = result_meta[order(result_meta$p.adjust.fdr),], file = paste0(wd,"all.genes.tsv"), sep = "\t", quote = FALSE)


# STEP 3. INFLUENCE AND SENSITIVITY ANALYSIS
# ===============================================================

#add 4 new variables about influence & sensitivity analysis:  

for (i in rownames(sig.genes.df)){
  print(i)
  #define studies for each function (not NA)
  estudios <- colnames(mat.logFC)[!mat.logFC.NA[i,]]
  
  #influence info 1: 
  #number of studies where the sign of the logOR is the same  of the global logOR:
  sig.genes.df[i, "infl.same.sign.logFC"] <- sum(sign(MA[[i]]$yi)== rep(sign(MA[[i]]$b),length(estudios)))
  
  #influence info 2: how many studies could be influencers?
  inf <- influence(MA[[i]])
  res <- paste(estudios[inf$is.infl], collapse = ",")  
  sig.genes.df[i, "infl.nstudies"] <- ifelse(res =="", "non", res)
  
  #sensivity analysis
  l1 <-as.data.frame(leave1out(MA[[i]]))
  rownames(l1) <- estudios
  
  #1. p.value about differences between all estimates from leave one out
  #   and global estimate)
  sig.genes.df[i, "sensi.global"] <-t.test(x= l1$estimate,
                                           mu=as.numeric(MA[[i]]$b))$p.value
  #2. number of  studies where pvalue > 0.05 
  # (we hope p-values < 0.05, significant estimates) 
  res2 <- paste(estudios[l1$pval > 0.05], collapse = ",")  
  sig.genes.df[i, "sensi.specific"] <- ifelse(res2 =="", "all.p.values < 0.05", res2)
}


## QUESTIONS TO ASSESS META-ANALYSIS FOR EACH FUNCTION:

#1. INFLUENCE STUDIES. How many logOR have the same sign to global logOR?
table(sig.genes.df$infl.same.sign.logFC)

#2. INFLUENCE STUDIES. How many functions including influence studies?
table(sig.genes.df$infl.nstudies=="non")

#3. SENSITIVITY. In global, are there many functions with differences in the estimate?
table(sig.genes.df$sensi.global < 0.05)

#4. SENSITIVITY.  How many functions including changes in the significance about 
# its new estimate  after leave1out? 
table(sig.genes.df$sensi.specific == "all.p.values < 0.05")


#save final results:
cat ("ID\t", file = paste0(wd,"sig.genes.df.txt"))
write.table(sig.genes.df, file = paste0(wd,"sig.genes.df.txt"), sep ="\t", quote = F, 
            append = TRUE, row.names = T)




# STEP 4. Visualization of significant genes
# ===============================================================

#select significant functions to visualize:
sig.results <- result_meta[result_meta[, "p.adjust.fdr"] < 0.05,]
#sig.results <- result_meta[result_meta[, "p.adjust.fdr"] < 0.2,]
sig.results
dim(sig.results)

dir.create(paste0(wd,"plots"), recursive = TRUE)
#setwd(paste0(wd,"plots"))

selMethod <- "DL"

for (i in 1:nrow(sig.results)){
  mygenes <- rownames(sig.results)[i]
  res <- rma(yi= mat.logFC[mygenes,], sei =mat.SE[mygenes,], method = "DL")
  
  #FOREST PLOT
  # png (filename = paste("FOREST_", mygenes,".png", sep =""), width = 960 , 
  #      height = 960, res = 200) 
  png (filename = paste0(wd,"plots/",gsub("-","_",mygenes),"_FOREST",".png"), width = 960 , 
       height = 960, res = 200) 
  forest(res, 
         slab = toupper(colnames(mat.logFC)), #Nombre de los estudios
         xlab="logFC", cex=0.7,
         mlab=paste(selMethod, "Model for All Studies", sep = " "), 
         border = "black", #Color del borde del rombo
         col = "red", #Color del rombo
         main = paste("\n", mygenes, sep=""))    
  text( 9,-3, "logFC [IC 95%]", pos=2, cex = 0.7)
  dev.off()
  
  #FUNNEL PLOT
  png (filename = paste0(wd,"plots/",gsub("-","_",mygenes),"_FUNNEL", ".png"), width = 960 , 
       height = 960, res = 200) 
  par(mfrow=c(2,2))
  funnel(res, main="Standard Error", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
         xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="seinv", main="Inverse Standard Error",
         back ="darkslategray1", xlab = paste("logFC (", mygenes, ")",sep =""))
  funnel(res, yaxis="vinv", main="Inverse Sampling Variance", 
         back ="darkslategray1",  xlab = paste("logFC (", mygenes, ")",sep =""))
  par(mfrow=c(1,1))
  dev.off()
  
  #INFLUENCE PLOTS 
  # That shows various diagnostic measures
  png (filename = paste0(wd,"plots/",gsub("-","_",mygenes), "_INFLUENCE", ".png"), width = 960 , 
       height = 960, res = 200) ##CAMBIAR
  inf <- influence(res)
  #plot(inf, plotfb = T)#"plotfb" is not a graphical parameter
  plot(inf)
  dev.off()
  
}

# STEP 5. Generating report
# ===============================================================
sig.genes.df <- sig.genes.df[order(sig.genes.df$p.adjust.fdr),]
save(sig.genes.df, result_meta, file = paste0(wd,metaanalysis_name, ".RData"))

# Function to create multiple tabs
make.tabs <- function(sig.genes.df){
  res <- NULL
  for(g in rownames(sig.genes.df)){
    file_name <- gsub("-","_", g)
    res <- c(res, '### ', g, '{-} \n',
             "**Statistics of ", g, " meta-analisys** \n",
             "```{r, fig.align='center'}", '\n',
             "kable(sig.genes.df['",g,"',1:11])", '\n',
             '```', '\n',
             "[Gene information](https://www.genecards.org/cgi-bin/carddisp.pl?gene=", g, ") \n\n",
             "**Forest plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', 'plots/', file_name, '_FOREST.png")\n',
             '```', '\n',
             "**Funnel plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', 'plots/', file_name, '_FUNNEL.png")\n',
             '```', '\n',
             "**Incluence plot** \n",
             "```{r, fig.align='center'}", '\n',
             'knitr::include_graphics("', 'plots/', file_name, '_INFLUENCE.png")\n',
             '```', '\n\n')
  }
  return(res)
}


# Create the Rmd to knit
ns <- nrow(sig.genes.df)
cat(
  '---
title: "Meta-analysis of genes [DRAFT]"
output:
  html_document:
    toc: false
    toc_float: false
    code_folding: hide
    number_sections: true
    theme: spacelab
---
## ', metaanalysis_name, ' {.tabset .tabset-pills -}
  
```{r, warning=F, message=F}
library(dplyr)
library(knitr)
load("', metaanalysis_name, '.RData")
```  \n
### Significant results {-}  \n',
  "```{r, fig.align='center'}", '\n',
  "kable(sig.genes.df[,1:11], caption='Statistics of ", ns, " significant genes')", '\n',
  '```', '\n\n',
  make.tabs(sig.genes.df), "\n\n",
  '### sessionInfo {-}  \n',
  "```{r, fig.align='center'}", '\n',
  "date()", "\n",
  "sessionInfo()", '\n',
  '```', '\n\n',
  sep = "",
  file = paste0(wd,metaanalysis_name, ".Rmd"))

# Render the Rmd created into html here
rmarkdown::render(paste0(wd,metaanalysis_name, ".Rmd"))
