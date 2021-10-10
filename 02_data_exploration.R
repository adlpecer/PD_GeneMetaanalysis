######################################
### 02_data_exploration.R
### Author: Adolfo López-Cerdán
### Contact: adlpecer@gmail.com
######################################

# The next plots are created from expression data:
# Barplot of phenotipic categories.
# Boxplot of expression values by subject.
# Dendogram of subjects according expression profiles
# PCA of expression matrix

# load libraries 
library(reshape2)
library(ggplot2)
library(dplyr)
library(Biobase)
library(ggplot2)
library(factoextra)
library(ggdendro)

# BOXPLOT FUNCTION
boxplot_ggplot <- function(data, groups, title = "Boxplot. Normalized expression", 
                           bottom = FALSE, save = NULL, width = 800, height = 600) {
  
  require(reshape2)
  
  box_data <- melt(data, varnames = c("gene", "sample"), as.is = TRUE)
  
  names(groups) <- colnames(data)
  box_data$Condition <- groups[(box_data$sample)]
  
  ## Boxplot
  g = ggplot(box_data, aes(x = sample, y = value, fill = Condition)) + 
    geom_boxplot(outlier.size = 0.001) + 
    xlab("Samples") + 
    ylab("Expression") + 
    ggtitle(title) + 
    labs(fill = "Groups") + 
    theme(axis.text.x = element_text(angle = 90, size = 0.5), 
          axis.text.y = element_text(size = 15), 
          title = element_text(size = 20), 
          text = element_text(size = 15))  
  
  if (bottom == TRUE) {
    # Plot legend
    g =  g + theme(legend.position ="bottom", legend.key.size = unit(0.5,"line"), legend.title=element_text(size=10))
  }
  
  if(!is.null(save)){
    png(filename = save, width = width, height = height)
    print(g)
    dev.off()
  }
  
  return (g)
  
}

# DENDOGRAM FUNCTION
treeclust_from_ggplot <- function(data, groups, title = "Clustering. Correlation distance", 
                                  subtitle = NULL, bottom = FALSE, dist = "cor", 
                                  save = NULL, width = 800, height = 600) {
  
  require(ggplot2)
  require(ggdendro)
  
  ## Create cluster
  if(dist == "cor"){
    correlacion <- cor(data)
    distancia <- as.dist((1 - correlacion) / 2)
  }else if(dist == "euclid"){
    distancia <- dist(t(data), method = "euclidean") # dist needs variables as columns
  }else{
    stop("Please specify an accepted distance. Options are 'cor' or 'euclid'")
  }
  
  cluster <- hclust(distancia)
  cluster$clase <- groups
  
  ##
  dendr <- ggdendro::dendro_data(cluster, type = "rectangle")
  ##
  clases <- as.character(cluster$clase)
  clust.df <- data.frame(label = cluster$labels, Condition = factor(clases))
  ##
  dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
  ##
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  g = ggplot() +
    geom_segment(data = ggdendro::segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_text(data = label(dendr), aes(x, y, label = label, hjust = 0, color = Condition), size = 2) +
    coord_flip() + scale_y_reverse(expand = c(0.2, 0)) +
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank(),
          # title = element_text(size = 20),
          text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    scale_colour_hue() +
    scale_fill_hue() 
  
  if(!is.null(subtitle)){
    # Añadir subtítulo
    g <- g + labs(title = title, subtitle = subtitle)
  }else{
    g <- g + ggtitle(title) + theme()
  }
  
  if (bottom == TRUE) {
    #leyenda bajo el gráfico
    g =  g + theme(legend.position ="bottom", legend.key.size = unit(0.5,"line"), legend.title=element_text(size=6))
  }
  
  if(!is.null(save)){
    png(filename = save, width = width, height = height)
    print(g)
    dev.off()
  }
  
  return (g)
}

# EXPLORATORY ANALYSIS FUNCTION
plotEDA <- function(study,l.studies){
  
  # Create plots directory
  wd <- "./Plots/"
  dir.create(wd, recursive = TRUE)
  
  # Extract expression set from list
  array <- l.studies[[study]]
  
  # BARPLOT
  dplot <- pData(array)
  agg <- count(dplot, sex, brain_region, disease, .drop = FALSE)
  agg_ord <- mutate(agg,
                    disease = reorder(disease, -n, sum),
                    sex = reorder(sex, -n, sum))
  ggplot(agg_ord, aes(x = disease, y = n, fill = sex)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(title =study) +
    geom_text(aes(label=n),
              vjust=1.6,
              position = position_dodge(0.9), size=4) +
    theme_minimal() +
    facet_wrap(~ brain_region) +
    theme(strip.text.x = element_text(size = 12),
          plot.background=element_rect(fill="white"))
  
  ggsave(paste0(wd,study,"_bar.png"),width = 12,height = 5)
  
  # Create groups from phenotypic variables
  dplot$groups <- with(dplot, paste(sex, brain_region, disease,sep = "_"))
  
  # BOXPLOT
  boxplot_ggplot(exprs(array),dplot$groups,title=paste0(study," - Boxplot"))
  ggsave(paste0(wd,study,"_box.png"),width = 12,height = 8)
  
  # DENDOGRAM
  treeclust_from_ggplot(exprs(array),dplot$groups)
  ggsave(paste0(wd,study,"_dendo.png"),width = 12,height = 10)
  
  # PCA
  res.pca <- prcomp(t(exprs(array)))
  fviz_pca_ind(res.pca,
               habillage = factor(dplot$groups),
               geom ="point",
               pointsize=4,
               invisible="quali"
  ) +
    theme(plot.background=element_rect(fill="white"))

  ggsave(paste0(wd,study,"_PCA.png"),width = 12,height = 8)
  
  
}

# START
load("./Output/arrays_st.RData")
for (s in names(arraysP)){
  plotEDA(s,arraysP)
}