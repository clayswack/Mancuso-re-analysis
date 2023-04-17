
#Load the packages
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("devtools")
library(devtools)
#devtools::install_github("jbisanz/qiime2R")

library(tidyverse)
library(vegan)
library(qiime2R)
  


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# To get these files you need to scp them from the cluster:
#
# first on  your laptop cd to the directory where you want to save them.
# Then use this code for our example dataset today:
# mkdir core-metrics-results/
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/* .
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/core-metrics-results/* core-metrics-results/.
##############################################


###Set your working directory
#setwd("path/to/ANSC516/ANSC-repo/ANSC516/data/moving-pictures")
setwd("~/OneDrive/OneDrive Docs/Purdue/2023 Bioinformatics class/Project analysis/R")

list.files()

if(!dir.exists("output"))
  dir.create("output")

#How to load a file into R
metadata2 <- read.delim("sample-metadata.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata2[1,]
metadata2[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata2 <- metadata2[-1,]

#Now the qiime2R method
metadata<-read_q2metadata("sample-metadata.txt")
str(metadata)
# levels(metadata$`body-site`)
# colnames(metadata)[3] <- "body.site"
# colnames(metadata)[8] <- "reported.antibiotic.usage"
# colnames(metadata)[9] <- "days.since.experiment.start"
# str(metadata)

row.names(metadata) <- metadata[,1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")
wUF <- read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")

body_colors <- c("Black", "Blue", "Green", "Gray")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Now we are going to make an ordination plot
# ggplot(bc_meta, aes(x=PC1, y=PC2, color=dilution)) +
#   geom_point() + #alpha controls transparency and helps when points are overlapping
#   theme_q2r() +
#   xlab("PC1 (32.27%)") +
#   ylab("PC2 (22.28%)") +
#   scale_color_manual(values=c("Blue", "Black", "Green", "Gray"), name = "dilution")

# Get the proportion of variance explained by each principle component
# PC1 = bc_PCoA$data$ProportionExplained[1]*100
# PC2 = bc_PCoA$data$ProportionExplained[2]*100
# Can also get from bray_curtis_emporer.qzv (look at it with qiime2 viewer)


# Now we are going to make our code a little more re-usable
my_column <- "frequency"
# my_column <- "dilution"

# For this one, we alternate using my_colum = frequency and my_column = frequency
# by correspondingly alternating facet_grid(~frequency) and facet_grid(~dilution)
# ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
#   geom_point() + #alpha controls transparency and helps when points are overlapping
#   theme_q2r() +
#   facet_grid(~frequency) +
#   xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
#   ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
#   scale_color_manual(values=body_colors, name = my_column)
# #ggsave(paste0("output/BC-basic_", my_column,".svg"), height=2, width=8, device="pdf") # save a PDF 3 inches by 4 inches

  
centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "frequency" # should be the same as the var indicated in my_column

# ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
#   geom_point() + #alpha controls transparency and helps when points are overlapping
#   geom_point(data=centroids, size = 3) +
#   theme_q2r() +
#   stat_ellipse(level = 0.95, type = "t") +
#   xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
#   ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
#   scale_color_manual(values=body_colors, name = my_column)
# #ggsave(paste0("output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


# This is the one that shows different dilutions with symbols and diff freqs with colors
# To make this plot we want my_column = frequency
# ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
#   geom_point(aes(shape= dilution)) + #alpha controls transparency and helps when points are overlapping
#   #geom_point(data=centroids, size = 3) +
#   theme_q2r() +
#   stat_ellipse(level = 0.95, type = "t") +
#   xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
#   ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
# #scale_color_manual(values=body_colors, name = my_column)
# # ggsave(paste0("output/BC-ellipse_", my_column,"-dilution.svg"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

## SAME thing but with weighted UniFrac

  
Wuni_PCoA<-read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

# ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
#   geom_point() + #alpha controls transparency and helps when points are overlapping
#   #geom_point(data=centroids, size = 3) +
#   theme_q2r() +
#   stat_ellipse(level = 0.95, type = "t") +
#   xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
#   ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
#   scale_color_manual(values=body_colors, name = "Body Site")
# #ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

# This is the one that shows different dilutions with symbols and diff freqs with colors
# To make this plot we want my_column = frequency
# ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
#   geom_point(aes(shape= dilution)) + #alpha controls transparency and helps when points are overlapping
#   #geom_point(data=centroids, size = 3) +
#   theme_q2r() +
#   stat_ellipse(level = 0.95, type = "t") +
#   xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
#   ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
#   #scale_color_manual(values=body_colors, name = "Body Site")
# # ggsave(paste0("output/Wuni-ellipse_", my_column,"-dilution.svg"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

bc_dist_mat<-read_qza("core-metrics-results/bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(bc_dm),metadata$SampleID),]
rownames(bc_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(bc_dm ~ frequency + dilution, data = metadata_sub)

write.table(PERMANOVA_out,"output/bray_curtis_Adonis_overall.csv",sep=",", row.names = TRUE)

######################################################################################
# WEIGHTED UNIFRAC
wu_dist_mat<-read_qza("core-metrics-results/weighted_unifrac_distance_matrix.qza")
wu_dm <- as.matrix(wu_dist_mat$data) 
rownames(wu_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(wu_dm),metadata$SampleID),]
rownames(wu_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(wu_dm ~ frequency + dilution, data = metadata_sub)

write.table(PERMANOVA_out,"output/weighted_unifrac_Adonis_overall.csv",sep=",", row.names = TRUE) 

######################################################################################
# UNWEIGHTED UNIFRAC
uw_dist_mat<-read_qza("core-metrics-results/unweighted_unifrac_distance_matrix.qza")
uw_dm <- as.matrix(uw_dist_mat$data) 
rownames(uw_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(uw_dm),metadata$SampleID),]
rownames(uw_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(uw_dm ~ frequency + dilution, data = metadata_sub)

write.table(PERMANOVA_out,"output/unweighted_unifrac_Adonis_overall.csv",sep=",", row.names = TRUE) 


uwUF <- read_qza("core-metrics-results/unweighted_unifrac_pcoa_results.qza")

uwUF_meta <- uwUF$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),uwUF_meta,mean)

# This is the one that shows different dilutions with symbols and diff freqs with colors
# To make this plot we want my_column = frequency
# ggplot(uwUF_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
#   geom_point(aes(shape= dilution)) + #alpha controls transparency and helps when points are overlapping
#   #geom_point(data=centroids, size = 3) +
#   theme_q2r() +
#   stat_ellipse(level = 0.95, type = "t") +
#   xlab(paste0("PC1 (", round(100*uwUF$data$ProportionExplained[1], digits = 2), "%)")) +
#   ylab(paste0("PC2 (", round(100*uwUF$data$ProportionExplained[2], digits = 2), "%)")) +
#   #scale_color_manual(values=body_colors, name = "Body Site")
#   # ggsave(paste0("output/uwUF-ellipse_", my_column,"-dilution.svg"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

######################################################################################
# JACCARD DISTANCE

jaccard_dist_mat<-read_qza("core-metrics-results/jaccard_distance_matrix.qza")
jaccard_dm <- as.matrix(jaccard_dist_mat$data) 
rownames(jaccard_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(jaccard_dm),metadata$SampleID),]
rownames(jaccard_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(jaccard_dm ~ frequency + dilution, data = metadata_sub)

write.table(PERMANOVA_out,"output/jaccard_Adonis_overall.csv",sep=",", row.names = TRUE) 


jaccard <- read_qza("core-metrics-results/jaccard_pcoa_results.qza")

jaccard_meta <- jaccard$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),jaccard_meta,mean)

# This is the one that shows different dilutions with symbols and diff freqs with colors
# To make this plot we want my_column = frequency
ggplot(jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= dilution)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*jaccard$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jaccard$data$ProportionExplained[2], digits = 2), "%)")) +
  #scale_color_manual(values=body_colors, name = "Body Site")
  ggsave(paste0("output/jaccard-ellipse_", my_column,"-dilution.svg"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

body.site_Pair <- pairwise.adonis2(bc_dm ~ body.site, data = metadata_sub)
write.table(body.site_Pair,"output/Body.site_Adonis_pairwise.csv",sep=",", row.names = TRUE) 


