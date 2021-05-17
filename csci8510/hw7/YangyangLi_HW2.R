## March, 2021
## Author: Sambhawa Priya, PhD candidate (BICB), Blekhman Lab. 
## Contact: priya030@umn.edu

## BICB 8510, Spring 2021
## HW2 
## Based on metadata and otu tables taken from https://github.com/LangilleLab

## Preparation: Create a directory/folder named "ML_HW2" at a relevant location on your computer,
## e.g. if you have a course directory for BICB_8510, create this directory there. 
## Next, place this Rscript in the directory ML_HW2. 
## Next, create a directory within ML_HW2 called "input" and place the input files
## i.e. metadata.txt and otu_table.txt in the "input" directory. 


##install libraries
install.packages("vegan") ## for distance computation for microbiome data
## vignette: https://cran.r-project.org/web/packages/vegan/vegan.pdf
##import libraries
library(vegan) 

## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

## Input metadata
## This will be located in directory "ML_HW2/input/", if you followed the instruction above. 
metadata <- read.table(paste0(current_dir,"/input/metadata.txt"),sep="\t",header=T,row.names=1,stringsAsFactors=TRUE, comment.char="")
## Input otu table. Also should be located at "ML_HW2/input/" directory
otu <- read.table(paste0(current_dir,"/input/otu_table.txt"),sep="\t",header=T,row.names=1,stringsAsFactors=FALSE, comment.char="")

## number of rows and columns in metadata and otu tables
dim(metadata) 
#[1] 40  2
dim(otu) 
# [1] 1000  40

## Number of inflamed and controls
table(metadata$state)
# control inflamed 
# 20       20 

## Transpose otu table to make samples as rows and otus as columns.
otu <- as.data.frame(t(otu))
dim(otu)
# [1]   40 1000

## identify common samples
common_samples <- intersect(rownames(metadata), rownames(otu))
common_samples

## Ensure order of samples in metadata and otu table are identical
otu <- otu[common_samples,]
metadata <- metadata[common_samples,]

## Explore data distribution
## Pick an OTU (OTU_77)
otu_random <- otu[,grep("OTU_77$", colnames(otu))]
hist(otu_random, br = 100)

###### Q1. How many otus were found in no samples (i,e, 0 samples)? Identify the names of those otus (10 points)
colnames(otu)[apply(otu==0, 2, all)]
# [1] "OTU_223" "OTU_284"

###### Q2. How many otus were found in only 1 sample (also called "singletons")? Identify the names of those otu (20 points) 
colnames(otu)[apply(otu!=0, 2, sum) == 1]
# [1] "OTU_33"  "OTU_38"  "OTU_44"  "OTU_55"  "OTU_86"  "OTU_312" "OTU_429" "OTU_432" "OTU_706" "OTU_918"
# [11] "OTU_940" "OTU_972"

###### Q3. Remove OTUs which are present in fewer than 10% of samples (20 points)
clean_otus <- colnames(otu)[apply(otu!=0, 2, sum)/40 >= 0.1]
otu <- otu[, clean_otus]

## Create PCoA plot -- similar to PCA, except can use non-euclidean distance/disimmilarity matrix 
# compute Bray-Curtis distance/dissimilarity for the otu table
## See "vegdist" in vegan for details: https://cran.r-project.org/web/packages/vegan/vegan.pdf
d.bray <- vegdist(otu) ## bray-curtis is default for vegan:vegdist()

## Perform non-metric multidimensional scaling for ordination (PCoA)
# See for details https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cmdscale.html
pc.bray <- cmdscale(d.bray,k=2, eig = T, add = T)
col.fill <- c("#3490DE","#A50F15") ##blue for control, red for inflammed.
## Create PCoA plot
plot(pc.bray$points[,1], pc.bray$points[,2], cex=3, pch=16,  
     col = col.fill[metadata$state], xlab="PC1", ylab = "PC2")

###### Q4. Add legend for the points colored by inflammation state (ie, inflammed or control). Show script and upload the plot. (10 points)
## Hint: Use legend() after the above command. Alternatively, if you're familar with ggplot(), you can also regenerate the plot
## above using ggplot including legend. 
tiff('pcoa.tiff', width=6, height = 4, units = 'in', res=200)
label <- factor(c('inflamed','control'))
plot(pc.bray$points[,1], pc.bray$points[,2], cex=3, pch=16,  
     col = col.fill[metadata$state], xlab="PC1", ylab = "PC2")
legend(x=-0.1, y=0.2, label, fill=col.fill[label])
dev.off()

## Extract Eigenvalues from the ordination (see vegan package for details on eigenvals())
pc.bc.eig <- eigenvals(pc.bray)
pc.bc.eig

####### Q5. Compute percent variation for eignevalues (20 points)
## Hint: percent_var = eigenvalue/sum(eigenvalues) * 100
percent_var = pc.bc.eig/sum(pc.bc.eig) * 100

####### Q6. Regenerate PCoA plot showing percentage variation along PC1 and PC2. Show script and upload the plot. (20 points)
## Hint: replace xlab and ylab in the plot() command above with first two percent variation computed. 
## These correspond to percent variation along PC1 and PC2.
tiff('pcoa_with_variance2.tiff', width=6, height = 4, units = 'in', res=200)
label <- factor(c('inflamed','control'))
plot(pc.bray$points[,1], pc.bray$points[,2], cex=3, pch=16,  
     col = col.fill[metadata$state], xlab=paste('PC1', percent_var[1], '%'), ylab = paste('PC2', percent_var[1], '%'))
legend(x=-0.1, y=0.2, label, fill=col.fill[label])
dev.off()

