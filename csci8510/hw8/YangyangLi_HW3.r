## Mar, 2021
## Author: Sambhawa Priya, PhD candidate (BICB), Blekhman Lab. 
## Contact: priya030@umn.edu

## This homework has been adapted based on the following papers:
## Code and data adapted from Zhou et al., 2019: https://www.frontiersin.org/articles/10.3389/fgene.2019.00579/full
## Dataset orginally published at:  
## Goodrich et al. 2014: http://dx.doi.org/10.1016/j.cell.2014.09.053
## Compiled by Duvallet et al.: https://www.ncbi.nlm.nih.gov/pubmed?Db=pubmed&Cmd=ShowDetailView&TermToSearch=29209090

## Preparation: Create a directory/folder named "ML_HW2" at a relevant location on your computer,
## e.g. if you have a course directory for BICB_8510, create this directory there. 
## Next, place this Rscript in the directory ML_HW1. 
## Next, create a directory within ML_HW1 called "input" and place the input files
## i.e. metadata.txt and otu_table.txt in the "input" directory. 

## Initialization
rm(list=ls()) ## Don't do this if other objects in workspace that you need later.

library(caret) ## Most popular R package for machine learning
library(randomForest)
library(pROC)
library(stringr)


## In Rstudio, find the path to the directory where the current script is located.
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)

############### Input metadata and otu table #################
## Goodrich et al. data
metadata <- read.table(paste0(current_dir,"/input/ob_goodrich.metadata.txt"),sep="\t",header=T,row.names=1, stringsAsFactors=TRUE, check.names = F)
dim(metadata)


#Keep samples with n_sample == 0 and disease-states "OB" and "H"
metadata <- metadata[(metadata$n_sample == 0 & metadata$DiseaseState %in% c("H","OB")),]
dim(metadata)


#select one sample per individual and one individual per twin pair
metadata <- metadata[!duplicated(metadata$familyid),]
dim(metadata)


## drop unused level from metadata
metadata <- droplevels(metadata)
dim(metadata)
## How many samples per disease state?
table(metadata$DiseaseState)
# #H  OB 
# 279 135
## This will be our classification problem!

## Please be patient. This will take sometime (1-2 mins).
otu <- read.table(paste0(current_dir,"/input/ob_goodrich.otu_table.100.denovo.rdp_assigned"),sep="\t",header=T,row.names=1, stringsAsFactors=TRUE, check.names = F)
dim(otu)

## identify common samples between metadata and otu table
common_samples <- intersect(rownames(metadata), colnames(otu))
length(common_samples)

## Ensure order of samples in metadata and otu table are identical
otu <- otu[,common_samples]
dim(otu) 

metadata <- metadata[common_samples,]
dim(metadata) 

all(rownames(metadata) == colnames(otu)) 
# should be TRUE

##################### Do Preprocessing (see specific instructions where provided) ########################

## Q1. Remove OTUs with fewer than 10 reads (10 points)
otu <- otu[which(rowSums(otu)>=10),]
dim(otu) 
# [1] 34490   414

## Q2. Remove OTUs which were present in fewer than 5% of samples (10 points)
otu <- otu[which(rowSums(otu>0) >= ncol(otu)*.05),]
dim(otu) 
# [1] 11225   414

## Binning
# collapse OTUs to genus level by summing their abundance counts
## Split taxa names by ";" into 8 parts
tax_table <- str_split_fixed(rownames(otu),";",n=8)
## Append genus name (6th column in genus) to the otu table
otu <- data.frame(genus=tax_table[,6],otu)
## Filter out otus where genus is not characterized. 
otu <- otu[!(otu$genus=="g__"),]
dim(otu)

## summarize otu table by genus label
otu <- aggregate(otu[,-1], by=list(otu$genus), sum)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
dim(otu)


#calculate relative abundance of each genus by dividing its value by the total reads per sample
otu <- sweep(otu,2,colSums(otu),"/")

## Prepare for training
x <- data.matrix(otu)
## transpose to make rows are samples and feature (i.e. genera) as columns. 
x <- t(x)

## relevel disease state to allow "OB" to be case for ML model downstream
levels(metadata$DiseaseState) #[1] [1] "H"  "OB"
metadata$DiseaseState <- relevel( metadata$DiseaseState, "OB")
levels(metadata$DiseaseState) #[1] "OB" "H" 


################ Train and Test ###############
## Q3. Split data by 90% training and 10% test, and report the output of training (best mtry),
## and the output of testing (confusion matrix, sensitivity, specificity, precision, AUC).  
## Also show the ROC curve. (30 points)
ratio <- 0.9 
set.seed(1000)
train_index <- createDataPartition(metadata$DiseaseState, ## outcome
                                   p = ratio, ## percentage of training samples
                                   list = FALSE) ## show subsamples as matrix, not list
x.train <- x[train_index,] 
y.train <- metadata$DiseaseState[train_index]
x.test <- x[-train_index,]
y.test <- metadata$DiseaseState[-train_index]

train_control <- trainControl(
  method = "cv",
  number = 5, ## also try 10
  summaryFunction=twoClassSummary, # computes area under the ROC curve
  classProbs = TRUE ## required for scoring models using ROC
)
set.seed(1000)
rf_train <- train( x = x.train, y = as.factor(y.train),
                   method='rf',
                   metric="ROC", ## default accuracy
                   trControl = train_control)


rf_train
# mtry  ROC        Sens        Spec     
# 2    0.6696120  0.04966667  0.9803137
# 50    0.6305593  0.19700000  0.9286275
# 99    0.6109578  0.23800000  0.9247843

## mtry: Number of variables randomly sampled as candidates at each split.

rf_train$resample
# ROC       Sens      Spec Resample
# 1 0.5522876 0.04166667 0.9411765    Fold1
# 2 0.7047059 0.04000000 0.9803922    Fold2
# 3 0.7075000 0.16666667 1.0000000    Fold5
# 4 0.6391667 0.00000000 0.9800000    Fold4
# 5 0.7444000 0.00000000 1.0000000    Fold3


rf_test <- predict(rf_train, x.test) 
rf_test

# compare predicted outcome and true outcome
conf_matrix <- confusionMatrix(rf_test, y.test)
conf_matrix$table
# Reference
# Prediction OB  H
# OB  1  2
# H  12 25

## Compute precision

conf_matrix$byClass
#          Sensitivity          Specificity       Pos Pred Value       Neg Pred Value            Precision               Recall 
#           0.07692308           0.92592593           0.33333333           0.67567568           0.33333333           0.07692308 
#                   F1           Prevalence       Detection Rate Detection Prevalence    Balanced Accuracy 
#           0.12500000           0.32500000           0.02500000           0.07500000           0.50142450

## Can also spit out probability instead of predicted class
rf_test <- predict(rf_train, x.test, type = "prob")
rf_test
rf_test <- rf_test[,1]

############ Plot ROC curve ############
## ROC curve
rf <- roc(y.test,rf_test) ## pROC package
auc <- rf$auc
auc
# Area under the curve: 0.5983

pdf(paste0(current_dir,"/0.9_ROC_Singh.pdf"))
plot(rf, col="blue",legacy.axes = T)
dev.off()


## Q4. Split data by 70% training and 30% test, and report the output of training (best mtry),
## and the output of testing (confusion matrix, sensitivity, specificity, precision, AUC).  
## Also show the ROC curve. Did the AUC increase or decrease compared to result from Q3? (30 points) 

ratio <- 0.7 
set.seed(1000)
train_index <- createDataPartition(metadata$DiseaseState, ## outcome
                                   p = ratio, ## percentage of training samples
                                   list = FALSE) ## show subsamples as matrix, not list
x.train <- x[train_index,] 
y.train <- metadata$DiseaseState[train_index]
x.test <- x[-train_index,]
y.test <- metadata$DiseaseState[-train_index]

train_control <- trainControl(
  method = "cv",
  number = 5, ## also try 10
  summaryFunction=twoClassSummary, # computes area under the ROC curve
  classProbs = TRUE ## required for scoring models using ROC
)
set.seed(1000)
rf_train <- train( x = x.train, y = as.factor(y.train),
                   method='rf',
                   metric="ROC", ## default accuracy
                   trControl = train_control)


rf_train
# mtry  ROC        Sens       Spec     
# 2    0.6858266  0.1157895  0.9693590
# 50    0.6804656  0.2842105  0.9080769
# 99    0.6578947  0.2526316  0.8932051

## mtry: Number of variables randomly sampled as candidates at each split.

rf_train$resample
# ROC       Sens      Spec Resample
# 1 0.5796053 0.10526316 0.9750000    Fold1
# 2 0.7489879 0.05263158 1.0000000    Fold2
# 3 0.7692308 0.21052632 1.0000000    Fold5
# 4 0.6093117 0.10526316 0.8974359    Fold4
# 5 0.7219973 0.10526316 0.9743590    Fold3


rf_test <- predict(rf_train, x.test) 
rf_test

# compare predicted outcome and true outcome
conf_matrix <- confusionMatrix(rf_test, y.test)
conf_matrix$table
# Reference
# Prediction OB  H
# OB  3  4
# H  37 79

## Compute precision

conf_matrix$byClass
#   Sensitivity          Specificity       Pos Pred Value       Neg Pred Value            Precision               Recall 
#   0.07500000           0.95180723           0.42857143           0.68103448           0.42857143           0.07500000 
#           F1           Prevalence       Detection Rate Detection Prevalence    Balanced Accuracy 
#   0.12765957           0.32520325           0.02439024           0.05691057           0.51340361 

## Can also spit out probability instead of predicted class
rf_test <- predict(rf_train, x.test, type = "prob")
rf_test
rf_test <- rf_test[,1]

############ Plot ROC curve ############
## ROC curve
rf <- roc(y.test,rf_test) ## pROC package
auc <- rf$auc
auc
# Area under the curve: 0.5509

pdf(paste0(current_dir,"/0.7_ROC_Singh.pdf"))
plot(rf, col="blue",legacy.axes = T)
dev.off()



######## AUC decreases compared toresult from Q3.



########## List feature importance in random forest ###########
## Q5. Show top 10 features ranked by their importance corresponding to Q4 output. (20 points) 
## Hint: Use varImp() as shown in class. Sort by importance and pick top-10 features


impVars <- varImp(rf_train)

ImpMeasure <- data.frame(impVars$importance)

ImpMeasure[order(ImpMeasure$Overall, decreasing = T),][1:10]
# [1] 100.00000  94.39560  83.37954  81.72939  81.20243  80.11882  77.37541  76.47546  74.60497  74.25717

row.names(ImpMeasure)[  order(ImpMeasure$Overall, decreasing = T) ][1:10]
# [1] "g__Sporobacter"      "g__Roseburia"        "g__Gemmiger"         "g__Oscillibacter"    "g__Clostridium_XlVb"
# [6] "g__Barnesiella"      "g__Dorea"            "g__Clostridium_IV"   "g__Blautia"          "g__Clostridium_XI"





