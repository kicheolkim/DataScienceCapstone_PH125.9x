###################################################################################################
##### Capstone - Choose Your Own project #####
###################################################################################################
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# BiocManager::install("tximport")
# BiocManager::install("DESeq2")
# BiocManager::install("GEOquery")

invisible(
  lapply(c("GEOquery","tidyverse","DESeq2","tximport","RColorBrewer",
           "pheatmap","ggbeeswarm","gridExtra","caret","pROC"), 
         library, character.only=TRUE)
)

rm(list = ls())
workDir <- "~/Desktop/Learning/edX-DataScience/9_Capstone/2_ChooseYourOwn/CYO_R"
setwd(workDir)



##### Loading and QC dataset ####
# The dataset (gene expression) has been downloaded from GEO (Acc. GSE137143; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137143)
# This dataset includes 3 different cell types (CD4+ T cell, CD8+ T cells, CD14+ monocytes). 
# Since I have interested in CD4+ T cells, I kept only CD4+ T cells and removed CD8+ & CD14+ cells.

### Metadata - downloading from GEO
gse <- getGEO("GSE137143", GSEMatrix=TRUE)
sTable <- pData(phenoData(gse[[1]]))[,c(1,48:54)]
colnames(sTable) <- make.names(colnames(sTable))
colnames(sTable) <- str_replace(colnames(sTable), ".ch1", "")
# patients used in the analysis are all treatment-naive patients
# disease status: SP (secondary progressive) and PP (primary progressive) considers PMS (progressive MS)
sTable <- sTable %>% rownames_to_column("geoID") %>% 
  mutate(disease.state2 = str_split(title, ", ", simplify = TRUE)[,1],
         sampleID = str_split(title, ", ", simplify = TRUE)[,2],
         cell.type2 = str_extract(cell.type, pattern = "CD\\d+"),
         age.at.exam = as.numeric(age.at.exam),
         edss = as.numeric(edss),
         disease.state2 = case_when(disease.state2 == "Healthy controls" ~ "HC", 
                                    disease.state2 == "Treatment naïve MS patients" ~ "MS"),
         disease.subtype2 = case_when(disease.subtypes %in% c("PP","SP") ~ "PMS",
                                      disease.subtypes %in% c("CIS") ~ "CIS",
                                      disease.subtypes %in% c("RR") ~ "RR",
                                      disease.state2 == "HC" ~ "HC")) %>% 
  mutate(disease.subtype2 = factor(disease.subtype2))
sTable <- sTable %>% filter(!is.na(disease.subtype2))
# binned age into 5 groups since age will be used as a covariate in differential expression analysis
sTable$age.at.exam[is.na(sTable$age.at.exam)] <- median(sTable$age.at.exam, na.rm = TRUE)
sTable$AgeAtExamGrp <- cut(sTable$age.at.exam, breaks=c(10, 30, 40, 50, 60, 90), 
                           labels = c("10_30","30_40","40_50","50_60","60_90"))
row.names(sTable) <- sTable$sampleID
rm(gse)

# select CD4 cells for further analysis
sTable.cd4 <- sTable %>% filter(cell.type2 %in% "CD4")



###### Count dataset - CD4 cells ####
## I've downloaded RSEM outputs from GEO (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE137143&format=file) and decompressed all files.
## Count matrix created using tximport.
# files <- list.files(path="GSE137143_RAW", pattern="*.genes.results.txt", recursive = TRUE, full.names = TRUE)
# names(files) <- unlist(lapply(strsplit(files, "/"), "[[",2))
# names(files) <- unlist(lapply(strsplit(names(files), "_"), "[[", 2))
# names(files) <- str_replace(names(files), ".genes.results.txt", "")
# files.cd4 <- files[sTable.cd4$sampleID]

### Loading RSEM gene counts into DESeq2 - CD4
# rsem.genes.cd4 <- tximport(files.cd4, type = "rsem", txIn = FALSE, txOut = FALSE)
# rsem.genes.cd4$length[rsem.genes.cd4$length == 0] <- 1
# save(rsem.genes.cd4, file = "rsem.genes.cd4.RData")

### Normalization of gene expression using DESeq2 and variance stabilization transformation with all samples
## confirm matching between 'column of sample table' and 'row of count table'
## loading count matrix from prepared data. R code for this matrix is in above.
load(paste0(workDir,"/rsem.genes.cd4.RData"))
mean(colnames(rsem.genes.cd4$counts) == sTable.cd4$sampleID)

## running DESeq2 and get normalized count matrix & vst transformed matrix
cds <- DESeqDataSetFromTximport(txi = rsem.genes.cd4, 
                                colData = sTable.cd4, 
                                design = ~ disease.subtype2 + gender + AgeAtExamGrp - 1)
cds <- cds[ rowSums(counts(cds)) > 1, ]
cds <- DESeq(cds)

rm(rsem.genes.cd4, files.cd4)



###### PCA plot to checking outlier samples ####
vst.cds <- vst(cds, blind = TRUE)
pcaData <- plotPCA(vst.cds, intgroup=c("disease.state2", "gender"), ntop=10000, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=disease.state2, shape=gender, label = name)) +
  geom_text(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

rm(vst.cds, percentVar, pcaData)



##### Differential expression analysis without outlier samples ####
### Normalization and differential expression analysis without outlier samples
## based on PCA plot, considered 2 samples as outliers: 75216a-CD4, 84917a-CD4
sTable.cd4 <- sTable.cd4 %>% filter(!sampleID %in% c("75216a-CD4", "84917a-CD4"))

## running DESEq2 without outlier samples
cds <- cds[, sTable.cd4$sampleID]   # select samples without outlier
cds <- cds[ rowSums(counts(cds)) > 1, ]   # filter low expressed genes (total raw count > 1)
cds <- DESeq(cds)
vst.cds <- vst(cds)



###### DEG (differentially expressed genes) ####
### get DESeq2 results
res_hc_cis <- results(cds, alpha=0.2, contrast=c("disease.subtype2","HC","CIS"))
res_hc_cis_df <- as.data.frame(res_hc_cis) %>% arrange(pvalue)

summary(res_hc_cis)



###### Heatmap ####
heatmap.expr <- heatmap.expr <- assay(normTransform(cds))
heatmap.expr <- heatmap.expr - rowMeans(heatmap.expr)

select <- c(res_hc_cis_df %>% rownames_to_column("gene") %>%  filter(baseMean > 1) %>% 
              arrange(pvalue) %>% top_n(100, -pvalue) %>% pull(gene))
length(select)

df <- data.frame(Subtype = sTable.cd4$disease.subtype2, Status = sTable.cd4$disease.state2)
row.names(df) <- sTable.cd4$sampleID
pheatmap(heatmap.expr[select,], color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),   
         scale = "row", clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", clustering_method = "ward.D2", 
         cluster_rows=TRUE, cluster_cols=TRUE, fontsize = 8, 
         show_rownames=FALSE, show_colnames = FALSE, annotation_col=df,
         main = "Heatmap with DEG top 100")  

rm(select, df, heatmap.expr)




##### Machine learning ####
###  

## select CIS and HC samples only
sTable.cd4.cis <- sTable.cd4 %>% filter(disease.subtype2 %in% c("HC","CIS")) %>% 
  mutate(disease.subtype2 = factor(disease.subtype2))

## variance stabilization transformed counts
count.vst <- assay(vst.cds)
count.vst <- count.vst[, sTable.cd4.cis$sampleID]
dim(count.vst)

## random split training (80%) and test (20%) sets
test_idx <- createDataPartition(sTable.cd4.cis$disease.subtype2, times=1, p=0.2, list=FALSE)     # random split
test_set <- t(count.vst[,test_idx])
train_set <- t(count.vst[,-test_idx])
test_set_meta <- sTable.cd4.cis[test_idx,]
train_set_meta <- sTable.cd4.cis[-test_idx,]

## function to retrieve modeling output
fit_output <- function(fit, test_set, test_set_meta, model_name, gene_n){
  confusionMatrix <- confusionMatrix(predict(fit, test_set), test_set_meta$disease.subtype2)
  confusionMatrix
  
  pred_for_roc <- as.data.frame(predict(fit, test_set, type="prob"))
  pred_for_roc$predict <- names(pred_for_roc)[1:2][apply(pred_for_roc[,1:2], 1, which.max)]
  pred_for_roc$observed <- test_set_meta$disease.subtype2
  
  roc_obj <- roc(pred_for_roc$observed, as.numeric(pred_for_roc$CIS))
  
  ml_res <- data.frame(Method = model_name, NumAnalyte = gene_n,
                       AUC = auc(roc_obj),
                       Accuracy = confusionMatrix$overall["Accuracy"],
                       Sensitivity = confusionMatrix$byClass["Sensitivity"],
                       Specificity = confusionMatrix$byClass["Specificity"],
                       row.names = NULL)
  
  return(ml_res)
}


###### 1. Find best method and number of analytes ####
### repeat 5 different modeling methods (lda, rf, glmnet, knn, svmlinear) with 11 different number of analytes (10 genes ~ 500 genes)
gene_n <- c(10, seq(50, 500, 50))

ml_res_num <- data.frame(Method = character(), NumAnalyte = numeric(), AUC = numeric(), 
                     Accuracy = numeric(), Sensitivity = numeric(), Specificity = numeric())
for (gene_n in gene_n){
  topGenes <- res_hc_cis_df %>% rownames_to_column("gene") %>% 
    filter(!is.na(padj)) %>% top_n(gene_n, -pvalue) %>% pull(gene)
  
  lda_fit <- train(x=train_set[,topGenes], y=train_set_meta$disease.subtype2, 
                   method = "lda")
  ml_res_num <- rbind(ml_res_num, fit_output(lda_fit, test_set, test_set_meta, "lda", gene_n))
  
  rf_fit <- train(x=train_set[,topGenes], y=train_set_meta$disease.subtype2, 
                  method="rf", ntree = 1000, importance = TRUE)
  ml_res_num <- rbind(ml_res_num, fit_output(rf_fit, test_set, test_set_meta, "rf", gene_n))
  
  glm_fit <- train(x=train_set[,topGenes], y=train_set_meta$disease.subtype2, 
                   method = "glmnet")
  ml_res_num <- rbind(ml_res_num, fit_output(glm_fit, test_set, test_set_meta, "glmnet", gene_n))
  
  knn_fit <- train(x=train_set[,topGenes], y=train_set_meta$disease.subtype2, 
                   method = "knn")
  ml_res_num <- rbind(ml_res_num, fit_output(knn_fit, test_set, test_set_meta, "knn", gene_n))
  
  svml_fit <- train(x=train_set[,topGenes], y=train_set_meta$disease.subtype2, 
                    method = "svmLinear2", probability = TRUE)
  ml_res_num <- rbind(ml_res_num, fit_output(svml_fit, test_set, test_set_meta, "svmLinear2", gene_n))
  
  rm(topGenes, lda_fit, rf_fit, glm_fit, knn_fit, svml_fit)
}

## creating plot - accuracy results from the test
ml_res_num %>% 
  ggplot(aes(NumAnalyte, Accuracy)) + geom_line(aes(color=Method)) + 
  geom_vline(xintercept = max(ml_res_num$NumAnalyte[which(ml_res_num$Accuracy == max(ml_res_num$Accuracy))]), color = "red") +
  geom_vline(xintercept = min(ml_res_num$NumAnalyte[which(ml_res_num$Accuracy == max(ml_res_num$Accuracy))]), color = "blue") +
  ggtitle(paste("Best method and number of analytes for accuracy"))

## creating plot - AUC results from the test
ml_res_num %>% 
  ggplot(aes(NumAnalyte, AUC)) + geom_line(aes(color=Method)) + 
  geom_vline(xintercept = max(ml_res_num$NumAnalyte[which(ml_res_num$AUC == max(ml_res_num$AUC))]), color = "red") +
  geom_vline(xintercept = min(ml_res_num$NumAnalyte[which(ml_res_num$AUC == max(ml_res_num$AUC))]), color = "blue") +
  ggtitle(paste("Best method and number of analytes for AUC"))



###### 2. rf: Random Forest ####
### Because dataset is small, random split cannot reflect heterogeneity of human samples. 
### Therefore, I used 5-fold cross validation in the modeling
train_control <- trainControl(method="cv", number=10, savePredictions = TRUE)

## based on previous analysis, I choose random forest for final modeling
ml_res_num_rf <- ml_res_num %>% filter(Method == "rf")
## select number of analytes that showed best AUC in previous test
gene_n_best <- min(ml_res_num_rf$NumAnalyte[which(ml_res_num_rf$AUC == max(ml_res_num_rf$AUC))])
## select 100 most significant genes, 100 number is based on previous observation (best performance in random forest)
topGenes <- res_hc_cis_df %>% rownames_to_column("gene") %>% 
  filter(!is.na(padj)) %>% top_n(gene_n_best, -pvalue) %>% pull(gene)

## random forest modeling
rf_fit <- train(x=t(count.vst)[,topGenes], y=sTable.cd4.cis$disease.subtype2, 
                method="rf", trControl=train_control,
                ntree = 3000, importance = TRUE)

## retrieve predictions from cross-validation
predictions <- rf_fit$pred %>% filter(mtry == rf_fit$bestTune$mtry)
predictions$sampleID <- sTable.cd4.cis$sampleID[predictions$rowIndex]

## result - getting confusion matrix
confusionMatrix <- confusionMatrix(predictions$pred, predictions$obs)
confusionMatrix

## retrieve prediction scores from cross-validation
pred_for_roc <- as.data.frame(rf_fit$finalModel$votes) %>% rownames_to_column("sampleID")
pred_for_roc$sampleID <- str_replace(pred_for_roc$sampleID, "\\.", "\\-")
pred_for_roc$sampleID <- str_replace(pred_for_roc$sampleID, "X", "")
pred_for_roc <- left_join(pred_for_roc, predictions %>% select(sampleID, pred, obs), by="sampleID")

## result - getting ROC curve
roc_obj <- roc(pred_for_roc$obs, as.numeric(pred_for_roc$CIS),
               # arguments for ci
               ci=TRUE, ci.alpha=0.9, stratified=FALSE,
               # arguments for plot
               plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
               print.auc=TRUE, show.thres=TRUE)

## summarized results (AUC, accuraci, sensitivity, specificity)
rf_result <- data.frame(Method = "RF", NumAnalyte = gene_n_best,
                     AUC = auc(roc_obj),
                     Accuracy = confusionMatrix$overall["Accuracy"],
                     Sensitivity = confusionMatrix$byClass["Sensitivity"],
                     Specificity = confusionMatrix$byClass["Specificity"],
                     row.names = NULL)

## retrieve importance scores
rf_imp <- varImp(rf_fit)
rf_imp <- as.data.frame(rf_imp$importance) %>% rownames_to_column("gene") %>% 
  select(-HC) %>% arrange(-CIS) %>% dplyr::rename(Importance = CIS)

## adding random forest results into DEG result table
res_hc_cis_df <- res_hc_cis_df %>% rownames_to_column("gene") %>% 
  left_join(., rf_imp, by="gene")

## comparison between DEG adjusted p-value and random forest importance score
res_hc_cis_df %>% filter(!is.na(Importance)) %>% 
  mutate(label = str_split(gene, "_", simplify = TRUE)[,2]) %>% 
  ggplot(aes(padj, Importance, label = label)) + geom_point() + geom_label(size = 2) + theme_bw()



###### Boxplot for genes with top 6 importance from random forest ####
## creating boxplot for genes with highest importance score in the random forest results
p <- list()
for (i in 1:6){
  geneCounts <- plotCounts(cds, gene = rf_imp$gene[i], intgroup = c("disease.subtype2"), returnData = TRUE)
  colnames(geneCounts) <- c("norm.count","status")
  geneCounts <- geneCounts %>% filter(status %in% c("CIS", "HC"))
  p[[i]] <- ggplot(geneCounts, aes(x = status, y = norm.count, color = status)) + geom_boxplot() + 
    scale_y_log10() +  geom_beeswarm(cex = 1.5) + 
    ggtitle(paste(str_split(rf_imp$gene[i], "_", simplify = TRUE)[2], "gene")) + 
    xlab("") + theme_bw() + theme(axis.title.y = element_blank(), legend.position = "none")
}

grid.arrange(grobs = p, nrow = 2, 
             left = textGrob("log10(DESeq2 normalized count)", rot = 90, vjust = 1))

