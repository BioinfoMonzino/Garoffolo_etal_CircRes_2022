#######################################################################################################
###      Code for differential expression analysis of:                                              ###
###      "Strain-sensitive fibrotic programming of human cardiac stromal cells can be reverted      ###
###       by interfering with YAP-dependent transcriptional activation"                             ###
#######################################################################################################

# Note: Rawcount matrix amd metadata were deposited @ GSE203358

#######################################################################################################

## Clean R global einvironment:
rm(list=ls())

## Load libraries:
library(EDASeq)
library(RUVSeq)
library(edgeR)
library(RColorBrewer)
library(DESeq2)
library(DaMiRseq)
library(Hmisc)

###############################################################
##################        Import Data        ##################
###############################################################
setwd("User_Directory") # set working directory
Raw_Counts <- read.delim("User_Directory/Raw_Counts.txt", stringsAsFactors=TRUE)
metadata <- read.delim("User_Directory/metadata.txt", stringsAsFactors=TRUE)


###############################################################
###############        Data Preprocessing       ###############
###############################################################
## Removed genes with 5 counts in at least 25% of samples (i.e. 6 samples):
idx <- rowSums(Raw_Counts>5) >= 6
Raw_Counts_filt <- Raw_Counts[idx,]


###############################################################
###     Identifiying latent variables by RUVSeq package    ####
###############################################################
## Use edgeR and RUVr method of RUVSeq for estimating latent variables,
## associated with unwanted variation (see RUVseq Manual)

## Step1: Assess the residuals of differential analysis, by edgeR:
set <- newSeqExpressionSet(Raw_Counts_filt, phenoData = data.frame(metadata), row.names=colnames(Raw_Counts_filt))
set <- betweenLaneNormalization(set, which="upper", offset=TRUE)
design <- model.matrix(~0+cell_line_ID+treatment, data=pData(set))
y <- DGEList(counts=counts(set), group=treatment)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmQLFit(y, design)
res <- residuals(fit, type="deviance")

## Step2: Use RUVr to estimate the factors of unwanted variation using residuals:
genes <- row.names(y)
# Testing different number of k:
set_k1 <- RUVr(set, cIdx=genes, k=1, res)
set_k2 <- RUVr(set, cIdx=genes, k=2, res)
set_k3 <- RUVr(set, cIdx=genes, k=3, res)
set_k8 <- RUVr(set, cIdx=genes, k=8, res)
set_k10 <- RUVr(set, cIdx=genes, k=10, res)

## Step3: Generate diagnostic plots to identify the optimal number (k) of latent variables
# plot RLE:
plotRLE(set@assayData$counts, outline=FALSE, col=colors[treatment],main="Raw Counts", las=2, cex.axis=0.7)
plotRLE(cpm(y$counts), outline=FALSE, col=colors[treatment],main="cpm", las=2, cex.axis=0.7)
plotRLE(set@assayData$normalizedCounts, outline=FALSE, col=colors[treatment], main="upperquartile", las=2, cex.axis=0.7)
plotRLE(set_k8@assayData$normalizedCounts, outline=FALSE, col=colors[treatment], main="W_8", las=2, cex.axis=0.7)
# plot PCA:
plotPCA(set@assayData$counts, outline=FALSE, col=colors[treatment],main="Raw Counts", cex=0.8)
plotPCA(cpm(y$counts), outline=FALSE, col=colors[treatment],main="cpm", cex=0.8)
plotPCA(set@assayData$normalizedCounts, outline=FALSE, col=colors[treatment], main="upperquartile", cex=0.8)
plotPCA(set_k8@assayData$normalizedCounts, outline=FALSE, col=colors[treatment], main="W_8", cex=0.8)


###############################################################
###########           Differential Analysis          ##########
###############################################################
## performing differential expression aanalysis, adjusting for latent variables (Ws):

design <- model.matrix(~0+treatment + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + cell_line_ID, data=pData(set_k8))
colnames(design) <- c("CTRL","TGFB","TGFB_VTP","VTP","W1","W2","W3","W4","W5","W6","W7","W8","M07","M078","M081","M089","M099")
y <- DGEList(counts=counts(set_k8), group=treatment)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y, design)

# set the contrast matrix:
my_contrasts <- makeContrasts(TGFBvsCTRL=TGFB-CTRL,
                              VTPvsCTRL=VTP-CTRL,
                              TGFB_VTPvsTGFB=TGFB_VTP-TGFB,
                              TGFB_VTPvsCTRL=TGFB_VTP-CTRL,
                              TGFBvsVTP = TGFB-VTP,
                              levels = design)
fit2 <- glmQLFit(y, design)

# chose a specific contrast:
qlf <- glmQLFTest(fit2, contrast = my_contrasts[, c("TGFBvsCTRL")])
DE.Results.table <- topTags(qlf, n = nrow(set))$table

# plot histogram of pvalue distribution for a specific contrast:
hist(DE.Results.table[,4], main = "TGFB vs CTRL", breaks = 20, col = "green", xlab = "pval")
# summary of DE genes by filtering:
sum(DE.Results.table$FDR < 0.05 & (DE.Results.table$logFC > 0.58 | DE.Results.table$logFC < (-0.58)), na.rm=TRUE)


###############################################################
################        Export Results        #################
###############################################################
write.table(DE.Results.table, "DE.Results.txt", sep="\t", quote=FALSE)


