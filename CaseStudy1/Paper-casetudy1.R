######################################################################
#
# Case study1: COmparison of OXIC breast cancer signatures in cell line
# data and in patient data 
#
# In-vitro: Normalized CCLERNAseq with 26 breast cancer cell lines
# Patient: METABRIC
# Subtyping: SCMOD2 function in genefu package
# Correlation of scores: Spearman
# Signature score computation: weighted average of genes in a
# gene signature (sig.score function in genefu package)
# Scale the cell line and patient datasets for fair comparison
# of signature scores
# 
######################################################################

library(biomaRt)
library(genefu)
library(xlsx)
library(org.Hs.eg.db)
library(readxl)
library(RColorBrewer)

setwd("E:/RadiationSigDB/CaseStudy1/")

load('OxicSignatures_RTResponse.RData')
# Piening oxic signature
Piening.brca1 <- Oxic_RTsignatures[[1]] 
# Speers oxic signature
speers.brca.sf2 <- Oxic_RTsignatures[[4]]

# load cell line data
load('CCLE_RNAseq-Breast.RData')

### Computation of signature score
# Piening signature
annot = data.frame(rownames(edata4))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(edata4)
data <- edata4
b <- Piening.brca1$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res1 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

# speers
b <- speers.brca.sf2$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res2 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

cor.test(res1$score,res2$score,method = "spearman")
# cor = 0.701 with p value = 0.0001007

###########################################################
###########################################################
###########################################################

# METABRIC

load('TNBC_METABRIC_RT.RData')
load('her2_METABRIC_RT.RData')
load('erher2.Highprolif_METABRIC_RT.RData')
load('erher2.Lowprolif_METABRIC_RT.RData')
metabric <- cbind(tnbc.metabric1,her2.metabric1,erher2.Highprolif.metabric1,erher2.Lowprolif.metabric1)
metabric1 <- metabric

annot = data.frame(rownames(metabric1))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(metabric1)
data <- metabric1
b <- Piening.brca1$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res1 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

# speers
b <- speers.brca.sf2$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res2 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

cor.test(res1$score,res2$score,method = "spearman")

# cor = 0.686 and p-value < 2.2e-16

#########################
# Based on subtypes

load('TNBC_METABRIC_RT.RData')

annot = data.frame(rownames(tnbc.metabric1))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(tnbc.metabric1)
data <- tnbc.metabric1
b <- Piening.brca1$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res1 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

# speers
b <- speers.brca.sf2$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res2 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

cor.test(res1$score,res2$score,method = "spearman")
# cor = 0.66 and p-value < 1.7e-08

# her2
load('her2_METABRIC_RT.RData')

annot = data.frame(rownames(her2.metabric1))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(her2.metabric1)
data <- her2.metabric1
b <- Piening.brca1$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res1 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

# speers
b <- speers.brca.sf2$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res2 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

cor.test(res1$score,res2$score,method = "spearman")
# cor = 0.664 and p-value = 0.001

# erhigh
load('erher2.Highprolif_METABRIC_RT.RData')

annot = data.frame(rownames(erher2.Highprolif.metabric1))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(erher2.Highprolif.metabric1)
data <- erher2.Highprolif.metabric1
b <- Piening.brca1$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res1 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

# speers
b <- speers.brca.sf2$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res2 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

cor.test(res1$score,res2$score,method = "spearman")
# cor = 0.258 and p-value < 0.06

# erlow
load('erher2.Lowprolif_METABRIC_RT.RData')

annot = data.frame(rownames(erher2.Lowprolif.metabric1))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(erher2.Lowprolif.metabric1)
data <- erher2.Lowprolif.metabric1
b <- Piening.brca1$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res1 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

# speers
b <- speers.brca.sf2$EntId
b1 <- as.numeric(+1)
x = data.frame(b,b,b1)
colnames(x) <- c("probe","EntrezGene.ID","coefficient")
x$probe <- as.character(x[,1])
x$EntrezGene.ID <- as.character(x[,2])
x <- x[complete.cases(x),]
rownames(x) <- x$EntrezGene.ID
res2 <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)

cor.test(res1$score,res2$score,method = "spearman")
# cor = 0.60 and p-value < 2.2e-16

############################################################################
############################################################################
############################################################################

# plot

a <- c(0.701, 0.686,0.67,0.664, 0.25, 0.60) 
pdf("Brsig-barplot.pdf")
par(mar=c(9, 6, 4, 2) + 0.1)
a1=barplot(a, col=brewer.pal(6,"Set2"), las=1, names.arg="",ylab="Correlation (Spearman)",cex.axis = 1.3,cex.lab=1.5,ylim = c(0,0.80))
text(a1[,1], -0.02,srt = 60, adj= 1, xpd = TRUE, labels = c("CCLE-RNAseq","METABRIC-All","TNBC","HER2","ERHigh","ERLow"), 
     cex=1.3,cex.lab=1.3,cex.axis=1.3)
dev.off()






