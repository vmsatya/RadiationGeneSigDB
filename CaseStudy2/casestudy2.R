##########################################################################
#
# Case study2: Comparison of 3 HYPOXIC Head and Neck cancer signatures 
# Toustrup_2011; Eustace_2013; Lendahl_2009 in theTCGA data and,
# Compare Pathways enriched using these signatures.
#
# Patient: TCGA HNSCC (99 treated with RT only)
# Correlation of scores using: Spearman
# Signature score computation: weighted average of genes in a
# gene signature (sig.score function in genefu package)
# Scale the patient datasets 
# 
# Pathway analysis: 
# Gene sets: GO terms
# Package: runGSAHyper in Piano package
#
##########################################################################

library(biomaRt)
library(genefu)
library(xlsx)
library(org.Hs.eg.db)
library(readxl)
library(RColorBrewer)
require(piano)
require(VennDiagram)
require(GSA)
library(corrplot)

setwd("E:/RadiationSigDB/CaseStudy2/")

# Load Hypoxic signatures from the "RadiationResponseSigs" folder
load("HypoxicSignatures_RTResponse.RData")
Toustrup <- Hypoxic_RTSignatures[[1]]
Eustace <- Hypoxic_RTSignatures[[2]]
Lendahl <- Hypoxic_RTSignatures[[3]]

# Load the patient data
load('HN-TCGA.RData')
edata1 <- HNTCGA[[1]]

g <- list(Toustrup$Gene,Eustace$Gene,Lendahl$Gene)

annot = data.frame(rownames(edata1))
colnames(annot) <- "EntrezGene.ID"
rownames(annot) <- rownames(edata1)
data <- edata1

scores <- list()
mat <- matrix(NA,length(g),ncol(edata1))
for(k in 1:length(g)){
  b <- g[[k]]
  b1 <- as.numeric(+1)
  x = data.frame(b,b,b1)
  colnames(x) <- c("probe","EntrezGene.ID","coefficient")
  x$probe <- as.character(x[,1])
  x$EntrezGene.ID <- as.character(x[,2])
  rownames(x) <- b
  res <- genefu::sig.score(x, t(data), annot,do.mapping=TRUE, signed=TRUE, verbose=TRUE)
  mat[k,] <- res$score
  scores[[k]] <- res$score
}
rownames(mat) <- c("Toustrup","Eustace","Lendahl")
colnames(mat) <- colnames(edata3)
c <- cor(t(mat))

Jpeg("Fig2A:CorPlot_Scores.jpg")
corrplot::corrplot(c, type="upper",col=brewer.pal(n=8, name="RdBu"),tl.col="black", tl.srt=45)
dev.off()

#Jpeg("BoxplotScores.jpg")
#boxplot(mat[1,],mat[2,],mat[3,],col=brewer.pal(3,"Set2"),names=c("Toustrup","Eustace","Lindahl"),ylab="Hypoxia Score",cex.lab=1.5,cex.axis=1.5)
#dev.off()

#draw.triple.venn(area1 = nrow(Toustrup), area2 = nrow(Eustace), area3 = nrow(Lendahl),
#                n12 = length(intersect(Toustrup$Gene,Eustace$Gene)),
#                n23 = length(intersect(Eustace$Gene,Lendahl$Gene)),
#                 n13 = length(intersect(Lendahl$Gene,Toustrup$Gene)),
#                n123 = length(Reduce(intersect,list(Toustrup$Gene,Eustace$Gene,Lendahl$Gene))),
#                 category = c("Toustrup", "Eustace", "Lendahl"), lty = "blank",
#                fill = brewer.pal(3,"Set2"))

################################################################################################
################################################################################################
################################################################################################

# PATHWAY ANALYSIS

gSets <- GSA.read.gmt("c5.bp.v6.2.symbols.gmt")
dfgSets <- as.data.frame(cbind(unlist(gSets$genesets))) ## genes
dfgSNames <- as.data.frame(cbind(unlist(gSets$geneset.names)))
listS <- lapply(gSets$genesets,length)
gTogs <- data.frame(dfgSets$V1,rep(dfgSNames$V1,listS))
names(gTogs) <- c("V1","V2")
gTogs <- gTogs[gTogs$V1!="",]

a <- loadGSC(gTogs)

res.Tou <- runGSAhyper(g[[1]], gsc=a,adjMethod="fdr")
b <- data.frame(res.Tou$resTab)
b1.Tou <- subset(b,b[,"Adjusted.p.value"] < 0.1) #fdr

res.Eus <- runGSAhyper(g[[2]], gsc=a,adjMethod="fdr")
b <- data.frame(res.Eus$resTab)
b1.Eus <- subset(b,b[,"Adjusted.p.value"] < 0.1) #fdr

res.Len <- runGSAhyper(g[[3]], gsc=a,adjMethod="fdr")
b <- data.frame(res.Len$resTab)
b1.Len <- subset(b,b[,"Adjusted.p.value"] < 0.1) #fdr

n12 = length(intersect(rownames(b1.Tou),rownames(b1.Eus)))
n13 = length(intersect(rownames(b1.Tou),rownames(b1.Len)))
n23 = length(intersect(rownames(b1.Eus),rownames(b1.Len)))
n123 = length(Reduce(intersect,list(rownames(b1.Tou),rownames(b1.Eus),rownames(b1.Len))))

#draw.triple.venn(area1 = nrow(b1.Tou), area2 = nrow(b1.Eus), area3 = nrow(b1.Len), n12, n23, n13,
#                 n123, category = c("Toustrup", "Eustace", "Lendahl"), lty = "blank",
#                 fill = c("skyblue", "pink1", "mediumorchid"),cex = 1.5,
#                 cat.cex = 2,cat.col=c("skyblue", "pink1", "mediumorchid"))


pdf("Fig2B-Pathways-HypoxiaSigs.pdf")
draw.triple.venn(area1 = nrow(b1.Tou), area2 = nrow(b1.Eus), area3 = nrow(b1.Len), n12, n23, n13,
                 n123, category = c("Toustrup", "Eustace", "Lendahl"), lty = "blank",
                 fill = brewer.pal(3,"Set1"),cex = 1.5,
                 cat.cex = 2,cat.col=brewer.pal(3,"Set1"))
dev.off()










