######################################################################
#
# This script is used to test the quality of gene signatures using
# sigQC method developed by Andrew Dhawan et al. 
#
######################################################################

library(sigQC) #load the necessary package

# CaseStudy1
p <- load('BrPatientSubtypes.RData')
q <- load('CCLE.RData')
r <- load('OxicBrSigs.RData')

#caseStudy2 data
s <- load('HypSigs.RData')
t <- load('TCGA.RData')

names(BrPatientSubtypes) <- c('TNBC','HER2','ERHigh','ERLow') #name the datasets

names(Oxic_BrSigs) <- c('Piening','Speers') #name the gene signatures

gene_sigs_list <- list() #create gene signatures list
gene_sigs_list[['Piening']] <- as.matrix(Oxic_BrSigs[['Piening']][,3])
gene_sigs_list[['Speers']] <- as.matrix(Oxic_BrSigs[['Speers']][,3])

data_list <- c(BrPatientSubtypes,list(CCLE=edata4)) #create the dataset list

#run sigQC
make_all_plots(gene_sigs_list = gene_sigs_list, mRNA_expr_matrix = data_list,showResults=F,out_dir='sigQC_case1_out',doNegativeControl=F )
 
#case study 2

names(g) <- c('Toustrup','Eustace','Lendahl') #name the gene signatures in their list

#run sigQC
make_all_plots(gene_sigs_list = g, mRNA_expr_matrix = list(TCGA=edata),showResults=F,out_dir='sigQC_case2_out',doNegativeControl=F )
