# This tutorial demonstrates the basic steps necessary to download raw DNA methylation
# data generated using Illumina's Infinium methylation array 450k from the Gene Expression
# Omnibus public repository, preprocess the data, perform differential methylation
# analysis, and generate plots to visualize the results.
# This tutorial assumed that R has been already installed.

# First, install package BiocManager which enables installing Bioconductor packages 
# if it has not been already installed 
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 1) Use package GEOquery to download raw files from Gene Expression Omnibus (GEO) repository
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
library(GEOquery)
# GEO series GSE47915 has 4 prostate tumor samples and 4 benign samples 
# Raw data will be downloaded as a compressed folder under a created folder carrying the name of the GEO series
getGEOSuppFiles("GSE47915")
# set the working directory to the newly created folder
setwd(paste(getwd(), "/GSE47915", sep=""))
# Decompress the download tar folder GSE47915_RAW.tar
untar("GSE47915_RAW.tar")
# delete the compressed file since no longer needed
file.remove("GSE47915_RAW.tar")

# There are 2 files per sample: One for the red channel and one for the green channel
# create the base names for raw files
idat.files <- list.files(pattern="idat")
idat.basenames <- sub(x=idat.files, pattern="_Grn.idat.gz", replacement="")
idat.basenames <- sub(x=idat.basenames, pattern="_Red.idat.gz", replacement="")
idat.basenames <- unique(idat.basenames)
# create group labels for samples
group.labels <- c("benign", "tumor", "tumor", "tumor", "tumor", "benign", "benign", "benign")

# 2) Use package minfi to read raw idat files into R
# Proper annotation and chip manafest packages are needed as well
if(!requireNamespace("minfi", quietly=TRUE)) BiocManager::install("minfi")
if(!requireNamespace("IlluminaHumanMethylation450kmanifest", quietly=TRUE)) BiocManager::install("IlluminaHumanMethylation450kmanifest")
if(!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly=TRUE)) BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Read raw files into one RGChannelSet object
RGset <- read.metharray(basenames=idat.basenames)
# Identify probes with poor quality based on detection p-values
detP <- detectionP(RGset, type = "m+u")
poorProbes <- which(rowSums(detP > 0.01) > 0)
# remove probes with poor quality
RGset <- RGset[-poorProbes,]
# Preprocess raw intensities using the Illumina Genome Analyzer method
Mset <- preprocessIllumina(RGset)
# Annotate probes
GMset <- mapToGenome(Mset)
# Calculate the methylation metric (Beta value)
rc <- ratioConvert(GMset)

# 3) Remove probes at genomic locations overlapping known SNPs with non-zero allele frequency
rc <- addSnpInfo(rc)
rc <- dropLociWithSnps(rc, snps=c("SBE","CpG"), maf=0)

# Get annotations for probe sites
rc.ann <- getAnnotation(object=rc, what="everything", lociNames=NULL, orderByLocation=FALSE, dropNonMapping=TRUE)
# Keep only UCSC genic annotations
rc.ann.short <- rc.ann[,c("chr","pos","strand","UCSC_RefGene_Name","UCSC_RefGene_Group")]

# 4) Perform differential methylation analysis using package limma
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
library(limma)
# get probe-level methylation ratio (Beta)
probeBeta <- getBeta(rc)
design <- model.matrix(~group.labels)
fit <- lmFit(log2(1+probeBeta), design=design)
fit <- eBayes(fit)
limma.results <- topTable(fit, n=nrow(fit))

# 5) Extract an annotated table of significant probes at stringent threshold levels
sig <- which((limma.results$adj.P.Val < 0.01) & (abs(limma.results$logFC) > 0.3))
# Extract significant probes from the limma results Table and add probe annotations
limma.sig.ann <- cbind(limma.results[sig,], rc.ann.short[rownames(limma.results[sig,]),])
# Create a dataframe object and save results 
df <- as.data.frame(limma.sig.ann)
write.csv(df, file="limma_significant_probes.csv")

# 6) Find significant probes within 200 base pairs upstream from transcription start sites (TSS)
TSS200.probes <- df[grep("TSS200", df[, "UCSC_RefGene_Group"]),]
# Identify the gene with most significant probes within 200 base pairs upstream from its TSS
names(which.max(table(TSS200.probes$UCSC_RefGene_Name)))
# Find how many significant probes overlapping the previous gene's TSS200 region
max(table(TSS200.probes$UCSC_RefGene_Name))
 
# 7) Create a volcano plot to highlight significant probes
color <- rep("gray", nrow(limma.results))
# highlight DE genes with red color
color[sig] <- "red"
plot(limma.results$logFC, -log10(limma.results$adj.P.Val), type="p", 
pch=16, col=color, xlab="log2FC", ylab="-log10(FDR)", main="Tumor vs. Benign")

# 8) Create PCA scatter plot using only the significant probes
pca <- prcomp(t(log2(1+probeBeta)[rownames(limma.results[sig,]),]), scale=T)
a <- pca$sdev
pc1 <- a[1]^2/sum(a^2)
pc2 <- a[2]^2/sum(a^2)
colors <- ifelse(group.labels=="benign", "green", "red")
plot(pca$x[,1:2], pch=16, cex=1.5, col=colors, main="2D PCA plot",
xlab=paste("PC1 (", round(100*pc1, 2), "%)", sep=""), 
ylab=paste("PC2 (", round(100*pc2, 2), "%)", sep=""))
legend("bottomleft", legend=c("Benign","Tumor"), col=c("green","red"), pch=c(16,17), pt.cex=1.5)


