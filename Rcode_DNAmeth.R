# This module demonstrates the basic steps necessary to download raw DNA methylation array
# data generated using Illumina's Infinium methylation array 450k from the Gene Expression
# Omnibus public repository, preprocess the data, perform differential methylation
# analysis, and generate plots to visualize the results.

# First, install package BiocManager which enables installing Bioconductor packages 
rm(list=ls())
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install Bioconductor packages needed for this module
# Package GEOquery is used to download raw files from the GEO repository
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
# Package minfi is used to read raw methylation data into R, and perform preprocessing steps
if(!requireNamespace("minfi", quietly=TRUE)) BiocManager::install("minfi")
# The chip (450k) manifest and annotation packages are needed to annotate probe identifiers with genomic features
if(!requireNamespace("IlluminaHumanMethylation450kmanifest", quietly=TRUE)) BiocManager::install("IlluminaHumanMethylation450kmanifest")
if(!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly=TRUE)) BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# Package limma is used to perform a differential methylation test
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")

# Load the needed packages
library(GEOquery)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)

# 1) Prepare files.
#    Use package GEOquery to download raw files of series GSE47915 from the GEO repository. 
#    This series has 4 prostate tumor samples and 4 benign samples. 
#    Raw data files will be downloaded as a compressed folder under a newly created folder carrying the name of the GEO series.

#    Object downloaded will be a data frame with rownames representing the full path of the resulting downloaded files. 
#    The records in this data frame represent information for each downloaded file
downloaded <- getGEOSuppFiles("GSE47915")

# Set the working directory to the newly created folder
setwd(paste(getwd(), "/GSE47915", sep=""))

# Decompress the download tar folder GSE47915_RAW.tar
untar("GSE47915_RAW.tar")
# Delete the downloaded compressed files since they are no longer needed
file.remove("GSE47915_RAW.tar")
file.remove("GSE47915_signal_intensities.txt.gz")

# There are 2 raw files (with .idat extension) per sample: One for the red channel and one for the green channel
# Create the base names for raw files
idat.files <- list.files(pattern="idat")
idat.basenames <- sub(x=idat.files, pattern="_Grn.idat.gz", replacement="")
idat.basenames <- sub(x=idat.basenames, pattern="_Red.idat.gz", replacement="")
idat.basenames <- unique(idat.basenames)
# Create group labels for samples
group.labels <- c("benign", "tumor", "tumor", "tumor", "tumor", "benign", "benign", "benign")
idat.basenames
group.labels

# 2) Read raw data.
#    Use package minfi to read raw idat files into R
#    Proper annotation and chip manafest packages are needed

#    Read raw files into one RGChannelSet data object
RGset <- read.metharray(basenames=idat.basenames)
RGset

# 3) preprocess raw data:
#    Identify probes with poor quality (detection p-value >0.01) based on detection p-values
detP <- detectionP(RGset, type = "m+u")
poorProbes <- which(rowSums(detP > 0.01) > 0)
#    remove probes with poor quality
RGset <- RGset[-poorProbes,]
#    Preprocess raw intensities using the Illumina Genome Analyzer method (one option among others)
Mset <- preprocessIllumina(RGset)
Mset
#    Annotate probes in the methylation set data object
GMset <- mapToGenome(Mset)
GMset
#    Calculate the methylation metric (Beta values) from the methylated and unmethylated probe values
rc <- ratioConvert(GMset)
rc

# 4) Remove additional probes.
#    Remove CpG probes at genomic locations overlapping known SNPs with non-zero allele frequency.
#    Checking the methylation status at locations with known SNPs is not useful.
rc <- addSnpInfo(rc)
rc <- dropLociWithSnps(rc, snps=c("SBE","CpG"), maf=0)

# 5) Annotate probe sites
#    Get all available annotation information for probe sites
rc.ann <- getAnnotation(object=rc, what="everything", lociNames=NULL, orderByLocation=FALSE, dropNonMapping=TRUE)
#    There are too many columns in the Table but not all this information is needed
colnames(rc.ann)
#    Keep only UCSC genic region annotations
rc.ann.short <- rc.ann[,c("chr","pos","strand","UCSC_RefGene_Name","UCSC_RefGene_Group")]
head(rc.ann.short, n=6)

# 6) Perform differential methylation (DM) analysis using package limma
#    Get probe-level methylation ratio values (Beta)
probeBeta <- getBeta(rc)
design <- model.matrix(~group.labels)
fit <- lmFit(log2(1+probeBeta), design=design)
fit <- eBayes(fit)
limma.results <- topTable(fit, n=nrow(fit))
head(limma.results, n=6)

# 7) Extract an annotated results Table of significant probes only and save it
#    Probes with significant DM are those that meet assigned thresholds (small adjusted p-value and large absolute fold-change) 
sig <- which((limma.results$adj.P.Val < 0.001) & (abs(limma.results$logFC) > 0.4))
print(paste("There are ", length(sig), " significant probes", sep=""), quote=FALSE)
#    Extract significant probes from the limma results Table and add probe annotations
limma.sig.ann <- cbind(limma.results[sig,], rc.ann.short[rownames(limma.results[sig,]),])
#    Create a dataframe object and save results as a csv file 
df <- as.data.frame(limma.sig.ann)
write.csv(df, file="limma_significant_probes.csv")
head(df, n=4)

# 8) Find significant probes within 200 base pairs upstream from transcription start sites (TSS) of genes according to the used annotations
TSS200.probes <- df[grep("TSS200", df[, "UCSC_RefGene_Group"]),]
#    Identify one gene with most significant probes within 200 base pairs upstream from its TSS
names(which.max(table(TSS200.probes$UCSC_RefGene_Name)))
#    Find how many significant probes overlapping the previous gene's TSS200 region
max(table(TSS200.probes$UCSC_RefGene_Name))
#    Show DM probes overlapping the APC gene
APC.DM <- df[grep("APC", df$UCSC_RefGene_Name),]
APC.DM
#    Show DM probes overlapping the TSS200 region of the APC gene
APC.DM[grep("TSS200", APC.DM$UCSC_RefGene_Group),]

# 9) Results visualization.
#    Create a volcano plot to highlight significant probes
color <- rep("gray", nrow(limma.results))
#    highlight DM probes with red color
color[sig] <- "red"
plot(limma.results$logFC, -log10(limma.results$adj.P.Val), type="p", cex.lab=1.5, cex.main=1.5,
pch=16, col=color, xlab="log2FC", ylab="-log10(FDR)", main="Tumor vs. Benign")

# 10) Generate PCA plot using only DM probes
pca <- prcomp(t(log2(1+probeBeta)[rownames(limma.results[sig,]),]), scale=T)
a <- pca$sdev
pc1 <- a[1]^2/sum(a^2)
pc2 <- a[2]^2/sum(a^2)
colors <- ifelse(group.labels=="benign", "green", "red")
plot(pca$x[,1:2], pch=ifelse(group.labels=="benign", 16, 17), cex=2.5, cex.lab=1.5, cex.main=1.5, col=colors, main="2D PCA plot",
xlab=paste("PC1 (", round(100*pc1, 2), "%)", sep=""), 
ylab=paste("PC2 (", round(100*pc2, 2), "%)", sep=""))
legend("bottomleft", legend=c("Benign","Tumor"), col=c("green","red"), pch=c(16,17), cex=1.5, pt.cex=2)
