#' Cancer Genomics workshop
#' ========================
#' 
#' Author: Tony Papenfuss
#' Date: Friday 1st May 2020
#'  
#' The code for this workshop is available here:
#' http://bioinf.wehi.edu.au/~papenfuss/cancer_genomics_workshop/workshop.R
#'
#' The html version of this document is available here: 
#' http://bioinf.wehi.edu.au/~papenfuss/cancer_genomics_workshop/workshop.html
#'
#'
#' In this practical class, we are going to perform copy number analysis on deep whole genome sequencing
#' from a patient with well-differentiated liposarcoma from this paper:
#' 
#' Garsed/Marshall/Corbin et al, The Architecture and Evolution of Cancer Neochromosomes, Cancer Cell 26:5, P653-667, 2014
#' https://www.cell.com/cancer-cell/fulltext/S1535-6108(14)00373-0
#' 
#' The data is tumour only sequencing, so we will use a simple approach to the copy number calling 
#' using the R/Bioconductor package QDNA-seq, 
#' 
#' https://www.bioconductor.org/packages/release/bioc/html/QDNAseq.html,
#' 
#' and explore the data in R and IGV. The raw data is publically available here:
#' 
#' https://www.ebi.ac.uk/ena/browser/view/PRJEB4696
#' 
#' First we install some dependencies.


# Download and install R packages
#install.packages("tidyverse")


# Download Bioconductor packages

# Install scripts for older versions of Bioconductor
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("QDNAseq")

# For Bioconductor 3.10
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10") # The current version (1/5/2020) is 3.11, but requires R 4.0.
BiocManager::install(c("QDNAseq"))


# Load R libraries
library(QDNAseq)
library(tidyverse)


#' By editing and uncommenting the following line you can change 
#' to whatever working directory you like.
#' 
setwd("~/Desktop/cancer_genomics_workshop")

#' I have already aligned the reads to the hg19 reference genome and 
#' summarised the reads into counts.
#' 
#' QDNA-seq uses sets of pre-defined bins of various sizes with precalculated 
#' GC content and mappability. These were selected using the following commands.

#  bins <- getBinAnnotations(binSize = 30, type = "SR50")
#  saveRDS(bins, file = "bins30K.dat")

#' binSize = 30 means that 30kb bins are used.
#' type = "SR50" means that mappability is defined for 50nt single end reads. 
#' This is not quite right, as the data was 100nt paired end. Generating 
#' this data for 100nt long reads is very slow, but would improve the analysis 
#' a bit.
#' 
#' First we load the saved version of this.

#' If you have downloaded all the data or created your own bin30K.dat file, you can load the
#' bin data from a local file:

# bins.filename = "file://bins30K.dat"

bins.filename = "http://bioinf.wehi.edu.au/~papenfuss/cancer_genomics_workshop/bins30K.dat"
bins = readRDS(gzcon(url(bins.filename)))

#' Binning reads into these bins is also pretty slow, so I have already done this.
#' It was done using the following commands.
#' 
#  filename = "/Users/papenfuss/projects/liposarcoma/data/ST059/ST059_NoIndex_bt2.bam"
#  readCounts = binReadCounts(bins, bamfiles = filename)
#  saveRDS(readCounts, file = "scratch/ST059_counts_30kb_bins.dat")
# 
#' We just load it from a binary dat file.

counts.filename = "http://bioinf.wehi.edu.au/~papenfuss/cancer_genomics_workshop/ST059_counts_30kb_bins.dat"
# counts.filename = "file://ST059_counts_30kb_bins.dat"
readCounts = readRDS(gzcon(url(counts.filename)))

#' To see what these look like, we can use the R command head

head(bins)
head(readCounts)

#' These are R objects. Data sits in slots (featureData, assayData, etc)
#' and can be accessed using the @ operator and then the $ operator.

head(readCounts@assayData$counts)

#' Here's another amazing copy number profile if you want to look.
#' 
#  readCounts <- binReadCounts(bins, bamfiles = filename.all))
#  saveRDS(readCounts, file = "scratch/T1000_counts_30kb_bins.dat")
# 
# readCounts = readRDS("scratch/T1000_counts_30kb_bins.dat")

#' Let's have a look at the GC bias in the count data.

GC.data = data.frame(
  gc = readCounts@featureData@data$gc,
  counts = readCounts@assayData$counts[, 1]
)
plot(counts~gc, GC.data, pch=20, cex=0.4, col="#00000044", log="y", main="GC bias")

#' This looks quite flat (excellent), but does show some banding which corresponds to 
#' different copy number states.
#'
#' Now lets have a look at the count data across the genome using the default
#' QDNA-seq plots.

plot(readCounts, logTransform=FALSE, main="Read counts")
highlightFilters(
  readCounts,
  logTransform = FALSE,
  residual = TRUE,
  blacklist = TRUE
)

#' The default plotting style is OK, but squished and difficult to read. 
#' It mainly hightlights high levels of amplification on chr12 and 4. What else can you make out?
#' Let's try to make something nicer using ggplot.

dat = data.frame(
  chrom = readCounts@featureData@data$chromosome,
  Position = readCounts@featureData@data$start/1e6,
  Counts = readCounts@assayData$counts[, 1]
)
ggplot(dat, aes(x=Position, y=Counts)) + 
  geom_point(size=0.5, alpha = 0.2) + 
  facet_wrap(~chrom, scales="free_x") + 
  scale_y_log10() +
  ggtitle("Read counts")

#' Noice, but the chromosomes are out of order (doh). You can fix this by making chrom an ordered factor.

chroms = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
dat$chrom = factor(dat$chrom, levels=chroms, ordered=TRUE)

ggplot(dat, aes(x=Position, y=Counts)) + 
  geom_point(size=0.5, alpha = 0.2) + 
  facet_wrap(~chrom, scales="free_x") + 
  scale_y_log10() +
  ggtitle("Read counts")
             
#' Coolio. What do you think is going on with he genome?
#'         
#' QDNA-seq allows us to mask outliers and regions in dodgy parts of the genome.
#' These are called blacklisted regions. We can also smooth and segment the data.

readCountsFiltered <-
  applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)

copyNumbers <- correctBins(readCountsFiltered)
plot(copyNumbers, main="copyNumbers")

copyNumbersNormalized <- normalizeBins(copyNumbers)
plot(copyNumbersNormalized, main="copyNumbersNormalized")

copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
plot(copyNumbersSmooth, main="copyNumbersSmooth")

#' We'll export these results to look at in IGV if time permits.

exportBins(copyNumbers, file = "ST059", format = "bed", type="copynumber")
exportBins(copyNumbers, file = "ST059.igv", format = "igv", type="copynumber")
exportBins(copyNumbers, file = "ST059.txt", type="copynumber")

#' We can also segment the data and call copy gain and loss.

copyNumbersSegmented <-
  segmentBins(copyNumbersSmooth, transformFun = "sqrt")
copyNumbersSegmented <-
  normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented, main="copyNumbersSegmented")

copyNumbersCalled <- callBins(copyNumbersSegmented)
plot(copyNumbersCalled, main="copyNumbersCalled")

cgh <- makeCgh(copyNumbersCalled)

#' Now let's take a closer look at chr12.

dat = data.frame(
  chrom = copyNumbers@featureData@data$chromosome,
  Position = copyNumbers@featureData@data$start/1e6,
  CopyNumber = copyNumbers@assayData$copynumber[, 1]
)

dat1 = subset(dat, chrom=="12")
ggplot(dat1, aes(x=Position, y=CopyNumber)) +
  geom_point(size=0.5, alpha = 0.2) + 
  scale_y_log10() +
  ylim(0.1, 75) + 
  ggtitle("Chr12")

#' ...and zooming in to some genes of interest

dat = data.frame(
  chrom = copyNumbersSmooth@featureData@data$chromosome,
  Position = copyNumbersSmooth@featureData@data$start/1e6,
  CopyNumber = copyNumbersSmooth@assayData$copynumber[, 1]
)
dat2 = subset(dat, chrom=="12")
ggplot(dat2, aes(x=Position, y=CopyNumber)) +
  geom_line(size=0.5, alpha = 0.2) + 
  scale_y_log10() +
  xlim(65, 75) +
  ylim(0.1, 75) + 
  ggtitle("Chr12")


#' The oncogenes MDM2, CDK4 and HMGA2 are located on chr12 and
#' frequently amplified in well-differentiated liposarcomas. 
#' The easiest way to check if any of these oncogenes are amplified 
#' in this tumour is to use IGV. If time permits, load the file ST059.igv 
#' into IGV and have a look at these genes.
