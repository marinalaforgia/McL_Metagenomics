library(tximportData)
library(tximport)

#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon

dir <- "/Users/Cassie/Desktop/Grass_cluster/salmon_gene"
list.files(dir)

samples <- read.table(file.path(dir, "GRASS_salmon_merge/file_names.txt"), header = TRUE)
samples

files <- file.path(dir, paste0(samples$FILE, ".quant"), "quant.sf")
names(files) <- paste0(samples$FILE)
all(file.exists(files))

#two columns transcript id and gene id - this is the same for us 
tx2gene <- read.delim2("combined_names.txt")

txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
# reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 
# summarizing abundance
# Error: vector memory exhausted (limit reached?)


library(wasabi)

dir <- "/Users/Cassie/Desktop/Grass_cluster/salmon_gene"
list.files(dir)

samples <- read.table(file.path(dir, "GRASS_salmon_merge/file_names.txt"), header = TRUE)
samples

files <- file.path(dir, paste0(samples$FILE, ".quant"))
names(files) <- paste0(samples$FILE)
all(file.exists(files))

#converts salmon files to kallisto format
prepare_fish_for_sleuth(files)

#try to convert kallisto format to a combined file
files_a5 <- file.path(dir, paste0(samples$FILE, ".quant"),"abundance.h5")
txi.kallisto <- tximport(files_a5, type = "kallisto", txOut = TRUE)

head(txi.kallisto$counts)
colnames(txi.kallisto$counts) <- samples$FILE
colnames(txi.kallisto$abundance) <- samples$FILE

saveRDS(txi.kallisto, "txi.kallisto.08042021.RDS")

txi.kallisto <- readRDS("txi.kallisto.08042021.RDS")

#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


colnames(txi.kallisto$counts) <- sub("_S.*", "", colnames(txi.kallisto$counts))
colnames(txi.kallisto$counts) <- sub("-", "_", colnames(txi.kallisto$counts))
colnames(txi.kallisto$abundance) <- sub("_S.*", "", colnames(txi.kallisto$abundance))
colnames(txi.kallisto$abundance) <- sub("-", "_", colnames(txi.kallisto$abundance))

metadata <- read.csv("McL_metag_metadata.csv")
metadata <- metadata[-c(1,3,5,6,8:10),] # only sequenced one water, one kit, and one mock
row.names(metadata) <- metadata$SampleID

library(DESeq2)

#need to make sampleTable (we probably have this) with metadata
#condition -> formula for deseq where forumula based on sampleTable e.g. watering_status + plant_comm_treatment + watering_status*plant_comm_treatment
ddsTxi <- DESeqDataSetFromTximport(txi.kallisto,
                                   colData = metadata,
                                   design = ~ New_Treatment)

#remove any rows with no counts
#dds <- ddsTxi[ rowSums(counts(ddsTxi)) > 10, ]
dds <- ddsTxi[ rowMeans(counts(ddsTxi)) > 1, ]

dds <- estimateSizeFactors(dds)

#run deseq diff abundance analysis
dds.test <- DESeq(dds)
saveRDS(dds.test, "GL_dds.deseq.obj.10142021.RDS")

#extract deseq normalized counts
normalized_counts <- counts(dds.test, normalized=TRUE)
saveRDS(normalized_counts, "GL_normcounts.10142021.RDS")


#get results
#this is where we would need contrasts / loops
#if we had more than two groups to compare
res <- results(dds.test) 

#get only results that are sig & have log fold change >2
res.sig <- subset(res, (padj < 0.01 & abs(log2FoldChange) >2))
