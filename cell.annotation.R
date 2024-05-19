setwd("~/xwj/T4/zs")
setwd("~")
setwd('xwj/T4/zs/')
library(Seurat)
library(Matrix)
### 1. Preparation ----
# Load unwanted genes, which are genes being induced during protoplasting in our case
pp.genes <- as.character(read.table("Protoplasting_DEgene_FC2_list.txt", header=F)$V1)

# load Data
load("jean.RData")
jean<-load("jean.RData")

head(ken)   
dim(ken$count)  #37336 x 2085

# Obtain the name of the cell whose treatment_id value is Control from meta
wt.cell <- rownames(ken$meta)[ken$meta$treatment_id == "Control"]

count <- ken$count[,wt.cell]
dim(count)  #37336  1076

#change column name of the matrix
colnames(count) <- gsub("-1_agg_HS_analysis","-2",colnames(count))

#load data
load("Shulse.RData")
load("Ryu3.RData")
load("Zheng5.RData")

colnames(Zheng) <- paste0(colnames(Zheng),"_5")

rownames(Zheng)
colnames(Zheng)
#dim(Zheng) #[1] 32833  7695

tariAnno <- merge(ensgene1,test1,by="ensembl_gene_id",all=TRUE)#####缺少ensgene1,test1

tariAnno$symbol[is.na(tariAnno$symbol)] <- tariAnno$external_gene_name[is.na(tariAnno$symbol)]

tariAnno <- tariAnno[,-3]

write.csv(tariAnno,file="tariAnno.csv")
#tariAnno
#dim(tariAnno)#[1] 27701     3


#Return the position of the gene name that matches tariAnno in the dataset Zheng (in tariAnno),
#Replace with NA_integer_ for unmatched ones
nm <- match(rownames(Zheng), tariAnno$symbol, nomatch = NA_integer_, incomparables = NULL)

#Take the gene from tariAnno that matches Zheng
ngn <- tariAnno[nm,]$ensembl_gene_id
#length(ngn) [1] 32833

#For positions where the median of ngn is NA, use the value of the corresponding position in Zheng's row name to represent it
ngn[is.na(ngn)] <- rownames(Zheng)[is.na(ngn)]

#Searching for gene IDs that match AT [0-9] G in ngn
ss <- grepl("AT[0-9]G",ngn)
#length(ss)[1] 32833

#take subset
Zheng <- Zheng[ss,]
#dim(Zheng)[1] 30544  7695


#GeneM belongs to Ryu3. RData
#In the Seurat object, counts stores the raw data and is a sparse matrix
counts<- geneM@assays$RNA@counts 
#dim(counts)[1] 21087  7522 


colnames(counts) <- gsub("(_3$)","_3_3",colnames(counts)) #Replace column names
counts3 <- counts[,grep("_3_3",colnames(counts))]


#samples
use.sample <- c("sc_12","sc_11","sc_30","sc_31","sc_37","sc_40","sc_51","sc_9_at","sc_10_at","sc_1","tnw1","tnw2","col0")
use.sample <- c("sc_20","sc_21","sc_25","sc_36","sc_52","sc_53")

#Process each sample and create a Seurat object
for (i in 1:length(use.sample)) {
  cat("---",i,"---")
  sample.name <- use.sample[i]
  #sample.name<-"sc_12"
  dir_splice <- paste0('counts/', sample.name,"/",sample.name,'_mtx/spliced_counts_filtered/')
  dir_unsplice <- paste0('counts/', sample.name,"/",sample.name,'_mtx/unspliced_counts_filtered/')
  # Read in the quality-filtered spliced and unspliced counts matrices
  #Read the quality filtered concatenated and non concatenated count matrices
  #So, the concatenated and un concatenated counting matrices have already been processed
  spliced <- readMM(paste0(dir_splice,"matrix.mtx")) # load raw mtx
  genes <- read.csv(paste0(dir_splice,"genes.tsv"), sep = '\t', header = F) # load genes
  rownames(spliced) <- genes[,1] # attach gene_ids
  colnames(spliced) <- read.csv(paste0(dir_splice,"barcodes.tsv"), sep = '\t', header = F)[,1] # attach barcodes
  
  unspliced <- readMM(paste0(dir_unsplice,"matrix.mtx")) # load raw mtx
  genes <- read.csv(paste0(dir_unsplice,"genes.tsv"), sep = '\t', header = F) # load genes
  rownames(unspliced) <- genes[,1] # attach gene_ids
  colnames(unspliced) <- read.csv(paste0(dir_unsplice,"barcodes.tsv"), sep = '\t', header = F)[,1] # attach barcodes
  # Combined spliced and unspliced transcripts
  combined <- spliced + unspliced
  # Create Seurat Object
  #geneM <- suppressWarnings(CreateSeuratObject(counts = combined, assay = "RNA", project = sample.name))
  geneMSeuratObject <- suppressWarnings(CreateSeuratObject(counts = combined, assay = "RNA", project = sample.name))
  #geneM[["spliced_RNA"]] <- CreateAssayObject(spliced)
  #geneM[["unspliced_RNA"]] <- CreateAssayObject(unspliced)
  DefaultAssay(geneMSeuratObject) <- "RNA"
}


# Calculate the correlation coefficient between each cell and each reference expression profile,
# And label the cell with the highest correlation coefficient
celltypeTest <- function(celltype,celltype_label,method){
  
  ncell <- length(celltype_label)
  celltype_stat <- suppressWarnings(sapply((ncell+1):ncol(celltype), function(i) sapply(1:ncell, function(j) cor.test(celltype[,i],celltype[,j],method = method)[c(3,4)])))
  celltype_cor <- celltype_stat[seq(2,nrow(celltype_stat),2),]
  celltype_pvalue <- celltype_stat[seq(1,nrow(celltype_stat)-1,2),]
  celltype_max <- sapply(1:(ncol(celltype)-ncell), function(i) max(as.numeric(celltype_cor[,i])))
  celltype_ident <- sapply(1:(ncol(celltype)-ncell), function(i) celltype_label[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  celltype_maxp <- sapply(1:(ncol(celltype)-ncell), function(i) as.numeric(celltype_pvalue[,i])[which(as.numeric(celltype_cor[,i])==max(as.numeric(celltype_cor[,i])))])
  names(celltype_max) <- celltype_ident
  result <- list()
  result$celltype_max <- celltype_max
  result$celltype_maxp <- celltype_maxp
  return(result)
}

# Create object
# Handle Shulse Content in the RData dataset
# Contains counting matrix count [10], and the original s data meta contains cell description information:
for(i in 1:length(Shulse$count)){
  cat(i)
  count <- Shulse$count[[i]]
  #count <- Shulse$count[[1]]
  counts <- as.data.frame(count)
  rownames(counts) <- counts[,1]
  counts <- counts[,-1]
  
  colnames(counts) <- paste0(colnames(counts),"-4-",i)
  
  geneM <- CreateSeuratObject(counts = Zheng, project = "Gene", min.cells = 3, min.features = 200)
  cat("--",dim(geneM),"--\n")
  # As an alternative to NormalizeData, FindVariableFeatures, and ScaleData workflows
  geneM <- SCTransform(geneM, variable.features.n = nrow(geneM), assay = "RNA", new.assay.name = "SCT", verbose = FALSE)
  
  # Prepare list of genes to be use in PCA, UMAP ... etc. Here we use "all the genes/features" except for mitochondrial genes, chloroplast genes and genes involved in protoplasting 
  use.genes <- rownames(geneM)[-c(grep(paste("ATMG", collapse = "|"),rownames(geneM)),grep(paste("ATCG", collapse = "|"),rownames(geneM)),sort(match(pp.genes, rownames(geneM))))]
  
  ### 2. Correlation-based annotation ----
  
  ## 2.1 cell annotation ----
  # Load reference expression profiles
  load(file="Root_bulk_arabidopsis_curated.RD") # time celltype
  #' time: Reference expression profiles for time
  #' celltype: Reference expression profiles for cell type
  #' test1<-load(file="Root_bulk_arabidopsis_curated.RD")
  #' "time"           "celltype"       "Long"           "Rad"            "merge.rownames"
  
  # Extract matrix of SCTransformed expression value 
  rc <- as.matrix(geneM@assays$SCT@data)
  # Merge the reference expression profile with the normalized expression matrix of our sample  
  
  merge.rownames <- function (x,y){
    dat <- merge(x = x, y = y, by = "row.names")
    rownames(dat) <- dat$Row.names
    dat <- dat[,-1]
    return(dat)
  }
  
  time <- Reduce(merge.rownames, list(time,rc))
  celltype <- Reduce(merge.rownames, list(celltype,rc))
  # Prepare customized label name (optional)
  time_label=c("Elongation", "Maturation", "Meristem")
  celltype_label=c("phloem & companion cells", "developing cortex", "hair cells", "matured cortex",
                   "matured endodermis", "non-hair cells", "columella", "phloem pole pericycle",
                   "matured xylem pole", "protophloem & metaphloem","developing xylem", "endodermis & QC cells", "LRC & non-hair cells","QC cells")
  
  
  cellanno <- celltypeTest(celltype,celltype_label,method = "pearson")
  
  timeanno <- celltypeTest(time,time_label,method = "pearson")
  
  # Store the annotation, correlation coefficient and the p-value in Seurat object
  geneM@meta.data$celltype.ID.P <- as.character(names(cellanno$celltype_max))
  geneM@meta.data$timezone.ID.P <- as.character(names(timeanno$celltype_max))
  geneM@meta.data$celltype.cor.P <- cellanno$celltype_max
  geneM@meta.data$timezone.cor.P <- timeanno$celltype_max
  geneM@meta.data$celltype.pvalue.P <- cellanno$celltype_maxp
  geneM@meta.data$timezone.pvalue.P <- timeanno$celltype_maxp
  
  # In case there is cell with insufficient information for annotation, label them as "unknown"
  geneM@meta.data$celltype.ID.P[which(geneM@meta.data$celltype.ID.P=='character(0)')]="unknown"
  geneM@meta.data$timezone.ID.P[which(geneM@meta.data$timezone.ID.P=='character(0)')]="unknown"
  
  
  saveRDS(geneM, file = paste0("Zhang.rds"))
  
}

#The Seurat objects created in these two for loops are saved in the/gene directory

## 2.2 cluster annotation (optional)----

# Run PCA
geneM <- RunPCA(geneM, verbose = FALSE, approx = FALSE, npcs = 50, features=use.genes)
# Run UMAP
suppressMessages(suppressWarnings(
  geneM <- RunUMAP(geneM, reduction = "pca", dims = 1:50, umap.method = "umap-learn", metric = "correlation")
))
# Find nearest neighbors
suppressMessages(suppressWarnings(
  geneM <- FindNeighbors(geneM, reduction = "pca",dims = 1:50)
))
# Find clusters, here we choose Leiden clustering algorithm with resolution 0.5. Parameter "algorithm": 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm
suppressMessages(suppressWarnings(
  geneM <- FindClusters(geneM, resolution = 0.5, algorithm = 1)
))

DimPlot(geneM,reduction = "umap")

# Extract integrated (batch-corrected and scaled) expression matrix
afm <- as.matrix(geneM@assays$SCT@data)

# Pool (average) expression values of each cluster
pooled <- matrix(nrow=nrow(afm), ncol = 0)
for (i in 0:(length(unique(geneM@meta.data$seurat_clusters))-1)) {
  m <- afm[,which(geneM@meta.data$seurat_clusters==i)]
  pooled <- cbind(pooled, rowSums(m)/ncol(m))
}
colnames(pooled) <- 0:(ncol(pooled)-1)
# Load reference expression profiles
load(file="Root_bulk_arabidopsis_curated.RD")

# Annation (see notebook "1-Correlation-Based-Annotation.ipynb")
merge.rownames <- function (x,y){
  dat <- merge(x = x, y = y, by = "row.names")
  rownames(dat) <- dat$Row.names
  dat <- dat[,-1]
  return(dat)
}

time <- Reduce(merge.rownames, list(time,pooled))
celltype <- Reduce(merge.rownames, list(celltype,pooled))

time_label=c("Elongation", "Maturation", "Meristem")
celltype_label=c("phloem & companion cells", "developing cortex", "hair cells", "matured cortex",
                 "matured endodermis", "non-hair cells", "columella", "phloem pole pericycle", 
                 "matured xylem pole", "protophloem & metaphloem","developing xylem", "endodermis & QC cells", "LRC & non-hair cells","QC cells")
cellanno <- celltypeTest(celltype,celltype_label,method = "spearman")
geneM@meta.data$seucluster.celltype.cor <- as.numeric(FoldNames(cellanno$celltype_max,geneM))
geneM@meta.data$seucluster.celltype.ID <- FoldNames(names(cellanno$celltype_max),geneM)
geneM@meta.data$seucluster.celltype.P <- as.numeric(FoldNames(cellanno$celltype_maxp,geneM))


timeanno <- celltypeTest(time,time_label,method = "spearman")
geneM@meta.data$seucluster.time.cor <- as.numeric(FoldNames(timeanno$celltype_max,geneM))
geneM@meta.data$seucluster.time.ID <- FoldNames(names(timeanno$celltype_max),geneM)
geneM@meta.data$seucluster.time.P <- as.numeric(FoldNames(timeanno$celltype_maxp,geneM))


Long <- Reduce(merge.rownames, list(Long,pooled))
Rad <- Reduce(merge.rownames, list(Rad,pooled)) 

Long_label=c("Columella", "Meri-1", "Meri-2", "Meri-3", "Meri-4", "Meri-5", "Meri-6", "Elong-7", "Elong-8", "Mat-9", "Mat-10", "Mat-11", "Mat-12")
Rad_label=c("QC", "Hair Cell", "Cortex", "Non-Hair Cell", "Xylem Pole Pericycle", "LRC", 
                 "Columella", "Phloem Pole Pericycle", "Mat.Xylem", "Meri.Xylem", "Phloem [S32]", "Endodermis", "Phloem [SUC2]")

Longanno <- celltypeTest(Long,Long_label,method = "spearman")
geneM@meta.data$seucluster.Long.cor <- as.numeric(FoldNames(Longanno$celltype_max,geneM))
geneM@meta.data$seucluster.Long.ID <- FoldNames(names(Longanno$celltype_max),geneM)
geneM@meta.data$seucluster.Long.P <- as.numeric(FoldNames(Longanno$celltype_maxp,geneM))


Radanno <- celltypeTest(Rad,Rad_label,method = "spearman")
geneM@meta.data$seucluster.Rad.cor <- as.numeric(FoldNames(Radanno$celltype_max,geneM))
geneM@meta.data$seucluster.Rad.ID <- FoldNames(names(Radanno$celltype_max),geneM)
geneM@meta.data$seucluster.Rad.P <- as.numeric(FoldNames(Radanno$celltype_maxp,geneM))

#/gene
saveRDS(geneM, file = "ken2.rds")



##2.3 bulk annotation ----
rc <- as.matrix(gene_har@assays$RNA@data)


celltype <- Reduce(merge.rownames, list(Rad,rc[,1:10000]))

time_label=c("Columella", "Meri-1", "Meri-2", "Meri-3", "Meri-4", "Meri-5", "Meri-6", "Elong-7", "Elong-8", "Mat-9", "Mat-10", "Mat-11", "Mat-12")


timeanno <- list()

time <- Reduce(merge.rownames, list(Long,rc[,1:30000]))
timeanno[[1]] <- celltypeTest(time,time_label,method = "pearson")

time <- Reduce(merge.rownames, list(Long,rc[,30001:80001]))
timeanno[[2]] <- celltypeTest(time,time_label,method = "pearson")

time <- Reduce(merge.rownames, list(Long,rc[,80002:128549]))
timeanno[[3]] <- celltypeTest(time,time_label,method = "pearson")

gene_har$Long.cor.P<- c(timeanno[[1]]$celltype_max,timeanno[[2]]$celltype_max,timeanno[[3]]$celltype_max)
gene_har$Long.ID.P<- c(names(timeanno[[1]]$celltype_max),names(timeanno[[2]]$celltype_max),names(timeanno[[3]]$celltype_max))
gene_har$Long.pvalue.P<- c(timeanno[[1]]$celltype_maxp,timeanno[[2]]$celltype_maxp,timeanno[[3]]$celltype_maxp)

saveRDS(gene_har@meta.data,"gene_meta.rds")

celltype_label=c("QC", "Hair Cell", "Cortex", "Non-Hair Cell", "Xylem Pole Pericycle", "LRC", 
                 "Columella", "Phloem Pole Pericycle", "Mat.Xylem", "Meri.Xylem", "Phloem [S32]", "Endodermis", "Phloem [SUC2]")

celltype_cor <- sapply(14:ncol(celltype), function(i) sapply(1:13, function(j) cor(celltype[,i],celltype[,j],method = "pearson")))
rownames(celltype_cor) <- celltype_label
celltype_max <- sapply(1:(ncol(celltype)-13), function(i) max(celltype_cor[,i]))
celltype_ident <- sapply(1:(ncol(celltype)-13), function(i) celltype_label[which(celltype_cor[,i]==max(celltype_cor[,i]))])
names(celltype_max) <- celltype_ident


seucluster.Rad.ID.P1 <- celltype_ident
seucluster.Rad.cor.P1 <- celltype_max

seucluster.Rad.ID.P <- c(seucluster.Rad.ID.P1,seucluster.Rad.ID.P2,seucluster.Rad.ID.P3,seucluster.Rad.ID.P4,seucluster.Rad.ID.P5,seucluster.Rad.ID.P6)
seucluster.Rad.cor.P <- c(seucluster.Rad.cor.P1,seucluster.Rad.cor.P2,seucluster.Rad.cor.P3,seucluster.Rad.cor.P4,seucluster.Rad.cor.P5,seucluster.Rad.cor.P6)


gene_har@meta.data$Rad.ID.P <- as.character(seucluster.Rad.ID.P)
gene_har@meta.data$Rad.cor.P <- as.numeric(seucluster.Rad.cor.P)




#gene_har=geneM


##2.4 Final plot----
load("~/xwj/opt/Data/T4/zs/COPILOT-master/supp_data/color_scheme_at.RData")
#[1] "celltypeorder"       "celltypeorder_ici"   "celltypepalette"     "celltypepalette_ici" "longorder"  "longpalette"
# Plot correlation-based annotation
color <- celltypepalette[sort(match(unique(gene_har$celltype.ID.P),celltypeorder))]
#DimPlot(gene_har, reduction = "umap", group.by = "celltype.ID.P", cols = color)+ggtitle("RNA-seq correlation annotation")
gene_har$celltype.ID.P_significance <- rep("> 0.01", nrow(gene_har@meta.data))
gene_har$celltype.ID.P_significance[which(gene_har$celltype.pvalue.P <= 0.01)] = "<= 0.01"
gene_har$celltype.cor.P_significance <- rep(">= 0.6", nrow(gene_har@meta.data))
gene_har$celltype.cor.P_significance[which(gene_har$celltype.cor.P < 0.6)] = "< 0.6"
DimPlot(gene_har, reduction = "combat_umap", group.by = "celltype.ID.P_significance", order = c("<= 0.01","> 0.01"),cols = c("#cccccc", "#ff4040"))
FeaturePlot(gene_har, reduction = "combat_umap", features = "celltype.cor.P")
DimPlot(gene_har, reduction = "combat_umap", group.by = "celltype.cor.P_significance", order = c(">= 0.6","< 0.6"),cols = c("#cccccc", "#008081"))
gene_har$Rad.ID.P <- factor(gene_har$Rad.ID.P, levels = radorder[sort(match(unique(gene_har$Rad.ID.P),radorder))])
color <- radpalette[sort(match(unique(gene_har$Rad.ID.P),radorder))]
DimPlot(gene_har, reduction = "combat_umap", group.by = "Rad.ID.P", cols = color)+ggtitle("Microarray correlation annotation")
gene_har$Rad.ID.P_significance <- rep("> 0.01", nrow(gene_har@meta.data))
gene_har$Rad.ID.P_significance[which(gene_har$Rad.pvalue.P <= 0.01)] = "<= 0.01"
gene_har$Rad.cor.P_significance <- rep(">= 0.6", nrow(gene_har@meta.data))
gene_har$Rad.cor.P_significance[which(gene_har$Rad.cor.P < 0.6)] = "< 0.6"
DimPlot(gene_har, reduction = "umap", group.by = "Rad.ID.P_significance", order = c("<= 0.01","> 0.01"),cols = c("#cccccc", "#ff4040"))
FeaturePlot(gene_har, reduction = "umap", features = "Rad.cor.P")
DimPlot(gene_har, reduction = "umap", group.by = "Rad.cor.P_significance", order = c(">= 0.6","< 0.6"),cols = c("#cccccc", "#008081"))




### 3. Integration ----

library(dplyr)
library(Seurat)
# Prepare a list of sample names that is going to be used for integration
use.sample <- c("sc_12","sc_11","sc_30","sc_31","sc_37","sc_40","sc_51","sc_9_at","sc_10_at","sc_1","tnw1","tnw2","col0")
# Read in Seurat objects and make a list out of them, here we put all the Seurat objects under the foler ./COPILOT_RDS 

read_seu <- function(dir,sample.name) { 
  seu <- readRDS(dir)
  
  # remove unused data to save some memory (optional)
  if(gsub(".rds","",sample.name) %in% use.sample){
    seu@assays$spliced_RNA <- NULL
    seu@assays$spliced_SCT <- NULL
    seu@assays$unspliced_RNA <- NULL
    seu@assays$unspliced_SCT <- NULL
  }
  return(seu)
  
}

list.filenames <- list.files(path = "gene/",pattern=".rds$")
list.filenames[1] <- "sc_12.rds"
list.filenames[9] <- "col0.rds"

rc.list <- list()

for (i in 1:length(list.filenames)){
  cat(i,"/",length(list.filenames),"\n")
  rc.list[[i]]<- read_seu(dir = paste0("gene/",list.filenames[i]),list.filenames[i])
}

names(rc.list) <- list.filenames %>% gsub(".rds","",.)

for (i in use.sample) {
  cat(i)
  rc.list[[i]] <- RenameCells(object = rc.list[[i]] , add.cell.id = i)
}
saveRDS(rc.list,file="rclist.rds")
rc.list <- readRDS("relist.rds")

# Select genes shared among samples that will be used for integration, for Arabidopsis, 25000 genes is about the upper limit of genes a single cell can have
subrs <- rc.list[1:2]
rc.features <- SelectIntegrationFeatures(object.list = subrs)
length(rc.features)
#12191
# Remove mitochondrial, chloroplast and protoplasting-induced genes from the shared gene list
rc.features <- rc.features[-c(grep("ATMG",rc.features),grep("ATCG",rc.features),sort(match(pp.genes, rc.features)))]
length(rc.features)
#10920
# Prepare for integration
subrs <- lapply(X = subrs, FUN = function(x) {
  cat("--")
  #DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
  x <- ScaleData(x, features = rc.features, verbose = FALSE)
  x <- RunPCA(x, features = rc.features,verbose = FALSE)
  #x@assays$SCT <- NULL
})


anchors <- FindIntegrationAnchors(object.list = subrs, reference = c("1"),reduction = "rpca",dims = 1:50)


rc.listp <- PrepSCTIntegration(object.list =subrs, anchor.features = rc.features, verbose = TRUE)

#reference_dataset <- which(names(rc.listp) == "sc_12")
reference_dataset <- "1"

anchors <- FindIntegrationAnchors(object.list = rc.listp,normalization.method = "SCT", anchor.features = rc.features)

integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT",verbose = FALSE)


anchors <- FindIntegrationAnchors(object.list = rc.listp, reference = reference_dataset,
                                  normalization.method = "SCT", reduction = "rpca", anchor.features = rc.features,dims = 1:20)

anchors <- FindIntegrationAnchors(object.list = rc.listp, normalization.method = "SCT", reference = reference_dataset,
                                           anchor.features = rc.features)

anchors <- FindIntegrationAnchors(object.list = subr, reference = reference_dataset,
                                  reduction = "rpca",dims = 1:20)

saveRDS(anchors,file="anchors.rds")
anchors <- readRDS("anchors.rds")



########## 3.2 harmony + Combat-------

genelist <- list()
for (i in 1:length(rc.list)) {
  genelist[[i]] <- rownames(rc.list[[i]])
}

gene_inter <- Reduce(intersect,genelist)
rm(genelist)
gc()

dim(rc.list[[2]])
ncg <- 0
for (i in 1:length(rc.list)) {
  ncg = ncg + dim(rc.list[[i]])[2]
}
gcounts <- matrix(nrow = length(gene_inter),ncol = ncg)
start <- 1
colns <- c()
for (i in 1:length(rc.list)) {
  end <- start + dim(rc.list[[i]])[2] -1
  xs <- GetAssayData(rc.list[[i]], slot = "counts", assay = "SCT")
  xs <- as.matrix(xs[gene_inter,])
  gcounts[,start:end]<- xs
  colns <-c(colns,colnames(rc.list[[i]]))
  start <- end +1
  cat(i,"/28:",start,":",end,"\n")
  gc()
}
rownames(gcounts)<- gene_inter
colnames(gcounts)<- colns
saveRDS(gcounts,file="gcounts.rds")
gcounts<-readRDS('gcounts.rds')
merger <- function(x,y){
  
  if(class(x) != "matrix"){
    xs <- GetAssayData(x, slot = "counts", assay = "SCT")
    xs <- xs[gene_inter,]
  }else{
    xs <- x
  }
  xs <- 
  
  cat("x:",dim(xs),"\t")
  ys <- GetAssayData(y, slot = "counts", assay = "SCT")
  ys <- ys[gene_inter,]
  cat("y:",dim(ys),"\n")
  dat <- cbind(as.matrix(xs),as.matrix(ys))
  return(dat)
}

merger.pca <- function(x,y){
  
  if(class(x) != "matrix"){
    
    pcax <- x@reductions[["pca"]]@cell.embeddings
    xs <- t(pcax)
  }else{
    xs <- x
  }
  
  cat("x:",dim(xs),"\t")
  pcay <- y@reductions[["pca"]]@cell.embeddings
  
  ys <- t(pcay)
  cat("y:",dim(ys),"\n")
  dat <- cbind(as.matrix(xs),as.matrix(ys))
  return(dat)
}

#gcounts <- Reduce(merger,rc.list)

gcounts<-readRDS('gcounts.rds')

gene_har <- CreateSeuratObject(counts = gcounts, project = "genes",min.cells = 5,min.features = 500)
#gene_har <- NormalizeData(gene_har, verbose = FALSE)
gene_har <- FindVariableFeatures(gene_har, nfeatures = nrow(gene_har), verbose = FALSE)
gene_har <- ScaleData(gene_har, verbose = FALSE,vars.to.regress = "sample")

gene_har<-readRDS('gene_har.rds')

samples.all <- c("sc_12_","sc_11_","sc_30_","sc_31_","sc_37_","sc_40_",
                 "sc_51_","sc_9_at_","sc_10_at_","sc_1_","tnw1_","tnw2_",
                 "col0_","[A-Z]-2","[A-Z]_3_1","[A-Z]_3_2","[A-Z]_3_3","[A-Z]-4-1","[A-Z]-4-2",
                 "[A-Z]-4-3","[A-Z]-4-4","[A-Z]-4-5","[A-Z]-4-6","[A-Z]-4-7","[A-Z]-4-8","[A-Z]-4-9","[A-Z]-4-10","[A-Z]_5")

samples3 <- c("6_sc12","6_sc11","6_sc30","6_sc31","6_sc37","6_sc40",
              "6_sc51","6_sc9","6_sc10","6_sc1","6_tnw1","6_tnw2",
              "6_col0","2_1","3_1","3_2","3_3","4_A","4_B","4_C",
              "4_D","4_E","4_F","4_G","4_H","4_I","4_J","5_1")


samples4 <- c("Shahan et al.","Shahan et al.","Shahan et al.","Shahan et al.","Shahan et al.","Shahan et al.",
              "Shahan et al.","Shahan et al.","Shahan et al.","Shahan et al.","Shahan et al.","Shahan et al.",
              "Shahan et al.","Jean-Baptiste et al.","Ryu et al.","Ryu et al.","Ryu et al.","Shulse et al.",
              "Shulse et al.","Shulse et al.","Shulse et al.","Shulse et al.","Shulse et al.","Shulse et al.",
              "Shulse et al.","Shulse et al.","Shulse et al.","Zhang et al.")

samplesn <- character(length = ncol(gene_har))
sources <- character(length = ncol(gene_har))


for (j in 1:length(samples.all)) {
  samplesn[grepl(samples.all[j],colnames(gene_har))] <- samples3[j]
  sources[grepl(samples.all[j], colnames(gene_har))] <- samples4[j]
}


gene_har@meta.data$sample <- samplesn
gene_har@meta.data$source <- sources


#-- celltype metadata
celltypes <- data.frame()
for (i in 1:length(rc.list)) {
  tmp<- rc.list[[i]]@meta.data[,c("celltype.ID.P","timezone.ID.P","celltype.cor.P","timezone.cor.P","celltype.pvalue.P","timezone.pvalue.P")]
  celltypes <- rbind(celltypes,tmp)
}
celltypes$ids <- rownames(celltypes)
gene_har@meta.data$ids <- rownames(gene_har@meta.data)
nw <- merge(gene_har@meta.data,celltypes,by="ids", all.x=T)
rownames(nw) <- nw$ids
nw <- nw[match(rownames(gene_har@meta.data),nw$ids),]
gene_har@meta.data <- nw
#--end

#
load("~/xwj/opt/Data/T4/zs/COPILOT-master/supp_data/color_scheme_at.RData")

celltype.color <- celltypepalette[sort(match(unique(gene_har$celltype.ID.P),celltypeorder))]

celltype.color <-c(celltype.color,"#ffffff")

time.color <- celltypepalette[sort(match(unique(gene_har$timezone.ID.P),long))]


gene_har <- RunPCA(gene_har,verbose = FALSE)
gene_har <- RunUMAP(gene_har, reduction = "pca",dims = 1:50,reduction.name = "combat_umap",metric = "correlation")

Idents(gene_har) <- "celltype.ID.P"
DimPlot(gene_har,reduction = "combat_umap")

DimPlot(gene_har,reduction = "combat_umap",order = c("Maturation","Elongation","Meristem"),cols = c("#DCEDC8", "#42B3D5", "#1A237E"))

library(SeuratWrappers)
library(harmony)

gene_har <- gene_har %>% 
  RunHarmony("sample", plot_convergence = FALSE)

gene_har <- RunUMAP(gene_har, reduction = "harmony",dims = 1:50,reduction.name = "harmony2_umap")
Idents(gene_har) <- "timezone.ID.P2"
DimPlot(gene_har,reduction = "harmony2_umap",cols = celltype.color)

saveRDS(gene_har,file="gene_har_com.rds")
gene_har <- readRDS("gene_har_com.rds")

library(sva)

edata <- as.matrix(gene_har@assays$RNA@scale.data)
batch <- factor(gene_har$sample)

#combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
saveRDS(combat_edata1,"combat_edata1.rds")

gene_har@assays$RNA@scale.data <- combat_edata1


### 4. ICI annotation ----
library(ICITools)
library(dplyr)
atlas_expression_data <- Seurat::GetAssayData(gene_har_com, slot = "data") %>% as_tibble(rownames = "Locus")

spec_table_filtered_rnaseq_only <- readr::read_csv("spec_table_filtered_rnaseq_only.csv")
spec_table_filtered <- readr::read_csv("spec_table.csv")

# Use both spec score tables to compute ICI scores for 100 cells of atlas
iter <- dim(atlas_expression_data)[2]/10
ici_filtereds <- data.frame()

#Single threaded operation
##First position: Create a new progress bar
pb <- txtProgressBar(style=3)

star_time <- Sys.time() 

for(i in 11962:iter){
    setTxtProgressBar(pb, i/iter)
    start <- (i-1)*10+1
    if(start == 1) start <- 2
    end <- start+9
    if(end>dim(atlas_expression_data)[2]) end=dim(atlas_expression_data)[2]
    atlas_sub <- atlas_expression_data[,c(1,start:end)]
    
    ici_rnaseq_ath1 <- ICITools::compute_ici_scores(expression_data = atlas_sub,
                                                    spec_table = spec_table_filtered,
                                                    min_spec_score = 0.15,
                                                    information_level = 50, sig = TRUE)
    
    ici_rnaseq_only <- ICITools::compute_ici_scores(expression_data = atlas_sub,
                                                    spec_table = spec_table_filtered_rnaseq_only,
                                                    min_spec_score = 0.15,
                                                    information_level = 50, sig = TRUE)
    
    ici <- bind_rows(ici_rnaseq_ath1 %>% mutate(analysis = "ATH1RNASeq"),
                     ici_rnaseq_only %>% mutate(analysis = "RNASeq"))
    
    ici_filtered <- ici %>% group_by(Cell) %>% filter(p_adj == min(p_adj)) %>% filter(ici_score_norm == max(ici_score_norm)) %>% filter(p_val == max(p_val)) %>% filter(ici_score == max(ici_score)) 
    
    ici_filtereds <- rbind(ici_filtereds,ici_filtered)
    
}
end_time <- Sys.time()  
close(pb)
run_time <- end_time - star_time   

ici_filtereds <- ici_filtereds[!duplicated(ici_filtereds$Cell),]

saveRDS(ici_filtereds,"ici_filtereds2.rds")


ici_filtereds <- readRDS("ici_filtereds2.rds")
# Parallel operations
library(foreach)
library(doParallel)
# Real physical cores in the computer
cores <- detectCores(logical=F)
cl <- makeCluster(2)
registerDoParallel(cl, cores=cores)
registerDoParallel(2)  # use multicore, set to the number of our cores

# split data by ourselves

chunk.size <- len/cores

system.time(
  ici_filtered2 <- foreach(i=247:250, .combine='rbind') %dopar%
    {
      start <- (i-1)*10+1
      end <- start+9
      if(end>dim(atlas_expression_data)[2]) end=dim(atlas_expression_data)[2]
      atlas_sub <- atlas_expression_data[,c(1,start:end)]
      
      ici_rnaseq_ath1 <- ICITools::compute_ici_scores(expression_data = atlas_sub,
                                                      spec_table = spec_table_filtered,
                                                      min_spec_score = 0.15,
                                                      information_level = 50, sig = TRUE)
      
      ici_rnaseq_only <- ICITools::compute_ici_scores(expression_data = atlas_sub,
                                                      spec_table = spec_table_filtered_rnaseq_only,
                                                      min_spec_score = 0.15,
                                                      information_level = 50, sig = TRUE)
      
      ici <- bind_rows(ici_rnaseq_ath1 %>% mutate(analysis = "ATH1RNASeq"),
                       ici_rnaseq_only %>% mutate(analysis = "RNASeq"))
      
      ici_filtered <- ici %>% group_by(Cell) %>% filter(p_adj == min(p_adj)) %>% filter(ici_score_norm == max(ici_score_norm)) %>% filter(p_val == max(p_val)) %>% filter(ici_score == max(ici_score)) 
      
      cat(i,"--\t")
      ici_filtered
    }
)

stopImplicitCluster()
stopCluster(cl)

ici_atlas <- readRDS("ici_filtereds2.rds")

#

length(colnames(gene_har))
length(ici_atlas$Cell)

idx <- match(colnames(gene_har),ici_atlas$Cell)

ici_atlas2 <- ici_atlas[idx,]

if(identical(colnames(gene_har),ici_atlas2$Cell)){
  rownames(ici_atlas2) <- ici_atlas2$Cell
  gene_har$ici_celltype <- ici_atlas2$Cell_Type;
  gene_har$ici_celltype_score <- ici_atlas2$ici_score;
  gene_har$ici_celltype_p_val <- ici_atlas2$p_val;
  gene_har$ici_celltype_score_norm <- ici_atlas2$ici_score_norm;
  gene_har$ici_celltype_p_adj <- ici_atlas2$p_adj;
}


# Merge labels
gene_har$ici_celltype[which(gene_har$ici_celltype == "Columella_Stem")]="Columella";
gene_har$ici_celltype[which(gene_har$ici_celltype == "Ground_Cell_Stem")]="Unknown";
gene_har$ici_celltype[which(gene_har$ici_celltype == "Epidermis_LRC_Stem")]="Lateral_Root_Cap";
gene_har$ici_celltype[which(gene_har$ici_celltype == "Xylem_Stem")]="Meristematic_Xylem";
gene_har$ici_celltype[is.na(gene_har$ici_celltype)] = "Unknown"


# Plot ICI annotation
library(ggplot2)
ici.order <- unique(gene_har$ici_celltype)
ici.palette <- c('#C7EA46','#FF9900','#B67721','#21B6A8','#0082C8','#FE7F9C','#FF00FF','#0000FF','#008081','#FF0800','#C0C0C0','#5AB953','#9400D3','#7EF9FF')
gene_har$ici_celltype <- factor(gene_har$ici_celltype, levels = ici.order[sort(match(unique(gene_har$ici_celltype),ici.order))]) 
color.ici <- ici.palette[sort(match(unique(gene_har$ici_celltype),ici.order))]
options(repr.plot.width=10, repr.plot.height=8)
DimPlot(gene_har, reduction = "combat_umap", group.by = "ici_celltype", cols = color.ici)+ggtitle("ici annotation")


#### 5. Markers annotation ----

#### 5.1 AUC -----
markerall <- read.csv(file="root_marker.csv")

xy.list <- split(markerall$markers, factor(markerall$celltype))


library(GSEABase)
all.sets <- lapply(names(xy.list), function(x) {
  GeneSet(xy.list[[x]], setName=x)        
})
all.sets <- GeneSetCollection(all.sets)

library(AUCell)
ncells <- dim(gene_har)[2]
iters <- ncells/1000
pb <- txtProgressBar(style=3)

star_time <- Sys.time() 
auc.P <- c() 
auc.id <- c()
for (i in 1:iters) {
  setTxtProgressBar(pb, i/iters)
  start <- (i-1)*1000 +1
  end<- start+999
  Wt <- as.matrix(gene_har@assays$RNA@data[,128001:128549])
  rankings <- AUCell_buildRankings(Wt,plotStats=FALSE, verbose=FALSE)
  cell.aucs <- AUCell_calcAUC(all.sets, rankings,aucMaxRank=nrow(rankings)*0.05,verbose=FALSE)
  results <- t(cell.aucs@assays@data$AUC)
 
  new.labels <- colnames(results)[max.col(results)]
  s <- apply(results, 1,max)
  auc.id <- c(auc.id,new.labels)
  auc.P <- c(auc.P,s)
}
end_time <- Sys.time()  
close(pb)
run_time <- end_time - star_time  

gene_har@meta.data <- gene_meta


identical(colnames(gene_har),names(auc.P))

gene_har@meta.data$auc.P <- as.numeric(auc.P)
gene_har@meta.data$auc.id <- as.character(auc.id)

saveRDS(gene_har@meta.data,file="gene_meta.rds")

### 5.2 --CellAssign-- (Destroy)----

exp <- gene_har@assays$RNA@data
markers <- as.character(unique(unlist(xy.list)))

exps <- exp[rownames(exp) %in% markers,]
rho <- marker_list_to_mat(xy.list, include_other = FALSE)
rho<- rho[rownames(rho) %in% rownames(exps),]

exps1 <- as.matrix(exps[,1:1000])

time1 <- Sys.time()  

result <- cellassign(t(exps1),
                     marker_gene_info = rho,
                     s = colSums(exps1),
                     learning_rate = 1e-2,
                     shrinkage = TRUE,
                     verbose = FALSE)
time2 <- Sys.time()  

cost <- time2 -time1

gene_har@meta.data$celltype.ID.P
gene_har@meta.data$auc.id

###6 Annotation Merge----
#1 Determine the most accurate classification through primary classification;
#2 Calculate the markers for each category (calculated by coefficient of variation)
#3 Reclassify using new markers as reference data.

#mix Lateral Root Cap and non hair cells
###6.1 Cell types----
##6.1.1 First level classification----

# For making consensus annotation, we only keep those with correlation coefficient > 0.6 and adjusted p-value/p-value < 0.01

# Correlation-based annotation
hc.rna.anno <- as.character(gene_har$celltype.ID.P)
hc.rna.anno[which(gene_har$celltype.pvalue.P > 0.01)] = NA
#hc.rna.anno[which(gene_har$celltype.cor.P < 0.6)] = NA


hc.ma.anno <- as.character(gene_har$Rad.ID.P)
#hc.ma.anno[which(gene_har$Rad.pvalue.P > 0.01)] = NA;
#hc.ma.anno[which(gene_har$Rad.cor.P < 0.6)] = NA;

# ICI annotation
hc.ici.anno <- as.character(gene_har$ici_celltype)
#hc.ici.anno[which(gene_har$ici_celltype_p_adj > 0.01)] = NA
# Marker-based annotation
hc.km.anno <- as.character(gene_har$auc.id)

hc.km.anno[which(hc.km.anno == "Hair Cells")]="Trichoblast"
hc.km.anno[which(hc.km.anno == "Non-hair Cells")]="Atrichoblast"
hc.km.anno[which(hc.km.anno == "Protophloem sieve elements")]="Phloem"

hc.ma.anno[which(hc.ma.anno == "Hair Cell")]="Trichoblast"
hc.ma.anno[which(hc.ma.anno == "LRC")]="Lateral Root Cap"
hc.ma.anno[which(hc.ma.anno == "Mat.Xylem")]="Xylem"
hc.ma.anno[which(hc.ma.anno == "Meri.Xylem")]="Xylem"
hc.ma.anno[which(hc.ma.anno == "Non-Hair Cell")]="Atrichoblast"
hc.ma.anno[which(hc.ma.anno == "Phloem [S32]")]="Phloem"
hc.ma.anno[which(hc.ma.anno == "Phloem [SUC2]")]="Phloem"
hc.ma.anno[which(hc.ma.anno == "QC")]="Quiescent Center"
hc.ma.anno[which(hc.ma.anno == "Phloem Pole Pericycle")]="Pericycle"
hc.ma.anno[which(hc.ma.anno == "Xylem Pole Pericycle")]="Pericycle"

hc.rna.anno[which(hc.rna.anno == "columella")]="Columella"
hc.rna.anno[which(hc.rna.anno == "developing cortex")]="Cortex"
hc.rna.anno[which(hc.rna.anno == "developing xylem")]="Xylem"
hc.rna.anno[which(hc.rna.anno == "endodermis & QC cells")]="Endodermis"
#hc.rna.anno[which(hc.rna.anno == "endodermis & QC cells")]="Quiescent Center"
hc.rna.anno[which(hc.rna.anno == "hair cells")]="Trichoblast"
hc.rna.anno[which(hc.rna.anno == "LRC & non-hair cells")]="Lateral Root Cap"
#hc.rna.anno[which(hc.rna.anno == "LRC & non-hair cells")]="Atrichoblast"
hc.rna.anno[which(hc.rna.anno == "matured cortex")]="Cortex"
hc.rna.anno[which(hc.rna.anno == "matured endodermis")]="Endodermis"
hc.rna.anno[which(hc.rna.anno == "matured xylem pole")]="Xylem"
hc.rna.anno[which(hc.rna.anno == "non-hair cells")]="Atrichoblast"
hc.rna.anno[which(hc.rna.anno == "phloem & companion cells")]="Phloem"
hc.rna.anno[which(hc.rna.anno == "phloem pole pericycle")]="Pericycle"
hc.rna.anno[which(hc.rna.anno == "protophloem & metaphloem")]="Phloem"
hc.rna.anno[which(hc.rna.anno == "QC cells")]="Quiescent Center"

hc.ici.anno[which(hc.ici.anno == "Mature_Cortex")]="Cortex"
hc.ici.anno[which(hc.ici.anno == "Phloem_Pole_Pericycle")]="Pericycle"
hc.ici.anno[which(hc.ici.anno == "Maturing_Xylem")]="Xylem"
hc.ici.anno[which(hc.ici.anno == "Developing_Cortex")]="Cortex"
hc.ici.anno[which(hc.ici.anno == "Phloem_Companion_Cells")]="Phloem"
hc.ici.anno[which(hc.ici.anno == "Meristematic_Xylem")]="Xylem"
hc.ici.anno[which(hc.ici.anno == "Quiescent_Center")]="Quiescent Center"
hc.ici.anno[which(hc.ici.anno == "Lateral_Root_Cap")]="Lateral Root Cap"
hc.ici.anno[which(hc.ici.anno == "Xylem_Pole_Pericycle")]="Pericycle"

# Merge annotations and metadata
cb.anno <- data.frame(hc.rna.anno=hc.rna.anno,hc.ma.anno=hc.ma.anno,hc.ici.anno=hc.ici.anno,hc.km.anno=hc.km.anno)
cb.anno$annotations <- "NA"

num.unique <- apply(cb.anno,1,function(x){length(unique(x[which(!is.na(x))]))})
num.not.NA <- apply(cb.anno,1,function(x){length(which(!is.na(x)))})
num.max <- apply(cb.anno,1,function(x){max(table(as.character(x)))})

ici.not.NA <- apply(cb.anno,1,function(x){length(which(!is.na(x[3])))})
km.not.NA <- apply(cb.anno,1,function(x){length(which(!is.na(x[4])))})
cb.anno$num.unique <- num.unique 
cb.anno$num.not.NA <- num.not.NA 
cb.anno$num.max <- num.max 
cb.anno$ici.not.NA <- ici.not.NA  
cb.anno$km.not.NA <- km.not.NA  
cb.anno$Cell_ID <- colnames(gene_har) 
cb.anno$order <- seq(1,ncol(gene_har)) 


# Voting, a cell needs at least two supports from different annotation methods in order to remain annotated. Marker-based annotation is given higher weight on final annotation decision
cb.anno$level <- 0

# Level 1: num.max >= 3

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv <- uniqv[!is.na(uniqv)]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

cb.anno.sub1 <- cb.anno %>% filter(num.max>=3)

cb.anno.sub1$annos <- apply(cb.anno.sub1[,1:4],1,getmode)

cb.anno[cb.anno.sub1$order,]$level <- 1
cb.anno[cb.anno.sub1$order,]$annotations <- cb.anno.sub1$annos
table(cb.anno$annotations)

# Level 2: num.max == 2: 2v2;2v1v1;2v0
cb.anno.sub2 <- cb.anno %>% filter(level==0 & num.max == 2)

cb.anno.sub2$annos <- apply(cb.anno.sub2[,1:4],1,getmode)

cb.anno.sub2[cb.anno.sub2$num.not.NA == 4 & cb.anno.sub2$num.unique ==2,]$annos <- cb.anno.sub2[cb.anno.sub2$num.not.NA == 4 & cb.anno.sub2$num.unique ==2,]$hc.km.anno

cb.anno[cb.anno.sub2$order,]$level <- 2
cb.anno[cb.anno.sub2$order,]$annotations <- cb.anno.sub2$annos
table(cb.anno$annotations)
 
cb.anno[cb.anno$hc.km.anno == "Procambium" & cb.anno$annotations=="Pericycle",]$annotations <- "Procambium"

cb.anno.sub3 <- cb.anno %>% filter(level==0 & num.unique >1) #6062

cb.anno[idx1 & cb.anno$annotations == "Atrichoblast",]

gene_har$le <- cb.anno$annotations
gene_har@meta.data[gene_har$le =="NA",]$le <- NA
DimPlot(gene_har, reduction = "combat_umap", group.by = "le")

gene_har$le0.01 <- cb.anno$annotations
gene_har@meta.data[gene_har$le0.01 =="NA",]$le0.01 <- NA
DimPlot(gene_har, reduction = "combat_umap", group.by = "le0.01")

# Level 3: num.max == 1

gene_har$cb.anno.1 <- cb.anno$annotations

# Plot initial consensus annotation



# Merge cell type consensus annotation with time zone correlation-based annotation
long.merge <- as.character(gene_har$Long.ID.P)
cb.anno.t <- data.frame(rna=as.character(gene_har$timezone.ID.P), ma=long.merge)
cb.anno.t$order <- seq(1,ncol(gene_har))

cb.anno.t$annotation <- apply(cb.anno.t,1,function(x){if(TRUE){x[2]}else{"Unknown"}})
cb.anno.t.sub <- cb.anno.t %>% filter(annotation!="Unknown")
cb.anno.m <- left_join(cb.anno.sub,cb.anno.t.sub,by="order")
cb.anno.m <- cb.anno.m %>% filter(annotations != "NA" & annotations != "Unknown" & annotation != "NA")

# Dissect each cell type based on microarray time zone correlation-based annotation 
# to have more grouping of cells/ to give finest resolution 
cb.anno.m$final.annotation <- paste(cb.anno.m$annotations, cb.anno.m$annotation, sep="-")
cb.anno.m$final.annotation[grep(paste(names(table(cb.anno.m$final.annotation)[which(table(cb.anno.m$final.annotation) <= 2)]),collapse = "|"),cb.anno.m$final.annotation)] = "Unknown"
cb.anno.m <- cb.anno.m %>% filter(final.annotation != "Unknown")

# 88 new reference expression profiles are going to be build
length(sort(table(cb.anno.m$final.annotation)))
gene_har$cb.anno <- rep("Unknown", ncol(gene_har))
gene_har$cb.anno[as.numeric(cb.anno.m$order)] <- as.character(cb.anno.m$final.annotation)

sorttypes <- as.character(unique(gene_har$cb.anno[which(gene_har$cb.anno!="Unknown")]))


###6.1.2 Refactoring classification-----
#Here, new markers are used as reference data for reclassification

# Extract integrated (batch-corrected) expression matrix

sorttypes <- as.character(names(table(gene_har$le)))
afm <- gene_har@assays$RNA@data
# Pool (average) expression values of each grouping
new_ref <- matrix(nrow=nrow(afm), ncol = 0)

for (i in 1:length(sorttypes)) {
  m <- afm[,which(gene_har$le==sorttypes[i])]
  new_ref <- cbind(new_ref, rowSums(m)/ncol(m))
}

colnames(new_ref) <- sorttypes
gene.var <- apply(new_ref,1,var)

# Select top 200 highly variable genes 
new_ref_sub <- new_ref[names(sort(gene.var,decreasing = TRUE)[1:200]),]

# Merge 88 newly-built reference with atlas

merge.rownames <- function (x,y){
  dat <- merge(x = x, y = y, by = "row.names")
  rownames(dat) <- dat$Row.names
  dat <- dat[,-1]
  return(dat)
}

nr_label=colnames(new_ref)

nranno <- list()
afm_sub <- afm[,colnames(afm) %in% colnames(gene_har)[is.na(gene_har$le)]]
nr <- Reduce(merge.rownames, list(new_ref_sub,afm_sub))
nr <- as.matrix(nr)

nranno <- celltypeTest(nr,nr_label,method = "pearson")


gene_har$final.ID <- gene_har$le
gene_har$final.ID[gene_har$ids %in% colnames(gene_har)[is.na(gene_har$le)]]<-as.character(names(nranno$celltype_max))

qc <- gene_har@meta.data[gene_har$final.ID == "Quiescent Center",c("order","celltype.ID.P2","Rad.ID.P2","ici_celltype2","auc.id")]

idxSCN<-qc[rowSums(qc=="Quiescent Center") <= 1,]$order

gene_har$final.ID2 <- as.character(gene_har$final.ID)

gene_har$final.ID2[idxSCN] <- as.character("Stem Cell Niche")


idx1 <- gene_har[["combind_umap"]]@cell.embeddings[,1] >6
idx2 <- gene_har$final.ID == "Atrichoblast"
gene_har@meta.data[idx1&idx2,c("order","celltype.ID.P2","Rad.ID.P2","ici_celltype2","auc.id")]

gene_har$ici_p_sig <- rep(">= 0.01", nrow(gene_har@meta.data))
gene_har$ici_p_sig[which(gene_har$ici_celltype_p_adj < 0.05)] <- "< 0.01"

DimPlot(gene_har, reduction = "combat_umap", group.by = "ici_p_sig")

saveRDS(gene_har@meta.data,file="gene_meta.rds")


###6.2 Main cell type branches----
hc.fina <- as.character(gene_har$final.ID2)

hc.fina[which(hc.fina == "Columella")]="Root cap"
hc.fina[which(hc.fina == "Lateral Root Cap")]="Root cap"


hc.fina[which(hc.fina == "Trichoblast")]="Epidermis"
hc.fina[which(hc.fina == "Atrichoblast")]="Epidermis"

hc.fina[which(hc.fina == "Cortex")]="Ground tissue"
hc.fina[which(hc.fina == "Endodermis")]="Ground tissue"

hc.fina[which(hc.fina == "Xylem")]="Stele"
hc.fina[which(hc.fina == "Phloem")]="Stele"
hc.fina[which(hc.fina == "Procambium")]="Stele"
hc.fina[which(hc.fina == "Pericycle")]="Stele"

gene_har$Main.tp <- hc.fina

DimPlot(gene_har,reduction = "combat_umap",group.by = "Main.tp")

DimPlot(gene_har,reduction = "combat_umap",group.by = "celltype.ID.P2")
DimPlot(gene_har,reduction = "combat_umap",group.by = "ici_celltype2")
DimPlot(gene_har,reduction = "combat_umap",group.by = "Rad.ID.P2")

gene_har@meta.data[gene_har$Main.tp == "Stele" & gene_har$Long.ID.P2 == "Maturation",]


ss <- gene_har@meta.data[gene_har$auc.id == "Trichoblast" & gene_har$ici_celltype2 == "Trichoblast" & gene_har$cb.anno.2 != "Trichoblast",]

ss[,c("celltype.ID.P2","ici_celltype2","Rad.ID.P2","cb.anno.2")]



###6.3 Time classification----

# Merge cell type consensus annotation with time zone correlation-based annotation
long.merge <- as.character(gene_har_com$Long.ID.P)
long.merge[grep("Meri",long.merge)]="Meristem";
long.merge[grep("Elong",long.merge)]="Elongation";
long.merge[grep("Mat",long.merge)]="Maturation";
long.merge[grep("Columella",long.merge)]="Maturation"
gene_har_com$Long.ID.P2 <- long.merge

time.anno <- data.frame(order=seq(1,ncol(gene_har_com)),time=gene_har_com$timezone.ID.P,Long=gene_har_com$Long.ID.P2,timeanno=NA)

time.anno[time.anno$time == time.anno$Long,]$timeanno <- time.anno[time.anno$time == time.anno$Long,]$time

gene_har_com$timelv1 <- time.anno$timeanno


gene_har_com@meta.data[(gene_har_com$final.ID2 == "Stem Cell Niche" | gene_har_com$final.ID2 =="Quiescent Center") & gene_har_com$timelv1 != "Meristem" & !is.na(gene_har_com$timelv1),]$timelv1 <- NA


sorttypes <- as.character(names(table(gene_har_com$timelv1)))

afm <- gene_har_com@assays$RNA@data
new_ref <- matrix(nrow=nrow(afm), ncol = 0)

for (i in 1:length(sorttypes)) {
  m <- afm[,which(gene_har_com$timelv1==sorttypes[i])]
  new_ref <- cbind(new_ref, rowSums(m)/ncol(m))
}

colnames(new_ref) <- sorttypes
gene.var <- apply(new_ref,1,var)

# Select top 200 highly variable genes 
new_ref_sub <- new_ref[names(sort(gene.var,decreasing = TRUE)[1:200]),]

# Merge 88 newly-built reference with atlas

merge.rownames <- function (x,y){
  dat <- merge(x = x, y = y, by = "row.names")
  rownames(dat) <- dat$Row.names
  dat <- dat[,-1]
  return(dat)
}

nr_label=colnames(new_ref)

afm_sub <- afm[,colnames(afm) %in% colnames(gene_har_com)[is.na(gene_har_com$timelv1)]]
nr <- Reduce(merge.rownames, list(new_ref_sub,afm_sub))
nr <- as.matrix(nr)
timeannos <- celltypeTest(nr,nr_label,method = "pearson")

gene_har_com$timelv1[gene_har_com$ids %in% colnames(gene_har_com)[is.na(gene_har_com$timelv1)]] <- as.character(names(timeannos$celltype_max))


gene_har_com@meta.data[(gene_har_com$final.ID2 == "Stem Cell Niche" | gene_har_com$final.ID2 =="Quiescent Center") & gene_har_com$timelv1 != "Meristem" & !is.na(gene_har_com$timelv1),]$timelv1 <- "Meristem"


saveRDS(gene_har_com,file="gene_har_com.rds")

saveRDS(gene_har_com@meta.data ,file= "gene_meta.rds")

### Figure 1 ----
#Gene_har_com.rds -- Results of Gene Integration Analysis
gene_har_com=readRDS('gene_har_com.rds')


gene_har_com$Source<-gene_har_com$source
gene_har_com$Branch.tp<-gene_har_com$Main.tp
gene_har_com$Cell.tp<-gene_har_com$final.ID2
gene_har_com$Develop.st<-gene_har_com$timelv1

#filter non pacs cells
figure1A <- DimPlot(gene_har_com, reduction = "combat_umap", group.by = "Source")
#figure1
order <- c("Epidermis",
            "Ground tissue", 
            "Root cap",
            "Stele",
            "Stem Cell Niche","Quiescent Center")
palette <- c("#7fb80e", 
              "#2a5caa", 
              "#843900", 
              "#f7b13f",
              "#f58f98","#d71345")


order <- c("Epidermis",
           "Ground tissue", 
           "Root cap",
           "Stele")
palette <- c("#7fb80e", 
             "#2a5caa", 
             "#843900", 
             "#f7b13f")

figure1B <- DimPlotOrder(gene_har_com,reduction = "combat_umap",group.by = "Branch.tp",order=order,palette=palette)
#DimPlot(gene_har_com,reduction = "combat_umap",group.by = "Main.tp",order=order,palette=palette)
DimPlotOrder(gene_har_com,reduction = "combat_umap",group.by = "Main.tp",order=order,palette=palette)
gene_har_com=readRDS("../gene_har_com.rds")


#figure3
order <- c("Atrichoblast", "Trichoblast",
           "Cortex", "Endodermis", 
           "Columella", "Lateral Root Cap",
           "Pericycle", "Phloem", "Xylem", "Procambium",
           "Stem Cell Niche","Quiescent Center")
palette <- c("#7fb80e", "#1d953f",
             "#2a5caa", "#a565ef",
             "#843900", "#69541b",
             "#f7b13f", "#f9e264", "#009db2", "#375830",
             "#f58f98","#d71345")

#Key(gene_har_com@reductions$combat_umap) <- "UMAP_"
figure1C <- DimPlotOrder(gene_har_com,reduction = "combat_umap",group.by = "Cell.tp",order=order,palette=palette)


order <- c("Maturation","Elongation","Meristem")
palette <- c("#d71345","#f391a9","#feeeed")
figure1D <- DimPlotOrder(gene_har_com,reduction = "combat_umap",group.by = "Develop.st",order=order,palette=palette)


library(cowplot)

Figure1 <- ggdraw() +     
  draw_plot(figure1A, 0, 0, 0.5, 0) + 
  draw_plot(figure1B, 0.5, 0.5, 0.5, 0.5) + 
  draw_plot(figure1C, 0.5, 0, 0.5, 0.5) + 
  draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0.5), c(1, 1, 0.5))  

plot_grid(figure1A,figure1B,figure1C,figure1D,nrow=2,ncol=2,labels=c('A','B','C','D'))
### Figure S1 ----

figureS1A <- DimPlot(gene_har_com, reduction = "combat_umap", group.by = "source")

# gene_har_com@meta.data <- gene_meta
# gene_har_com@meta.data$le <- as.character(gene_har_com@meta.data$le)
# gene_har_com@meta.data$le[!is.na(gene_har_com@meta.data$le)] <- "marked"
# gene_har_com@meta.data$le[is.na(gene_har_com@meta.data$le)] <- "unmarked"

order <- c("marked", "unmarked")
palette <- c("#7fb80e",rgb(252,157,154, maxColorValue = 255))
figureS1B <- DimPlotOrder(gene_har_com,reduction = "combat_umap",group.by = "le",order=order,palette=palette)


order <- c("Maturation","Elongation","Meristem")
palette <- c("#d71345","#f391a9","#feeeed")
figureS1C <- DimPlotOrder(gene_har_com,reduction = "combat_umap",group.by = "timelv1",order=order,palette=palette)

FigureS1 <- ggdraw() +     
  draw_plot(figureS1A, 0, 0.5, 0.5, 0.5) +  
  draw_plot(figureS1B, 0.5, 0.5, 0.5, 0.5) +  
  draw_plot(figureS1C, 0, 0, 0.5, 0.5) + #  
  draw_plot_label(c("A", "B", "C"), c(0, 0.5, 0), c(1, 1, 0.5))  


### END

### >> DE markers----

gene_har<-readRDS('gene_har.rds')
gene_har_com<-readRDS('gene_har_com.rds')

# Idents(gene_har) <- "final.ID2"
Idents(gene_har_com) <- "final.ID2"
#gmarkers <- FindAllMarkers(gene_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gmarkers <- FindAllMarkers(gene_har_com, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

gmarkers1<-readRDS('gmarkers.rds')

saveRDS(gmarkers,"gmarkers.rds")

obj_har<-readRDS('obj_har.rds')

Idents(obj_har) <- "final.ID2"
pamarkers <- FindAllMarkers(obj_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pamarkers,file="pamarkers.rds")
saveRDS(obj_har,file="obj_har.rds")

Idents(obj_har) <- "Main.tp"
mainmarkers <- FindAllMarkers(obj_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#pamarkers %>% group_by(cluster) %>% top_n(n=2,wt=avg_logFC)



Idents(ratio_har) <- "final.ID2"
ramarkers <- FindAllMarkers(ratio_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(ramarkers,file="ramarkers.rds")

Idents(mura_har) <- "final.ID2"
mumarkers <- FindAllMarkers(mura_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mumarkers %>% group_by(cluster) %>% top_n(n=2,wt=avg_logFC)



saveRDS(mumarkers,file="mumarkers.rds")

Idents(mura_har) <- "Main.tp"
mu_main <- FindAllMarkers(mura_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(mu_main,file="mu_main.rds")





### overlapping pamarkers 3.3.1----
pamarkers2<-readRDS('pamarkers2.rds')
pamarkers=readRDS('pamarkers.rds')
obj_meta=readRDS('obj_meta.rds')

tairsem <- read.csv('genes.tsv', header=F, sep="\t")

colnames(tairsem) <- c("genes","symbol")

pamarkers$genes <- gsub("-.*","",gsub("PAC","",pamarkers$gene))

pamarkers_ <- merge(pamarkers,tairsem,by="genes",all.x = T)

inte_pa <- list()
for(i in names(table(pamarkers2$cluster.x))){
  sp <- markerall[markerall$celltype == i,]$markers
  mp <- pamarkers2[pamarkers2$cluster.x == i,]$genes
  inte_pa[[i]] <- intersect(sp,mp)
  #inte_pa[[i]] <- pamarkers2[pamarkers2$genes %in% is,]$`Gene Name`
}

library('readr')
marker_metadata <- read_csv("COPILOT-master/supp_data/marker_metadata.csv")

s <- intersect(marker_metadata$Marker, pamarkers2$symbol)

pamarkers2[pamarkers2$symbol %in% s,]


mumarkers<-readRDS('mumarkers.rds')
mumarkers$genes <- gsub("-.*","",gsub("PAC","",mumarkers$gene))
inte_mu <- list()

#markerall?
for(i in names(table(mumarkers$cluster))){
  sp <- markerall[markerall$celltype == i,]$markers
  mp <- mumarkers[mumarkers$cluster == i,]$genes
  inte_mu[[i]] <- intersect(sp,mp)
}



##*****filter markers********####
filter.ge <- filter.ratio(gene_har,group.by = "final.ID2",delist = gmarkers,pt.1=0.2,pt.2=0.2)
# 3051    7
table(filter.ge$cluster)
# Atrichoblast      Trichoblast           Cortex       Endodermis        Columella Lateral Root Cap        Pericycle           Phloem 
# 391              399              127              273              172              126              227              107 
# Xylem       Procambium  Stem Cell Niche Quiescent Center 
# 334               26              283              586

readr::write_csv(filter.ge,file = "filter.ge.csv")

filter.pa <- filter.ratio(obj_har,group.by = "final.ID2",delist = pamarkers,pt.1=0.2,pt.2=0.2)


table(filter.pa$cluster)

#readr::write_csv(filter.pa,file = "filter.pa.csv")

pamarkers<-readRDS('pamarkers.rds')
filter.main <- filter.ratio(obj_har,group.by = "Main.tp",delist = pamarkers,pt.1=0.2,pt.2=0.2)


mainmarkers<-readRDS('mainmarkers.rds')
filter.main <- filter.ratio(obj_har,group.by = "Main.tp",delist = mainmarkers,pt.1=0.2,pt.2=0.2)




filter.ra <- filter.ratio(ratio_har,group.by = "final.ID2",delist = ramarkers,pt.1=0.2,pt.2=0.2)

table(filter.ra$cluster)

#readr::write_csv(filter.ra,file = "filter.ra.csv")

#filter.mu <- filter.ratio(mura_har,group.by = "final.ID2",delist = mumarkers,pt.1=0.2,pt.2=0.2)
#table(filter.mu$cluster)
#readr::write_csv(filter.mu,file = "filter.mu.csv")


#### Overlapping pa vs ge 3.3.2----
gene_har<-readRDS('gene_har.rds')
gmarkers<-readRDS('gmarkers.rds')

NDEG <- rownames(gene_har)[!(rownames(gene_har) %in% gmarkers$gene)]

### nge vs pa
pamarkers<-readRDS('pamarkers.rds')

pa_ges <- gsub("PAC","",gsub("-.*","",pamarkers$gene))
table(pamarkers[pa_ges %in% NDEG,]$cluster)


#----(/..OK../)
filter.pa11<-read.csv('filter.pa.csv')
#filter.pa<-readRDS('filter.pa.rds')
fpa_ges11 <- gsub("PAC","",gsub("-.*","",filter.pa11$gene))
table(filter.pa11[fpa_ges11 %in% NDEG,]$cluster)

### nge vs ra
filter.ra<-read.csv('filter.ra.csv')
ramarkers<-readRDS('ramarkers.rds')

ra_ges <- gsub("PAC","",gsub("-.*","",ramarkers$gene))
table(ramarkers[ra_ges %in% NDEG,]$cluster)


fra_ges <- gsub("PAC","",gsub("-.*","",filter.ra$gene))
table(filter.ra[fra_ges %in% NDEG,]$cluster)


### nge vs mu

mu_ges <- gsub("PAC","",gsub("-.*","",mumarkers$gene))
table(mumarkers[mu_ges %in% NDEG,]$cluster)

muges <-mumarkers[mu_ges %in% NDEG,]$gene
length(unique(gsub("-.*","",gsub("PAC","",muges))))

#fmu_ges <- gsub("PAC","",gsub("-.*","",filter.mu$gene))
#table(filter.mu[fmu_ges %in% NDEG,]$cluster)


### Plot
gene_har<-readRDS('gene_har.rds')
gene_har_com<-readRDS('gene_har_com.rds')
obj_har<-readRDS('obj_har.rds')
NGPA3<-read.csv('NGPA.csv')
NGPAG3<-NGPA3$gene

s11 <- FeaturePlot(obj_har, features = NGPA3$gene[286],reduction="harmony_umap")

NGPAG3<-gsub("-.*","",gsub("PAC","",NGPAG3))
#s22 <- FeaturePlot(gene_har, features = NGPAG[286],reduction="combat_umap")+ggtitle("Non-differential gene: AT2G25710")
s22 <- FeaturePlot(gene_har_com, features =NGPAG3[286],reduction="combat_umap")
#s2$data$AT3G10610 <- as.numeric(gene_har@assays$RNA@scale.data["AT3G10610",])
s11$data$cell <- obj_har@meta.data$ids
s22$data$cell <- rownames(s22$data)
ns1 <- left_join(s22$data,s11$data,by="cell")
s33 <- s22
s33$data$AT5G64816 <- ns$PACAT5G64816.1
s33$labels$x <- "UMAP1"
s33$labels$y <- "UMAP2"

ids<-gsub(".*_","",gsub("-.*","",obj_har@meta.data$ids))


s22$labels$title <- "Non-differential gene: AT5G64816"
s33$labels$title <- "differential APA: AT5G64816-1"

plot_grid(s22,s33,labels=c("A","B"))


v1 <- VlnPlot(obj_har, features =  NGPA$gene[286], slot = "data", log = TRUE,pt.size=0)

v2 <- VlnPlot(gene_har, features =  NGPAG[286], slot = "data", log = TRUE,pt.size=0)

which.min(apply(gene_har@assays$RNA@data[NGPAG,],1,sd))

plot_grid(v2,v1,ncol=1)

DimPlot(gene_har,reduction = "combat_umap")




#### GO----

library(clusterProfiler)
library("org.At.tair.db")

cps <- names(table(gmarkers$cluster))

syms <- gmarkers$gene[gmarkers$p_val_adj<0.05 & gmarkers$cluster ==cps[1]]

y2 = enrichGO(syms, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")

enrichplot::dotplot(y2,showCategory=5)

b1 <- barplot(y2, showCategory=10,title="Gene") 
b2 <- barplot(y3, showCategory=10,title="APA") 

plot_grid(b1,b2,ncol=2)



syms <- filter.pam$gene[filter.pam$p_val_adj<0.05 & filter.pam$cluster ==cps[1]]
syms <- gsub("PAC","",gsub("-.*","",syms))
y3 = enrichGO(syms, OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")


library(ggplot2)
library(dplyr)
library(stringr)


dot_df <- arrange(y2@result,pvalue)[1:5,]

## plot
library(forcats) ## for reordering the factor
ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")

p <- ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")

p + facet_grid(.~type)




#### Hair vs Non-Hair 3.3.3-----

gs <- gsub("-.*","",gsub("PAC","",rownames(mura_har)))

kr <- data.frame("gene"=gs,"PAC"=rownames(mura_har))
kp <- utr3_anno[!duplicated(utr3_anno$gene),][,c("strand","gene")]

kr2 <- merge(kr,kp,by="gene",all.x=T)
pacs <- c()
for(i in levels(factor(kr2$gene))){
  strand <- kr2[kr2$gene == i,"strand"][1]
  PAC <- kr2[kr2$gene == i,"PAC"]
  nums <- gsub("PAC.*-","",PAC)
  
  if(strand == "+"){
    pacs <- c(pacs,PAC[which.min(nums)])
  }else{
    pacs <- c(pacs,PAC[which.max(nums)])
  }
}

single_har <- mura_har[rownames(mura_har) %in% pacs,]
single_har <- NormalizeData(single_har, verbose = FALSE)
single_har <- FindVariableFeatures(single_har, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
single_har <- ScaleData(single_har, verbose = FALSE)

single_har <- RunPCA(single_har,verbose = FALSE)
single_har <- RunUMAP(single_har, reduction = "pca", dims = 1:10,reduction.name = "combind_umap")
Idents(single_har) <- "samples"
DimPlot(single_har,reduction = "combind_umap")
library(harmony)
single_har <- single_har %>% 
  RunHarmony("samples", plot_convergence = FALSE)

single_har <- RunUMAP(single_har, reduction = "harmony",dims = 1:ncol(single_har[["harmony"]]),reduction.name = "harmony_umap")
DimPlot(single_har,label = T)
Idents(single_har) <- "final.ID2"

Atri_markers <- FindMarkers(object = single_har, ident.1 = 'Atrichoblast', ident.2 = "Trichoblast")

simarkers <- FindAllMarkers(object = single_har)


saveRDS(single_har,file="single_har.rds")


#### Trichoblast vs Atrichoblast
mura_har<-readRDS('mura_har.rds')

Idents(mura_har) <- "final.ID2"


#Identify the difference in DEAPA ratio markers between Trichoblast and Atrichoblast
Atri_markers <- FindMarkers(object = mura_har, ident.1 = 'Atrichoblast', ident.2 = "Trichoblast")

test<- FindMarkers(object = mura_har, ident.1 = 'Atrichoblast', ident.2 = "Trichoblast")


Asgene <- gsub("PAC","",gsub("-.*","",rownames(Atri_markers)))

pas <- rownames(obj_har)
gs <- gsub("PAC","",gsub("-.*","",rownames(obj_har)))
dtls <- data.frame(PAC_id = rownames(obj_har), gene=gs)


new_3utr2<-readRDS('new_3utr2.rds')

nt <- new_3utr2[,c("PAC_id","dis2")]
nt$PAC_id <- gsub("_","-",nt$PAC_id)

pac_len <- merge(dtls,nt,by="PAC_id")

sel_pac <- pac_len[pac_len$gene %in% Asgene,]$PAC_id

istype1 <- obj_har$final.ID2 == "Atrichoblast"
ns1 <- rowSums(obj_har@assays$RNA@counts[sel_pac,istype1]) / sum(istype1)

istype2 <- obj_har$final.ID2 == "Trichoblast"
ns2 <- rowSums(obj_har@assays$RNA@counts[sel_pac,istype2]) / sum(istype2)

sub_dtls <- pac_len[pac_len$PAC_id %in% names(ns1),]
if(identical(sub_dtls$PAC_id,names(ns1))){
  sub_dtls$b1 <- ns1
}
if(identical(sub_dtls$PAC_id,names(ns2))){
  sub_dtls$b2 <- ns2
}
sub_dtls$c1 <- sub_dtls$dis2 * sub_dtls$b1
sub_dtls$c2 <- sub_dtls$dis2 * sub_dtls$b2

g1 <- tapply(sub_dtls$c1, INDEX=sub_dtls$gene, FUN=sum) / tapply(sub_dtls$b1, INDEX=sub_dtls$gene, FUN=sum)
g2 <- tapply(sub_dtls$c2, INDEX=sub_dtls$gene, FUN=sum) / tapply(sub_dtls$b2, INDEX=sub_dtls$gene, FUN=sum)

dislen <- g1-g2

lengthingA <- dislen[dislen>0] 
shortingA <- dislen[dislen<0] 

ns1 <- obj_har@assays$RNA@data[sel_pac,istype1]

sub_dtls <- pac_len[pac_len$PAC_id %in% rownames(ns1),]

if(identical(sub_dtls$PAC_id,rownames(ns1))){
  type_dis1 <- sub_dtls$dis2 * ns1
  type_dis1 <- as.data.table(type_dis1)
  type_dis1$gene <- sub_dtls$gene
  
  type_w1 <- type_dis1[,lapply(.SD, sum),by=gene]
}

ns2 <- obj_har@assays$RNA@data[sel_pac,istype2]

sub_dtls <- pac_len[pac_len$PAC_id %in% rownames(ns2),]

if(identical(sub_dtls$PAC_id,rownames(ns2))){
  type_dis2 <- sub_dtls$dis2 * ns2
  type_dis2 <- as.data.table(type_dis2)
  type_dis2$gene <- sub_dtls$gene
  
  type_w2 <- type_dis2[,lapply(.SD, sum),by=gene]
}

result <- data.frame()
for (n in 1:nrow(type_w1)) {
  type_w1s <- as.numeric(type_w1[n,-1])
  type_w2s <- as.numeric(type_w2[n,-1])
  
  gene_n <- data.frame(group=c(rep(1,length(type_w1s)),rep(2,length(type_w2s))),value= c(type_w1s,type_w2s))
  
  p_value <- wilcox.test(value~group, gene_n)$p.value
  
  result <- rbind(result,data.frame(gene=type_w1$gene[n],pvalue=p_value))
}

res.001 <- result[result$pvalue <0.01,]

lenA <- lengthingA[names(lengthingA) %in% res.001$gene]
shoA <- shortingA[names(shortingA) %in% res.001$gene]
length(lenA)
length(shoA)

ATswitch <- list(len=lenA,sho=shoA)

saveRDS(ATswitch,file="ATswitch_single.rds")

library(clusterProfiler)

ATswitch<-readRDS('ATswitch.rds')
ATswitch<-readRDS('ATswitch_single.rds')


AlenGO <- enrichGO(names(ATswitch$len), OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")
AshoGO <- enrichGO(names(ATswitch$sho), OrgDb = "org.At.tair.db", keyType="TAIR", ont="BP")


b1 <- barplot(AlenGO, showCategory=20,title="Atrichoblast APA switch lengthening")
b2 <- barplot(AshoGO, showCategory=20,title="Atrichoblast APA switch shorting")

plot_grid(b1,b2,ncol=2)

write.csv(AlenGO@result,file="Atrichoblast_lengthening.csv")
write.csv(AshoGO@result,file="Atrichoblast_shorting.csv")

#### Stele
table(mura_har$final.ID2)
stele.list <- c("Pericycle","Phloem","Xylem","Procambium")
stele.markers <- list()
for (i in 1:4) {
  for (j in (i+1):4) {
    cat(i,":",j,"-->")
    if(j<5){
      sname <- paste0(stele.list[i],".",stele.list[j])
      stele.markers[[sname]] <- FindMarkers(object = mura_har, ident.1 = stele.list[i], ident.2 = stele.list[j])
    
    }
  }
}
stele.markers.neg <- list()

for (i in 1:4) {
  for (j in (i+1):4) {
    cat(i,":",j,"-->")
    if(j<5){
      sname <- paste0(stele.list[j],".",stele.list[i])
      stele.markers.neg[[sname]] <- FindMarkers(object = mura_har, ident.1 = stele.list[j], ident.2 = stele.list[i])
      
    }
  }
}
stele.lengthing<- list()
for (i in 1:length(stele.markers.neg)) {
  lis <- unlist(strsplit(names(stele.markers)[i],"[.]"))
  sname1 <- paste0(lis[1],"_",names(stele.markers)[i])
  sname2 <- paste0(lis[2],"_",names(stele.markers)[i])
  l1 <-  rownames(stele.markers[[i]][stele.markers[[i]]$avg_logFC>0,])
  l2 <-  rownames(stele.markers[[i]][stele.markers[[i]]$avg_logFC<0,])
  l3 <-  rownames(stele.markers.neg[[i]][stele.markers.neg[[i]]$avg_logFC>0,])
  l4 <-  rownames(stele.markers.neg[[i]][stele.markers.neg[[i]]$avg_logFC<0,])
  stele.lengthing[[sname1]] <- intersect(l1,l4)
  stele.lengthing[[sname2]] <- intersect(l2,l3)
  
}





