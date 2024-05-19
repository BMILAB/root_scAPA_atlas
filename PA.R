library(data.table)
library(Seurat)
library(movAPA)
#load Ryu.3.wt
setwd("~/xwj/opt/Data/T4/zs/PAC/")
files <- list.files(pattern = "RDA")

source1 <- c("6_col0","4_A","4_B","4_C","4_D","4_E","4_F","4_G","4_H","4_I","4_J",
             "6_sc1","6_sc10","6_sc11","6_sc12","6_sc20","6_sc21","6_sc25","6_sc30",
             "6_sc31","6_sc36","6_sc37","6_sc40","6_sc51","6_sc52","6_sc53","6_sc9",
             "3_1","3_2","3_3","5_1","6_tnw1","6_tnw2","2_1")


for (i in 1:length(files)) {
  cat(i)
  load(files[i])
}

ls(pattern = "_APA")

#write.csv(annos,file = "Ryu_3wt.csv")
as.scPAClist<- function(obj,name){
  if(class(obj)[1]=="PACdataset"){
    obj@anno$tagc <- rowSums(obj@counts)
    obj@anno$expnum <- rowSums(obj@counts>0)
    obj@anno$taga <- obj@anno$tagc/obj@anno$expnum
    
    scPAClist <- list(
      "name" = name,
      "counts" = as.data.table(obj@counts,keep.rownames = "PA_id"),
      "colData" = as.data.table(obj@colData,keep.rownames = "PA_id"),
      "anno" = as.data.table(obj@anno,keep.rownames = "PA_id"),
      "supp" = as.data.table(obj@supp,keep.rownames = "PA_id")
    )
    
  }
  return(scPAClist)
}

objList <- list()

objList[[1]] <- as.scPAClist(col0_APA,source1[1])
objList[[2]] <- as.scPAClist(repA_APA,source1[2])
objList[[3]] <- as.scPAClist(repB_APA,source1[3])
objList[[4]] <- as.scPAClist(repC_APA,source1[4])
objList[[5]] <- as.scPAClist(repD_APA,source1[5])
objList[[6]] <- as.scPAClist(repE_APA,source1[6])
objList[[7]] <- as.scPAClist(repF_APA,source1[7])
objList[[8]] <- as.scPAClist(repG_APA,source1[8])
objList[[9]] <- as.scPAClist(repH_APA,source1[9])
objList[[10]] <- as.scPAClist(repI_APA,source1[10])
objList[[11]] <- as.scPAClist(repJ_APA,source1[11])
objList[[12]] <- as.scPAClist(sc_1_APA,source1[12]) 
objList[[13]] <- as.scPAClist(sc_10_at_APA,source1[13]) 
objList[[14]] <- as.scPAClist(sc_11_APA,source1[14]) 
objList[[15]] <- as.scPAClist(sc_12_APA,source1[15]) 
objList[[16]] <- as.scPAClist(sc_30_APA,source1[19]) 
objList[[17]] <- as.scPAClist(sc_31_APA,source1[20]) 
objList[[18]] <- as.scPAClist(sc_37_APA,source1[22]) 
objList[[19]] <- as.scPAClist(sc_40_APA,source1[23]) 
objList[[20]] <- as.scPAClist(sc_51_APA,source1[24]) 
objList[[21]] <- as.scPAClist(sc_9_at_APA,source1[27]) 
objList[[22]] <- as.scPAClist(scPACdsrep13UTR,source1[28])
objList[[23]] <- as.scPAClist(scPACdsRep23UTR,source1[29])
objList[[24]] <- as.scPAClist(scPACdsrep33UTR,source1[30])
objList[[25]] <- as.scPAClist(SRP182008_APA,source1[31])
objList[[26]] <- as.scPAClist(tnw1_APA,source1[32])
objList[[27]] <- as.scPAClist(tnw2_APA,source1[33])
objList[[28]] <- as.scPAClist(whole_root_Control_2APA,source1[34])


rm(list=ls(pattern="_APA"))
rm(list=ls(pattern="3UTR"))


genelist <- c()
combind_anno <- data.table()
for(i in 1:length(objList)){
  
  genelist <- c(genelist,list(objList[[i]]$anno$gene))
  objList[[i]]$anno[,batch_id := i]
  combind_anno <- rbind(combind_anno,objList[[i]]$anno[,.(PA_id,chr,strand,coord,start,end,gene,ftr,ftr_start,ftr_end,gene_stop_codon,batch_id,tagc,expnum,taga)])
}

saveRDS(combind_anno,"combind_anno.rds")


#排序
setkey(combind_anno,gene,start)
utr3_anno <- combind_anno[ftr=="3UTR"]


combind_PAC <- function(annos){
  PAC_id <- c()
  len <- length(table(factor(annos$gene)))
  pb <- txtProgressBar(style=3)
  star_time <- Sys.time()  
  
  for(i in levels(factor(annos$gene))){
    dip <- which(levels(factor(annos$gene)) == i)
    setTxtProgressBar(pb, dip/len)
    
    tmp_anno <- annos[annos$gene == i,]
    ns <- tmp_anno[1,]$start
    ne <- tmp_anno[1,]$end
    kk <- 1
    n_id <- paste0("PAC",i,"_",kk)
    PAC_id <- c(PAC_id,n_id)
    if(nrow(tmp_anno) > 1){
      for (j in 2:nrow(tmp_anno)) {
        if(tmp_anno[j,]$start <= ne){
          #overlap
          #set PAC_id
          n_id <- paste0("PAC",i,"_",kk)
          PAC_id <- c(PAC_id,n_id)
          #reset:ne
          if(tmp_anno[j,]$end > ne){
            ne <- tmp_anno[j,]$end 
          }
        }else{
          # no overlap
          #set PAC_id
          kk <- kk+1
          n_id <- paste0("PAC",i,"_",kk)
          PAC_id <- c(PAC_id,n_id)
          #reset:ne,ns
          ne <- tmp_anno[j,]$end 
          ns <- tmp_anno[j,]$start 
        }
      }
    }
  }
  end_time <- Sys.time() 
  close(pb)
  run_time <- end_time - star_time
  
  return(PAC_id)
}

#This step involves merging PACs from different samples, which takes a long time (mainly due to double-layer nested for loops)
PAC_id <- combind_PAC(utr3_anno)

utr3_anno$PAC_id <- PAC_id

saveRDS(utr3_anno,"utr3_anno.rds")


utr3_anno[,.SD(),by=PAC_id]

utr3_anno$PAC_id
utr3_anno[1:5]

mergeRange <- function(anno){
  anno<-utr3_anno[1:5,]
  if(as.character(anno[1]$strand) == "+"){
    tmp <- anno[which.max(anno$coord),]
    anno$
  }
}

utr3_anno %>% group_by(PAC_id) 

dtlists<- split(utr3_anno,utr3_anno$batch_id)



### Sequence info ----
library(movAPA)
data(PACds)
library("BSgenome.Oryza.ENSEMBL.IRGSP1")
bsgenome <- BSgenome.Oryza.ENSEMBL.IRGSP1
faFiles=faFromPACds(PACds, bsgenome, what='updn', fapre='updn', up=-100, dn=100, byGrp='ftr')
faFiles=c("updn.3UTR.fa", "updn.Ext_3UTR.fa", "updn.intergenic.fa", "updn.intron.fa" )
## plot single nucleotide profile for a fa file
plotATCGforFAfile (faFiles="updn.3UTR.fa", ofreq=TRUE, opdf=TRUE, refPos=301)
## plot multiple fa files
fafiles=faFromPACds(PACds, bsgenome, what='updn', fapre='400nt', up=-300, dn=100, byGrp='ftr', chrCheck=FALSE)
plotATCGforFAfile (faFiles=fafiles, ofreq=TRUE, opdf=TRUE, refPos=301)
## plot multiple fa files into one PDF
fafiles=faFromPACds(PACds, bsgenome, what='updn', fapre='400nt', up=-300, dn=100, byGrp='ftr', chrCheck=FALSE)
plotATCGforFAfile (faFiles=fafiles, ofreq=TRUE, opdf=TRUE, refPos=301, mergePlots=TRUE, filepre='allplots')

sample_utr3 <- utr3_anno[,.SD[sample(.N, min(200,.N))],by = batch_id]



  
ncounts <- list()

#This step is to merge under a single sample, with the aim of merging the counts of some special loci in a single sample
PAC_counts <- list()
for (i in 1:length(objList)) {
  setkey(objList[[i]]$counts,PA_id)
  setkey(dtlists[[i]],PA_id)
  ncounts[[i]]<- objList[[i]]$counts[dtlists[[i]][,.(PA_id,PAC_id)]]
  tmp <- ncounts[[i]][,-c("PA_id")]
  PAC_counts[[i]] <- tmp[,lapply(.SD, sum),by=PAC_id]
}

PAC_counts
for(i in 1:length(PAC_counts)){
  PAC_counts[[1]]
}

saveRDS(PAC_counts,"PAC_counts.rds")


PAC_counts <- readRDS("PAC_counts.rds")

nrow(ss[grepl("sc_12",rownames(ss)),])

#Merge by PAC_id
library(Seurat)

obj.list <- list()

for (i in 1:length(PAC_counts)) {
  cat(i,"\n")
  tmp <- as.data.frame(PAC_counts[[i]])
  rownames(tmp) <- PAC_counts[[i]]$PAC_id
  tmp <- tmp[,-1]
  colnames(tmp) <- paste0(colnames(tmp),"_",i)
  obj <- CreateSeuratObject(counts = tmp, project = paste0("OBJ",i), min.cells = 3, min.features = 200)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nrow(tmp), verbose = FALSE)
  obj.list <- c(obj.list,obj)
}


ggplot(obj.list[[1]]@meta.data, aes(x=nFeature_RNA)) + geom_bar(stat="bin")
ggplot(obj.list[[1]]@meta.data, aes(x=nCount_RNA)) + geom_bar(stat="bin")

ggplot(obj.list[[2]]@meta.data, aes(x=nFeature_RNA)) + geom_bar(stat="bin")
ggplot(obj.list[[2]]@meta.data, aes(x=nCount_RNA)) + geom_bar(stat="bin")

ggplot(obj.list[[3]]@meta.data, aes(x=nFeature_RNA)) + geom_bar(stat="bin")
ggplot(obj.list[[3]]@meta.data, aes(x=nCount_RNA)) + geom_bar(stat="bin")

### harmony integrate DE----
library(SeuratWrappers)
library(Seurat)
library(harmony)
Com_PAcount <- data.table()
table(gene_har@meta.data$sample)
# 2_1    3_1    3_2    3_3    4_A    4_B    4_C    4_D    4_E    4_F    4_G    4_H    4_I    4_J    5_1 6_col0 
# 1076   4406   1473   1643   2169   1672     39    971     30    409    388    148   2009   1622   7690   6779 
# 6_sc1 6_sc10 6_sc11 6_sc12 6_sc30 6_sc31 6_sc37 6_sc40 6_sc51  6_sc9 6_tnw1 6_tnw2 
# 10590  10040   9913  11236  11118   9443   6156   7719   6736   3926   5083   4065 

source1 <- c("6_col0","4_A","4_B","4_C","4_D","4_E","4_F","4_G","4_H","4_I","4_J",
             "6_sc1","6_sc10","6_sc11","6_sc12","6_sc20","6_sc21","6_sc25","6_sc30",
             "6_sc31","6_sc36","6_sc37","6_sc40","6_sc51","6_sc52","6_sc53","6_sc9",
             "3_1","3_2","3_3","5_1","6_tnw1","6_tnw2","2_1")


samples3 <- c("6_sc12","6_sc11","6_sc30","6_sc31","6_sc37","6_sc40",
              "6_sc51","6_sc9","6_sc10","6_sc1","6_tnw1","6_tnw2",
              "6_col0","2_1","3_1","3_2","3_3","4_A","4_B","4_C",
              "4_D","4_E","4_F","4_G","4_H","4_I","4_J","5_1")

names(PAC_counts) <- source1

ys <- source1 %in% samples3
relist <- source1[ys]

PAC_S <- PAC_counts[ys]

for (i in 1:length(PAC_S)) {
  
  ns <- names(PAC_S[i])
  cs <- colnames(PAC_S[[i]])
  df<- gene_har@meta.data
  rs <- gsub("[^A-Z]", "", rownames(df[df$sample == ns,]))
  res <- cs %in% rs
  res[1] <- T
  if((sum(res)-1) != length(rs)){cat(i,"-",ns,":",sum(res)-1,"-",length(rs),"\n")}
  PAC_S[[i]] <- PAC_S[[i]][,res,with=F]
}

l<-0
for (i in 1:length(PAC_S)) {
l <- l+ncol(PAC_S[[i]])
}
#l:128162
#ncol(gene_har) 128549

saveRDS(PAC_S,"PAC_S.rds")
PAC_S <- readRDS("~/xwj/opt/Data/T4/zs/PAC/PAC_S.rds")
APAsum <- 0
for (i in 1:length(PAC_S)) {
  cat(names(PAC_S[i]))
  cat(dim(PAC_S[[i]]))
  APAsum <- nrow(PAC_S[[i]]) + APAsum
}



for (i in 1:length(PAC_S)) {
  colnames(PAC_S[[i]]) <- paste0(colnames(PAC_S[[i]]),"_",i)
  colnames(PAC_S[[i]])[1] <- "PAC_id"
}

Com_PAcount <- merge(PAC_S[[1]],PAC_S[[2]],by="PAC_id",all=T)

for (j in 3:10) {
  cat('-',j)
  Com_PAcount <- merge(Com_PAcount,PAC_S[[j]],by="PAC_id",all=T)
  PAC_S[[j]] <- data.table()
  gc()
}

saveRDS(Com_PAcount,"Com_PAcount.rds")
Com_PAcount <- readRDS("~/zs/opt/Data/T4/zs/PAC/Com_PAcount.rds")
#as(as.matrix(Com_PAcount), "dgCMatrix")

metaref <- gene_har@meta.data
dd <- gsub("[A-Z]","",colnames(obj_har))
table(dd)

dim(Com_PAcount)
#[1]  29784 128135
Com_PAcount[is.na(Com_PAcount)] <- 0

Com_PAc <- as(as.matrix(Com_PAcount[,-1]), "dgCMatrix")
dim(Com_PAc)
#[1]  29784 128134
#Com_PAc <- as.data.frame(Com_PAcount)
rownames(Com_PAc) <- Com_PAcount$PAC_id
saveRDS(Com_PAc,"Com_PAc.rds")

saveRDS(obj_har,"obj_har.rds")
saveRDS(obj_har,"obj_har_sim.rds")

obj_har <- readRDS("~/zs/opt/Data/T4/zs/PAC/obj_har.rds")

obj_har <- CreateSeuratObject(counts = Com_PAc, project = "harmony")
obj_har <- NormalizeData(obj_har, verbose = FALSE)
obj_har <- FindVariableFeatures(obj_har, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj_har <- ScaleData(obj_har, verbose = FALSE)

#combat
library(sva)

edata <- as.matrix(obj_har@assays$RNA@scale.data)
batch <- factor(obj_har$samples)

#combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

combat_var = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
saveRDS(combat_var,"combat_var.rds")

obj_har@assays$RNA@scale.data <- combat_var

obj_har <- RunPCA(obj_har,verbose = FALSE)
obj_har <- RunUMAP(obj_har, reduction = "pca", dims = 1:30,reduction.name = "combat_umap")
Idents(obj_har) <- "samples"
DimPlot(obj_har,reduction = "combat_umap")


###### Annotation ####
samples <- gsub("[A-Z]*_","",rownames(obj_har@meta.data))
obj_har@meta.data$samples <- samples

library(stringr)
gene_har@meta.data$ids2 <- str_extract(gene_har@meta.data$ids, "[A-Z]+")
gene_har@meta.data$os <- 0

for (i in 1:length(relist)) {
  gene_har@meta.data[gene_har@meta.data$sample == relist[i],]$os <- i
}

gene_har@meta.data$ids3 <- paste(gene_har@meta.data$ids2,"_",gene_har@meta.data$os,sep="")

obj_har@meta.data$ids3 <- rownames(obj_har@meta.data)
newmeta <- merge(obj_har@meta.data,gene_har@meta.data,by="ids3",all=F,sort=F)
rownames(newmeta) <- newmeta$ids3
obj_har@meta.data <- newmeta

qc <- obj_har@meta.data[obj_har$final.ID == "Quiescent Center",c("order","celltype.ID.P2","Rad.ID.P2","ici_celltype2","auc.id")]

idxSCN<-qc[rowSums(qc=="Quiescent Center") <= 1,]$order

obj_har$final.ID2 <- as.character(obj_har$final.ID)

obj_har$final.ID2[idxSCN] <- as.character("Stem Cell Niche")

obj_har@meta.data
obj_har$timelv1 <- as.character(obj_har$timelv1)
obj_har$timelv1 <- as.character(gene_har$timelv1[match(obj_har$order,gene_har$order)])

saveRDS(obj_har@meta.data,file="obj_meta.rds")


obj_har <- RunPCA(obj_har,verbose = FALSE)
obj_har <- RunUMAP(obj_har, reduction = "pca", dims = 1:10,reduction.name = "combind_umap")
Idents(obj_har) <- "samples"
DimPlot(obj_har,reduction = "combind_umap",group.by = "source.x")

obj_har <- obj_har %>% 
  RunHarmony("samples", plot_convergence = FALSE)

obj_har <- RunUMAP(obj_har, reduction = "harmony",dims = 1:ncol(obj_har[["harmony"]]),reduction.name = "harmony_umap")

Idents(obj_har) <- "timezone.ID.P.x"
DimPlot(obj_har,reduction = "harmony_umap")

Idents(obj_har) <- "samples"
DimPlot(obj_har,reduction = "harmony_umap")

Idents(obj_har) <- "final.ID2"
obj_har$final.ID2 <- factor(obj_har$final.ID2, levels = order[sort(match(unique(obj_har$final.ID2),order))]) 
color <- palette[sort(match(unique(obj_har$final.ID2),order))]
DimPlot(obj_har,reduction = "harmony_umap",cols = color)

##DE PA-----
Idents(obj_har) <- "final.ID2"
pamarkers <- FindAllMarkers(obj_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(pamarkers,file="pamarkers.rds")

Idents(obj_har) <- "Main.tp"
mainmarkers <- FindAllMarkers(obj_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

##filter DEPAs
filter.pa <- filter.ratio(obj_har,group.by = "final.ID2",delist = pamarkers,pt.1=0.2,pt.2=0.2)
table(filter.pa$cluster)


### APA ratio DE----

#PA ratio + harmony
APAratio <- Cal.APAratio(Com_PAc)
saveRDS(APAratio,"APAratio.rds")
ratio_har <- CreateSeuratObject(counts = APAratio, project = "ratio")
ratio_har <- NormalizeData(ratio_har, verbose = FALSE)
ratio_har <- FindVariableFeatures(ratio_har, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
ratio_har <- ScaleData(ratio_har, verbose = FALSE)


if(identical(rownames(ratio_har@meta.data),rownames(obj_har@meta.data))){
  ratio_har@meta.data <- cbind(ratio_har@meta.data,obj_har@meta.data[,c(-1,-2,-3,-4)])
}

ratio_har <- RunPCA(ratio_har,verbose = FALSE)
ratio_har <- RunUMAP(ratio_har, reduction = "pca", dims = 1:10,reduction.name = "combind_umap")
Idents(ratio_har) <- "samples"
DimPlot(ratio_har,reduction = "combind_umap")

ratio_har <- ratio_har %>% 
  RunHarmony("samples", plot_convergence = FALSE)

ratio_har <- RunUMAP(ratio_har, reduction = "harmony",dims = 1:ncol(ratio_har[["harmony"]]),reduction.name = "harmony_umap")

DimPlot(ratio_har,reduction = "harmony_umap",label = T)

saveRDS(ratio_har,file="ratio_har.rds")


b1<- ggplot(ratio_har@meta.data, aes(x=nFeature_RNA)) + geom_bar(stat="bin")+ggtitle("PAC>0 Number per cell")
b2<- ggplot(ratio_har@meta.data, aes(x=nCount_RNA)) + geom_bar(stat="bin")+ggtitle("PAC count SUM per cell")

plot_grid(b1,b2,labels=c("A","B"))

##DE PA ratio
Idents(ratio_har) <- "final.ID2"
ramarkers <- FindAllMarkers(ratio_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#saveRDS(ramarkers,file="ramarkers.rds")

#topratios <- ratiomarkers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
#write.csv(ratiomarkers,file="PAratioDE.csv")
#FeaturePlot(ratio_har, features = topratios$gene,reduction = "harmony_umap")
save(ratio_har,file="ratio_har.RData")

filter.ra <- filter.ratio(ratio_har,group.by = "final.ID2",delist = ramarkers,pt.1=0.2,pt.2=0.2)
table(filter.ra$cluster)

### APA ratio multiple ----
ls <- gsub("_.*","",rownames(APAratio))
per <- unique(ls[duplicated(ls)])
APAratio_mu <- APAratio[ls %in% per,]

mura_har <- CreateSeuratObject(counts = APAratio_mu, project = "ratio_muli")
mura_har <- NormalizeData(mura_har, verbose = FALSE)
mura_har <- FindVariableFeatures(mura_har, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
mura_har <- ScaleData(mura_har, verbose = FALSE)

if(identical(rownames(mura_har@meta.data),rownames(obj_har@meta.data))){
  mura_har@meta.data <- cbind(mura_har@meta.data,obj_har@meta.data[,c(-1,-2,-3,-4)])
}

mura_har <- RunPCA(mura_har,verbose = FALSE)
mura_har <- RunUMAP(mura_har, reduction = "pca", dims = 1:10,reduction.name = "combind_umap")
Idents(mura_har) <- "samples"
DimPlot(mura_har,reduction = "combind_umap")

mura_har <- mura_har %>% 
  RunHarmony("samples", plot_convergence = FALSE)

mura_har <- RunUMAP(mura_har, reduction = "harmony",dims = 1:ncol(ratio_har[["harmony"]]),reduction.name = "harmony_umap")

DimPlot(mura_har,reduction = "harmony_umap",label = T)

saveRDS(mura_har,file="mura_har.rds")

Idents(mura_har) <- "final.ID2"
mumarkers <- FindAllMarkers(mura_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(mumarkers,file="mumarkers.rds")

sum(table(mumarkers$cluster))


### APA switching ----
muliPAC <- rownames(APA_ratio)[!grepl("_1",rownames(APA_ratio))]
muligenes <- gsub("_.*","",muliPAC)
swAPA <- APA_ratio[gsub("_.*","",rownames(APA_ratio)) %in% muligenes,]

sw_har <- CreateSeuratObject(counts = swAPA, project = "ratio", min.cells = 3, min.features = 200)
sw_har <- NormalizeData(sw_har, verbose = FALSE)
sw_har <- FindVariableFeatures(sw_har, selection.method = "vst", nfeatures = nrow(sw_har), verbose = FALSE)
sw_har <- ScaleData(sw_har, verbose = FALSE)

sw_har@meta.data$barcode <- rownames(sw_har@meta.data)

sw_har@meta.data <- left_join(sw_har@meta.data, subanno, by = "barcode")
rownames(sw_har@meta.data) <- sw_har@meta.data$barcode

sw_har <- RunPCA(sw_har,verbose = FALSE)
sw_har <- RunUMAP(sw_har, reduction = "pca", dims = 1:10,reduction.name = "combind_umap") 
Idents(sw_har) <- "sample"
DimPlot(sw_har,reduction = "combind_umap")

sw_har <- sw_har %>% 
  RunHarmony("sample", plot_convergence = FALSE)

sw_har <- RunUMAP(sw_har, reduction = "harmony",dims = 1:ncol(sw_har[["harmony"]]),reduction.name = "harmony_umap")

DimPlot(sw_har,reduction = "harmony_umap")

Idents(sw_har) <- "celltype"
swmarkers <- FindAllMarkers(sw_har, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topsw <- swmarkers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

write.csv(swmarkers,file="PAswDE.csv")

FeaturePlot(sw_har, features = topsw$gene,reduction = "harmony_umap")
save(sw_har,file="sw_har.RData")

### gene DE----
library(Seurat)

load("F:/NewProject/scAPAmarker/reference/3-Ryu 2019/count/Ryu.RData")
gcount <- as.data.frame(Ryu$count)
rownames(gcount) <- gcount[,1]
gcount <- gcount[,-1]
gcount <- gcount[,grepl("WT.*",colnames(gcount))]
dim(gcount)

gs <- data.table(x=colnames(gcount))
gs <- gs[, c("c1", "c2","c3") := tstrsplit(x, "_", fixed=TRUE)][]

gsn <- gsub("WT","",paste0(gs$c3,"_",gs$c2))

colnames(gcount) <- gsn

#transfer to ensembel
library(data.table)
ensgene1 <- read.csv("F:/NewProject/scAPAmarker/reference/3-Ryu 2019/count/Sample_WT-WERGFP/filtered_gene_bc_matrices/TAIR10/genes.tsv",sep="\t",header = F)
ensgene2 <- read.csv("F:/NewProject/scAPAmarker/reference/3-Ryu 2019/count/Sample_WT-WERGFP_2/filtered_gene_bc_matrices/TAIR10/genes.tsv",sep="\t",header = F)
ensgene3 <- read.csv("F:/NewProject/scAPAmarker/reference/3-Ryu 2019/count/Sample_WT-WERGFP_3/filtered_gene_bc_matrices/TAIR10/genes.tsv",sep="\t",header = F)

nm <- match(rownames(gcount), ensgene1$V2, nomatch = NA_integer_, incomparables = NULL)
nm <- nm[!is.na(nm)]
ngn <- ensgene1[nm,]$V1
gcounts <- gcount[rownames(gcount) %in% ensgene1$V2,]
rownames(gcounts) <- ngn

geneM <- CreateSeuratObject(counts = gcounts, project = "Gene", min.cells = 3, min.features = 200)
geneM <- NormalizeData(geneM, verbose = FALSE)
geneM <- FindVariableFeatures(geneM, selection.method = "vst", nfeatures = nrow(geneM), verbose = FALSE)
geneM <- ScaleData(geneM, verbose = FALSE)

geneM@meta.data$barcode <- rownames(geneM@meta.data)

geneM@meta.data <- left_join(geneM@meta.data, subanno, by = "barcode")
rownames(geneM@meta.data) <- geneM@meta.data$barcode

geneM <- RunPCA(geneM,verbose = FALSE)
geneM <- RunUMAP(geneM, reduction = "pca", dims = 1:10,reduction.name = "combind_umap") 
Idents(geneM) <- "sample"
s1 <- DimPlot(geneM,reduction = "combind_umap")

geneM <- geneM %>% RunHarmony("sample", plot_convergence = FALSE)

geneM <- RunUMAP(geneM, reduction = "harmony",dims = 1:ncol(geneM[["harmony"]]),reduction.name = "harmony_umap")

s2 <- DimPlot(geneM,reduction = "harmony_umap")
plot_grid(s1,s2,labels=c("A","B"))
Idents(geneM) <- "celltype"
DimPlot(geneM,reduction = "harmony_umap",label = T)

d1<- ggplot(geneM@meta.data, aes(x=nFeature_RNA)) + geom_bar(stat="bin")+ggtitle("PAC>0 Number per cell")
d2<- ggplot(geneM@meta.data, aes(x=nCount_RNA)) + geom_bar(stat="bin")+ggtitle("PAC count SUM per cell")

plot_grid(d1,d2,labels=c("A","B"))

geneMmarkers <- FindAllMarkers(geneM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topgeneM <- geneMmarkers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

write.csv(geneMmarkers,file="geneDE.csv")

FeaturePlot(geneM, features = topgeneM$gene,reduction = "harmony_umap")
save(geneM,file="geneM.RData")

#### Hair vs Non-Hair-----

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


### 绘图 ----

Atrichoblast_lengthening <- read.csv("D:/project/rootmarker/paper/Atrichoblast_lengthening.csv")
Atrichoblast_shorting <- read.csv("D:/project/rootmarker/paper/Atrichoblast_shorting.csv")

library(clusterProfiler)

res <- enrichplot::dotplot(Atrichoblast_lengthening,orderBy = "x")
res$data$Description

num.rep <- 10
go_L <- Atrichoblast_lengthening[1:num.rep,]
go_S <- Atrichoblast_shorting[1:num.rep,]


go_enrich_df<-data.frame(ID=c(go_L$ID, go_S$ID),
                         Description=c(as.character(go_L$Description), 
                                       as.character(go_S$Description)),
                         GeneNumber=c(go_L$Count, go_S$Count),
                         type=factor(c(rep("Lengthening", num.rep), 
                                       rep("Shorting", num.rep)),
                                     levels=c("Lengthening", 
                                              "Shorting")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=25){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > n_char))
  {
    if (nchar(x) > n_char) x <- substr(x, 1, n_char)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}
go_enrich_df$Description <- factor(go_enrich_df$Description)
labels=(sapply(levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
               shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))

go_enrich_df$Description <- factor(go_enrich_df$Description)
CPCOLS <- c("#8DA1CB", "#ffa020", "#66C3A5","#FD8D62","#00a020")


ggplot(data=go_enrich_df,aes(x=number, y=GeneNumber,fill=type))+
  ggplot2::geom_bar(stat="identity",position="dodge")+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust=1,size=10),plot.margin = unit(c(1,1,2,4),"lines"),
        panel.grid =element_blank())+
  scale_x_discrete(labels=labels)+
  scale_fill_manual(values = CPCOLS)+
  facet_wrap(~type,scales = "free_x",ncol=5)+xlab(NULL)+ylab("Gene counts")+NoLegend()

library(eoffice)
topptx(filename = "D:/project/rootmarker/paper/ALS.pptx", width = 7.5,height = 4)



### overlapping ----

PAratioDE <- read.csv("PAratioDE.csv")
PAratioDE$gene <- gsub("PAC","",PAratioDE$gene)
PAratioDE$gene <- gsub("-.*","",PAratioDE$gene)

PAcountDE <- read.csv("PAcountDE.csv")
PAcountDE <- pmarkers
PAcountDE$gene <- gsub("PAC","",PAcountDE$gene)
PAcountDE$gene <- gsub("-.*","",PAcountDE$gene)



geneDE <- read.csv("geneDE.csv")
allgenes <- rownames(geneM)
allgenes <- rownames(gene_sub)

PAswDE <-  read.csv("PAswDE.csv")
PAswDE$gene <- gsub("PAC","",PAswDE$gene)
PAswDE$gene <- gsub("-.*","",PAswDE$gene)

celltypes <- levels(factor(PAcountDE$cluster))
REAGenes <- allgenes[!(allgenes %in% geneDE$gene)]

overlap_genes <- list()
overlap_num <- data.frame()


for(i in 1:length(celltypes)){
  tmpratio <- subset(PAratioDE,cluster==celltypes[i])
  tmpcount <- subset(PAcountDE,cluster==celltypes[i])
  tmpsw <- subset(PAswDE,cluster==celltypes[i])
  tmpgene <- subset(geneDE,cluster==celltypes[i])
  ReGenes <- allgenes[!(allgenes %in% tmpgene$gene)]
  
  GC <- intersect(tmpgene$gene,tmpcount$gene)
  GR <- intersect(tmpgene$gene,tmpratio$gene)
  GS <- intersect(tmpgene$gene,tmpsw$gene)
  RC <- intersect(tmpratio$gene,tmpcount$gene)
  RS <- intersect(tmpratio$gene,tmpsw$gene)
  CS <- intersect(tmpsw$gene,tmpcount$gene)
  
  NGC <- intersect(ReGenes,tmpcount$gene)
  NGR <- intersect(ReGenes,tmpratio$gene)
  NGS <- intersect(ReGenes,tmpsw$gene)
  
  NAGC <- intersect(REAGenes,tmpcount$gene)
  NAGR <- intersect(REAGenes,tmpratio$gene)
  NAGS <- intersect(REAGenes,tmpsw$gene)
  
  overlap_genes[[celltypes[i]]] <- list(GC=GC,GR=GR,RC=RC,RS =RS,CS = CS,NGC=NGC,NGR=NGR,NGS=NGS,NAGC=NAGC,NAGR=NAGR,NAGS=NAGS)
  
  overlap_num <- rbind(overlap_num,data.frame(celltype=celltypes[i],
                                              GC=length(unique(GC)),
                                              GR=length(unique(GR)),
                                              RC=length(unique(RC)),
                                              RS=length(unique(RS)),
                                              CS=length(unique(CS)),
                                              NGC=length(unique(NGC)),
                                              NGR=length(unique(NGR)),
                                              NGS=length(unique(NGS)),
                                              NAGC=length(unique(NAGC)),
                                              NAGR=length(unique(NAGR)),
                                              NAGS=length(unique(NAGS))))
}

for(i in 1:length(celltypes)){
  tmpcount <- subset(PAcountDE,cluster==celltypes[i])
  tmpgene <- subset(geneDE,cluster==celltypes[i])
  ReGenes <- allgenes[!(allgenes %in% tmpgene$gene)]
  
  GC <- intersect(tmpgene$gene,tmpcount$gene)
  NGC <- intersect(ReGenes,tmpcount$gene)
  NAGC <- intersect(REAGenes,tmpcount$gene)
  
  overlap_genes[[celltypes[i]]] <- list(GC=GC,NGC=NGC,NAGC=NAGC)
  
  overlap_num <- rbind(overlap_num,data.frame(celltype=celltypes[i],
                                              GC=length(unique(GC)),
                                              NGC=length(unique(NGC)),
                                              NAGC=length(unique(NAGC))
                                              ))
}


tmpcox1 <- overlap_genes$Stele$RC


ss <- PAratioDE %>% group_by(cluster) %>% arrange(desc(avg_logFC))
ss <- ss[ss$cluster == "Stele",]
ns <- data.table(gene = tmpcox1,id=match(tmpcox1,ss$gene))
nss <- arrange(ns,id)$gene


tmpcox2 <- paste0("PAC",nss,"-1")

s1 <- FeaturePlot(geneM, features = nss,reduction="harmony_umap")
s1 <- FeaturePlot(obj_har, features = tmpcox2[4],reduction="harmony_umap")

s2 <- FeaturePlot(ratio_har, features = tmpcox2[4],reduction="harmony_umap")

SA <- data.frame()

for (i in 1:length(s1)) {
  SA <- rbind(SA,data.frame(gene = sd(s1[[i]]$data[,4]),PA = sd(s2[[i]]$data[,4])))
  
}
nn <- which.max(SA[,2]-SA[,1])


nn
s3 <- FeaturePlot(geneM, features = nss[8],reduction="harmony_umap",label = T)
s3s <- FeaturePlot(geneM, features = nss[1],reduction="harmony_umap",label = T)+ggtitle("PACAT2G25710-1")
s4 <- FeaturePlot(ratio_har, features = tmpcox2[8],reduction="harmony_umap",label = T)

s3s$data[,4] <- s4$data[,4]

plot_grid(s3,s3s,labels=c("A","B"))

VlnPlot(geneM, features = c("AT2G25710"))
VlnPlot(ratio_har, features = c("PACAT2G25710-1"))


########Correlation analysis of expression levels-------
library(ggplot2)

gene_count <- gene_har@assays$RNA@counts
dim(gene_count)
gene_count_sub <- gene_count[,colnames(gene_count) %in% obj_har@meta.data$ids]

count1 <- rowSums(gene_count_sub)
rm(gene_count)
rm(gene_count_sub)
gc()

pa_count <- obj_har@assays$RNA@counts
count2 <- rowSums(pa_count)
rm(pa_count)
gc()
newn <- gsub("-[0-9]","",gsub("PAC","",names(count2)))
dtpa <- data.table::data.table(pa=count2,gene=newn)
dtpa2 <- dtpa[,lapply(.SD, sum),by=gene]
identical(dtpa2$gene,names(count1))
ings <- intersect(names(count1),dtpa2$gene)
counts1 <- count1[ings]
idx <- match(ings,dtpa2$gene)
dtpa3 <- dtpa2[idx,]
identical(dtpa3$gene,names(counts1))
dtpa3$ge <- counts1

pa.ge.cor <- ggplot(data=dtpa3, aes(x=pa, y=ge)) + geom_point(color="blue")+
  stat_smooth(method="lm",se=FALSE,size=0.5)+ labs(x="scPA",y="scGE")#se=FALSE：不添加置信区间
  geom_rug()

  
###Time development sequence----
  mura_har@assays$RNA@scale.data <- matrix()
  
  rownames(mura_har@assays$RNA@counts)
  mura_har@assays$RNA@counts[1:3,][,mura_har@assays$RNA@counts[2,] >0]
  ppui.mtx <- mura_har@assays$RNA@counts[grepl("-1",rownames(mura_har)),]
  saveRDS(ppui.mtx,"ppui.mtx.rds")
  
  obj_har_sim@meta.data$nPPUI <- colSums(ppui.mtx>0)
  obj_har_sim@meta.data$PPUIsum <- colSums(ppui.mtx)
  obj_har_sim@meta.data$PPUIavg <- obj_har_sim@meta.data$PPUIsum / obj_har_sim@meta.data$nPPUI
  obj_har_sim@meta.data$nPPUImax <- colSums(ppui.mtx == 1)
  
  DimPlotOrder(obj_har_sim,reduction = "harmony_umap",group.by = "nPPUI",order=order,palette=palette)
  DimPlot(obj_har_sim, reduction = "harmony_umap", group.by = "nPPUI",raster=FALSE)+NoLegend()
  DimPlot(obj_har_sim, reduction = "harmony_umap", group.by = "nPPUI",raster=FALSE)+NoLegend()

  p1 <- FeaturePlot(object = obj_har_sim, reduction = "harmony_umap", features = 'nPPUI')
  p1$data$UMAP_1 <- obj_har_sim@reductions$harmony_umap@cell.embeddings[,1]
  p1$data$UMAP_2 <- obj_har_sim@reductions$harmony_umap@cell.embeddings[,2]
  p1  
  
  p2 <- FeaturePlot(object = obj_har_sim, reduction = "harmony_umap", features = 'PPUIsum')
  p2$data$UMAP_1 <- obj_har_sim@reductions$harmony_umap@cell.embeddings[,1]
  p2$data$UMAP_2 <- obj_har_sim@reductions$harmony_umap@cell.embeddings[,2]
  p2  
  
  p3 <- FeaturePlot(object = obj_har_sim, reduction = "harmony_umap", features = 'nPPUImax',raster=FALSE)
  p3$data$UMAP_1 <- obj_har_sim@reductions$harmony_umap@cell.embeddings[,1]
  p3$data$UMAP_2 <- obj_har_sim@reductions$harmony_umap@cell.embeddings[,2]
  p3
  

  dt <- data.frame(RUD_avg = obj_har_sim@meta.data$PPUIsum,barcode = rownames(obj_har_sim@meta.data))
  
  layer <- obj_har_sim@reductions$harmony_umap@cell.embeddings
  layer <- as.data.frame(layer)
  p <- ggplot(data = dt,aes(x = as.numeric(layer$UMAP_1),y = as.numeric(layer$UMAP_2)))+
    geom_point(aes(color=RUD_avg),size=0.5,shape=16)+
    scale_color_distiller(palette = "Spectral")+
    labs(x = "UMAP_1", y = "UMAP_2",color="RUD_avg")+
    theme_bw()+theme(panel.grid =element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     legend.text = element_text(size=12),
                     legend.title  = element_blank(),
                     legend.position = "right")
  
  
#Check correlation coefficient
  obj_har_sim@meta.data$time.num <- as.numeric(as.factor(obj_har_sim@meta.data$timelv1))
  
  library(rstatix) 
  library(tidyr)
  
  cor.test(obj_har_sim@meta.data$nPPUI,obj_har_sim@meta.data$time.num,method = "pearson")
  
  library(ggpubr)
  obj_har_sim@meta.data$timelv1
  
  select_df <- obj_har_sim@meta.data[,c("Main.tp","final.ID2","nPPUI","PPUIsum","nPPUImax","time.num","timelv1")]
  
  select_df$timelv1 <- factor(select_df$timelv1,levels = c("Meristem","Elongation","Maturation"))
  
  subcelltype <- unique(select_df$final.ID2)
  select_df1 <- select_df[select_df$final.ID2 == subcelltype[1],]
  
  table(select_df$final.ID2)
  table(select_df$Main.tp)
  
  ggboxplot(select_df, x = "Main.tp", y = "nPPUI")
  
  select_df1 %>% anova_test(nPPUI ~ time.num)
  table(select_df[select_df$PPUIsum<300,]$Main.tp)
  table(select_df[select_df$PPUIsum>300 & select_df$PPUIsum<600,]$Main.tp)
  table(select_df[select_df$PPUIsum>600,]$Main.tp)
  
  select_df$nPPUI
  

  d <- density(select_df$nPPUI) # returns the density data
  plot(d) # plots the results
  
  library(ggplot2)
  # Basic violin plot
  select_df$Main.tp <- as.character(select_df$Main.tp)
  select_df[select_df$Main.tp %in% c("Quiescent Center","Stem Cell Niche"),]$Main.tp <- "Stem Cell"
  select_df$Main.tp <- factor(select_df$Main.tp, levels= c("Stem Cell","Root cap","Stele",
                                                           "Ground tissue","Epidermis"),
                              labels = c("Stem Cell","Root cap","Stele",
                                         "Ground tissue","Epidermis"))
  
  
  p2 <- ggplot(select_df, aes(x=Main.tp, y=nPPUI, fill=Main.tp)) + 
    geom_boxplot(alpha=0.3) +
    theme(legend.position="none") +
    scale_fill_brewer(palette="Dark2")+theme_bw()+
    labs(x = "Tissue", y = "PPUI Sum",color="Main.tp")+
    theme_bw()+theme(panel.grid =element_blank(),
                     plot.title = element_text(hjust = 0.5),
                     legend.text = element_text(size=12),
                     legend.title  = element_blank(),
                     legend.position = "right")
  
  ##Combining gene expression and APA profiles to refine cell identity ------
  
  BiocManager::install("AUCell")
  
  library(GSEABase)
  library(AUCell)
  ## QC
  # find markers
  
  # label 1 
  ###markers based on pa count ----
  #pamarkers=readRDS("G:/plant roots scRNA-seq论文整理/data/BXY-20240322/filter.pa.rds")
  pacounts <- pa_har@assays$RNA@data
  dim(pacounts)#29784 128134
  
  pa_qc <- pacounts[,pa_har$final.ID2 == "Quiescent Center"| pa_har$final.ID2 == "Stem Cell Niche"]
  dim(pa_qc)#29784  1483
  
  markers_pa <- pamarkers[pamarkers$cluster == "Stem Cell Niche" | pamarkers$cluster == "Quiescent Center",]
  table(markers_pa$cluster)
  #QC 540;STC 363
  pa.list <- split(markers_pa$gene, factor(markers_pa$cluster))
  
  ##构建marker基因集
  library(GSEABase)
  all.sets <- lapply(names(pa.list), function(x) {
    GeneSet(pa.list[[x]], setName=x)        
  })
  #Use GeneSetCollection to construct a collection of gene sets from GeneSet arguments, or a list of GeneSets.
  all.sets <- GeneSetCollection(all.sets)
  
  
  library(AUCell)
  Wt <- as.matrix(pa_qc)
  #Construct gene expression ranking in each cell
  rankings <- AUCell_buildRankings(Wt,plotStats=FALSE, verbose=FALSE)
  #Calculates the 'AUC' for each gene-set in each cell
  #Calculate the activity level of different gene sets in each cell, that is, use the area under the curve (AUC) to calculate whether the input gene sets are enriched in the expressed genes of each cell
  cell.aucs <- AUCell_calcAUC(all.sets, rankings,aucMaxRank=nrow(rankings)*0.05,verbose=FALSE)
  results <- t(cell.aucs@assays@data$AUC)
  #Identify the cell type corresponding to the more active markers gene set in each cell, and re annotate the cell type accordingly
  new.labels1 <- colnames(results)[max.col(results)]
  
  ##label2
  ####markers based on APA ratio ----
  mucounts <- mura_har@assays$RNA@data
  mu_qc <- mucounts[,mura_har$final.ID2 == "Quiescent Center"| mura_har$final.ID2 == "Stem Cell Niche"]
  dim(mu_qc)#11211  1483
  
  markers_mu <- mumarkers[mumarkers$cluster == "Stem Cell Niche" | mumarkers$cluster == "Quiescent Center",]
  table(markers_mu$cluster)
  
  
  mu.list <- split(markers_mu$gene, factor(markers_mu$cluster))
  
  
  all.sets <- lapply(names(mu.list), function(x) {
    GeneSet(mu.list[[x]], setName=x)        
  })
  all.sets <- GeneSetCollection(all.sets)
  
  library(AUCell)
  
  Wt <- as.matrix(mu_qc)
  rankings <- AUCell_buildRankings(Wt,plotStats=FALSE, verbose=FALSE)
  cell.aucs <- AUCell_calcAUC(all.sets, rankings,aucMaxRank=nrow(rankings)*0.05,verbose=FALSE)
  results <- t(cell.aucs@assays@data$AUC)
  new.labels2 <- colnames(results)[max.col(results)]
  
  ##New.labels1: According to the QC and SCN annotated by pa markers; New.labels2: QC, SCN according to APA ratio markers annotations
  label.qc <- data.frame(pa=new.labels1,mu=new.labels2)
  label.qc$ge <- mura_har$final.ID2[mura_har$final.ID2 == "Quiescent Center"| mura_har$final.ID2 == "Stem Cell Niche"]
  #label.qc[label.qc$pa == label.qc$mu,]$final <- label.qc[label.qc$pa == label.qc$mu,]$pa
  
  ##function:getmode ----
  #Calculate the element with the highest number of repetitions in the vector
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv <- uniqv[!is.na(uniqv)]
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  
  label.qc$final <- apply(label.qc,1,getmode)
  
  table(label.qc$final)
 
 
  ##combind gene expression and APA profiles to annotate unknown cells ----
  
  # label 1 
  pacounts <- pa_har@assays$RNA@data
  table(is.na(pa_meta$le))
  #FALSE   TRUE 
  #122077   6057 
  
  
  pa_qc <- pacounts[,is.na(pa_har$le)]
  #29784  6057
  
  markers_pa <- pamarkers
  pa.list <- split(markers_pa$gene, factor(markers_pa$cluster))
  
  #Constructing pa markets sets in different cell types
  library(GSEABase)
  all.sets <- lapply(names(pa.list), function(x) {
    GeneSet(pa.list[[x]], setName=x)        
  })
  all.sets <- GeneSetCollection(all.sets)
  
  library(AUCell)
  
  Wt <- as.matrix(pa_qc)
  #Calculate the ranking of gene expression in each cell
  rankings <- AUCell_buildRankings(Wt,plotStats=FALSE, verbose=FALSE)
  #Calculate the activity level of different gene sets in each cell, that is, use the area under the curve (AUC) to calculate whether the input gene sets are enriched in the expressed genes of each cell
  cell.aucs <- AUCell_calcAUC(all.sets, rankings,aucMaxRank=nrow(rankings)*0.05,verbose=FALSE)
  results <- t(cell.aucs@assays@data$AUC)
  head(results)
  new.labels1 <- colnames(results)[max.col(results)]
  
  ##label 2
  mucounts <- mura_har@assays$RNA@data
  table(is.na(mura_har$le))
  
  mu_qc <- mucounts[,is.na(mura_har$le)]# 11211  6057
  
  markers_mu <- mumarkers
  mu.list <- split(markers_mu$gene, factor(markers_mu$cluster))
  
  #Constructing APA ratio markers sets for different cell types
  library(GSEABase)
  all.sets <- lapply(names(mu.list), function(x) {
    GeneSet(mu.list[[x]], setName=x)        
  })
  all.sets <- GeneSetCollection(all.sets)
  
  library(AUCell)
  
  Wt <- as.matrix(mu_qc)
  rankings <- AUCell_buildRankings(Wt,plotStats=FALSE, verbose=FALSE)
  cell.aucs <- AUCell_calcAUC(all.sets, rankings,aucMaxRank=nrow(rankings)*0.05,verbose=FALSE)
  results <- t(cell.aucs@assays@data$AUC)
  head(results)
  new.labels2 <- colnames(results)[max.col(results)]
  
  #New.labels1: According to the cell type annotated by pa markers; New.labels2: Cell types annotated based on APA ratio markers
  label.qc <- data.frame(pa=new.labels1,mu=new.labels2)
  label.qc$ge <- mura_har$final.ID2[is.na(mura_har$le)]
  #label.qc<- label.qc[,c(3,1,2)]
  label.qc$newfinal <- apply(label.qc,1,getmode)
  
  table(label.qc$newfinal)
  table(label.qc$ge)  