
#' @description  Calculating the correlation coefficient of each cell to each reference expression profile and annotate the cell as the label 
#' that it has the highest correlation coefficient with.  
#' @param method  One of "pearson", "kendall", or "spearman"
#' 
#' Calculate the correlation coefficient between each cell and each reference expression profile,
#' And label the cell with the highest correlation coefficient
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


FoldNames <- function(ids,obj){
  names(ids) <- levels(obj@active.ident)
  obj <- RenameIdents(object = obj, ids)
  return(as.character(obj@active.ident))
}


SYMtoEN.tari <- function(names){
  
  library(biomaRt)
  
  mart = useMart(biomart="plants_mart",host="plants.ensembl.org")
  ensembl <- useDataset("athaliana_eg_gene",mart= mart)
  
  res <- getBM(attributes=c("external_gene_name",'ensembl_gene_id'),
                 filters=c("external_gene_name"),
                 values = names,
                 mart = ensembl)
  return(res)
}

DimPlotOrder <- function(object,reduction,group.by,order=NULL,palette=NULL,size=NULL){
  groups <- object@meta.data[,group.by]
  if(!is.null(order)){
    object@meta.data[,group.by] <- factor(groups, levels = order[sort(match(unique(groups),order))]) 
  }
  if(!is.null(palette)){
    color <- palette[sort(match(unique(groups),order))]
    figure <- DimPlot(object, reduction = reduction, group.by = group.by, cols = color,pt.size=size,raster=FALSE)
  }else{
    figure <- DimPlot(object, reduction = reduction, group.by = group.by,pt.size=size,raster=FALSE)
  }
  return(figure)
}


##get APA ratio matrix
ratio <- function(x){
  if(sum(x)!=0){
    round(x/sum(x),3)
  }else{
    x
  }
}

Cal.APAratio <- function(x){
  require(data.table)
  cat("Calculate the APA ratio")
  pb <- txtProgressBar(style=3)
  ncols <- ncol(x)
  iter <- floor(ncols/1000)+1
  star_time <- Sys.time()
  APA_ratio <- matrix(nrow = nrow(x))
  for (i in 1:iter) {
    setTxtProgressBar(pb, i/iter)
    start <- (i-1)*1000+1
    end <- start+999
    if(end>ncols) end=ncols
    PAC.mtx <- as.data.table(x[,start:end],keep.rownames = "PAC_id")
    PAC.mtx[,Genes:= gsub("_.*","",PAC_id)]
    PAC.mtx <- PAC.mtx[,-1]
    ratio.tmp <- PAC.mtx[,lapply(.SD,ratio),by=Genes]
    ratio.tmp<- ratio.tmp[,-1]
    ratio.tmp <- as(ratio.tmp,"Matrix")
    APA_ratio <- cbind(APA_ratio,ratio.tmp)
  }
  end_time <- Sys.time()
  close(pb)
  cat(end_time - star_time,"\n")
  rownames(APA_ratio) <- rownames(x)
  APA_ratio <- APA_ratio[,-1]
  return(APA_ratio)
}

## filter by expressed ratio
filter.ratio <- function(object,group.by,delist,pt.1=0.15,pt.2=0.15){
  gbys<- object@meta.data[,group.by]
  li <- names(table(gbys))
  count <- object@assays$RNA@counts
  res <- c()
  pb <- txtProgressBar(style=3)
  for(i in 1: length(li)){
    setTxtProgressBar(pb, i/length(li))
    markers_sub <- delist[delist$cluster == li[i],]
    g1 <- colnames(object)[ gbys %in% li[i]]
    g2 <- colnames(object)[!(gbys %in% li[i])]
    p1 <- (rowSums(count[markers_sub$gene,g1]>0) / length(g1)) > pt.1
    p2 <- (rowSums(count[markers_sub$gene,g2]>0) / length(g2)) < pt.2
    ys <- markers_sub$gene %in% names(p1[p1 & p2])
    res <- c(res,ys)
  }
  close(pb)
  return(delist[res,])
}


# main celltype

mainAnno <- function(gene_har){
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
  return(gene_har)
}


Overlapping<- function(x,y){
  cty <- unique(x$cluster)
  interes <- list()
  for(i in 1:length(cty)){
    ref <- x[x$cluster == as.character(cty[i]),]$gene
    con <- y[y$cluster == as.character(cty[i]),]$gene
    ref <- gsub("PA:","",ref)
    ref <- gsub(":.*","",ref)
    con <- gsub("PA:","",con)
    con <- gsub(":.*","",con)
    int <- intersect(ref,con)
    k <- as.character(cty[i])
    interes[[k]] <- int
  }
  return(interes)
}



getmode <- function(v) {
  uniqv <- unique(v)
  uniqv <- uniqv[!is.na(uniqv)]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}