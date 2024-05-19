library(Seurat)
library(movAPA)
library(ggplot2)
library(cowplot)
library(stringr)
library(eoffice)
library(dplyr)
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



BiocManager::install("org.At.tair.db") 
library(org.At.tair.db)
library(clusterProfiler)
library(dplyr)
library(openxlsx)
#pa.marker=read.xlsx("G:/plant roots scRNA-seq论文整理/TableSuppl(20240315).xlsx",sheet = 2,startRow = 2,colNames = TRUE)
pa.marker=read.csv("filter.pamarkers.csv")
pa.marker$gene=substring(pa.marker$PAC.id,4)
pa.marker$gene=gsub('-[0-9]',"",pa.marker$gene)
keytypes(org.At.tair.db)
symbol.gene<- bitr(pa.marker$gene, fromType = "TAIR", 
            toType = c("SYMBOL", "TAIR"), 
            OrgDb = org.At.tair.db,
            )

#symbol.gene=symbol.gene[!duplicated(symbol.gene$TAIR),]

loc=match(symbol.gene$TAIR,pa.marker$gene)
symbol.gene$celltype=pa.marker$cluster[loc]
symbol.gene$PAC.id=pa.marker$PAC.id[loc]

#gene count ----
#12191gene;128549cells
gene_meta=readRDS("gene_meta.rds")
gene_har=readRDS("gene_har_com.rds")


id=match(gene_har@meta.data$ids,gene_meta$ids)
gene_har@meta.data$final.ID2=gene_meta$final.ID2[id]
gene_har@meta.data$Main.tp=gene_meta$Main.tp[id]
gene_har@meta.data$timelv1=gene_meta$timelv1[id]



##******Figure 1 **********####
#Idents(gene_har) <- "final.ID2"
#gene_har$final.ID2 <- factor(gene_har$final.ID2, levels = order[sort(match(unique(gene_har$final.ID2),order))]) 
#color <- palette[sort(match(unique(gene_har$final.ID2),order))]
#DimPlot(gene_har,reduction = "combat_umap",group.by = "final.ID2",cols=color,raster=FALSE)

####source
p1=DimPlotOrder(gene_har,reduction = "combat_umap",group.by = "source")+labs(title = "")


####tissue
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

p2=DimPlotOrder(gene_har,reduction = "combat_umap",group.by = "Main.tp",order=order,palette=palette)+labs(title = "")

####celltype
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
p3=DimPlotOrder(gene_har,reduction = "combat_umap",group.by = "final.ID2",order=order,palette=palette)+labs(title = "")

#####time
order <- c("Maturation","Elongation","Meristem")
palette <- c("#d71345","#f391a9","#feeeed")
p4=DimPlotOrder(gene_har,reduction = "combat_umap",group.by = "timelv1",order=order,palette=palette)+labs(title = "")

library(cowplot)
#plot_grid(p1,p2,p3,p4,labels=c("A","B","C","D"),ncol = 2)
png(filename = "fig1_gene.count_unap.png", width = 12, height = 8,units = "in",res=600)
plot_grid(p1,p2,p3,p4,ncol = 2)
dev.off()



###PA count ----
#29784pa;128134cells
pa_har=readRDS("obj_har.rds")
colnames(pa_har@meta.data)

id=match(pa_har@meta.data$ids3,gene_meta$ids3)
pa_har@meta.data$source.x=gene_meta$source[id]
pa_har@meta.data$ids=gene_meta$ids[id]
pa_har@meta.data$le=gene_meta$le[id]
pa_har@meta.data$final.ID2=gene_meta$final.ID2[id]
pa_har@meta.data$Main.tp=gene_meta$Main.tp[id]
pa_har@meta.data$timelv1=gene_meta$timelv1[id]



##******Figure 2 *******####
####source
png( filename = "fig2_source.pa.png", width = 7, height = 6.5,units = "in",res=600)
p1=DimPlotOrder(pa_har,reduction = "harmony_umap",group.by = "source.x")+labs(title = "")
dev.off()

####tissue
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

png( filename = "fig2_source.tisue.png", width = 7, height = 6.5,units = "in",res=600)
DimPlotOrder(pa_har,reduction = "harmony_umap",group.by = "Main.tp",order=order,palette=palette)+labs(title = "")
dev.off()

####celltype
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

png( filename = "fig2_celltype.pa.png", width = 7, height = 6.5,units = "in",res=600)
p2=DimPlotOrder(pa_har,reduction = "harmony_umap",group.by = "final.ID2",order=order,palette=palette)+labs(title = "")
#fig.ct.ra <- DimPlotOrder(ratio_har,reduction = "harmony_umap",group.by = "final.ID2",order=order,palette=palette)
dev.off()

####time
order <- c("Maturation","Elongation","Meristem")
palette <- c("#d71345","#f391a9","#feeeed")

png( filename = "fig2_time.pa.png", width = 7, height = 6.5,units = "in",res=600)
p3=DimPlotOrder(pa_har,reduction = "harmony_umap",group.by = "timelv1",order=order,palette=palette)+labs(title = "")
dev.off()

png(filename = "fig2_pa.count_unap.png", width = 12, height = 8,units = "in",res=600)
plot_grid(p1,p2,p3,ncol = 2)
dev.off()

##pa count and gene count correlation ----

library(ggplot2)
library(psych)
library(ggpubr)

#obj_har<-readRDS('obj_har.rds')

gene_count <- gene_har@assays$RNA@counts
dim(gene_count)
gene_count_sub <- gene_count[,colnames(gene_count) %in% pa_har@meta.data$ids]

count1 <- rowSums(gene_count_sub)
rm(gene_count)
rm(gene_count_sub)
gc()

pa_count <- pa_har@assays$RNA@counts
count2 <- rowSums(pa_count)
rm(pa_count)
gc()

newn <- gsub("-[0-9]","",gsub("PAC","",names(count2)))
head(newn)
dtpa <- data.table::data.table(pa=count2,gene=newn)

dtpa2 <- dtpa[,lapply(.SD, sum),by=gene]

identical(dtpa2$gene,names(count1))
ings <- intersect(names(count1),dtpa2$gene)#12127genes
counts1 <- count1[ings]
idx <- match(ings,dtpa2$gene)
dtpa3 <- dtpa2[idx,]

identical(dtpa3$gene,names(counts1))
dtpa3$ge <- counts1

##计算相关性
cor = corr.test(dtpa3$pa,dtpa3$ge,  
                use = "pairwise", 
                method="pearson", 
                adjust = "fdr"   
)
cor$r 
cor$p 

cor1=ggplot(data=dtpa3, aes(x=pa, y=ge)) + geom_point(color="blue")+
  stat_smooth(method="lm",se=FALSE,size=0.5)+ labs(x="scPA",y="scGE")
topptx(cor1, filename = "fig2_cor1.pa.ge.pptx", width = 6, height = 5,units = "in")

# Add correlation coefficient: stat_cor()

ggscatter(dtpa3, x = "pa", y = "ge",
          add = "reg.line", conf.int = FALSE, 
          color = "#00AFBB",
          add.params = list(color = "#00AFBB",
                        cor.coef = TRUE    ))+
  stat_cor(method = "pearson")+
  labs(x="scPA",y="scGE")
topptx(cor2, filename = "fig2_cor2.pa.ge.pptx", width = 6, height = 5,units = "in")


###APA ratio ----
#29784pa ratio;128134cellls
ratio_har=readRDS("ratio_har.rds")
colnames(ratio_har@meta.data)
ratio.gene=gsub("PAC","",gsub("-.*","",row.names(ratio_har@assays$RNA@counts)))
ratio.gene=as.data.frame(table(ratio.gene))
table(ratio.gene$Freq==1)#18573 single pa gene;5082 APA gene


### APA ratio multiple ----
#11211 APA ratio;128134cells
mura_har=readRDS("mura_har.rds")
colnames(mura_har@meta.data)

id=match(row.names(mura_har@meta.data),gene_meta$ids3)
mura_har@meta.data$source.x=gene_meta$source[id]
mura_har@meta.data$ids=gene_meta$ids[id]
mura_har@meta.data$le=gene_meta$le[id]
mura_har@meta.data$final.ID2=gene_meta$final.ID2[id]
mura_har@meta.data$Main.tp=gene_meta$Main.tp[id]
mura_har@meta.data$timelv1=gene_meta$timelv1[id]
#saveRDS(mura_har,file="mura_har.rds")


mura_har@meta.data$final.ID2 <- factor(mura_har@meta.data$final.ID2, levels= c("Atrichoblast", "Trichoblast",
                                                                           "Cortex", "Endodermis", 
                                                                           "Columella", "Lateral Root Cap",
                                                                           "Pericycle", "Phloem", "Xylem", "Procambium",
                                                                           "Stem Cell Niche","Quiescent Center"),
                                     labels = c("Atrichoblast", "Trichoblast",
                                                "Cortex", "Endodermis", 
                                                "Columella", "Lateral Root Cap",
                                                "Pericycle", "Phloem", "Xylem", "Procambium",
                                                "Stem Cell Niche","Quiescent Center"))



##******Figure 3 APA markers plot********####
ms  <-c("WER1","COBL9","AAP1","CASP1","ACO4","TGG4","SOT18","APL","LAC11","PG3","AtHSBP","RGF3")

vs <-list()
Idents(pa_har)="final.ID2"
pa_har@meta.data$final.ID2 <- factor(pa_har@meta.data$final.ID2, levels= c("Atrichoblast", "Trichoblast",
                                                                           "Cortex", "Endodermis", 
                                                                           "Columella", "Lateral Root Cap",
                                                                           "Pericycle", "Phloem", "Xylem", "Procambium",
                                                                           "Stem Cell Niche","Quiescent Center"),
                            labels = c("Atrichoblast", "Trichoblast",
                                       "Cortex", "Endodermis", 
                                       "Columella", "Lateral Root Cap",
                                       "Pericycle", "Phloem", "Xylem", "Procambium",
                                       "Stem Cell Niche","Quiescent Center"))

for (i in 1:length(ms)) {
  k1 <- unique(symbol.gene[symbol.gene$SYMBOL %in% ms[i],]$PAC.id)
  ids <- str_sub(k1,str_locate(k1,"\\-.*")[1],str_locate(k1,"\\-.*")[2])
  y.title <- paste0(ms[i],":PAC",ids)
  vs[[i]] <- VlnPlot(pa_har, features =  k1, log = TRUE,pt.size=0)+
    labs(title = NULL) +ylab(y.title)+NoLegend()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
}


p1 <- plot_grid(vs[[1]],vs[[2]],vs[[3]],vs[[4]],vs[[5]],vs[[6]],ncol=1)

p2 <- plot_grid(vs[[7]],vs[[8]],vs[[9]],vs[[10]],vs[[11]],vs[[12]],ncol=1)

topptx(p1, filename = "pa.maikers1.pptx", width = 7.5, height = 6.5,units = "in")
topptx(p2, filename = "pa.maikers2.2.pptx", width = 7.5, height = 6.5,units = "in")

ncell <- as.data.frame(table(pa_har$final.ID2))

nmarkers <- as.data.frame(table(pa.marker$cluster))

p1 <- ggplot(ncell, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity")+ ylab("nCell")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + NoLegend()


nmarkers$Var1<- factor(nmarkers$Var1,levels = c("Atrichoblast", "Trichoblast",
                                                "Cortex", "Endodermis", 
                                                "Columella", "Lateral Root Cap",
                                                "Pericycle", "Phloem", "Xylem", "Procambium",
                                                "Stem Cell Niche","Quiescent Center"))

p2 <- ggplot(nmarkers, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity")+ ylab("nMarkers")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + NoLegend()


p=plot_grid(p2,p1,ncol=1)
topptx(p, filename = "ncells+markers2.pptx", width = 7.5, height = 6.5,units = "in")



##feature plot -----
pa_meta=readRDS("obj_meta.rds")
info.id=intersect(pa_meta$ids,gene_meta$ids)
gene_har=subset(gene_har, cells = info.id)

gene_har@meta.data$final.ID2 <- factor(gene_har@meta.data$final.ID2, levels= c("Atrichoblast", "Trichoblast",
                                                                               "Cortex", "Endodermis", 
                                                                               "Columella", "Lateral Root Cap",
                                                                               "Pericycle", "Phloem", "Xylem", "Procambium",
                                                                               "Stem Cell Niche","Quiescent Center"),
                                       labels = c("Atrichoblast", "Trichoblast",
                                                  "Cortex", "Endodermis", 
                                                  "Columella", "Lateral Root Cap",
                                                  "Pericycle", "Phloem", "Xylem", "Procambium",
                                                  "Stem Cell Niche","Quiescent Center"))


#NGPA$gene[286]:PACAT5G64816-1
Idents(gene_har)="final.ID2"
Idents(pa_har)="final.ID2"
Idents(ratio_har)="final.ID2"
Idents(mura_har)="final.ID2"

s1 <- FeaturePlot(pa_har, features ="PACAT5G64816-1" ,reduction="harmony_umap",raster=FALSE,slot = "data")
s1.1<-FeaturePlot(mura_har, features ="PACAT5G64816-1" ,reduction="combind_umap",raster=FALSE,slot = "data")
s2 <- FeaturePlot(gene_har, features = "AT5G64816",reduction="combat_umap",raster=FALSE,slot = "data")+
  ggtitle("Non-differential gene: AT5G64816")+NoLegend()
 
s3<- FeaturePlot(gene_har, features = "AT5G64816",reduction="combat_umap",raster=FALSE,slot = "data",max.cutoff = 3)+
  ggtitle("differential APA: AT5G64816-1")+NoLegend()

s3$data[,4] <- s1.1$data[,4]
#s3$data[,4] <- s1$data[,4]
png(filename = "fig3_ge+ratio.feature_unap.png", width = 5.5, height =9 ,units = "in",res=600)
plot_grid(s2,s3,nrow = 2,ncol=1)
  
dev.off()



###vlnplot

v1 <- VlnPlot(pa_har, features =  "PACAT5G64816-1", slot = "data", log = TRUE,pt.size=0)+labs(x="",y="log(pA Count)")+NoLegend()
v1.1 <- VlnPlot(mura_har, features =  "PACAT5G64816-1", slot = "data", log = TRUE,pt.size=0)+labs(x="",y="log(APA Ratio)")+NoLegend()

v2 <- VlnPlot(gene_har, features =  "AT5G64816", slot = "data", log = TRUE,pt.size=0)+labs(x="",y="log(Gene Count)")+NoLegend()

#which.min(apply(gene_har@assays$RNA@data[NGPAG,],1,sd))

p=plot_grid(v2,v1.1,ncol=1,nrow = 2)

#p=plot_grid(v2,v1,ncol=1,nrow=2)

topptx(p, filename = "fig3_vlnplot_ge+pa.count2.pptx", width = 4.5, height = 9,units = "in")


####Hair vs Non-Hair ----------------
library(clusterProfiler)
library("org.At.tair.db")

##绘图
Atrichoblast_lengthening <- read.csv("Atrichoblast_lengthening.csv")

Atrichoblast_shorting <- read.csv("Atrichoblast_shorting.csv")



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



##PPUI of main tissue ----
#ppui=readRDS("ppui.mtx.rds")
pa_har@meta.data$nPPUI <- colSums(ppui>0)
pa_har@meta.data$PPUIsum <- colSums(ppui)
pa_har@meta.data$PPUIavg <- pa_har@meta.data$PPUIsum / pa_har@meta.data$nPPUI
pa_har@meta.data$nPPUImax <- colSums(ppui == 1)


DimPlot(pa_har, reduction = "harmony_umap", group.by = "nPPUI",raster=FALSE)+NoLegend()
DimPlot(pa_har, reduction = "harmony_umap", group.by = "nPPUI",raster=FALSE)+NoLegend()

p1 <- FeaturePlot(object = pa_har, reduction = "harmony_umap", features = 'nPPUI',raster=FALSE)
p1$data$UMAP_1 <- pa_har@reductions$harmony_umap@cell.embeddings[,1]
p1$data$UMAP_2 <- pa_har@reductions$harmony_umap@cell.embeddings[,2]
p1  

p2 <- FeaturePlot(object = pa_har, reduction = "harmony_umap", features = 'PPUIsum')
p2$data$UMAP_1 <- pa_har@reductions$harmony_umap@cell.embeddings[,1]
p2$data$UMAP_2 <- pa_har@reductions$harmony_umap@cell.embeddings[,2]
p2  

p3 <- FeaturePlot(object = pa_har, reduction = "harmony_umap", features = 'nPPUImax',raster=FALSE)
p3$data$UMAP_1 <- pa_har@reductions$harmony_umap@cell.embeddings[,1]
p3$data$UMAP_2 <- pa_har@reductions$harmony_umap@cell.embeddings[,2]
p3

p4<- FeaturePlot(object = pa_har, reduction = "harmony_umap", features = 'PPUIavg',raster=FALSE)
p4$data$UMAP_1 <- pa_har@reductions$harmony_umap@cell.embeddings[,1]
p4$data$UMAP_2 <- pa_har@reductions$harmony_umap@cell.embeddings[,2]
p4



dt <- data.frame(PPUIsum = pa_har@meta.data$PPUIsum,barcode = rownames(pa_har@meta.data))

layer <- pa_har@reductions$harmony_umap@cell.embeddings
layer <- as.data.frame(layer)
library(RColorBrewer)

png(filename = "fig3_ppui.sum_umap.png", width = 7, height = 6.5,units = "in", res=600)
ggplot(data = dt,aes(x = as.numeric(layer$UMAP_1),y = as.numeric(layer$UMAP_2)))+
  geom_point(aes(color=PPUIsum),size=0.5,shape=16)+
  scale_color_distiller(palette = "Spectral")+
  labs(x = "UMAP_1", y = "UMAP_2",color="PPUIsum")+
  theme_bw()+theme(panel.grid =element_blank(),
                   plot.title = element_text(hjust = 0.5),
                   legend.text = element_text(size=12),
                   legend.title  = element_blank(),
                   legend.position = "right")
dev.off()

pa_har@meta.data$time.num <- as.numeric(as.factor(pa_har@meta.data$timelv1))

library(rstatix) 
library(tidyr)

cor.test(pa_har@meta.data$nPPUI,pa_har@meta.data$time.num,method = "pearson")

library(ggpubr)
pa_har@meta.data$timelv1

select_df <- pa_har@meta.data[,c("Main.tp","final.ID2","nPPUI","PPUIsum","nPPUImax","time.num","timelv1")]

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


p= ggplot(select_df, aes(x=Main.tp, y=nPPUI, fill=Main.tp)) + 
  geom_boxplot(alpha=0.3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Dark2")+theme_bw()+
  labs(x = "Tissue", y = "PPUI Sum",color="Main.tp")+
  theme_bw()+theme(panel.grid =element_blank(),
                   plot.title = element_text(hjust = 0.5),
                   legend.text = element_text(size=12),
                   legend.title  = element_blank(),
                   legend.position = "right")
topptx(p,filename = "fig3_ppui.sum_boxplot.pptx", width = 7, height = 6.5,units = "in")
#saveRDS(select_df,file = "BXY-20240322/select_df.rds")

##******Fingure 4 Combining gene expression and APA profiles to refine cell identity*******#####
BiocManager::install("AUCell")

library(GSEABase)
library(AUCell)
## QC
# find markers

# label 1 
###markers based on pa count ----
#pamarkers=readRDS("BXY-20240322/filter.pa.rds")
pacounts <- pa_har@assays$RNA@data
dim(pacounts)#29784 128134

pa_qc <- pacounts[,pa_har$final.ID2 == "Quiescent Center"| pa_har$final.ID2 == "Stem Cell Niche"]
dim(pa_qc)#29784  1483

markers_pa <- pamarkers[pamarkers$cluster == "Stem Cell Niche" | pamarkers$cluster == "Quiescent Center",]
table(markers_pa$cluster)
#QC 540;STC 363
pa.list <- split(markers_pa$gene, factor(markers_pa$cluster))

library(GSEABase)
all.sets <- lapply(names(pa.list), function(x) {
  GeneSet(pa.list[[x]], setName=x)        
})
#Use GeneSetCollection to construct a collection of gene sets from GeneSet arguments, or a list of GeneSets.
all.sets <- GeneSetCollection(all.sets)


library(AUCell)
Wt <- as.matrix(pa_qc)
rankings <- AUCell_buildRankings(Wt,plotStats=FALSE, verbose=FALSE)
#Calculates the 'AUC' for each gene-set in each cell
cell.aucs <- AUCell_calcAUC(all.sets, rankings,aucMaxRank=nrow(rankings)*0.05,verbose=FALSE)
results <- t(cell.aucs@assays@data$AUC)
new.labels1 <- colnames(results)[max.col(results)]

##label2
####markers based on APA ratio ----
mucounts <- mura_har@assays$RNA@data
mu_qc <- mucounts[,mura_har$final.ID2 == "Quiescent Center"| mura_har$final.ID2 == "Stem Cell Niche"]
dim(mu_qc)#11211  1483

#mumarkers=readRDS("BXY-20240322/mumarkers.rds")#3648
markers_mu <- mumarkers[mumarkers$cluster == "Stem Cell Niche" | mumarkers$cluster == "Quiescent Center",]
table(markers_mu$cluster)
#SCN 214;QC 333

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

label.qc <- data.frame(pa=new.labels1,mu=new.labels2)
label.qc$ge <- mura_har$final.ID2[mura_har$final.ID2 == "Quiescent Center"| mura_har$final.ID2 == "Stem Cell Niche"]

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv <- uniqv[!is.na(uniqv)]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


label.qc$final <- apply(label.qc,1,getmode)

table(label.qc$final)
#Quiescent Center  Stem Cell Niche 
#564              919 
saveRDS(label.qc,file="BXY-20240322/label.qc.rds")

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

#构建不同细胞类型中的pa markes集
library(GSEABase)
all.sets <- lapply(names(pa.list), function(x) {
  GeneSet(pa.list[[x]], setName=x)        
})
all.sets <- GeneSetCollection(all.sets)

library(AUCell)

Wt <- as.matrix(pa_qc)
rankings <- AUCell_buildRankings(Wt,plotStats=FALSE, verbose=FALSE)
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

label.qc <- data.frame(pa=new.labels1,mu=new.labels2)
label.qc$ge <- mura_har$final.ID2[is.na(mura_har$le)]
#label.qc<- label.qc[,c(3,1,2)]
label.qc$newfinal <- apply(label.qc,1,getmode)

table(label.qc$newfinal)
table(label.qc$ge)

order <- c("Atrichoblast", "Trichoblast",
           "Cortex", "Endodermis", 
           "Columella", "Lateral Root Cap",
           "Pericycle", "Phloem", "Xylem", "Procambium",
           "Stem Cell Niche","Quiescent Center","Marked")

palette <- c("#7fb80e", "#1d953f",
             "#2a5caa", "#a565ef",
             "#843900", "#69541b",
             "#f7b13f", "#f9e264", "#009db2", "#375830",
             "#f58f98","#d71345","#e7dad2")

pa_har$le2g <- "Marked"
pa_har$le2g[is.na(pa_har$le)] <- as.character(pa_har$final.ID2[is.na(pa_har$le)])

png(filename = "fig4_filter.gene.markerd_umap.png", width = 7, height = 6.5,units = "in", res=600)
DimPlotOrder(pa_har,reduction = "harmony_umap",group.by = "le2g",order=order,palette=palette)
dev.off()

pa_har$le2p <- "Marked"
pa_har$le2p[is.na(pa_har$le)] <- as.character(label.qc$newfinal)

png(filename = "fig4_filter.pa.markerd_umap.png", width = 7, height = 6.5,units = "in", res=600)
DimPlotOrder(pa_har,reduction = "harmony_umap",group.by = "le2p",order=order,palette=palette)
dev.off()

##绘制冲积图 ----
library(ggalluvial)
library(ggplot2)

table(label.qc$newfinal)
table(label.qc$ge)

label.qc$ge <- factor(label.qc$ge, levels= c("Atrichoblast", "Trichoblast",
                                                                               "Cortex", "Endodermis", 
                                                                               "Columella", "Lateral Root Cap",
                                                                               "Pericycle", "Phloem", "Xylem", "Procambium",
                                                                               "Stem Cell Niche","Quiescent Center"),
                                       labels = c("Atrichoblast", "Trichoblast",
                                                  "Cortex", "Endodermis", 
                                                  "Columella", "Lateral Root Cap",
                                                  "Pericycle", "Phloem", "Xylem", "Procambium",
                                                  "Stem Cell Niche","Quiescent Center"))

label.qc$newfinal <- factor(label.qc$newfinal, levels= c("Atrichoblast", "Trichoblast",
                                             "Cortex", "Endodermis", 
                                             "Columella", "Lateral Root Cap",
                                             "Pericycle", "Phloem", "Xylem", "Procambium",
                                             "Stem Cell Niche","Quiescent Center"),
                      labels = c("Atrichoblast", "Trichoblast",
                                 "Cortex", "Endodermis", 
                                 "Columella", "Lateral Root Cap",
                                 "Pericycle", "Phloem", "Xylem", "Procambium",
                                 "Stem Cell Niche","Quiescent Center"))

png(filename = "fig4_filter.Alluvial.map2.1.png", width = 10.5, height = 8,units = "in",res=600)
ggplot(data = label.qc,aes(axis1 = ge, axis2 = newfinal)) +
  geom_alluvium(aes(fill = newfinal)) +
  geom_stratum() + geom_text(stat = "stratum", aes(label = after_stat(stratum)))+
  theme_minimal() +
  ggtitle("Alluvial maps based on original unmarked cells and annotated information corrected by pA and APA ratio markers")+
  scale_x_discrete(limits = c("raw marked cell", "new marked cell"), expand = c(.2, .05))
dev.off()
 




