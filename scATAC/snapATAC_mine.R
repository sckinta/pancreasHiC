library(SnapATAC)
library(GenomicRanges)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
dir="/mnt/isilon/sfgi/suc1/analyses/grant/scATAC/pancreaticCells"

############################ snap object ###############################
### createSnap
snaptool_data=file.path(dir,"snaptools")

snap_files=dir(snaptool_data, pattern="\\.snap$", full.names=T) # "HPAP040", "HPAP045", "HPAP054", "HPAP066"
snap_files=snap_files[grepl("HPAP040|HPAP045|HPAP054", snap_files)] # "HPAP040", "HPAP045", "HPAP054"

x.sp = createSnap(
    file=snap_files,
    sample=gsub(".snap","",basename(snap_files)),
    num.cores=1
) #  17530 barcodes (17940 barcodes when have HPAP066)

x.sp@metaData %>% 
as_tibble() %>% 
bind_cols(sample=x.sp@sample) %>% 
group_by(sample) %>% 
count()
# sample              n
# <chr>           <int>
# 1 FGC2061_HPAP040  4029
# 2 FGC2061_HPAP045  4402
# 3 FGC2147_HPAP045  1691
# 4 FGC2147_HPAP054  3624
# 5 FGC2154_HPAP054  3784
# 6 FGC2185_HPAP066   410

## decide to remove FGC2185_HPAP066 due to 
## 1) too low number of cells
## 2) not included in scRNA
## 3) addBmatToSnap report error "bins does not match between snap files, please regenerate the cell-by-bin matrix by snaptools"

### filter sp by single cell fragment summary for each barcode derive from cell ranger
load(file.path(dir,"downloads/selected_donors/HPAP_scATAC_barcode_cellranger_4donors.RData")) # sc_summary_comb
sc_summary_comb = sc_summary_comb %>% as_tibble() %>% 
filter(sample %in% gsub(".snap","",basename(snap_files)))

sc_summary_comb = sc_summary_comb %>% 
mutate(logUMI=log10(sc_summary_comb$passed_filters + 1)) %>% 
mutate(promoter_ratio=(promoter_region_fragments+1) / (passed_filters + 1)) %>% 
mutate(MT_ratio=mitochondrial/(mitochondrial + passed_filters))

# filter passed_filters and mito ratio
sc_summary_comb_filtered = sc_summary_comb %>% 
filter(passed_filters > 1000) %>% 
filter(MT_ratio <= 0.1)

# filter logUMI and promoter_ratio
p <- sc_summary_comb_filtered %>% 
separate(sample,c("run","donor"), sep="_", remove=F) %>% 
ggplot(aes(x=logUMI, y=promoter_ratio)) +
facet_wrap(~donor) +
geom_point() +
geom_hline(yintercept=0.1, linetype=3, color="red") +
geom_hline(yintercept=0.6, linetype=3, color="red") +
geom_vline(xintercept=3, linetype=3, color="red") +
geom_vline(xintercept=6, linetype=3, color="red") +
theme_bw()
ggsave(file.path(dir,"plots","QC_logUMI_vs_promoterRatio.png"), p, width=9, height=3)

sc_summary_comb_filtered = sc_summary_comb_filtered %>% 
filter(logUMI >= 3 & logUMI <= 6) %>%
filter(promoter_ratio >= 0.1 & promoter_ratio <= 0.6)


# change barcode for filter temperarily (only used when there are multiple samples)
sc_summary_comb_filtered = sc_summary_comb_filtered %>% 
mutate(barcodeNew=glue::glue("{sample}_{barcode}"))
x.sp@barcode = glue::glue("{x.sp@sample}_{x.sp@barcode}")
x.sp@metaData$barcode = x.sp@barcode

# filter with sc_summary_comb_filtered
x.sp = x.sp[which(x.sp@barcode %in% sc_summary_comb_filtered$barcodeNew),]

# add all sc_summary_comb_filtered information to x.sp@metaData
x.sp@metaData = x.sp@metaData %>% as_tibble() %>% 
dplyr::rename(barcodeNew = barcode) %>% 
left_join(
        sc_summary_comb_filtered %>% 
        separate(sample, c("run","donor"), sep="_", remove=F)
        
) %>% as.data.frame() # TN is the passed_filters value from cellRanger; logUMI=log10(TN)

x.sp@metaData %>% 
as_tibble() %>% 
group_by(donor) %>% 
summarise(
        cell_n=n_distinct(barcodeNew),
        runs = paste(unique(run), collapse=",")
) %>% ungroup()
# donor   cell_n runs           
# <chr>    <int> <chr>          
# 1 HPAP040   2653 FGC2061        
# 2 HPAP045   3308 FGC2061,FGC2147
# 3 HPAP054   6512 FGC2147,FGC2154

save(x.sp, file=file.path(dir,"Rdata","sp_filtered.Rdata"))
# 12473 (12722 with HPAP066) cells

########################### snap object with bmat #############################
# need to make barcode consistent with .snap file ones (only used when there are multiple samples)
x.sp@barcode = sapply(x.sp@barcode, function(s){strsplit(s,"_")[[1]][3]})
x.sp@metaData$barcode = x.sp@barcode

# add bmat
lapply(snap_files, function(file){showBinSizes(file)}) # 1000  5000 10000
x.sp = addBmatToSnap(x.sp, bin.size=1000, num.cores=1) # 3137206 bins

# bmat binarization
x.sp = makeBinary(x.sp, mat="bmat") # 3137206 bins, bin coords saved in x.sp@feature

save(x.sp, file=file.path(dir,"Rdata","sp_filtered.Rdata"))

# bin filtering - black_list
black_list = vroom::vroom("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/encode_blacklists/hg19_encode_blacklists.bed.gz", col_names=c("chr","start","end"))
black_list.gr = GRanges(
    black_list$chr, IRanges(black_list$start, black_list$end)
)
idy = queryHits(
    findOverlaps(x.sp@feature, black_list.gr)
)
if(length(idy) > 0){
    x.sp = x.sp[,-idy, mat="bmat"];
}
# 3122175 bins

# bin filtering - non-standard chr
seqlevels(x.sp@feature)[!grepl("_|chrM", seqlevels(x.sp@feature))]
chr.exclude = seqlevels(x.sp@feature)[grepl("_|chrM", seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0){
    x.sp = x.sp[,-idy, mat="bmat"]
}
# 3081370 bins

# bin filtering - remove the top 5% bins that overlap with invariant features such as the house keeping gene promoters.
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
png(file.path(dir,"plots","QC_logBinCoverage.png"))
hist(
    bin.cov[bin.cov > 0], 
    xlab="log10(bin cov)", 
    main="log10(Bin Cov)", 
    col="lightblue", 
    xlim=c(0, 5)
)
dev.off()
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]
# 12473 cells x 2639549 bins
save(x.sp, file=file.path(dir,"Rdata","sp_filtered.Rdata"))

# cell filtering - remove any cells of bin coverage less than 1,000 (optional)
idx = which(Matrix::rowSums(x.sp@bmat) > 1000)
x.sp = x.sp[idx,]
# 10391 cells x 2639549 bins
save(x.sp, file=file.path(dir,"Rdata","sp_filtered2.Rdata"))

############################## clustering #########################
load(file.path(dir,"Rdata","sp_filtered2.Rdata"))
# # Dimentionality reduction - sampling 10K as landmarks based on the total read  distribution
# row.covs.dens = density(
#     x = x.sp@metaData[,"logUMI"], 
#     bw = 'nrd', adjust = 1
# )
# sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = x.sp@metaData[,"logUMI"])$y + .Machine$double.eps)
# set.seed(1)
# idx.landmark.ds = sort(sample(x = seq(nrow(x.sp)), size = 10000, prob = sampling_prob))
# x.landmark.sp = x.sp[idx.landmark.ds,]
# x.query.sp = x.sp[-idx.landmark.ds,]
# 
# # Dimentionality reduction - Run diffusion maps on the landmark cells
# x.landmark.sp_s1 = runDiffusionMaps(
#     obj= x.landmark.sp,
#     input.mat="bmat", 
#     num.eigs=50
# )
# # @ dmat  : num [1:10000, 1:50] 0.007117 -0.016846 0.003809 0.005648 0.000871 ...
# #   .. .. ..@ sdev  : num [1:50] 0.2442 0.1073 0.0798 0.0564 0.0421 ...
# save(x.landmark.sp_s1, file=file.path(dir,"Rdata","tmp_s1_smat_cluster.Rdata"))
# 
# rm(x.landmark.sp_s1)
# x.landmark.sp_s2 = runDiffusionMaps(
#     obj= x.landmark.sp,
#     input.mat="bmat", 
#     num.eigs=50
# )
# str(x.landmark.sp_s2)
# # @ dmat  : num [1:10000, 1:50] 0.007112 -0.016735 0.003801 0.005634 0.000886 ...
# #   .. .. ..@ sdev  : num [1:50] 0.2449 0.1075 0.0797 0.0565 0.0421 ...
# save(x.landmark.sp_s2, file=file.path(dir,"Rdata","tmp_s2_umap_cluster.Rdata"))
# x.landmark.sp@metaData$landmark = 1
# 
# # Dimentionality reduction - Project query cells to landmarks.
# x.query.sp = runDiffusionMapsExtension(
#     obj1=x.landmark.sp, 
#     obj2=x.query.sp,
#     input.mat="bmat"
# )
# x.query.sp@metaData$landmark = 0
# 
# # Dimentionality reduction - Combine landmark and query cells.
# x.sp = snapRbind(x.landmark.sp, x.query.sp)
# x.sp = x.sp[order(x.sp@metaData[,"sample"])] #IMPORTANT: To merge snap objects, all the matrix (bmat, gmat, pmat) and metaData must be of the same number of columns between snap objects.

# Dimentionality reduction (update @smat) without landmark way (since we only have 10391 cells, wont increase much computational time by landmark 10K cells)
set.seed(1)
x.sp = runDiffusionMaps(
    obj= x.sp,
    input.mat="bmat", 
    num.eigs=50
)
str(x.sp)
# @ smat    :Formal class 'dim.reduct' [package "SnapATAC"] with 5 slots
#   .. .. ..@ imat  : chr(0) 
#   .. .. ..@ dmat  : num [1:10391, 1:50] -0.007022 0.016292 -0.003721 -0.005542 -0.000886 ...
#   .. .. ..@ sdev  : num [1:50] 0.2453 0.1086 0.0791 0.0566 0.0426 ...
save(x.sp, file=file.path(dir,"Rdata","sp_filtered2_diffussionMapS1.Rdata"))

# Determine eigs number for cluster: select the number of eigen vectors that the scatter plot starts looking like a blob
png(file.path(dir,"plots","QC_DimReduction_PWeigs.png"))
plotDimReductPW(
    obj=x.sp, 
    eigs.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=NULL, 
    pdf.height=7, 
    pdf.width=7
)
dev.off()
ndim = 14 # select number of eigen as 14 based on QC_DimReduction_PWeigs.png

# remove run_donor (sample) batch effect
load(file.path(dir,"Rdata","sp_filtered2_diffussionMapS1.Rdata"))
x.sp = runHarmony(
    obj=x.sp,
    eigs.dim=1:ndim,
    meta_data=x.sp@sample
)
# @ smat    :Formal class 'dim.reduct' [package "SnapATAC"] with 5 slots
#   .. .. ..@ imat  : chr(0) 
#   .. .. ..@ dmat  : num [1:10391, 1:14] -0.00646 0.007516 -0.003113 -0.00489 0.000626 ...
#   .. .. ..@ sdev  : num [1:14] 0.2453 0.1086 0.0791 0.0566 0.0426 ...

# Run cluster - runKNN() to build graph, update @graph
x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:ndim, # A vector of the dimensions to use in construction of the KNN graph.
    k=15 # K for the k-nearest neighbor algorithm. default is 15
)
save(x.sp, file=file.path(dir,"Rdata","sp_filtered2_rmBatchKnn.Rdata"))

# Run cluster - cluster graph using leiden algorithm (opiton-1, recommanded), update @cluster
# library(leiden)
x.sp.leiden=runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="leiden",
    seed.use=10,
    resolution=0.7 # A numeric value that indicates the resolution for louvain clustering (default is 1)
)

x.sp@metaData$leiden_cluster <-x.sp.leiden@cluster # 12 clusters

# Run cluster - cluster graph using igraph algorithm (opiton-2), update @cluster
x.sp.igraph = runCluster(
        obj=x.sp, 
        tmp.folder=tempdir(), 
        louvain.lib="R-igraph", 
        seed.use=10
)
x.sp@metaData$igraph_cluster <- x.sp.igraph@cluster # 17 clusters

# visualize clusters - generate umap first (`runViz`), update @umap (based on @graph not the cluster)
library(umap)
x.sp = runViz(
    obj=x.sp, 
    tmp.folder=tempdir(),
    dims=2,
    eigs.dims=1:ndim, 
    method="umap",
    seed.use=10
)

save(x.sp, file=file.path(dir,"Rdata","sp_filtered2_cluster.Rdata"))

# visualize clusters - plotting with plotViz(point.color=) for category value in color and plotFeatureSingle(feature.value=) for numeric value in color.
load(file.path(dir,"Rdata","sp_filtered2_cluster.Rdata"))

## cluster leiden vs. igraph
pdf(file.path(dir,"plots","QC_scATAC_Cluster.pdf"), width=20, height=10)
par(mfrow = c(1,2))
plotViz(
    obj= x.sp,
    method="umap", 
    main="leiden cluster",
    point.color=x.sp@metaData$leiden_cluster, 
    point.size=0.2, 
    point.shape=19, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    # down.sample=10000,
    legend.add=FALSE
)
plotViz(
    obj= x.sp,
    method="umap", 
    main="igraph cluster",
    point.color=x.sp@metaData$igraph_cluster, 
    point.size=0.2, 
    point.shape=19, 
    text.add=TRUE,
    text.size=1,
    text.color="black",
    # down.sample=10000,
    legend.add=FALSE
)
dev.off()

## other metrics (logUMI, sample, donor)
pdf(file.path(dir,"plots","QC_scATAC_clusterMetrics.pdf"), width=14, height=14)
par(mfrow=c(2,2))
plotFeatureSingle(
    obj=x.sp,
    feature.value=x.sp@metaData[,"logUMI"],
    method="umap", 
    main="Read Depth",
    point.size=0.2, 
    point.shape=19, 
    # down.sample=10000,
    quantiles=c(0.01, 0.99)
)
plotViz(
    obj= x.sp,
    method="umap", 
    main="Sample",
    point.color=x.sp@metaData[,"sample"], 
    point.size=0.2, 
    point.shape=19,
    text.add=F,
    # text.size=1,
    # text.color="black"
    # down.sample=10000,
   legend.add=T
)
plotViz(
    obj= x.sp,
    method="umap", 
    main="Run",
    point.color=x.sp@metaData[,"run"], 
    point.size=0.2, 
    point.shape=19,
    text.add=F,
    # text.size=1,
    # text.color="black"
    # down.sample=10000,
   legend.add=T
)
plotViz(
    obj= x.sp,
    method="umap", 
    main="Donor",
    point.color=x.sp@metaData[,"donor"], 
    point.size=0.2, 
    point.shape=19,
    text.add=F,
    # text.size=1,
    # text.color="black",
    # down.sample=10000,
    legend.add=T
)
dev.off()

################### assign celltype  (Long's method) ################
load(file.path(dir,"Rdata","sp_filtered2_cluster.Rdata"))
#### based on accessibility across marker gene body and some promoter

### get gene grange coordinates (using gencode.v19)
genes.gr = read_delim(
        "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene.bed", 
        col_names=c("chr","start","end","gene","dot","strand"),
        delim="\t",
        comment="#"
) %>% 
separate(gene,c("name","gene_id"), sep="\\+") %>% 
dplyr::select(-dot) %>% 
makeGRangesFromDataFrame(
        starts.in.df.are.0based=T,
        keep.extra.columns=T
)

### selected gene gr (genes.sel.gr 210 genes)
cellType_geneMarkers = read_delim("/mnt/isilon/sfgi/suc1/analyses/grant/scRNA/pancreaticCells/downloads/pancreasGenes.txt", delim="\t", col_names="gene") %>% distinct()

genes.sel.gr = genes.gr[which(genes.gr$name %in% cellType_geneMarkers$gene)]
genome(genes.sel.gr) = "hg19"

### selected region gr 
region.sel.gr = makeGRangesFromDataFrame(
        data.frame(
                chr=c("chr11","chr2","chr2", "chr3", "chr17"),
                start=c(2177480, 163009884,163006548,187387437,42021533),
                end=c(2178750,163011028,163011131, 187389301,42021913),
                strand=c("+", "+", "+", "+","+"),
                name=c("INS_IGF2", "GCG_UP", "GCG_Prom", "SST_add", "PPY_add"),
                stringsAsFactors=F
        ),
        keep.extra.columns=T
)
genome(region.sel.gr) = "hg19"

### make @gmat (matrix row is cell, y is feature, value is accessibility)
sp.ga = x.sp
sp.ga = createGmatFromMat(
    obj=sp.ga,
    input.mat="bmat",
    genes=c(genes.sel.gr,region.sel.gr),
    do.par=TRUE,
    num.cores=10
)

sp.ga = scaleCountMatrix(
    obj=sp.ga,
    cov=sp.ga@metaData$passed_filters + 1,
    mat="gmat",
    method = "RPM"
) # normalize the cell-by-gene matrix

sp.ga = runMagic(
    obj=sp.ga,
    input.mat="gmat",
    step.size=3
) # smooth the cell-by-gene matrix

save(sp.ga, file=file.path(dir,"Rdata","spga_filtered2_clusterCelltype.Rdata"))

### marker gene accessibility plot
marker2cell = c(
        "CPA1"="acinar", "PRSS1"="acinar", "REG1A"="acinar", 
        "GCG"="alpha", "IRX1"="alpha", "IRX2"="alpha", "ARX"="alpha",
        "INS"="beta","NKX6-1"="beta", "MAFA"="beta", "IAPP"="beta","IGF2"="beta",
        "SST"="delta",
        "KRT19"="ductal",
        "VWF"="endo", "ESAM"="endo", "CD93"="endo",
        "GHRL"="eps",
        "COL1A1"="mes",
        "PPY"="pp", "SLC6A4"="pp",
        "SDS"="immune", "TPSAB1"="immune"
)
marker2cell = marker2cell[names(marker2cell) %in% colnames(sp.ga@gmat)]

### marker gene accessibility plot (igraph heatmap)
igraph_cluster = as.matrix(sp.ga@gmat)[,names(marker2cell)] %>% 
as.data.frame() %>% 
as_tibble() %>% 
mutate(cluster=sp.ga@metaData$igraph_cluster) %>% 
gather(key="gene", value="acc", -cluster) %>% 
group_by(gene, cluster) %>% 
# summarise(med_acc = median(acc)) %>% 
summarise(med_acc = mean(acc)) %>% 
ungroup() %>% 
left_join(
        tibble(
                gene=names(marker2cell),
                cellType=marker2cell
        )
) %>% 
spread(key=cluster, value=med_acc) %>% 
arrange(cellType)

igraph_cluster_m = igraph_cluster %>% 
select(-cellType) %>% 
as.data.frame() %>% 
column_to_rownames("gene") %>% 
as.matrix()

igraph_cluster_scale_m = t(apply(igraph_cluster_m,1,scale))
colnames(igraph_cluster_scale_m) = colnames(igraph_cluster_m)
idx = sapply(
        1:nrow(igraph_cluster_scale_m), 
        function(x){!all(is.na(igraph_cluster_scale_m[x,]))}
) # remove genes without variability across clusters
igraph_cluster_scale_m = igraph_cluster_scale_m[idx,]
igraph_cluster_m = igraph_cluster_m[idx,]
igraph_cluster = igraph_cluster[idx,]

library(ComplexHeatmap)
library(circlize)
cellType_annotation <- HeatmapAnnotation(
    cellType = igraph_cluster %>% select(gene, cellType) %>% deframe(),
    # col = list(donor=deframe(donor_colors[,c(1,2)])),
    annotation_legend_param = list(
            title = "cellType",
            at =unique(igraph_cluster$cellType),
            labels = unique(igraph_cluster$cellType)
    ),
    which = "row"
)
h <- Heatmap(
        igraph_cluster_scale_m, 
        name = "marker gene scaled median access", 
        col = colorRamp2(c(-2,0,2), c("blue","white" ,"red"),space="RGB"),
        cluster_columns = T, 
        clustering_distance_columns=function(x) as.dist(1-cor(t(x))),
        clustering_method_columns="average",
        show_column_dend= T,
        show_column_names = T, 
        cluster_rows=F,
        left_annotation=cellType_annotation,
)
X11.options(colortype="pseudo.cube")
pdf(file.path(dir,"plots","markerGene_igraphCluster_meanAcc_heatmap.pdf"))
draw(h)
dev.off()
X11.options(reset = T)

cellType_assignment = list(
        acinar=7,
        alpha=c(6,13,12,17,10,8,16), # ?4
        beta=c(2,9,1),
        delta=11,
        ductal=14, 
        endo=5, # c(5,15)
        eps=4, # c(6,13,12,17,10,8,16)
        mes=15, # c(5,15)
        pp=3
)

lapply(
        unique(marker2cell),
        function(celltype){
                genes = names(marker2cell[marker2cell==celltype])
                p <- as.matrix(sp.ga@gmat)[,genes] %>% 
                as.data.frame() %>% 
                as_tibble() %>% 
                mutate(igraph_cluster=sp.ga@metaData$igraph_cluster) %>% 
                mutate(leiden_cluster=sp.ga@metaData$leiden_cluster) %>% 
                gather(key="gene", value="acc", -igraph_cluster, -leiden_cluster) %>% 
                ggplot(aes(x=igraph_cluster, y=log2(acc))) +
                facet_wrap(~gene) +
                geom_violin() +
                geom_boxplot(width=0.1, outlier.shape=NA) +
                theme_bw()
                filename=paste0("violin.",celltype,".geneAcc.pdf")
                ggsave(file.path(dir, "plots", filename), p, height=5, width=5*length(genes))
        }
)

lapply(
        unique(marker2cell),
        function(celltype){
                genes = names(marker2cell[marker2cell==celltype])
                filename=paste0("umap.",celltype,".geneAcc.pdf")
                pdf(file.path(dir,"plots",filename), height=5, width=5*length(genes))
                par(mfrow=c(1,length(genes)))
                for (i in 1:length(genes)){
                        plotFeatureSingle(
                            obj=sp.ga,
                            feature.value=sp.ga@gmat[, genes[i]],
                            method="umap", 
                            main=genes[i],
                            point.size=0.1, 
                            point.shape=19, 
                            # down.sample=10000,
                            quantiles=c(0, 1)
                        )
                }
                dev.off()
        }
)


###################### transfer scRNA-seq cell cluster to x.sp #######################
load(file.path(dir,"Rdata","sp_filtered2_cluster.Rdata")) # x.sp
load("/mnt/isilon/sfgi/suc1/analyses/grant/scRNA/pancreaticCells/downloads/scRNAseq_SeuratData.RData") # panc_all

#### determine scRNA-seq cellType variable genes
library(Seurat)
unique(panc_all@active.ident) # Levels: acinar alpha beta delta ductal endo eps mes pp

# wilcox method geneMarkers from RNA-seq
cellType_geneMarkers1 = FindAllMarkers(
        panc_all, 
        test.use = "wilcox", # default
        only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
)# using default wilcox method

cellType_geneMarkers1 = cellType_geneMarkers1 %>%
as_tibble() %>% 
filter(p_val_adj < 0.05) %>% 
arrange(desc(avg_log2FC)) %>% 
group_by(cluster) %>% 
mutate(rank = row_number()) %>% 
ungroup() %>% 
arrange(cluster) # rank feature within cluster

cellType_geneMarkers1 %>% 
group_by(cluster) %>%
top_n(n = 3, wt = avg_log2FC) %>% 
ungroup() %>% 
as.data.frame()

# using known pancreasGenes guide the feature selections (can skip if no list was provided)
cellType_geneMarkers1 %>% 
filter(gene %in% c("ARX","KRT19","VWF","ESAM","CD93","GHRL","COL1A1"))

cellType_geneMarkers1 = cellType_geneMarkers1 %>% 
semi_join(
        read_delim("/mnt/isilon/sfgi/suc1/analyses/grant/scRNA/pancreaticCells/downloads/pancreasGenes.txt", delim="\t", col_names="gene")
) 

cellType_geneMarkers1 %>% 
distinct(cluster, gene) %>% 
group_by(cluster) %>% 
summarise(gene_n=n_distinct(gene))
# cluster gene_n
# <fct>    <int>
# 1 acinar      31
# 2 alpha       31
# 3 beta        15
# 4 delta       18
# 5 ductal      39
# 6 endo        21
# 7 eps         11
# 8 mes         42
# 9 pp          13

### get gene grange coordinates (using gencode.v19)
genes.gr = read_delim(
        "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene.bed", 
        col_names=c("chr","start","end","gene","dot","strand"),
        delim="\t",
        comment="#"
) %>% 
separate(gene,c("name","gene_id"), sep="\\+") %>% 
dplyr::select(-dot) %>% 
makeGRangesFromDataFrame(
        starts.in.df.are.0based=T,
        keep.extra.columns=T
)

### selected gene gr (genes.sel.gr 146 genes)
genes.sel.gr = genes.gr[which(genes.gr$name %in% cellType_geneMarkers1$gene)]
genome(genes.sel.gr) = "hg19"

marker2cell[names(marker2cell) %in% genes.sel.gr$name] # eps (GHRL) does not have representative

### add gmat (accessibility) with different geneset
sp.ga2 = x.sp
sp.ga2 = createGmatFromMat(
    obj=sp.ga2, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=10
)

save(sp.ga2, file=file.path(dir,"Rdata","spga2_filtered2_clusterCelltype.Rdata"))

### integrate atac snap with scRNA
ndim=14 # pre-determined dim used for cluster (see above)
se.ga2 = snapToSeurat(
    obj=sp.ga2, 
    eigs.dims=1:ndim, 
    norm=TRUE,
    scale=TRUE
) # make snap atac to seurat format (gene x 10392 cell), two assay data (@assay$ATAC, @assay$ACTIVITY), ACTIVITY contains 146 genes

transfer.anchors = FindTransferAnchors(
    reference = panc_all, 
    query = se.ga2, 
    features = genes.sel.gr$name, 
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    reduction = "cca"
) # find shared anchors (gene_name) between scRNA and atac => AnchorSet object

celltype.predictions = TransferData(
    anchorset = transfer.anchors, 
    refdata = "cellType", # name of meta.data in RNA for cluster
    reference = panc_all,
    weight.reduction = se.ga2[["SnapATAC"]],
    dims = 1:ndim
) # Transfer RNA-seq categorical data (cellType) based on ATAC reduction object through anchors

# unique(celltype.predictions$predicted.id)
sp.ga2@metaData$predicted.id = celltype.predictions$predicted.id # transfered RNA-seq celltype
sp.ga2@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max)
sp.ga2@metaData %>% as_tibble()

# cluster assignment
table(
        sp.ga2@metaData %>% 
        as_tibble() %>% 
        dplyr::select(igraph_cluster, predicted.id)
)

table(
        sp.ga2@metaData %>% 
        as_tibble() %>% 
        filter(predict.max.score > 0.5) %>% 
        dplyr::select(igraph_cluster, predicted.id)
)

atacCluster2rnaPredictedCellType = table(
        sp.ga2@metaData %>% 
        as_tibble() %>% 
        filter(predict.max.score > 0.5) %>% 
        dplyr::select(igraph_cluster, predicted.id)
) %>% as.data.frame() %>% as_tibble()

atacCluster2rnaPredictedCellType = atacCluster2rnaPredictedCellType %>% 
left_join(
        atacCluster2rnaPredictedCellType %>% 
        group_by(igraph_cluster) %>% 
        summarise(total_cell = sum(Freq))
) %>% 
mutate(ratio = Freq/total_cell) %>% 
group_by(igraph_cluster) %>% 
dplyr::slice(which.max(ratio)) %>% 
ungroup() %>% 
arrange(predicted.id)
# igraph_cluster predicted.id  Freq total_cell ratio
# <fct>          <fct>        <int>      <int> <dbl>
# 1 7              acinar        1055       1089 0.969
# 2 3              alpha          234        234 1    
# 3 4              alpha          332        430 0.772
# 4 6              alpha          461        462 0.998
# 5 8              alpha          481        481 1    
# 6 10             alpha          226        226 1    
# 7 12             alpha          644        644 1    
# 8 13             alpha         1051       1054 0.997
# 9 15             alpha            5          9 0.556
# 10 16             alpha          344        345 0.997
# 11 17             alpha         1255       1255 1    
# 12 1              beta           457        465 0.983
# 13 2              beta           309        309 1    
# 14 9              beta          1524       1542 0.988
# 15 11             delta          172        296 0.581
# 16 14             ductal         567        570 0.995
# 17 5              mes            147        248 0.593

save(sp.ga2, se.ga2, transfer.anchors, atacCluster2rnaPredictedCellType, file=file.path(dir,"Rdata","spga2_filtered2_clusterCelltype.Rdata"))

## replace accessibility count with gene expression in @gmat
refdata = GetAssayData(
    object = panc_all, 
    assay = "RNA", 
    slot = "data"
)

imputation = TransferData(
    anchorset = transfer.anchors, 
    refdata = refdata, 
    weight.reduction = se.ga2[["SnapATAC"]], 
    dims = 1:ndim
) # Transfer RNA-seq gene expression assay matrix based on ATAC reduction object through anchors


sp.ge2 = sp.ga2
sp.ge2@gmat = t(imputation@data)
save(sp.ge2, file=file.path(dir,"Rdata","spge2_filtered2_clusterCelltype.Rdata"))

sp.ge = sp.ga
sp.ge@gmat = t(imputation@data)
save(sp.ge, file=file.path(dir,"Rdata","spge_filtered2_clusterCelltype.Rdata"))

########## integrate scRNA with scATAC using known marker genes as anchor ##############
library(Seurat)
load(file.path(dir,"Rdata","spga_filtered2_clusterCelltype.Rdata"))
load("/mnt/isilon/sfgi/suc1/analyses/grant/scRNA/pancreaticCells/downloads/scRNAseq_SeuratData.RData") # panc_all

genes.gr = read_delim(
        "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene.bed", 
        col_names=c("chr","start","end","gene","dot","strand"),
        delim="\t",
        comment="#"
) %>% 
separate(gene,c("name","gene_id"), sep="\\+") %>% 
dplyr::select(-dot) %>% 
makeGRangesFromDataFrame(
        starts.in.df.are.0based=T,
        keep.extra.columns=T
)

### selected gene gr (genes.sel.gr 21 genes)
marker2cell = c(
        "CPA1"="acinar", "PRSS1"="acinar", "REG1A"="acinar", 
        "GCG"="alpha", "IRX2"="alpha", "ARX"="alpha",
        "INS"="beta", "MAFA"="beta", "IAPP"="beta","IGF2"="beta",
        "SST"="delta",
        "KRT19"="ductal",
        "VWF"="endo", "ESAM"="endo", "CD93"="endo",
        "GHRL"="eps",
        "COL1A1"="mes",
        "PPY"="pp"
)

variable.genes = VariableFeatures(object = panc_all)
marker2cell[names(marker2cell) %in% variable.genes]

genes.sel.gr = genes.gr[which(genes.gr$name %in% variable.genes)]
genome(genes.sel.gr) = "hg19"

### find anchorset
sp.ga3 = createGmatFromMat(
    obj=sp.ga, 
    input.mat="bmat",
    genes=genes.sel.gr,
    do.par=TRUE,
    num.cores=10
)

# integrate atac snap with scRNA
ndim=14 # pre-determined dim used for cluster (see above)
se.ga3 = snapToSeurat(
    obj=sp.ga3, 
    eigs.dims=1:ndim, 
    norm=TRUE,
    scale=TRUE
) # make snap atac to seurat format (gene x 10392 cell), two assay data (@assay$ATAC, @assay$ACTIVITY), ACTIVITY contains 146 genes

transfer.anchors = FindTransferAnchors(
    reference = panc_all, 
    query = se.ga3, 
    features = genes.sel.gr$name, 
    reference.assay = "RNA", 
    query.assay = "ACTIVITY", 
    reduction = "cca"
) # 990 anchors

### transfer expression data
refdata = GetAssayData(
    object = panc_all, 
    assay = "RNA", 
    slot = "data"
)

imputation = TransferData(
    anchorset = transfer.anchors, 
    refdata = refdata, 
    weight.reduction = se.ga3[["SnapATAC"]], 
    dims = 1:ndim
) # Transfer RNA-seq gene expression assay matrix based on ATAC reduction object through anchors


sp.ge = sp.ga
sp.ge@gmat = t(imputation@data)
# 


# plot transfer gene expression in sp.ge
igraph_cluster = as.matrix(sp.ge@gmat)[,names(marker2cell)] %>% 
as.data.frame() %>% 
as_tibble() %>% 
mutate(cluster=sp.ge@metaData$igraph_cluster) %>% 
gather(key="gene", value="exp", -cluster) %>% 
group_by(gene, cluster) %>% 
# summarise(med_acc = median(acc)) %>% 
summarise(med_exp = mean(exp)) %>% 
ungroup() %>% 
left_join(
        tibble(
                gene=names(marker2cell),
                cellType=marker2cell
        )
) %>% 
spread(key=cluster, value=med_exp) %>% 
arrange(cellType)

igraph_cluster_m = igraph_cluster %>% 
select(-cellType) %>% 
as.data.frame() %>% 
column_to_rownames("gene") %>% 
as.matrix()

igraph_cluster_scale_m = t(apply(igraph_cluster_m,1,scale))
colnames(igraph_cluster_scale_m) = colnames(igraph_cluster_m)
idx = sapply(
        1:nrow(igraph_cluster_scale_m), 
        function(x){!all(is.na(igraph_cluster_scale_m[x,]))}
) # remove genes without variability across clusters
igraph_cluster_scale_m = igraph_cluster_scale_m[idx,]
igraph_cluster_m = igraph_cluster_m[idx,]
igraph_cluster = igraph_cluster[idx,]

cellType_annotation <- HeatmapAnnotation(
    cellType = igraph_cluster %>% select(gene, cellType) %>% deframe(),
    # col = list(donor=deframe(donor_colors[,c(1,2)])),
    annotation_legend_param = list(
            title = "cellType",
            at =unique(igraph_cluster$cellType),
            labels = unique(igraph_cluster$cellType)
    ),
    which = "row"
)
h <- Heatmap(
        igraph_cluster_scale_m, 
        name = "marker gene scaled mean expression", 
        col = colorRamp2(c(-2,0,2), c("blue","white" ,"red"),space="RGB"),
        cluster_columns = T, 
        clustering_distance_columns=function(x) as.dist(1-cor(t(x))),
        clustering_method_columns="average",
        show_column_dend= T,
        show_column_names = T, 
        cluster_rows=F,
        left_annotation=cellType_annotation,
)
X11.options(colortype="pseudo.cube")
pdf(file.path(dir,"plots","markerGene_igraphCluster_meanExp_heatmap_2000feature.pdf"))
draw(h)
dev.off()
X11.options(reset = T)

cellType_assignment = list(
        acinar=7,
        alpha=c(6,13,12,17,10,8,16), # ?4
        beta=c(2,9,1),
        delta=11, # SST
        ductal=14, # KRT19
        endo=5, # c(5,15) # VMF,ESAM, CD9
        eps=4, # c(6,13,12,17,10,8,16),
        mes=15, # c(5,15) #  COL1A1
        pp=3
)

### assign cell type
cellType_assignment = list(
        acinar=7,
        alpha=c(6,13,12,17,10,8,16), # ?4
        beta=c(2,9,1),
        delta=11, # SST
        ductal=14, # KRT19
        endo=5, # c(5,15) # VMF,ESAM, CD9
        eps=4, # c(6,13,12,17,10,8,16),
        mes=15, # c(5,15) #  COL1A1
        pp=3
)

sp.ga@metaData = sp.ga@metaData %>% 
as_tibble() %>% 
left_join(
        enframe(cellType_assignment, name="cellType", value="igraph_cluster") %>% 
        unnest(igraph_cluster) %>% 
        mutate(igraph_cluster=factor(igraph_cluster))
)
sp.ga@cluster = factor(sp.ga@metaData$cellType)
# acinar  alpha   beta  delta ductal   endo    eps    mes     pp 
#   1120   4491   2327    747    572    277    531     92    234 
  
sp.ge@metaData = sp.ge@metaData %>% 
as_tibble() %>% 
left_join(
        enframe(cellType_assignment, name="cellType", value="igraph_cluster") %>% 
        unnest(igraph_cluster) %>% 
        mutate(igraph_cluster=factor(igraph_cluster))
)
sp.ge@cluster = factor(sp.ge@metaData$cellType)

save(sp.ga, file=file.path(dir,"Rdata","spga_filtered2_clusterCelltype.Rdata"))
save(sp.ge, file=file.path(dir,"Rdata","spge_filtered2_clusterCelltype.Rdata"))
