library(chromVAR) # BiocManager::install("chromVAR"); devtools::install_github("sckinta/chromVAR")
library(SummarizedExperiment)
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19) # BiocManager::install("BSgenome.Hsapiens.UCSC.hg19"), v1.4.3
library(JASPAR2020) # BiocManager::install("JASPAR2020")
library(motifmatchr) # BiocManager::install("motifmatchr")
library(ComplexHeatmap)
library(circlize)

dir="/mnt/isilon/sfgi/suc1/analyses/grant/scATAC/pancreaticCells"

load(file.path(dir,"Rdata","chromVAR_withallCellPmat.Rdata"))
# fragment_counts_filtered, 
# jaspar_motifs_anno, motif_ix, motif_pos, 
# dev, 
# variability, 



### create peaks.gr
peaks.df = read_delim(
        file.path(dir,"peaks2","pancreas_9cells.bed"), delim="\t",
        col_names=c("chr","start","end","cells")
) %>% 
mutate(id=row_number()) # 200270 peaks

peaks.gr = makeGRangesFromDataFrame(
        peaks.df %>% select(-cells), 
        keep.extra.columns=T,
        starts.in.df.are.0based=T
)
### dev_zscore
dev_zscore = deviationScores(dev) %>% 
as.data.frame() %>% 
rownames_to_column("motif") %>% 
as_tibble() %>% 
gather(key="cell", value="dev_zscore", -motif) %>% 
left_join(
        colData(dev) %>% 
        as.data.frame() %>% 
        rownames_to_column("cell") %>% 
        as_tibble() %>% 
        select(cell, cellType)
) %>% 
group_by(cellType, motif) %>% 
summarise(dev_zscore = mean(dev_zscore)) %>% 
ungroup() %>% 
arrange(motif)

dev_zscore = dev_zscore %>% 
left_join(jaspar_motifs_anno)

### variability
variability = variability %>% 
mutate(name=toupper(name)) %>% 
dplyr::rename(TF=name) %>% 
filter(p_value_adj < 0.05) %>% 
arrange(desc(variability))

# #### sample correlation
# sample_cor = getSampleCorrelation(dev)
# any(sample_cor, )
# library(pheatmap)
# pdf(file.path(dir,"plots","chromVAR_sample_correlation.pdf"))
# pheatmap(
#          sample_cor, 
#          annotation_row = colData(dev), 
#          clustering_distance_rows = as.dist(1-sample_cor), 
#          clustering_distance_cols = as.dist(1-sample_cor)
# )
# dev.off()
# 
# 
# #### Cell similarity tSNT
# tsne_results = deviationsTsne(dev, threshold = 1.5, perplexity = 10, shiny = FALSE)
# 
# tsne_results = tsne_results %>% 
# as.data.frame() %>% 
# rownames_to_column("cell") %>% 
# as_tibble() %>% 
# left_join(
#         colData(dev) %>% 
#         as.data.frame() %>% 
#         rownames_to_column("cell") %>% 
#         as_tibble() %>% 
#         select(cell, cellType)
# )
# p <- tsne_results %>% 
# ggplot(aes(x=V1, y=V2, color=cellType)) +
# geom_point() +
# theme_bw() +
# xlab("tSNE dim1") + xlab("tSNE dim2")
# ggsave(file.path(dir,"plots","chromVAR_cellSimilarity_tSNE.pdf"),p)

#### Cell similarity correlation
dev_zscore_m = dev_zscore %>% 
select(cellType, motif, dev_zscore) %>% 
spread(key=cellType, value=dev_zscore) %>% 
as.data.frame() %>% 
column_to_rownames("motif")

res1 <- cor(dev_zscore_m)
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(dev_zscore_m) # matrix of the p-value of the correlation
library(corrplot)
library(RColorBrewer)
pdf(file.path(dir,"plots","chromVAR_cellSimilarity_correlation.pdf"))
col <- colorRampPalette(c("blue", "white","red"))
corrplot(res1, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         p.mat = p.mat, sig.level = 0.01, insig = "blank", # Combine with significance
         diag=FALSE # hide correlation coefficient on the principal diagonal
         )
dev.off()

#### Differential accessibility and variability
dev_colData = colData(dev) %>% 
as.data.frame() %>% 
rownames_to_column("cell") %>% 
as_tibble() %>% 
mutate(tmp=1) %>% 
spread(key=cellType, value=tmp) %>% 
mutate_if(function(x){any(is.na(x))}, function(x){factor(ifelse(is.na(x),0,1))})
dev_colData_cellList = lapply(
        unique(colData(dev)$cellType),
        function(cell){
                deframe(dev_colData %>% select(cell,{cell}))
        }
)
names(dev_colData_cellList) = unique(colData(dev)$cellType)

dev@colData@listData = c(dev@colData@listData, dev_colData_cellList)

diff_acc = lapply(
        unique(colData(dev)$cellType),
        function(cell){
                diff_acc = differentialDeviations(dev, groups=cell, alternative="less")
                diff_acc = diff_acc %>% 
                rownames_to_column("motif") %>% 
                as_tibble() %>% 
                mutate(cell=cell) %>% 
                left_join(jaspar_motifs_anno) %>% 
                arrange(p_value_adjusted)
                diff_acc
        }
)
diff_acc = do.call("bind_rows", diff_acc)
diff_acc$p_value_adjusted = p.adjust(diff_acc$p_value, method = "BH")

# calculate zscore_change
cal_zscore_diff <- function(cell){
        df = dev_colData %>% 
        select(cell, detection={cell}) %>% 
        left_join(
                deviationScores(dev) %>% 
                as.data.frame() %>% 
                rownames_to_column("motif") %>% 
                as_tibble() %>% 
                gather(key="cell", value="dev_zscore", -motif)
        ) %>% 
        group_by(detection, motif) %>% 
        summarise(dev_zscore = mean(dev_zscore)) %>% 
        ungroup() %>% 
        mutate(detection=ifelse(detection==0,"no","yes")) %>% 
        spread(key=detection, value=dev_zscore) %>% 
        mutate(zscore_diff=yes-no) %>% 
        mutate(cell=cell)
        df
}
zscore_diff = do.call(
        "bind_rows",
        lapply(
                unique(colData(dev)$cellType), 
                cal_zscore_diff
        )
)

diff_acc = left_join(
        diff_acc,
        zscore_diff
) 

save(diff_acc, dev_zscore, file=file.path(dir, "Rdata", "tmp.Rdata"))

diff_acc %>% 
filter(zscore_diff >= 0.5, p_value_adjusted < 1e-50) %>% 
distinct(motif)

diff_acc %>% 
filter(zscore_diff >= 0.5, p_value_adjusted < 0.05) %>% 
distinct(motif)

diff_acc %>% 
filter(zscore_diff >= 0.5, p_value_adjusted < 0.05) %>% 
group_by(cell) %>% 
summarise(n_distinct(motif))

diff_acc %>% 
filter(zscore_diff >= 2, p_value_adjusted < 0.05) %>% 
filter(TF=="ARX") %>% 
select(cell, p_value_adjusted,yes, zscore_diff)

dev_zscore_df = dev_zscore %>% 
semi_join(
        diff_acc %>% 
        filter(p_value_adjusted < 0.05) %>% 
        distinct(motif)
) %>% 
spread(key=cellType, value=dev_zscore) %>% 
left_join(
        jaspar_motifs_anno
)

dev_zscore_df_m = dev_zscore_df %>% 
select_if(is.numeric) %>% 
as.matrix()
rownames(dev_zscore_df_m) = glue::glue("{dev_zscore_df$TF}|{dev_zscore_df$motif}")

dev_zscore_df_scale_m = t(apply(dev_zscore_df_m, 1, scale))
colnames(dev_zscore_df_scale_m) = colnames(dev_zscore_df_m)
dev_zscore_df_scale_m = dev_zscore_df_scale_m[,c("pp","acinar","ductal","endo","mes","alpha","eps","beta","delta")]

idx=sapply(
        rownames(dev_zscore_df_scale_m),
        function(name){
                TF_name=strsplit(name,"\\|")[[1]][1]
                if(TF_name %in% marker_TFs){
                        detect=T
                }else{
                        if (any(sapply(marker_TFs, function(x){grepl(TF_name, x)})) & any(sapply(marker_TFs, function(x){grepl(x,TF_name)}))){
                                detect=T
                        }else{
                                detect=F
                        }
                }
        }
)
rownames(dev_zscore_df_scale_m)[idx]

marker_TF1 = list(
        beta=c("PDX1","ISL1","NFATC","MAFB","FOXA2","RFX1","MEIS","NEUROD1","PAX6","MAFA","CREB"),
        alpha=c("NEUROD1","MEIS","RFX1","MAFB","ISL1","CREB","PAX6","ARX","FOXA2","GATA6"),
        ductal=c("SMAD","STAT3","TEAD","FOXP1","FOXA2","GATA6","GATA4")
)

marker_TF2 = c(
        "NEUROD1","INSM1","ISL1","NKX2-2","PAX6","PDX1","ETV1","MEIS2",
        "EGR4","MAFA","NKX6-1","SIX3","OLIG1","MAFB","IRX1","IRX2","ARX","POU6F2",
        "FEV","HHEX","POU3F1","NEUROD3","SIX2","ESR1","RXRG"
)

marker_TFs=unique(c(unlist(marker_TF1, use.names=F), marker_TF2))

ha_rownames = rowAnnotation(
        foo = anno_mark(
                at = which(idx), 
                labels = sapply(
                        rownames(dev_zscore_df_scale_m)[idx], function(x){strsplit(x,"\\|")[[1]][1]}
                )
        )
)

h <- Heatmap(
        dev_zscore_df_scale_m, 
        name = "scaled deviation z-score", 
        col = colorRamp2(c(-2,0,2), c("blue","white" ,"red"),space="RGB"),
        cluster_columns = F, 
        show_column_dend= F,
        show_column_names = T, 
        cluster_rows=T,
        clustering_distance_rows=function(x) as.dist(1-cor(t(x))),
        clustering_method_rows="complete",
        show_row_names = F, 
        use_raster=T,
        right_annotation = ha_rownames
)
X11.options(colortype="pseudo.cube")
pdf(file.path(dir,"plots","chromVAR_motif_devzscore_cluster.pdf"))
# set.seed(123)
p = draw(h)
dev.off()
X11.options(reset = T)

write_csv(
        dev_zscore %>% 
        spread(key=cellType,dev_zscore),
        file=file.path(dir,"plots","chromVAR_motif_devzscore_cluster.csv")
)

save(dev, dev_zscore, diff_acc, jaspar_motifs_anno,umap_zscore, file=file.path(dir,"Rdata","chromVAR_dev_plotData.Rdata"))

#### plot dev zscore on umap
library(SnapATAC)
load(file.path(dir,"Rdata","spga_filtered2_withallCellPmat.Rdata"))
umap = bind_cols(
        sp.ga@metaData %>% 
        select(barcodeNew) %>% 
        dplyr::rename(cell=barcodeNew),
        sp.ga@umap %>% 
        as_tibble() %>% 
        dplyr::rename(umap1=`umap-1`,umap2=`umap-2`)
)

umap_zscore = left_join(
        umap, 
        deviationScores(dev) %>% 
        as.data.frame() %>% 
        rownames_to_column("motif") %>% 
        as_tibble() %>% 
        gather(key="cell", value="dev_zscore", -motif)
) %>% filter(!is.na(dev_zscore))
rm(sp.ga)
# TF_name="ARX"
# jaspar_motifs_anno %>% filter(grepl("ETS",TF))
plot_TF_zscore_umap <- function(TF_name){
        p <- umap_zscore %>% 
        semi_join(
                jaspar_motifs_anno %>% 
                filter(TF==TF_name)
        ) %>% 
        mutate(dev_zscore = ifelse(dev_zscore > 3, 3, dev_zscore)) %>% 
        mutate(dev_zscore = ifelse(dev_zscore < -3, -3, dev_zscore)) %>% 
        ggplot(aes(x=umap1, y=umap2)) +
        geom_point(aes(color=dev_zscore), size = 0.2, alpha=0.2) +
        scale_color_gradient2(mid = "lightgray", low = "blue", high = "red", breaks=seq(-3,3)) +
        theme_bw() +
        theme(
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank()
        )
        filename=paste0('chromVAR_zscore_umap_',TF_name,'.pdf')
        ggsave(file.path(dir, 'plots', filename), p, width=7.5)
}
plot_TF_zscore_umap("MAFB")
plot_TF_zscore_umap("ARX")
plot_TF_zscore_umap("PDX1")
plot_TF_zscore_umap("ETS1")
plot_TF_zscore_umap("GATA4")
plot_TF_zscore_umap("ONECUT1")
plot_TF_zscore_umap("NR5A2")
plot_TF_zscore_umap("FOXA1")
plot_TF_zscore_umap("FOXA2")
plot_TF_zscore_umap("ISL1")
plot_TF_zscore_umap("NKX6-1")
plot_TF_zscore_umap("PAX6")
plot_TF_zscore_umap("PAX4")

jaspar_motifs_anno %>% filter(grepl("NKX6-1",TF))