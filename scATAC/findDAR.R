library(SnapATAC)
library(tidyverse)
library(GenomicRanges)
library(SummarizedExperiment)

dir="/mnt/isilon/sfgi/suc1/analyses/grant/scATAC/pancreaticCells"

##########################  using all peaks ############################
# read_peaks
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

#### add pmat
load(file.path(dir,"Rdata","spga_filtered2_clusterCelltype.Rdata")) # sp.ga
### create cell by peak matrix
sp.ga@file=file.path(dir, "snaptools", basename(sp.ga@file))
sp.ga = addPmatToSnap(sp.ga) # see commands_scATAC.sh to "add cell-by-peak matrix"
sp.ga = makeBinary(sp.ga, mat="pmat")
sp.ga # old Long's peak 323844, my peak is 200270
save(sp.ga, file=file.path(dir,"Rdata","spga_filtered2_withallCellPmat.Rdata"))


#### fragment_counts SummarizedExperiment object for chromVAR
load(file.path(dir,"Rdata","spga_filtered2_withallCellPmat.Rdata"))
data = sp.ga@pmat
rownames(data) = sp.ga@metaData$barcodeNew # rownames must be unique
fragment_counts = SummarizedExperiment(
        assays = list(counts = t(data)), # peak x cell sparseMatrix
        rowRanges = sp.ga@peak,
        colData = DataFrame(cellType=sp.ga@metaData$cellType, depth=Matrix::rowSums(data))
) # matrix is peak x single cell
save(fragment_counts, file=file.path(dir,"Rdata/chromVAR_fragment_counts.Rdata"))

### findDAR
load(file.path(dir,"Rdata","spga_filtered2_withallCellPmat.Rdata")) # sp.ga

clusters = as.character(unique(sp.ga@cluster))

DAR_df = lapply(
        clusters,
        function(cluster){
                DARs = findDAR(
                    obj=sp.ga,
                    input.mat="pmat",
                    cluster.pos=cluster,
                    cluster.neg.method="knn",
                    test.method="exactTest",
                    bcv=0.4, #0.4 for human, 0.1 for mouse
                    seed.use=10
                )

                DARs$FDR = p.adjust(DARs$PValue, method="BH")

                bind_cols(
                        peaks.df,
                        DARs %>% as_tibble() %>% 
                        mutate(cluster=cluster)
                ) %>% arrange(FDR)
        }
)
DAR_df = do.call("bind_rows", DAR_df)
save(DAR_df, file=file.path(dir,"Rdata","DAR_df_allCells.Rdata"))

####################### using 3 cell peaks #####################
#### read_peaks
peaks.df = read_delim(
        file.path(dir,"peaks2","pancreas_3cells.bed"), delim="\t",
        col_names=c("chr","start","end","cells")
) %>% 
mutate(id=row_number()) %>%
mutate(id=gsub("^","OCR_",id)) # 166,708 peaks

#### add pmat
load(file.path(dir,"Rdata","spga_filtered2_clusterCelltype.Rdata")) # sp.ga
sp.ga@file=file.path(dir, "snaptools_3cells", basename(sp.ga@file))
sp.ga = sp.ga[sp.ga@cluster %in% c("acinar", "alpha", "beta"),] # 7938 cells
sp.ga = addPmatToSnap(sp.ga) # see commands_scATAC.sh to "add cell-by-peak matrix"
sp.ga = makeBinary(sp.ga, mat="pmat")
sp.ga # 
save(sp.ga, file=file.path(dir,"Rdata","spga_filtered2_with3CellPmat.Rdata"))

#### fragment_counts SummarizedExperiment object for chromVAR
load(file.path(dir,"Rdata","spga_filtered2_with3CellPmat.Rdata"))
data = sp.ga@pmat
rownames(data) = sp.ga@metaData$barcodeNew # rownames must be unique
fragment_counts = SummarizedExperiment(
        assays = list(counts = t(data)), # peak x cell sparseMatrix
        rowRanges = sp.ga@peak,
        colData = DataFrame(cellType=sp.ga@metaData$cellType, depth=Matrix::rowSums(data))
) # matrix is peak x single cell
save(fragment_counts, file=file.path(dir,"Rdata/chromVAR_fragment_counts.3cells.Rdata"))

### findDAR
load(file.path(dir,"Rdata","spga_filtered2_with3CellPmat.Rdata")) # sp.ga

clusters = as.character(unique(sp.ga@cluster))

DAR_df = lapply(
        clusters,
        function(cluster){
                DARs = findDAR(
                    obj=sp.ga,
                    input.mat="pmat",
                    cluster.pos=cluster,
                    cluster.neg.method="knn",
                    test.method="exactTest",
                    bcv=0.4, #0.4 for human, 0.1 for mouse
                    seed.use=10
                )

                DARs$FDR = p.adjust(DARs$PValue, method="BH")

                bind_cols(
                        peaks.df,
                        DARs %>% as_tibble() %>% 
                        mutate(cluster=cluster)
                ) %>% arrange(FDR)
        }
)
DAR_df = do.call("bind_rows", DAR_df)
save(DAR_df, file=file.path(dir,"Rdata","DAR_df_3Cells.Rdata"))
