# 4.0.2
# https://bioconductor.org/packages/release/bioc/vignettes/multiHiCcompare/inst/doc/multiHiCcompare.html
# BiocManager::install("multiHiCcompare")
# BiocManager::install("S4Vectors")
# BiocManager::install("BiocParallel")
library(dplyr)
library(tidyr)
library(ggplot2)
library(multiHiCcompare)
library(BiocParallel)
library(edgeR)


args = commandArgs(trailingOnly=TRUE)

dir=args[1]
out_prefix=args[2]
numCores=args[3]

# dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells/data/4000/chr1"
# out_prefix="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells/data/4000/chr1/HiCcompare"
# numCores=1

if(numCores!=1){
        register(SnowParam(workers = numCores), default = TRUE)
        paral=T
}else{
        paral=F
}


files=dir(dir, pattern="*.txt", full.names=T)

data = lapply(
        files,
        function(file){
                df <- read.table(file, header = FALSE)
                # The chromosome number should be entered as just the number
                # Chromosomes such as X, Y, etc. should be entered as 23, 24, etc
                df$V1 = gsub("chr","",df$V1)
                df$V4 = gsub("chr","",df$V4)
                df = df %>% as_tibble() %>% 
                filter(V1!="M", V4!="M") %>% 
                mutate(V1=case_when(
                        V1=="X" ~ "23", 
                        V1=="Y" ~ "24",
                        T ~ V1
                )) %>% 
                mutate(V4=case_when(
                        V4=="X" ~ "23", 
                        V4=="Y" ~ "24",
                        T ~ V4
                ))
                df %>% select(chr=V1, region1=V2, region2=V5, IF=V7)
                # sparse <- HiCcompare::cooler2sparse(df) # it will report error if there is only one chr
                # sparse # region1, region2 are start point (0-based)
        }
)

names(data) = sapply(
        basename(files),
        function(s){strsplit(s,"\\.")[[1]][1]}
)

groups = as.integer(
        factor(sapply(
        names(data),
        function(s){strsplit(s,"_")[[1]][1]}
))
)

groups = as.factor(groups-1)

# make hicexp
hicexp1 <- make_hicexp(data_list=data, 
                       groups = groups, 
                       zero.p = 0.8, A.min = 5, filter = TRUE,
                       remove.regions = hg19_cyto)

# The zero.p option allows for filtering by the proportion of zero IFs for an interaction. The A.min allows for filtering by a minimum average IF value. These options can be used together or individually to filter your data. Filtering is important to remove interactions with lots of 0 IFs and low average expression.

# normalize
hicexp1 <- fastlo(hicexp1, verbose = FALSE, parallel = paral)
# hicexp1 <- cyclic_loess(hicexp1, verbose = FALSE, parallel = FALSE)
normIF = hic_table(hicexp1)

# design
design <- model.matrix(~0+groups)
rownames(design) <- names(data)

# contrast (change)
my.contrasts <- makeContrasts(
        Alpha_vs_Acinar=groups1-groups0,
        Beta_vs_Acinar=groups2-groups0,
        Beta_vs_Alpha=groups2-groups1,
        levels=design
)



# GLM test
results = lapply(
        colnames(my.contrasts),
        function(contrast_name){
                tmp <- hic_glm(
                        hicexp1, design = design, contrast = my.contrasts[,contrast_name], 
                        method = "QLFTest", p.method = "fdr", parallel = paral
                )
                results = results(tmp)
                results
        }
)
names(results) = colnames(my.contrasts)

save(normIF, results, file=paste0(out_prefix,".data.Rdata"))
