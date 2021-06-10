# devtools::install_github('xuranw/MuSiC')
library(MuSiC)
library(Biobase) # for ExpressionSet
library(tidyverse)
library(Seurat)
dir="/mnt/isilon/sfgi/suc1/analyses/grant/scRNA/pancreaticCells"

####################### building ExpressionSet #######################
### sc
load(file.path(dir,"downloads/scRNAseq_SeuratData.RData"))

sc = as.matrix(panc_all@assays$RNA@counts)

pData = panc_all@meta.data %>% as_tibble() %>% 
mutate(cell=colnames(panc_all)) %>% 
select(samples=cell, clusters=cellType, donor) %>% 
mutate(cell=samples) %>% 
as.data.frame() %>% 
column_to_rownames("cell")

sc.eset = ExpressionSet(
        assayData=sc,
        phenoData=new("AnnotatedDataFrame", data=pData)
)

### bulk seq
# minimalSet <- ExpressionSet(assayData=exprs)
bulk_dir="/mnt/isilon/sfgi/suc1/analyses/grant/rnaSeq/pancreaticCells/HTseq/nostrand"

bulkseq = read_delim(file.path(bulk_dir, "htseq_counts.txt"), delim="\t")

selected_samples=tibble(sample=colnames(bulkseq)) %>% 
filter(grepl("HTSeq", sample)) %>% 
separate(sample, c("cell","mut","ind","dump"), remove=F) %>% 
filter(ind %in% c("HPAP040","HPAP045","HPAP054")) %>% 
pull(sample)

bulkseq = bulkseq[,c("gene", selected_samples)]

gene_type = read_delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt", delim="\t", comment="#")

bulkseq = bulkseq %>% 
dplyr::rename(gene_id=gene) %>% 
left_join(
        gene_type %>% 
        distinct(gene_name, gene_id)
) %>% filter(!is.na(gene_name))

# filter genes that not match by maximize the total count
bulkseq %>% 
filter(gene_name %in% rownames(sc)) %>% 
distinct(gene_name) # 30,295 rows, 29,906 unique gene_names

bulkseq = bulkseq %>% 
filter(gene_name %in% rownames(sc))

bulkseq = bulkseq %>% 
mutate(total=rowSums(select(bulkseq, -gene_name, -gene_id))) %>% 
group_by(gene_name) %>% 
dplyr::slice(which.max(total)) %>% 
ungroup() %>% select(-total)

bulk_pData = tibble(sample=colnames(select(bulkseq, -gene_name, -gene_id))) %>% 
separate(sample, c("cell","mut","ind", "dump"), sep="_", remove=F) %>% 
select(-dump) %>% 
as.data.frame() %>% 
column_to_rownames("sample")

bulk.eset = ExpressionSet(
        assayData=bulkseq %>% 
        select(-gene_id) %>% 
        as.data.frame() %>% 
        column_to_rownames("gene_name") %>% 
        as.matrix(),
        phenoData=new("AnnotatedDataFrame", data=bulk_pData)
)

save(bulk.eset, sc.eset, file=file.path(dir,"Rdata","MuSiC.Rdata"))

#################### Bulk Tissue Cell Type Estimation ###############
library(xbioc)

head(pData(sc.eset))

# using all genes
total_estProp = music_prop(
        bulk.eset = bulk.eset, 
        sc.eset = sc.eset, 
        clusters = 'clusters',samples = 'samples', 
        verbose = F
)
total_estProp$Est.prop.weighted # estimated proportions
save(bulk.eset, sc.eset, total_estProp, file=file.path(dir,"Rdata","MuSiC.Rdata"))

# using selected marker genes
marker2cell = c(
        "CPA1"="acinar", "PRSS1"="acinar", "REG1A"="acinar", 
        "GCG"="alpha", "IRX1"="alpha", "IRX2"="alpha", "ARX"="alpha",
        "INS"="beta","NKX6-1"="beta", "MAFA"="beta", "IAPP"="beta","IGF2"="beta",
        "SST"="delta",
        "KRT19"="ductal",
        "VWF"="endo", "ESAM"="endo", "CD93"="endo",
        "GHRL"="eps",
        "COL1A1"="mes",
        "PPY"="pp", "SLC6A4"="pp"
)

marker_estProp = music_prop(
        bulk.eset = bulk.eset, 
        sc.eset = sc.eset, 
        clusters = 'clusters',samples = 'samples', 
        markers=names(marker2cell),
        verbose = F
)

marker_estProp$Est.prop.weighted
save(bulk.eset, sc.eset, total_estProp, marker_estProp, file=file.path(dir,"Rdata","MuSiC.Rdata"))

###################### plot estPop ######################
# marker_estProp
estProp = marker_estProp$Est.prop.weighted %>% 
as.data.frame() %>% 
rownames_to_column("sample") %>% 
as_tibble() %>% 
gather(key="cellType", value="prop", -sample) %>% 
arrange(sample) %>% 
separate(sample, c("cell","mut","ind","dump"), sep="_") %>% 
dplyr::select(-mut,  -dump) %>% 
mutate(sample=glue::glue("{cell}_{ind}"))

p <- ggplot(estProp) +
facet_wrap(~cell) +
geom_col(aes(x=ind, y=prop, fill=cellType)) +
xlab("donor") +
ylab("proportions of single cell cluster") +
coord_flip() +
theme_bw() +
theme(legend.position="bottom")

ggsave(file.path(dir,"plots","bulk_deconvultion_estProp_marker.pdf"), p, height=3)
write_csv(estProp, file=file.path(dir,"plots","bulk_deconvultion_estProp_marker.csv"))


# total_estProp
estProp = total_estProp$Est.prop.weighted %>% 
as.data.frame() %>% 
rownames_to_column("sample") %>% 
as_tibble() %>% 
gather(key="cellType", value="prop", -sample) %>% 
arrange(sample) %>% 
separate(sample, c("cell","mut","ind","dump"), sep="_") %>% 
dplyr::select(-mut,  -dump) %>% 
mutate(sample=glue::glue("{cell}_{ind}"))

estProp %>% 
filter(tolower(cell) == cellType)
p <- ggplot(estProp) +
facet_wrap(~cell) +
geom_col(aes(x=ind, y=prop, fill=cellType)) +
xlab("donor") +
ylab("proportions of single cell cluster") +
coord_flip() +
theme_bw() +
theme(legend.position="bottom")

ggsave(file.path(dir,"plots","bulk_deconvultion_estProp_total.pdf"), p, height=3)
write_csv(estProp, file=file.path(dir,"plots","bulk_deconvultion_estProp_total.csv"))





