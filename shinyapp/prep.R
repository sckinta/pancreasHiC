library(tidyverse)
library(DFbedtools)
library(GenomicRanges)

dir="/mnt/isilon/sfgi/suc1/analyses/grant/disease/pancreaticCells_ABC"
dir.create(dir,"tables")

######### reduce gene_id
gene_type = read_delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt", delim="\t", comment="#")

gene_type = gene_type %>%
filter(gene_type %in% c("protein_coding","noncoding_RNA_long")) %>%
select(-biotype)

prom_bed = prom_bed %>%
semi_join(
        gene_type %>%
        mutate(pro_anno=glue::glue("{gene_name}+{transcript_id}"))
)

reduce_gene_prom <- function(gene_id){
        df = prom_bed %>%
        filter(gene_id==!!gene_id) %>%
        makeGRangesFromDataFrame(
                seqnames.field="pro_chr",start.field="pro_start",end.field="pro_end",strand.field="strand",
                starts.in.df.are.0based=T
        )
        df = reduce(df) %>%
        as.data.frame() %>%
        as_tibble() %>%
        dplyr::rename(pro_chr=seqnames, pro_start=start, pro_end=end, strand=strand) %>%
        mutate(pro_start=pro_start-1) %>%
        mutate(gene_id=gene_id)
        df
}

prom_bed_reduced = do.call(
        "bind_rows",
        lapply(
                prom_bed %>% distinct(gene_id) %>% pull(gene_id),
                reduce_gene_prom
        )
)

gene_type_coord = left_join(
        gene_type %>%
        select(-transcript_id) %>%
        distinct(),
        prom_bed_reduced %>%
        select(-width) %>%
        mutate(strand=as.character(strand))
)

gz <- gzfile(file.path("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencodeV19_gene_type.protein_lncRNA.txt.gz"), "w")

write.table(
        gene_type_coord %>%
        mutate_at(c("pro_start","pro_end"), function(x){gsub(" ","",format(x,scientific=F))}),
        gz,
        sep="\t", row.names=F, quote=F
)


######### add gene promoter annotation to refOCR
load(file.path(dir,"Rdata","refOCR_cell.Rdata"))

gene_type_coord = gene_type_coord %>%
select(chr=pro_chr,start=pro_start,end=pro_end,strand,gene_id, gene_name)

ocr2promoter = overlap_df(
        refOCR_cell %>%
        distinct(chr,start,end,id),
        gene_type_coord,
        df1_0base=T, df2_0base=T, minoverlap=1L
)

ocr2promoter = bind_cols(
        ocr2promoter$overlap_df1,
        ocr2promoter$overlap_df2 %>%
        select(gene_id, gene_name)
)

refOCR_cell = refOCR_cell %>%
left_join(
        ocr2promoter %>%
        group_by(chr,start,end,id) %>%
        summarise(gene_ids=paste(gene_id, collapse="|")) %>%
        ungroup()
)

gz <- gzfile(file.path(dir,"tables","refOCR_cell_simple2.txt.gz"), "w")
write.table(
        refOCR_cell %>%
        mutate_at(c("start","end","id"), function(x){gsub(" ","",format(x,scientific=F))}),
        file=gz,
        sep="\t", row.names=F, quote=F
)
close(gz)
