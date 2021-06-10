library(tidyverse)
library(parseIbed)

dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells"
juicer_dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/pancreaticCells"

resols=c("1K","2000","4000")
comparisons=c("Alpha_vs_Acinar","Beta_vs_Acinar","Beta_vs_Alpha")

#################### sigDE_consensus_loop ###########################
extract_sigDE_consensus_loop <- function(resol, comp){
        load(file.path(dir,"data",resol,paste0(comp,".sigDE.Rdata")))
        if(resol=="1K"){resol="1000"}
        consensus_loop=vroom::vroom(
                file.path(juicer_dir,"consensus_loops/mustache_fithicFDR1e6",paste0("consensus_loops.res",resol,".bedpe")),
                col_names=c("chr_a","start_a","end_a","chr_b","start_b","end_b")
        )

        df = sigDE %>% as_tibble() %>% 
        dplyr::rename(chr_a=chr, start_a=region1,start_b=region2) %>% 
        mutate(chr_a=gsub("^","chr",chr_a)) %>% 
        left_join(consensus_loop) %>% 
        filter(!is.na(end_a))

        if (nrow(df)==0){
                df=df %>% 
                mutate(resol=as.character()) %>% 
                mutate(comp=as.character())
        }else{
                df=df %>% 
                mutate(resol=resol) %>% 
                mutate(comp=comp)
        }
        df
}

grid_set=expand_grid(resol=resols, comp=comparisons)

sigDE_consensus_loop = do.call(
        "bind_rows",
        mapply(
                extract_sigDE_consensus_loop,
                resol=grid_set$resol,
                comp=grid_set$comp,
                SIMPLIFY=F
        )
)

# summary
sigDE_consensus_loop_summary = grid_set %>% 
left_join(
        sigDE_consensus_loop %>% 
        group_by(resol, comp) %>% 
        count(name="loopN") %>% 
        ungroup()
) %>% 
mutate(loopN=ifelse(is.na(loopN),0,loopN))

#################### annotate sigDE_consensus_loop ###########################
data(prom_bed)
ocr_bed = vroom::vroom(
        "/mnt/isilon/sfgi/suc1/analyses/grant/scATAC/pancreaticCells/peaks/consensus_3cells.bed",
        col_names=c("chr","start","end","id")
) %>% mutate(id=paste0("OCR_",row_number()))

sigDE_consensus_loop_anno = annotate_bedpe2geneOCR(sigDE_consensus_loop, ocr_bed, prom_bed)

sigDE_consensus_loop_anno = sigDE_consensus_loop_anno %>%
left_join(
        prom_bed %>% 
        select(anno_a=pro_anno, gene_a=gene_id)
) %>% 
left_join(
        prom_bed %>% 
        select(anno_b=pro_anno, gene_b=gene_id)
)

# summary
sigDE_consensus_loop_summary = sigDE_consensus_loop_summary %>% 
left_join(
        sigDE_consensus_loop_anno %>%
        distinct(chr_a, start_a, end_a, chr_b, start_b, end_b, resol, comp) %>% 
        group_by(resol,comp) %>%
        count(name="anno_loopN") %>%
        ungroup()
) %>% 
mutate(anno_loopN=ifelse(is.na(anno_loopN),0,anno_loopN))

sigDE_consensus_loop_summary = sigDE_consensus_loop_summary %>% 
left_join(
        bind_rows(
                sigDE_consensus_loop_anno %>%
                select(ocr=ocr_a, gene=gene_b, resol, comp),
                sigDE_consensus_loop_anno %>%
                select(ocr=ocr_b, gene=gene_a, resol, comp)
        ) %>% 
        filter(!is.na(ocr), !is.na(gene)) %>% 
        distinct() %>% 
        group_by(resol,comp) %>%
        summarise(
                geneOCR_N=n_distinct(row_number()),
                gene_N=n_distinct(gene),
                OCR_N=n_distinct(ocr)
        ) %>% 
        ungroup()
) %>% 
mutate_at(c("geneOCR_N","gene_N","OCR_N"), function(x){ifelse(is.na(x),0,x)}) %>% 
arrange(comp)

# add original global differential
sigDE_consensus_loop_summary = do.call(
        "bind_rows",
        lapply(
                resols,
                function(resol){
                        read_csv(
                                file.path(dir,"data",resol,"sigDE_summary.csv")
                        ) %>% 
                        dplyr::rename(comp=comparison) %>% 
                        mutate(resol=resol)
                }
        )
) %>% 
left_join(
        sigDE_consensus_loop_summary
) %>% 
arrange(comp)

write_csv(
        sigDE_consensus_loop_summary,
        file=file.path(dir,"tables","sigDE_consensus_loop_summary.csv")
)

sigDE_geneOCR_summary = bind_rows(
        sigDE_consensus_loop_anno %>%
        select(ocr=ocr_a, gene=gene_b, comp),
        sigDE_consensus_loop_anno %>%
        select(ocr=ocr_b, gene=gene_a, comp)
) %>% 
filter(!is.na(ocr), !is.na(gene)) %>% 
distinct() %>% 
group_by(comp) %>%
summarise(
        geneOCR_N=n_distinct(row_number()),
        gene_N=n_distinct(gene),
        OCR_N=n_distinct(ocr)
) %>% 
ungroup()

write_csv(
        sigDE_geneOCR_summary,
        file=file.path(dir,"tables","sigDE_geneOCR_summary.csv")
)

save(sigDE_geneOCR_summary, sigDE_consensus_loop_summary, sigDE_consensus_loop, sigDE_consensus_loop_anno, file=file.path(dir,"Rdata","sigDE_consensus_loop.Rdata"))



