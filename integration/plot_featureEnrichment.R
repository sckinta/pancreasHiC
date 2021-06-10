library(tidyverse)

dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/pancreaticCells/featureEnrichment"

## plot featureEnrichment
chroHMM_annotation = read_delim("/mnt/isilon/sfgi/suc1/analyses/grant/chipSeq/pancreaticCells/downloads/processed/emission_annotation.txt", delim="\t", col_names=c("feature","color","annotation"), col_types="ccc")

files=dir(dir,pattern="summary\\.tsv", full.names=T)
enrichmentSummary = do.call(
        "bind_rows",
        lapply(
                files,
                function(file){
                        read_delim(file, delim="\t") %>% 
                        mutate(file=gsub(".featureEnrichment.tsv.summary.tsv","",basename(file))) %>% 
                        separate(file, c("cell","resol"), sep="\\.") %>% 
                        left_join(chroHMM_annotation) %>% 
                        mutate(annotation=ifelse(is.na(annotation),"open chromatin", annotation)) %>% 
                        mutate(color=ifelse(feature=="E8", "255,160,0", color)) %>% 
                        mutate(color=ifelse(is.na(color),"0,191,255",color)) %>% 
                        mutate(color=sapply(strsplit(color, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))) %>% mutate(annotation=glue::glue("{annotation} ({feature})"))
                }
        )
)

colors=enrichmentSummary %>% 
mutate(feature=factor(feature, levels=c(paste0("E",1:10),"atac"))) %>% 
arrange(feature) %>% 
mutate(color=factor(color, levels=unique(color))) %>% 
select(color) %>% 
distinct(color) %>% 
mutate(color=as.character(color)) %>% 
pull(color)

p <- enrichmentSummary %>% 
mutate(feature=factor(feature, levels=c(paste0("E",1:10),"atac"))) %>% 
mutate(cell=gsub("_2reps","",cell)) %>% 
arrange(feature) %>% 
mutate(annotation=factor(annotation, levels=unique(annotation))) %>% 
ggplot(aes(x = annotation)) +
facet_wrap(~cell) + 
geom_col(aes(y=FC, fill=annotation), color="black") +
geom_errorbar(aes(ymax=FC_up, ymin=FC_down), width=0.5) +
geom_hline(yintercept=1, linetype=3) +
theme_bw() +
xlab("") + ylab("enrichment fold change") +
coord_flip() +
scale_fill_manual(values=colors) +
theme(legend.position="none")
ggsave(file.path(dir,"featureEnrichment_chroHMM_atac.pdf"), p, height=3, width=10)

### write table 
files=dir(dir,pattern="\\.featureEnrichment\\.tsv$", full.names=T)
LoopEnd_feature_Overlap = lapply(
        files,
        function(file){
                read_delim(file, delim="\t") %>% 
                mutate(file=gsub(".featureEnrichment.tsv","",basename(file))) %>% 
                separate(file, c("cell","resol"), sep="\\.") %>% 
                left_join(chroHMM_annotation) %>% 
                mutate(annotation=ifelse(is.na(annotation),"open chromatin", annotation)) %>% 
                mutate(color=ifelse(feature=="E8", "255,160,0", color)) %>% 
                mutate(color=ifelse(is.na(color),"0,191,255",color)) %>% 
                select(-resol)
        }
)
LoopEnd_feature_Overlap = do.call(
        "bind_rows",
        LoopEnd_feature_Overlap
)
write_csv(
        LoopEnd_feature_Overlap,
        file = file.path(dir,"featureEnrichment_chroHMM_atac.csv")
)

## plot chroHMM overlap loopEnd count
p <- LoopEnd_feature_Overlap %>% 
filter(sample == "sig", feature!="atac") %>% 
distinct(loopEnd_n, feature, feature_loopEnd_n, cell, annotation) %>% 
mutate(cell=gsub("_2reps","", cell)) %>% 
mutate(feature=factor(feature, levels=c(paste0("E",1:10)))) %>% 
arrange(feature) %>% 
mutate(annotation=factor(annotation, levels=unique(annotation))) %>% 
ggplot(aes(x=cell, y=feature_loopEnd_n)) +
geom_col(aes(fill=annotation), color="black", size=0.5) + 
theme_bw() +
scale_fill_manual(values=colors[1:(length(colors)-1)])
ggsave(file.path(dir,"LoopEnd_feature_Overlap_chroHMM_barplot.pdf"), p)

