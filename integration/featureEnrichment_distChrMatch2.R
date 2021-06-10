library(tidyverse)
# devtools::install_github("sckinta/myRpackages", ref = "master", subdir = "DFbedtools")
library(DFbedtools)

args = commandArgs(trailingOnly=T)

bedpe_file=args[1]
feature_bed_list=args[2]
dump_dir=args[3]
out_file=args[4]

# bedpe_file="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/pancreaticCells/consensus_loops/mustache_fithicFDR1e6/Acinar_2reps.res4000.bedpe"
# feature_bed_list="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/pancreaticCells/featureEnrichment/Acinar.featureList.txt"
# dump_dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells/data/4000"
# out_file="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/pancreaticCells/featureEnrichment/Acinar_2reps.res4000.atac.enrichment.tsv"

# read bedpe and feature_bed
bedpe = read_delim(
        bedpe_file, delim="\t", 
        col_names=c("chr_a","start_a","end_a","chr_b","start_b","end_b")
)

featureList_files = read_delim(
        feature_bed_list, delim="\t",
        col_names=c("feature","file")
) %>% deframe()

featureList = lapply(
        featureList_files, 
        function(file){
                read_delim(file, delim="\t", col_names=c("chr","start","end"))
        }
)

resol=bedpe %>% mutate(resol=end_a-start_a) %>% distinct(resol) %>% pull(resol)

sampleN=100

### sampling using normIF.Rdata
# sampleN=100
# chr="chr22"
sample_nosig_per_chr <- function(chr, sampleN=100){
        load(file.path(dump_dir,chr,"normIF.Rdata"))

        signif = normIF %>% as_tibble() %>% 
        mutate(chr=gsub("^","chr",chr)) %>% 
        semi_join(
                bedpe %>% 
                dplyr::rename(chr=chr_a, region1=start_a, region2=start_b)
        ) %>% select(chr, region1, region2, D)

        nosignif = normIF %>% as_tibble() %>% 
        mutate(chr=gsub("^","chr",chr)) %>% 
        anti_join(signif) %>% 
        select(chr, region1, region2, D)

        rm(normIF)

        D_counts = signif %>% 
        group_by(D) %>% 
        count() %>% 
        ungroup()

        matched_nosignif_list = vector(mode = "list", length = sampleN)
        for (i in 1:sampleN){
                set.seed(i)
                matched_nosignif_list[[i]]=do.call(
                        "bind_rows",
                        lapply(
                                D_counts$D,
                                function(myD){
                                        df = nosignif %>% filter(D==myD) 
                                        if (nrow(df) > 0){
                                                df %>% 
                                                sample_n(D_counts %>% filter(D==myD) %>% pull(n), replace=T)
                                        }else{
                                                nosignif %>% 
                                                sample_n(D_counts %>% filter(D==myD) %>% pull(n), replace=T)
                                        }
                                }
                        )
                )
        }
        
        rm(nosignif)
        
        c(list(signif), matched_nosignif_list)
}

chrs=bedpe %>% distinct(chr_a) %>% pull(chr_a) 
# chrs=c("chr22","chrX")

int_list = lapply(
        chrs,
        sample_nosig_per_chr
)

names(int_list) = chrs

int_list = lapply(
        seq(1,sampleN+1),
        function(i){
                do.call(
                        "bind_rows",
                        lapply(
                                chrs,
                                function(chr){
                                        int_list[[chr]][[i]]
                                }
                        )
                )
        }
)


### overlap with feature_bed
loop_df = int_list[[1]]

overlap_with_feature <- function(loop_df, feature){
        loopEnd_df = bind_rows(
                loop_df %>% 
                select(chr,start=region1),
                loop_df %>% 
                select(chr,start=region2)
        ) %>% mutate(end=start+resol) %>% 
        distinct()

        loopEnd2feature = overlap_df(loopEnd_df, feature)

        tibble(
                loopEnd_n = nrow(loopEnd_df),
                feature_loopEnd_n=loopEnd2feature$overlap_df1 %>% 
                distinct() %>% count() %>% pull(n)
        ) %>% mutate(feature_loopEnd_ratio = feature_loopEnd_n/loopEnd_n)
}

results = do.call(
        "bind_rows",
        lapply(
                names(featureList),
                function(feature){
                        do.call(
                                "bind_rows",
                                lapply(
                                        int_list,
                                        overlap_with_feature,
                                        feature=featureList[[feature]]
                                )
                        ) %>% mutate(sample=c("sig",rep("nonsig",sampleN))) %>% 
                        mutate(feature=feature)
                }
        )
)


write.table(
        results,
        file=out_file,
        sep="\t", row.names=F, quote=F
)


summary = results %>% group_by(sample, feature) %>% 
summarise(feature_loopEnd_ratio = mean(feature_loopEnd_ratio)) %>% 
ungroup() %>% 
spread(key=sample, value=feature_loopEnd_ratio) %>% 
left_join(
        results %>% 
        filter(sample=="nonsig") %>% 
        group_by(feature) %>% 
        summarise(
                mean=mean(feature_loopEnd_ratio),
                sd = sd(feature_loopEnd_ratio)
        ) %>% 
        mutate(
                nonsig_up=mean+sd,
                nonsig_down=mean-sd,
        )
) %>% 
mutate(
        FC=sig/nonsig,
        FC_up=sig/nonsig_down,
        FC_down=sig/nonsig_up
) %>% select(feature, contains("FC"))

write.table(
        summary,
        file=paste0(out_file,".summary.tsv"),
        sep="\t", row.names=F, quote=F
)
