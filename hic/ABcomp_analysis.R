# R/4.0.2
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DFbedtools)
# matt's ABcompartment folder: /mnt/isilon/sfgi/pahlm/analyses/grant/captureC/Promoterome/ABCompartments/homer

dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells"
juicer_parent_dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer"
conditions=c("Acinar_2reps", "Alpha_2reps", "Beta_2reps")

##################### write PC1 and PC2 ##################
files=dir(file.path(juicer_parent_dir,conditions,"ABcompartments"), pattern="40K.cis.vecs.tsv", full.names=T)

## PC1
PC1=lapply(
        files,
        function(file){
                read_delim(file, delim="\t") %>% 
                select(chrom,start,end,E1)
        }
)
names(PC1) = gsub(".eigen.40K.cis.vecs.tsv","",basename(files))
lapply(PC1, function(df){df %>% filter(!is.na(E1)) %>% count() %>% pull(n)})

PC1_tmp = do.call(
        "bind_cols",
        lapply(
                PC1,
                function(df){
                        df %>% select(E1)
                }
        )
)
colnames(PC1_tmp)=names(PC1)

PC1_df = bind_cols(
        PC1[[1]] %>% select(chrom,start,end),
        PC1_tmp
)
# PC1_df = PC1_df %>% na.omit()
rm(PC1_tmp, PC1)

## PC2
PC2=lapply(
        files,
        function(file){
                read_delim(file, delim="\t") %>% 
                select(chrom,start,end,E2)
        }
)
names(PC2) = gsub(".eigen.40K.cis.vecs.tsv","",basename(files))
lapply(PC2, function(df){df %>% filter(!is.na(E2)) %>% count() %>% pull(n)})

PC2_tmp = do.call(
        "bind_cols",
        lapply(
                PC2,
                function(df){
                        df %>% select(E2)
                }
        )
)
colnames(PC2_tmp)=names(PC2)

PC2_df = bind_cols(
        PC2[[1]] %>% select(chrom,start,end),
        PC2_tmp
)
# PC2_df = PC2_df %>% na.omit()
rm(PC2_tmp, PC2)

save(PC1_df, PC2_df, file=file.path(dir,"Rdata/PCA.Rdata"))

#################### plot chr PC1 #######################
dir.create(file.path(dir,"plots/AB_PCA"))
plot_chr_PC1 <- function(chr){
        tmp = PC1_df %>% filter(chrom==chr) %>% 
        select(-contains("_2reps")) %>% 
        gather(key="sample", value="PC1", -chrom, -start, -end) %>% 
        mutate(direction=ifelse(PC1 > 0, "A", "B")) %>% 
        mutate(direction=ifelse(is.na(direction),"B",direction))

        p <- ggplot(tmp, aes(x=start, y=PC1, fill=direction)) +
        facet_wrap(~sample, ncol=1, strip.position = "right") + 
        geom_col() +
        scale_fill_manual(values=c("black","grey")) +
        theme_classic() +
        # theme(
        #         axis.title.x=element_blank(),
        #         axis.text.x=element_blank(),
        #         axis.ticks.x=element_blank(),
        #         legend.position="none"
        # ) +
        ylim(c(-1,1))
        filename=paste0("sample_",chr,"_PC1.pdf")
        ggsave(file.path(dir,"plots/AB_PCA",filename))
}


lapply(
        paste0("chr",c(1:22)),
        plot_chr_PC1
)

lapply(
        paste0("chr",c("X","Y")),
        plot_chr_PC1
)
plot_chr_PC1("chr1")

################## plot jaccard heatmap #####################
tmp = PC1_df %>% 
select(-contains("_2reps")) %>% 
na.omit() %>% 
filter(chrom %in% paste0("chr",c(1:22)))

tmp_AB = tmp %>% 
gather(key="sample",value="PC1",-chrom, -start,-end) %>% 
mutate(Acomp=ifelse(PC1>0,1,0)) %>% 
select(-PC1) %>% 
spread(key=sample, value=Acomp)
tmp_AB_m=tmp_AB %>% select(-chrom, -start, -end) %>% as.matrix()
tmp_AB_dist = philentropy::distance(t(tmp_AB_m), method = "jaccard", use.row.names = TRUE)

pdf(file.path(dir,"plots/PC1_jaccard_heatmap.pdf"))
gplots::heatmap.2(1-tmp_AB_dist,
          distfun = function(x) as.dist(1-x),
          hclust=function(x) hclust(x,method="complete"),
          Rowv = TRUE, Colv= T,
          col=RColorBrewer::brewer.pal(9,"Blues"),
          symm = T,
          margins = c(12,12),
          trace="none",
          density.info="none",
          key.title = NA,
          keysize=1,
          key.xlab="jaccard index"
 )
dev.off()

library(corrplot)
library(RColorBrewer)
pdf(file.path(dir,"plots/PC1_jaccard_heatmap_corrplot.pdf"))
corrplot(1-tmp_AB_dist,
               is.corr=FALSE, type="upper",
               p.mat = 1-tmp_AB_dist, insig = "p-value", sig.level = -1,
               order="hclust",
               cl.lim=c(0.5, 1),
               col=colorRampPalette(c("blue","white","red"))(10)
)
dev.off()

#################### A/B switch #####################
tmp = PC1_df %>% 
select(-contains("_2reps")) %>% 
na.omit() %>% 
filter(chrom %in% paste0("chr",c(1:22)))

# filter both replicates same direction
tmp = tmp %>% 
filter(Acinar_1*Acinar_2 > 0, Alpha_1*Alpha_2 > 0, Beta_1*Beta_2 > 0)

# condition Acomp
tmp = tmp %>% 
gather(key="sample",value="PC1",-chrom,-start, -end) %>% 
mutate(condition=sapply(sample, function(x){strsplit(x,"_")[[1]][1]})) %>% 
group_by(chrom, start, end, condition) %>% 
summarise(PC1=mean(PC1)) %>% 
ungroup() %>% 
mutate(Acomp=ifelse(PC1>0,1,-1)) %>% 
select(-PC1) %>% 
spread(key=condition, value=Acomp)

# AB switch
tmp=tmp %>% 
mutate(
        Beta_vs_Alpha=case_when(
                (Beta*Alpha==1&Beta>0&Alpha>0) ~ "A-A",
                (Beta*Alpha==1&Beta<0&Alpha<0) ~ "B-B",
                (Beta*Alpha==-1&Beta>0&Alpha<0) ~ "A-B",
                (Beta*Alpha==-1&Beta<0&Alpha>0) ~ "B-A",
        ),
        Beta_vs_Acinar=case_when(
                (Beta*Acinar==1&Beta>0&Acinar>0) ~ "A-A",
                (Beta*Acinar==1&Beta<0&Acinar<0) ~ "B-B",
                (Beta*Acinar==-1&Beta>0&Acinar<0) ~ "A-B",
                (Beta*Acinar==-1&Beta<0&Acinar>0) ~ "B-A",
        ),
        Alpha_vs_Acinar=case_when(
                (Alpha*Acinar==1&Alpha>0&Acinar>0) ~ "A-A",
                (Alpha*Acinar==1&Alpha<0&Acinar<0) ~ "B-B",
                (Alpha*Acinar==-1&Alpha>0&Acinar<0) ~ "A-B",
                (Alpha*Acinar==-1&Alpha<0&Acinar>0) ~ "B-A",
        )
) %>% 
select(chrom,start,end, contains("_vs_"))

ABswitch = tmp
save(ABswitch, file=file.path(dir,"Rdata","ABswitch.Rdata"))

# tmp = tmp %>% 
# gather(key="compare",value="ABcomp",contains("_vs_"))
# 
# tmp_sum = tmp %>% 
# group_by(ABcomp,compare) %>% 
# count() %>% 
# ungroup() %>% 
# arrange(compare) %>% 
# mutate(ABcomp=factor(ABcomp,levels=c("A-A","B-B","A-B","B-A")))
# 
# tmp_sum = tmp_sum %>% 
# left_join(
#         tmp_sum %>% 
#         group_by(compare) %>% 
#         summarise(total=sum(n)) %>% 
#         ungroup()
# ) %>% 
# mutate(ratio=n/total) %>% 
# mutate(ratio=scales::percent(ratio))
# 
# p <- ggplot(tmp_sum, aes(x=ABcomp, y=n)) +
# facet_wrap(~compare) +
# geom_col(aes(fill=ABcomp)) +
# geom_text(aes(label=ratio), nudge_y=10) +
# theme_bw()+
# theme(legend.position="none") +
# xlab("A/B compartment switch") +
# ylab("number of 40kbp bins")
# 
# ggsave(file.path(dir,"plots","ABcomp_switch_boxplot.pdf"), p, width=9)

p <- ABswitch %>% 
gather(key="comparison", value="ABswitch", -chrom, -start, -end) %>% 
group_by(comparison, ABswitch) %>% 
count() %>% ungroup() %>% 
left_join(
        ABswitch %>% 
        gather(key="comparison", value="ABswitch", -chrom, -start, -end) %>% 
        group_by(comparison) %>% 
        count(name="total") %>%  ungroup()
) %>% 
mutate(ratio=n/total) %>% 
mutate(ABswitch=gsub("-","->",ABswitch)) %>% 
mutate(ABswitch=factor(ABswitch, levels=c("A->B","B->A","B->B","A->A"))) %>% 
ggplot(aes(x=comparison, y=ratio)) +
geom_col(aes(fill=ABswitch)) +
theme_bw() +
xlab("") +
coord_flip() +
theme(legend.position="bottom")
ggsave(file.path(dir,"plots","ABswitch_ratio_barplot.pdf"), p, height=4)

################### ABswitch and gene expression (bulk) #################
# https://www.nature.com/articles/nature14222
load(file.path(dir,"Rdata","ABswitch.Rdata"))
load(file.path(dir,"data/rnaseq/bulkRNA-seq/tpm.Rdata"))

tss = read_delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.nonpseudo.TSS.bed", delim="\t", comment="#", col_names=c("chr","start","end","tss_anno", "dot","strand")) %>% select(-dot)

gene_type=read_delim(
        "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt", delim="\t", comment="#"
) %>% 
filter(gene_type!="pseudogene")

gene_tss = tss %>% 
left_join(
        gene_type %>% 
        mutate(tss_anno=glue::glue("{gene_name}+{transcript_id}"))
) %>% select(-tss_anno, -transcript_id, -biotype) %>% 
distinct() %>% 
mutate(gene_id=sapply(gene_id, function(x){strsplit(x,"\\.")[[1]][1]}))

ABswitch2gene = overlap_df(
        ABswitch,
        gene_tss,
        df1_chr="chrom", df2_chr="chr",
        df1_0base=F,df2_0base=F, minoverlap=1L
)

ABswitch2gene = bind_cols(
        ABswitch2gene$overlap_df1,
        ABswitch2gene$overlap_df2 %>% 
        select(gene_id, gene_name)
) %>% distinct()


# contrast="Beta_vs_Alpha"
add_exp_log2FC <- function(contrast){
        tmp = ABswitch2gene %>% 
        select(gene_id, gene_name, ABswitch=!!sym(contrast)) %>% 
        distinct()

        # filter genes with more two ABswitch category
        tmp = semi_join(
                tmp,
                tmp %>% 
                group_by(gene_id, gene_name) %>% 
                summarise(type_n=n_distinct(ABswitch)) %>% 
                ungroup() %>% 
                filter(type_n == 1)
        )

        # add gene expression log2FC
        conditions_tmp=strsplit(contrast,"_")[[1]][c(1,3)]
        tmp = left_join(
                tmp, 
                tpm %>% filter(condition %in% conditions_tmp) %>% 
                group_by(gene_name, gene_id, condition) %>% 
                summarise(tpm=mean(tpm)) %>% 
                ungroup() %>% 
                spread(key=condition, value=tpm) %>% 
                mutate(log2FC=(!!sym(conditions_tmp[1])+1)/(!!sym(conditions_tmp[2])+1)) %>% 
                mutate(log2FC=log2(log2FC))
        )

        tmp = tmp %>% na.omit()

        tmp %>% mutate(contrast=contrast)
}

# since 
tmp = add_exp_log2FC("Beta_vs_Alpha")
tmp = tmp %>% 
mutate(ABswitch=case_when(
        ABswitch %in% c("A-A","B-B") ~ "Stable",
        ABswitch=="A-B" ~ "A-B",
        ABswitch=="B-A" ~ "B-A"
))
tmp = tmp %>% mutate(ABswitch=factor(ABswitch, levels=c("Stable","A-B","B-A")))

tmp_sum=tmp %>% group_by(ABswitch) %>% summarise(gene_n=n_distinct(gene_id))
p <- ggplot(tmp, aes(x=ABswitch, y=log2FC)) +
geom_violin(aes(fill=ABswitch)) +
geom_boxplot(outlier.shape=NA, width=0.1, fill="white") +
geom_text(aes(x=ABswitch,y=5,label=gene_n), data=tmp_sum) +
theme_bw() +
theme(legend.position="none")

ggsave(file.path(dir,"plots","ABswitch_geneExpLog2FC.Beta_vs_Alpha.pdf"), p, width=4)

tmp_list=split(tmp$log2FC,tmp$ABswitch)
pvalues=mapply(
        function(x,y){
                wilcox.test(tmp_list[[x]], tmp_list[[y]], alternative = "two.sided")$p.value
        },
        x=c("Stable","Stable"),
        y=c("A-B","B-A"),
        SIMPLIFY=T
)
tibble(
        x=c("Stable","Stable"),
        y=c("A-B","B-A"),
        pvalues=pvalues
)

# x      y      pvalues
# <chr>  <chr>    <dbl>
# 1 Stable A-B   5.25e-72
# 2 Stable B-A   3.67e-31

################### ABswitch and gene expression (single cell) #################
load(file.path(dir,"Rdata","ABswitch.Rdata"))
load("/mnt/isilon/sfgi/suc1/analyses/grant/scRNA/pancreaticCells/Rdata/avgExp.Rdata")

tss = read_delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.nonpseudo.TSS.bed", delim="\t", comment="#", col_names=c("chr","start","end","tss_anno", "dot","strand")) %>% select(-dot)

gene_type=read_delim(
        "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt", delim="\t", comment="#"
) %>% 
filter(gene_type!="pseudogene")

gene_tss = tss %>% 
left_join(
        gene_type %>% 
        mutate(tss_anno=glue::glue("{gene_name}+{transcript_id}"))
) %>% select(-tss_anno, -transcript_id, -biotype) %>% 
distinct() %>% 
mutate(gene_id=sapply(gene_id, function(x){strsplit(x,"\\.")[[1]][1]}))

ABswitch2gene = overlap_df(
        ABswitch,
        gene_tss,
        df1_chr="chrom", df2_chr="chr",
        df1_0base=F,df2_0base=F, minoverlap=1L
)

ABswitch2gene = bind_cols(
        ABswitch2gene$overlap_df1,
        ABswitch2gene$overlap_df2 %>% 
        select(gene_id, gene_name)
) %>% distinct()

ABswitch2gene = ABswitch2gene %>% 
gather(key="comp", value="ABswitch", -chrom, -start, -end, -gene_name, -gene_id)

avgExp = avgExp %>% ungroup() %>% 
dplyr::rename(gene_name = feature) %>% 
filter(cellType %in% c("acinar","alpha","beta")) %>% 
mutate(cellType=stringr::str_to_title(cellType)) 

expressedGenes = avgExp %>% 
group_by(gene_name) %>% 
dplyr::slice(which.max(expression)) %>% 
ungroup() %>% 
filter(expression > 0)
quantile(expressedGenes$expression,c(0.1,0.25,0.5,0.75,0.9))
# 10%          25%          50%          75%          90% 
# 0.0004313547 0.0019780898 0.0234087877 0.1281709839 0.3454756200

ABswitch2gene_exp = ABswitch2gene %>% 
semi_join(expressedGenes) %>% 
left_join(
        avgExp %>% 
        semi_join(expressedGenes %>% select(gene_name)) %>% 
        spread(key=cellType, value=expression) %>% 
        mutate(
                Beta_vs_Alpha=log2((Beta+0.001)/(Alpha+0.001)),
                Beta_vs_Acinar=log2((Beta+0.001)/(Acinar+0.001)),
                Alpha_vs_Acinar=log2((Alpha+0.001)/(Acinar+0.001))
        ) %>% 
        select(gene_name, contains("_vs_")) %>% 
        gather(key="comp", value="log2FC", -gene_name)
)

cal_wilcox_p <- function(df, AB_comp1, AB_comp2, cell_comp){
        val1 = df %>% 
        filter(comp==cell_comp, ABswitch==AB_comp1) %>% 
        pull(log2FC)
        val2 = df %>% 
        filter(comp==cell_comp, ABswitch==AB_comp2) %>% 
        pull(log2FC)
        test = wilcox.test(val1,val2)
        tibble(
                cell_comp=cell_comp,
                AB_comp1=AB_comp1,
                AB_comp2=AB_comp2,
                median1=median(val1),
                median2=median(val2),
                n1=length(val1),
                n2=length(val2),
                pvalue=test$p.value
        )
}

### all ABswitch seperate
p <- ABswitch2gene_exp %>% 
mutate(ABswitch=gsub("-","->",ABswitch)) %>% 
mutate(ABswitch=factor(ABswitch, levels=c("A->B","B->A","B->B","A->A"))) %>% 
mutate(comp=gsub("_vs_","->",comp)) %>% 
ggplot(aes(x=ABswitch, y=log2FC)) +
facet_wrap(~comp) +
stat_boxplot(geom = "errorbar",outlier.shape=NA, width = 0.2) + 
geom_boxplot(aes(fill=ABswitch),outlier.shape=NA, notch = TRUE) +
theme_bw() +
ylim(c(-5,5))
ggsave(file.path(dir,"plots","ABswitch_geneExpLog2FC.scRNA.pdf"), p, height=5)

comp_grid = expand_grid(
        AB_comp1=ABswitch2gene_exp %>% distinct(ABswitch) %>% pull(ABswitch),
        AB_comp2=ABswitch2gene_exp %>% distinct(ABswitch) %>% pull(ABswitch)
) %>% filter(AB_comp1 > AB_comp2)

out = do.call(
        "bind_rows",
        lapply(
                ABswitch2gene_exp %>% 
                distinct(comp) %>% 
                pull(comp),
                function(cell){
                        do.call(
                                "bind_rows",
                                mapply(
                                        function(A, B){
                                                cal_wilcox_p(
                                                        df = ABswitch2gene_exp, 
                                                        AB_comp1=A, AB_comp2=B, 
                                                        cell_comp = cell)
                                        },
                                        A=comp_grid$AB_comp1, 
                                        B=comp_grid$AB_comp2,
                                        SIMPLIFY=F
                                )
                                
                        )
                }
        )
)
write_csv(
        out,
        file=file.path(dir,"plots","ABswitch_geneExpLog2FC.scRNA.csv")
)

### combine AB
ggplot_build(p)$data[[2]] %>% 
as_tibble() %>% distinct(fill)
# 1 #F8766D
# 2 #7CAE00
# 3 #00BFC4
# 4 #C77CFF

p2 <- ABswitch2gene_exp %>% 
mutate(ABswitch=ifelse(ABswitch %in% c("B-B","A-A"),"Stable",ABswitch)) %>% 
mutate(ABswitch=gsub("-","->",ABswitch)) %>% 
mutate(ABswitch=factor(ABswitch, levels=c("Stable","A->B","B->A"))) %>% 
mutate(comp=gsub("_vs_","->",comp)) %>% 
ggplot(aes(x=ABswitch, y=-log2FC)) +
facet_wrap(~comp) +
stat_boxplot(geom = "errorbar",outlier.shape=NA, width = 0.2) + 
geom_boxplot(aes(fill=ABswitch),outlier.shape=NA, notch = TRUE) +
theme_bw() +
ylim(c(-5,5)) +
scale_fill_manual(values=c("grey","#F8766D","#7CAE00"))
ggsave(file.path(dir,"plots","ABswitch_geneExpLog2FC.scRNA2.pdf"), p2, height=5)

comp_grid = expand_grid(
        AB_comp1=ABswitch2gene_exp %>% 
        mutate(ABswitch=ifelse(ABswitch %in% c("B-B","A-A"),"Stable",ABswitch)) %>% 
        distinct(ABswitch) %>% 
        pull(ABswitch),
        AB_comp2=ABswitch2gene_exp %>% 
        mutate(ABswitch=ifelse(ABswitch %in% c("B-B","A-A"),"Stable",ABswitch)) %>% 
        distinct(ABswitch) %>% 
        pull(ABswitch)
) %>% filter(AB_comp1 > AB_comp2)

out = do.call(
        "bind_rows",
        lapply(
                ABswitch2gene_exp %>% 
                distinct(comp) %>% 
                pull(comp),
                function(cell){
                        do.call(
                                "bind_rows",
                                mapply(
                                        function(A, B){
                                                cal_wilcox_p(
                                                        df = ABswitch2gene_exp %>% 
                                                        mutate(ABswitch=ifelse(ABswitch %in% c("B-B","A-A"),"Stable",ABswitch)), 
                                                        AB_comp1=A, AB_comp2=B, 
                                                        cell_comp = cell)
                                        },
                                        A=comp_grid$AB_comp1, 
                                        B=comp_grid$AB_comp2,
                                        SIMPLIFY=F
                                )
                                
                        )
                }
        )
)
write_csv(
        out,
        file=file.path(dir,"plots","ABswitch_geneExpLog2FC.scRNA2.csv")
)

save(ABswitch2gene_exp, ABswitch2gene, expressedGenes, file=file.path(dir,"Rdata","ABswitch2gene.Rdata"))

##################### ABswitch functional enrichment ###############################
load(file.path(dir,"Rdata/PCA.Rdata"))
ABcompartment = PC1_df %>% 
select(-contains("_2reps")) %>% 
na.omit() %>% 
filter(chrom %in% paste0("chr",c(1:22))) %>% 
filter(Acinar_1*Acinar_2 > 0, Alpha_1*Alpha_2 > 0, Beta_1*Beta_2 > 0) %>% 
gather(key="sample",value="PC1",-chrom,-start, -end) %>% 
mutate(condition=sapply(sample, function(x){strsplit(x,"_")[[1]][1]})) %>% 
group_by(chrom, start, end, condition) %>% 
summarise(PC1=mean(PC1)) %>% 
ungroup() %>% 
mutate(compartment=ifelse(PC1>0,"A","B"))

load("/mnt/isilon/sfgi/suc1/analyses/grant/scRNA/pancreaticCells/Rdata/avgExp.Rdata")

tss = read_delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.nonpseudo.TSS.bed", delim="\t", comment="#", col_names=c("chr","start","end","tss_anno", "dot","strand")) %>% select(-dot)

gene_type=read_delim(
        "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt", delim="\t", comment="#"
) %>% 
filter(gene_type!="pseudogene")

gene_tss = tss %>% 
left_join(
        gene_type %>% 
        mutate(tss_anno=glue::glue("{gene_name}+{transcript_id}"))
) %>% select(-tss_anno, -transcript_id, -biotype) %>% 
distinct() %>% 
mutate(gene_id=sapply(gene_id, function(x){strsplit(x,"\\.")[[1]][1]}))

ABcompartment2gene = overlap_df(
        ABcompartment,
        gene_tss,
        df1_chr="chrom", df2_chr="chr",
        df1_0base=F,df2_0base=F, minoverlap=1L
)

ABcompartment2gene = bind_cols(
        ABcompartment2gene$overlap_df1,
        ABcompartment2gene$overlap_df2 %>% 
        select(gene_id, gene_name)
) %>% distinct()

# filter unbiougous genes
unbiougousGenes=ABcompartment2gene %>% 
group_by(gene_id, gene_name, condition) %>% 
summarise(compart_n=n_distinct(compartment)) %>% 
ungroup() %>% 
filter(compart_n>1) %>% 
distinct(gene_name) # unbiougous genes (genes assigned to two different compartments in any given cell type)

ABcompartment2gene = anti_join(
        ABcompartment2gene,
        unbiougousGenes
) %>% arrange(gene_name, gene_id) %>% 
distinct(gene_id, gene_name, condition, compartment)

# filter genes with more than 1 gene_id
ABcompartment2gene = ABcompartment2gene %>% 
anti_join(
        ABcompartment2gene %>% 
        distinct(gene_name, gene_id) %>% 
        group_by(gene_name) %>% 
        summarise(gene_id_n=n_distinct(gene_id)) %>% 
        filter(gene_id_n > 1) # remove genes with double gene_id 
) %>% spread(key=condition, value=compartment)

# filter genes with no expression in given cells
load(file.path(dir,"Rdata","ABswitch2gene.Rdata")) # ABswitch2gene_exp, ABswitch2gene, expressedGenes
ABcompartment2gene = semi_join(ABcompartment2gene, expressedGenes)
save(ABcompartment2gene, file=file.path(dir,"Rdata","ABcompartment2gene.Rdata"))

ABcompartment2gene %>% 
filter(gene_name %in% c("INS", "MAFA", "GCG", "MAFB", "ARX","CPA1", "PRSS1"))

compartment_cellspecific_genes=bind_rows(
        bind_rows(
                ABcompartment2gene %>% 
                filter(Acinar=="A", Alpha=="B", Beta=="B"),
                ABcompartment2gene %>% 
                filter(Acinar=="B", Alpha=="A", Beta=="A")
        ) %>% select(gene_name, compartment=Acinar) %>% 
        mutate(cell="Acinar"),
        bind_rows(
                ABcompartment2gene %>% 
                filter(Acinar=="B", Alpha=="A", Beta=="B"),
                ABcompartment2gene %>% 
                filter(Acinar=="A", Alpha=="B", Beta=="A")
        ) %>% select(gene_name, compartment=Alpha) %>% 
        mutate(cell="Alpha"),
        bind_rows(
                ABcompartment2gene %>% 
                filter(Acinar=="B", Alpha=="B", Beta=="A"),
                ABcompartment2gene %>% 
                filter(Acinar=="A", Alpha=="A", Beta=="B")
        ) %>% select(gene_name, compartment=Beta) %>% 
        mutate(cell="Beta")
)
compartment_cellspecific_genes %>% 
group_by(cell) %>% 
count()
# cell       n
# <chr>  <int>
# 1 Acinar  1761
# 2 Alpha   2066
# 3 Beta     980
save(ABcompartment2gene, compartment_cellspecific_genes, file=file.path(dir,"Rdata","ABcompartment2gene.Rdata"))

# write_csv(
#         ABcompartment2gene %>% 
#         filter(Acinar=="A", Alpha=="B", Beta=="B"),
#         file=file.path(dir,"plots","tmp_Acinar_ABgenes.csv")
# )
# 
# write_csv(
#         ABcompartment2gene %>% 
#         filter(Acinar=="B", Alpha=="A", Beta=="B"),
#         file=file.path(dir,"plots","tmp_Alpha_ABgenes.csv")
# )
# 
# write_csv(
#         ABcompartment2gene %>% 
#         filter(Acinar=="B", Alpha=="B", Beta=="A"),
#         file=file.path(dir,"plots","tmp_Beta_ABgenes.csv")
# )



### phyper
phyper_enrichment <- function(query,pathway,N){
        q = length(intersect(query, pathway))
        m = length(pathway)
        n = N-m
        k = length(query)
        p_value=phyper(q-1, m, n, k, lower.tail = FALSE)
        tibble(
                pathway_gene_n=m, 
                expected=m/(m+n)*k, 
                expected_ratio=m/(m+n), 
                observed=q, 
                observed_ratio=q/m, 
                p_value=p_value,
                genes=paste(intersect(query, pathway), collapse=",")
        )
}

## phyper reactome
load("/mnt/isilon/sfgi/suc1/msigdb/c2.cp.reactome.v7.1.symbols.Rdata")
pathway=db[[1]]
query = compartment_cellspecific_genes %>% 
filter(cell=="Beta") %>% pull(gene_name)
query = intersect(query, unique(unlist(db)))


# using all genes in db as universe
N=length(unique(unlist(db)))
enrich_results=lapply(db, phyper_enrichment, query=query, N=N)
enrich_results= do.call("bind_rows", enrich_results) %>% 
mutate(pathway=names(db)) %>% 
select(pathway, pathway_gene_n, expected, expected_ratio, observed, observed_ratio, p_value, genes)
enrich_results = enrich_results %>% 
mutate(fdr=p.adjust(enrich_results$p_value, method="fdr"))

enrich_results %>% 
select(pathway, p_value, fdr) %>% 
arrange(fdr)

enrich_results = enrich_results %>% 
filter(fdr < 0.05) %>% 
select(pathway, fdr, pathway_gene_n, observed, genes) %>% 
arrange(fdr) %>% 
mutate(pathway = gsub("REACTOME_","",pathway)) %>% 
mutate(pathway=tolower(pathway))

## phyper GO
load("/mnt/isilon/sfgi/suc1/msigdb/c5.all.v7.1.symbols.Rdata")
pathway=db[[1]]
query = compartment_cellspecific_genes %>% 
filter(cell=="Acinar", compartment=="A") %>% pull(gene_name)
query = intersect(query, unique(unlist(db)))

# using all genes in db as universe
N=length(unique(unlist(db)))
enrich_results=lapply(db, phyper_enrichment, query=query, N=N)
enrich_results= do.call("bind_rows", enrich_results) %>% 
mutate(pathway=names(db)) %>% 
select(pathway, pathway_gene_n, expected, expected_ratio, observed, observed_ratio, p_value, genes)
enrich_results = enrich_results %>% 
mutate(fdr=p.adjust(enrich_results$p_value, method="fdr"))
enrich_results %>% 
select(pathway, p_value, fdr) %>% 
arrange(fdr)

enrich_results = enrich_results %>% 
filter(fdr < 0.05) %>% 
select(pathway, fdr, pathway_gene_n, observed, genes) %>% 
arrange(fdr) %>% 
mutate(pathway = gsub("GO_","",pathway)) %>% 
mutate(pathway=tolower(pathway))

write_csv(enrich_results, file=file.path(dir,"plots","GOEnrichment_Acinar_AcompartmentSpecific.csv"))

p <- enrich_results %>% 
mutate(ratio=observed/pathway_gene_n) %>% 
mutate(log10FDR=-log10(fdr)) %>% 
head(10) %>% 
ggplot(aes(reorder(pathway,log10FDR),log10FDR)) +
geom_col(aes(alpha=ratio)) +
theme_bw() +
xlab("") + ylab("-log10FDR") +
coord_flip()
ggsave(file.path(dir,"plots","GOEnrichment_Acinar_AcompartmentSpecific.pdf"), p, height=3)

## endocrine
pathway=db[[1]]
query = compartment_cellspecific_genes %>% 
filter(cell=="Acinar", compartment=="B") %>% pull(gene_name)
query = intersect(query, unique(unlist(db)))

# using all genes in db as universe
N=length(unique(unlist(db)))
enrich_results=lapply(db, phyper_enrichment, query=query, N=N)
enrich_results= do.call("bind_rows", enrich_results) %>% 
mutate(pathway=names(db)) %>% 
select(pathway, pathway_gene_n, expected, expected_ratio, observed, observed_ratio, p_value, genes)
enrich_results = enrich_results %>% 
mutate(fdr=p.adjust(enrich_results$p_value, method="fdr"))
enrich_results %>% 
select(pathway, p_value, fdr) %>% 
arrange(fdr)

enrich_results = enrich_results %>% 
filter(fdr < 0.1) %>% 
select(pathway, fdr, pathway_gene_n, observed, genes) %>% 
arrange(fdr) %>% 
mutate(pathway = gsub("GO_","",pathway)) %>% 
mutate(pathway=tolower(pathway))

write_csv(enrich_results, file=file.path(dir,"plots","GOEnrichment_Acinar_BcompartmentSpecific.csv"))

p <- enrich_results %>% 
mutate(ratio=observed/pathway_gene_n) %>% 
mutate(log10FDR=-log10(fdr)) %>% 
head(10) %>% 
ggplot(aes(reorder(pathway,log10FDR),log10FDR)) +
geom_col(aes(alpha=ratio)) +
theme_bw() +
xlab("") + ylab("-log10FDR") +
coord_flip()
ggsave(file.path(dir,"plots","GOEnrichment_Acinar_BcompartmentSpecific.pdf"), p, height=2)

################### ABswitch and peak accessibility (single cell) #################
load(file.path(dir,"Rdata","ABswitch.Rdata"))
load("/mnt/isilon/sfgi/suc1/analyses/grant/scRNA/pancreaticCells/Rdata/avgExp.Rdata")

tss = read_delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.nonpseudo.TSS.bed", delim="\t", comment="#", col_names=c("chr","start","end","tss_anno", "dot","strand")) %>% select(-dot)

gene_type=read_delim(
        "/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt", delim="\t", comment="#"
) %>% 
filter(gene_type!="pseudogene")

gene_tss = tss %>% 
left_join(
        gene_type %>% 
        mutate(tss_anno=glue::glue("{gene_name}+{transcript_id}"))
) %>% select(-tss_anno, -transcript_id, -biotype) %>% 
distinct() %>% 
mutate(gene_id=sapply(gene_id, function(x){strsplit(x,"\\.")[[1]][1]}))

