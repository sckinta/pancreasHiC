library(chromVAR) # BiocManager::install("chromVAR"); devtools::install_github("sckinta/chromVAR")
library(SummarizedExperiment)
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19) # BiocManager::install("BSgenome.Hsapiens.UCSC.hg19"), v1.4.3
library(JASPAR2020) # BiocManager::install("JASPAR2020")
library(motifmatchr) # BiocManager::install("motifmatchr")

dir="/mnt/isilon/sfgi/suc1/analyses/grant/scATAC/pancreaticCells"

load(file.path(dir,"Rdata/chromVAR_fragment_counts.Rdata"))

### add GC bias
fragment_counts = addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg19)

### filter cell
fragment_counts_filtered = filterSamples(
        fragment_counts, 
        min_depth = 1500, 
        min_in_peaks = 0.15, shiny = FALSE
) # 200270 peaks x 9736 cells

### filter peak
fragment_counts_filtered = filterPeaks(
        fragment_counts_filtered, 
        non_overlapping = TRUE
) # 200268 peaks x 9736 cells

### get motif matches
jaspar_motifs = TFBSTools::getMatrixSet(
        JASPAR2020::JASPAR2020, 
        opts=list(tax_group = "vertebrates", collection = "CORE")
)

jaspar_motifs_anno = do.call(
        "bind_rows",
        lapply(
                names(jaspar_motifs),
                function(motif){
                        TF = jaspar_motifs[[motif]]@name
                        TF_uniprot_acc=jaspar_motifs[[motif]]@tags$acc
                        TF_family = jaspar_motifs[[motif]]@tags$family
                        tibble(
                                motif=motif,
                                TF=TF,
                                TF_uniprot_acc=paste(TF_uniprot_acc, collapse="|"),
                                TF_family=paste(TF_family, collapse="|")
                        )
                }
        )
)

jaspar_motifs_anno = jaspar_motifs_anno %>% mutate(TF=toupper(TF))


motif_ix = motifmatchr::matchMotifs(
        jaspar_motifs, 
        fragment_counts_filtered, 
        genome = BSgenome.Hsapiens.UCSC.hg19
) # 200268 peak x 746 motif

motif_pos = motifmatchr::matchMotifs(
        jaspar_motifs, 
        fragment_counts_filtered, 
        genome = BSgenome.Hsapiens.UCSC.hg19,
        out = "positions"
)

motif_pos_list = lapply(
        names(motif_pos),
        function(motif){
                TF = jaspar_motifs[[motif]]@name
                motif_pos[[motif]] %>% 
                as.data.frame() %>% 
                as_tibble() %>% 
                mutate(motif=motif) %>% 
                mutate(TF=TF)
        }
)

motif_pos = do.call("bind_rows", motif_pos_list)

motif_pos = motif_pos %>% mutate(TF=toupper(TF))

save(
        fragment_counts_filtered, 
        jaspar_motifs_anno, motif_ix, motif_pos,
        file=file.path(dir,"Rdata","chromVAR_withallCellPmat.Rdata")
)

### calculate dev for each motif
load(file.path(dir,"Rdata","chromVAR_withallCellPmat.Rdata"))
dev = computeDeviations(
        object = fragment_counts_filtered, 
        annotations = motif_ix
) # 746 motifs x 9736 cells

save(
        fragment_counts_filtered, 
        jaspar_motifs_anno, motif_ix, motif_pos, 
        dev, 
        file=file.path(dir,"Rdata","chromVAR_withallCellPmat.Rdata")
)

######### calculate Variability
variability = computeVariability(dev)
pdf(file.path(dir,"plots","chromVAR_motif_variability.pdf"))
plotVariability(variability, use_plotly = FALSE)
dev.off()

variability = variability %>% 
rownames_to_column("motif") %>% 
as_tibble()

save(
        fragment_counts_filtered, 
        jaspar_motifs_anno, motif_ix, motif_pos, 
        dev, 
        variability, 
        file=file.path(dir,"Rdata","chromVAR_withallCellPmat.Rdata")
)

# variability %>% 
# # filter(p_value_adj < 0.05) %>% 
# filter(grepl("SCRT",name))
