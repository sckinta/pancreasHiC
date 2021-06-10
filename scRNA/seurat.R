library(data.table)
library(dplyr)


donors <- c("HPAP040", "HPAP045", "HPAP054","HPAP066")
#load pancreatic markers
pancreasGenes <- unique(read.table("pancreasGenes.txt")[,1])                                               

#scRNA-seq data
StudyInfo <- fread("/project/ngsc_data/PI_INVESTIGATIONS/Kaestner_Klaus/KKMG-HPPAP/AAA-StudyInfo.txt", head=T, sep="\t")
StudyInfo <- StudyInfo[StudyInfo$ASSAY_name == "scRNA-Seq",]
StudyInfo <- StudyInfo[grepl("I", StudyInfo$RULA_Status) == FALSE, ]
StudyInfo <- StudyInfo[grepl("HPAP", StudyInfo$SAMP_name),]
#StudyInfo <- StudyInfo[StudyInfo$SAMP_name %in% donors,]
#StudyInfo <- StudyInfo[StudyInfo$SAMP_id != "83624", ]
sid2name        <- unique(StudyInfo$SAMP_sid)
sid2name        <- gsub(" ", "", sid2name)
sid2name        <- gsub("-beta", "", sid2name)
sid2name        <- gsub("_scRNASeq", "", sid2name)
sid2name        <- gsub("HPAP_", "HPAP", sid2name)
sid2name[nchar(sid2name) == 8] <- gsub("00", "0", sid2name)[which(nchar(sid2name) == 8)]
names(sid2name) <- unique(StudyInfo$SAMP_id)

samp2name <- sid2name[which(sid2name %in% donors)]
run_info  <- c("FGC2063", "FGC2063", "FGC2146")
names(run_info) <- names(samp2name)
samp2name <- paste(run_info, samp2name, sep="_")
names(samp2name) <- names(run_info)

scRNAseq_cr_dir <- "/project/kaestnerlab/gaolong/HPAP/ChromInteraction/scRNA_cellranger/"
cr_summary_df    <- c()
for(samp in names(samp2name)){
	samp_name <- samp2name[samp]
	print(samp_name)
    cr_folder  <- paste(scRNAseq_cr_dir, samp_name, sep="")
    cr_out     <- paste(cr_folder, "/outs/", sep="")
    cr_summary <- paste(cr_out, "metrics_summary.csv", sep="")
    if(dir.exists(cr_folder) == FALSE){
        next
    }
    cr_df      <- read.csv(cr_summary, colClasses = 'character')
    if(is.null(cr_summary_df)){
        cr_summary_df <- cr_df
    }
    else{
        shared_col    <- intersect(colnames(cr_summary_df), colnames(cr_df))
        cr_summary_df <- rbind(cr_summary_df[,shared_col], cr_df[,shared_col])
    }
}
cr_summary_df <- cbind(samp2name, cr_summary_df)                                                             
cr_summary_df <- select(cr_summary_df, samp2name, Estimated.Number.of.Cells)                                 
cr_summary_df$samp2name <- as.character(cr_summary_df$samp2name)                                             
cr_summary_df$Estimated.Number.of.Cells <- as.numeric(gsub(",", "", cr_summary_df$Estimated.Number.of.Cells))


library(Seurat)
for(i in 1:dim(cr_summary_df)[1]){
     sid <- row.names(cr_summary_df)[i]
     sname <- cr_summary_df[i, 1]
     print(sname)
     matrix_dir    <- paste(scRNAseq_cr_dir, sname, "/outs/filtered_feature_bc_matrix/", sep="")
     panc_temp.data <- Read10X(data.dir = matrix_dir)
     panc_temp <- CreateSeuratObject(counts = panc_temp.data, project = sname)
     panc_temp <- RenameCells(object=panc_temp, new.names=paste(sname, row.names(panc_temp@meta.data), sep="_"))
     if(i == 1){
         panc <- panc_temp
     }
     else{
         panc <- merge(panc, panc_temp, project = "scRNAseq")
     }
}

panc[["percent.mt"]] <- PercentageFeatureSet(panc, pattern = "^MT-")              

panc_st <- panc                                                                   
panc_st <- subset(panc_st, subset = nFeature_RNA > 200 & percent.mt < 10)         
rcounts <- GetAssayData(panc_st, assay = "RNA")                                   
rcounts <- rcounts[row.names(rcounts) %in% pancreasGenes, ]                       
csd <- apply(rcounts, 2, sd)                                                      
if(length(which(csd == 0)) > 0){                                                  
    rcounts <- rcounts[,-which(csd == 0)]                                         
}                                                                                 
panc_st <- CreateSeuratObject(counts = rcounts, project = "scRNAseq")
donor_info <- unlist(lapply(row.names(panc_st@meta.data), function(x) strsplit(x, "_")[[1]][2]))
panc_st@meta.data$donor <- donor_info 
#panc_st <- SCTransform(object = panc_st, verbose = TRUE)                          
panc_st <- NormalizeData(object = panc_st)                                       
panc_st <- FindVariableFeatures(object = panc_st)                                
panc_st <- ScaleData(object = panc_st)                                           
panc_st <- RunPCA(object = panc_st)                                               
panc_st <- RunUMAP(object = panc_st, dims=1:50, verbose = FALSE)                  
panc_st <- RunTSNE(object = panc_st)                                              
panc_st <- FindNeighbors(object = panc_st, dims=1:2, reduction = "umap")          
panc_st <- FindClusters(object = panc_st, reduction.type = "umap", resolution=0.5)


library(plotly)

xaxis_layout <- list(title="UMAP 1",
				  autotick = TRUE,
				  zeroline = FALSE,
				  showline = TRUE,
				  ticks = "outside",
				  tick0 = 0,
				  dtick = 0.25,
				  ticklen = 5,
				  tickwidth = 1,
				  showgrid = FALSE)
yaxis_layout <- list(title="UMAP 2",
				  autotick = TRUE,
				  zeroline = FALSE,
				  showline = TRUE,
				  ticks = "outside",
				  tick0 = 0,
				  dtick = 0.25,
				  ticklen = 5,
				  tickwidth = 1,
				  showgrid = FALSE)
																												
#umap
seurat_umap <- panc_st@reductions$umap@cell.embeddings
seurat_clt  <- panc_st@active.ident
norm_data  <- GetAssayData(panc_st)
signal <- as.numeric(norm_data["SST",])
clt_info <- as.character(panc_st@active.ident)
panc_st@active.ident <- as.factor(clt_info)

seurat_umap <- panc_st@reductions$umap@cell.embeddings
seurat_clt  <- panc_st@active.ident
norm_data  <- GetAssayData(panc_st)
signal <- as.numeric(norm_data["SST",])
seurat_umap_df <- data.frame(row.names(seurat_umap), seurat_umap, seurat_clt, signal, panc_st@meta.data$donor)
colnames(seurat_umap_df) <- c("cell", "X", "Y", "cluster", "signal", "donor")                                     
                                                                                                                    
palette <- colorRampPalette(c("white", "red"))
p <- plot_ly(seurat_umap_df, x=~X, y=~Y, color=~cluster)
p
																												
panc_all <- panc                                                           
panc_all <- subset(panc_all, subset = nFeature_RNA > 200 & percent.mt < 10)
panc_all <- NormalizeData(object = panc_all)       
panc_all <- FindVariableFeatures(object = panc_all)
panc_all <- ScaleData(object = panc_all)           

PlotFeature <- function(obj, obj_all, method="umap", gene){                               
    res <-  slot(slot(obj, "reductions")[[method]], "cell.embeddings")                    
    clt  <- slot(obj, "active.ident")                                                     
    norm_data  <- GetAssayData(obj_all)                                                   
    expression <- norm_data[gene,]                                                        
    seurat_umap_df <- data.frame(row.names(res), res, clt, expression)                    
    colnames(seurat_umap_df) <- c("cell", "X", "Y", "cluster", "expression")              
    palette <- colorRampPalette(c("lightgrey", "blue"))                                   
    p <- plot_ly(seurat_umap_df, x=~X, y=~Y, color=~expression, colors = palette(50)) %>% 
        layout(xaxis = xaxis_layout, yaxis = yaxis_layout)                                
    p                                                                                     
}                                                                                         

geneSet <- pancreasGenes
for(i in 49:length(geneSet)){
	sym <- geneSet[i]
    print(sym)
    tryCatch({p <- PlotFeature(panc_st, panc_all, gene=sym);
    orca(p, paste0("../Figures/scRNAseq_", sym, ".pdf"), height=600, width=600)}, error=function(e){})
}



#plot by cell type
alpha_clt <- c("1", "2", "3","8", "9", "14", "15", "16", "17", "20", "21")
beta_clt  <- c("0", "5", "7", "10", "13", "19", "22", "27")
delta_clt <- c("30")
acinar_clt <- c("4", "6", "11", "24")
duct_clt  <- c("12")
endo_clt  <- c("25")
mes_clt  <- c("18")
pp_clt   <- c("26")
epsilon_clt <- c("23")
clt <- c(alpha_clt, beta_clt, delta_clt, acinar_clt, duct_clt, endo_clt, mes_clt, pp_clt, epsilon_clt)
other_clt <- setdiff(as.character(1:30), clt)
clt <- c(clt, other_clt)

clt2cell <- c(rep("alpha",length(alpha_clt)), 
			  rep("beta", length(beta_clt)),
			  rep("delta", length(delta_clt)),
			  rep("acinar",length(acinar_clt)),
			  rep("ductal",length(duct_clt)), 
			  rep("endo", length(endo_clt)),
			  rep("mes",length(mes_clt)),
			  rep("pp", length(pp_clt)),
			  rep("epsilon", length(epsilon_clt)),
			  rep("other", length(other_clt)))
names(clt2cell) <- clt

seurat_umap <- panc_st@reductions$umap@cell.embeddings                                                        
seurat_clt  <- clt2cell[as.character(panc_st@active.ident)]
                                                                                                              
palette <- colorRampPalette(c("white", "red"))                                                                
p <- plot_ly(seurat_umap_df, x=~X, y=~Y, color=~cluster)                                                      
p                                                                                                             
