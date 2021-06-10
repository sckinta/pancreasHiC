# 4.0.2
library(magrittr)
library(multiHiCcompare)
library(BiocParallel)
library(edgeR)
library(data.table)


args = commandArgs(trailingOnly=TRUE)
dir=args[1]
cutoff=as.numeric(args[2]) # apply to FDR after re-adjust p.value

# dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/compare/naiveT_stim/data/4000"
# cutoff=0.05
data_files=dir(dir, pattern="HiCcompare\\.data\\.Rdata", recursive=T)
# data_file=data_files[1]

read_HiCcompare <- function(data_file){
        samples=dir(dirname(file.path(dir,data_file)), pattern="*.txt")
        samples=sapply(samples, function(l){strsplit(l,"\\.")[[1]][1]}, USE.NAMES=F)

        load(file.path(dir,data_file))
        colnames(normIF) = c("chr","region1","region2","D", samples)
        save(normIF, file=file.path(dir,dirname(data_file),"normIF.Rdata"))
        results
}

results = lapply(
        data_files,
        read_HiCcompare
)

comparisons = names(results[[1]])

# convert comp list
results = lapply(
        comparisons,
        function(comp){
                rbindlist(
                        lapply(
                                results,
                                function(l){
                                        l[[comp]]
                                }
                        )
                )
        }
)
names(results) = comparisons 

# recal p.adjust
results = lapply(
        results,
        function(df){
                df %>% 
                .[,p.adj:=p.adjust(df$p.value, method="fdr")]
                df
        }
)

names(results) = comparisons 

## write Rdata and summary
sigDE_summary = lapply(
        comparisons,
        function(comp){
                filename=paste0(comp,".sigDE.Rdata")
                sigDE = results[[comp]] %>% .[p.adj < cutoff]
                save(
                        sigDE,
                        file=file.path(dir,filename)
                )
                filename=paste0(comp,".results.Rdata")
                result=results[[comp]]
                save(
                        result,
                        file=file.path(dir,filename)
                )
                data.table(
                        comparison=comp,
                        totalInput=nrow(result),
                        sigDE=nrow(sigDE)
                )
        }
)

sigDE_summary = rbindlist(sigDE_summary)
write.csv(sigDE_summary, file=file.path(dir,"sigDE_summary.csv"), row.names=F)




