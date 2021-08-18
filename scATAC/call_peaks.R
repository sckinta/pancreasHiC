library(SnapATAC)
library(tidyverse)
library(parallel)

dir="/mnt/isilon/sfgi/suc1/analyses/grant/scATAC/pancreaticCells"

load(file.path(dir,"Rdata","spga_filtered2_clusterCelltype.Rdata")) # sp.ga

cellType = levels(sp.ga@cluster)

run_macs2 <- function(cell){
        print(paste("Run MACS2 on", cell));
        dir.create(file.path(dir,"peaks2",cell), showWarnings=F)
        peaks = runMACS(
                obj=sp.ga[which(sp.ga@cluster==cell),], 
                output.prefix=file.path(dir,"peaks2",cell, cell),
                path.to.snaptools="/mnt/isilon/sfgi/programs/miniconda3/envs/py37/bin/snaptools",
                path.to.macs="/mnt/isilon/sfgi/programs/miniconda3/envs/py27/bin/macs2",
                gsize="hs", # mm, hs, etc
                buffer.size=500, 
                num.cores=1,
                macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
                tmp.folder=file.path(dir,"peaks2",cell)
        );
	       peaks
}


peaks.ls <- mclapply(cellType, run_macs2, mc.cores=12)

save(peaks.ls, cellType, file=file.path(dir,"Rdata","peaks.Rdata"))
