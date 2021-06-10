dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/Acinar_2reps"
condition="Acinar_2reps"
cd $dir
mkdir -p $dir/scripts
# link cool files
mkdir -p $dir/cool
cd $dir/cool
files=($(ls /mnt/isilon/sfgi/suc1/analyses/grant/hiC/hicup/Acinar_1/matrix/*.cool))
for file in ${files[@]};do
        ln -s $file
done
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/hiC/hicup/Acinar_1/matrix/Acinar_1.8000.cool

files=($(ls /mnt/isilon/sfgi/suc1/analyses/grant/hiC/hicup/Acinar_2/matrix/*.cool))
for file in ${files[@]};do
        ln -s $file
done
ln -s /mnt/isilon/sfgi/suc1/analyses/grant/hiC/hicup/Acinar_2/matrix/Acinar_2.8000.cool


# merge cool
mkdir -p $dir/cool/merge
condition="Acinar_2reps"
resolutions=("2.5M" "1M" "500K" "250K" "100K" "50K" "40K" "25K" "10K" "5K" "4000" "2500" "2000" "1500" "1K" "500" "8000")
for resol in ${resolutions[@]}; do
        files=($(ls $dir/cool/*.$resol.cool))
        echo "conda activate py37
cooler merge $dir/cool/merge/$condition.$resol.cool ${files[@]}
cooler balance $dir/cool/merge/$condition.$resol.cool
" > $dir/scripts/cooler_merge.$resol.sh
        sed -i '1i#!/bin/bash' $dir/scripts/cooler_merge.$resol.sh
        cd $dir/scripts
        # sbatch --mem 8G cooler_merge.$resol.sh
done

cd $dir/cool/merge/
ln -s Acinar_2reps.1K.cool Acinar_2reps.1000.cool

# A/B compartments - HiCexplore hicPCA
mkdir -p $dir/hicPCA
cd $dir/hicPCA
sed 1d /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins_gc.40K.bed | awk 'BEGIN{OFS="\t"; FS="\t"}{print $1,$2,$3,NR,".",$4}' > /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins_gc.40K.6col.bed

files=($(ls $dir/cool/*.40K.cool) $(ls $dir/cool/merge/*.40K.cool))
for file in ${files[@]}; do
        prefix=$(basename $file | cut -d "." -f 1)
        echo "conda activate py37
hicPCA --matrix $file \
--format bedgraph \
--outputFileName $dir/hicPCA/$prefix.PC1.bedGraph $dir/hicPCA/$prefix.PC2.bedGraph \
--numberOfEigenvectors 2 \
--extraTrack /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins_gc.40K.6col.bed \
--pearsonMatrix $dir/hicPCA/$prefix.pearson.cool \
--obsexpMatrix $dir/hicPCA/$prefix.obsexp.cool
" > $prefix.hicPCA.sh
# qsub -cwd -l h_vmem=32G Acinar_2reps.hicPCA.sh
done


# A/B compartments - cooltools call-compartments
mkdir -p $dir/ABcompartments
cd $dir/ABcompartments

# https://github.com/open2c/cooltools/issues/116

files=($(ls $dir/cool/*.40K.cool) $(ls $dir/cool/merge/*.40K.cool))
for file in ${files[@]}; do
        prefix=$(basename $file | cut -d "." -f 1)
        echo "conda activate py37
cd $dir/ABcompartments
cooltools call-compartments --bigwig --reference-track /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins_gencodeV19Count.40K.bed -o $dir/ABcompartments/$prefix.eigen.40K $file 
" > $dir/scripts/eigen.$prefix.sh
        cd $dir/scripts
        qsub -cwd eigen.$prefix.sh
        cd $dir/ABcompartments
done

# create bigwig and bedGraph
files=($(ls $dir/ABcompartments/*.eigen.40K.cis.vecs.tsv))
for file in ${files[@]}; do
        prefix=$(basename $file | sed "s/\.eigen\.40K\.cis\.vecs\.tsv//")
        sed 1d $file | awk 'BEGIN{OFS="\t"; FS="\t"}{if(!$7){$7=0} print $1,$2,$3,$7}' > $dir/ABcompartments/$prefix.bedGraph
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($4<0){$4=0}; print $1,$2,$3,$4}' $dir/ABcompartments/$prefix.bedGraph > $dir/ABcompartments/$prefix.A.bedGraph
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($4>0){$4=0}else{$4=-$4}; print $1,$2,$3,$4}' $dir/ABcompartments/$prefix.bedGraph > $dir/ABcompartments/$prefix.B.bedGraph
        bedGraphToBigWig $dir/ABcompartments/$prefix.A.bedGraph /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $dir/ABcompartments/$prefix.A.bw
        bedGraphToBigWig $dir/ABcompartments/$prefix.B.bedGraph /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $dir/ABcompartments/$prefix.B.bw
        rm $prefix*bedGraph
done

# plot A/B compartment and create comp.bed
files=($(ls $dir/ABcompartments/*.eigen.40K.cis.vecs.tsv))
for file in ${files[@]}; do
        prefix=$(basename $file | sed "s/\.eigen\.40K\.cis\.vecs\.tsv//")
        echo "module load R/3.3.2
Rscript /home/suc1/hic/scripts/R/plot_ABcomp.R 1 $file
Rscript /home/suc1/hic/scripts/R/plot_ABcomp.R 2 $file
" > $dir/scripts/$prefix.plot_ABcomp.sh
        cd $dir/scripts/
        qsub -cwd -l h_vmem=4G $prefix.plot_ABcomp.sh
done

# TAD boundaries - cooltools diamond-insulation
mkdir -p $dir/insulation_TADs
### The window size ({}) has to be a multiple of the bin size {}
files=($(ls $dir/cool/*.10K.cool))
for file in ${files[@]}; do
        prefix=$(basename $file | cut -d "." -f 1)
        echo "conda activate py37
cd $dir/insulation_TADs
cooltools diamond-insulation $file 500000 > $dir/insulation_TADs/$prefix.bin10K.window500K.tsv
" > $dir/scripts/insulation.$prefix.sh
        cd $dir/scripts
        # qsub -cwd insulation.$prefix.sh
        cd $dir/insulation_TADs
done

echo "res:10000" > $dir/hitad_TADs/meta_10k.txt
files=($(ls $dir/cool/*.10K.cool))
for i in `seq 1 ${#files[@]}`; do
        j=$((i-1))
        echo " rep$i:${files[$j]}"
done >>  $dir/hitad_TADs/meta_10k.txt

echo "conda activate py37
cd $dir/hitad_TADs
hitad -p 4 -O $dir/hitad_TADs/hiTADs_10K.bed -d $dir/hitad_TADs/meta_10k.txt --logFile $dir/hitad_TADs/hitad_10k.log
" > $dir/scripts/hitad_10k.sh

files=($(ls $dir/cool/*.10K.cool))
for file in ${files[@]}; do
        prefix=$(basename $file | cut -d "." -f 1)
echo "output-DI -O $dir/hitad_TADs/$prefix.10K.DI.bedGraph -p $file"
done > $dir/scripts/hitad_10k.sh

cd $dir/scripts
qsub -cwd -pe smp 4 hitad_10k.sh
cd $dir/hitad_TADs

# # plot
# files=($(ls $dir/cool/*.40K.cool))
# for file in ${files[@]}; do
#         prefix=$(basename $file | cut -d "." -f 1)
# echo "conda activate py37
# cd $dir/hitad_TADs
# tad-plot -O $dir/hitad_TADs/$prefix.chr22.40K.png -p $file -T $dir/hitad_TADs/hiTADs_40K.bed -C chr22
# " > $dir/scripts/tad_plot.40k.$prefix.sh
# done

## create DI tracks
files=($(ls $dir/hitad_TADs/*.10K.DI.bedGraph))
for file in ${files[@]}; do
        prefix=$(basename $file | sed "s/\.10K\.DI\.bedGraph//")
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($4>0){n=0}else{n=-$4}; print $1,$2,$3,n}' $file > $dir/hitad_TADs/$prefix.10K.negDI.bedGraph
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($4>0){n=$4}else{n=0}; print $1,$2,$3,n}'  $file > $dir/hitad_TADs/$prefix.10K.posDI.bedGraph
        bedGraphToBigWig $dir/hitad_TADs/$prefix.10K.DI.bedGraph /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $prefix.10K.DI.bw
        bedGraphToBigWig $dir/hitad_TADs/$prefix.10K.negDI.bedGraph /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $dir/hitad_TADs/$prefix.10K.negDI.bw
        bedGraphToBigWig $dir/hitad_TADs/$prefix.10K.posDI.bedGraph /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $dir/hitad_TADs/$prefix.10K.posDI.bw
        rm $dir/hitad_TADs/$prefix.10K.negDI.bedGraph $dir/hitad_TADs/$prefix.10K.posDI.bedGraph
done

#################### interaction loops (mustache) ################################
mkdir -p $dir/mustache_ICE_loops
cd $dir/mustache_ICE_loops
condition="Acinar_2reps"
resolutions=("500" "1K" "1500" "2000" "2500" "4000" "5K")
for resolution in ${resolutions[@]}; do
        res=$(echo $resolution | awk '{print tolower($0)}')
        if [[ $res =~ "k" || $res =~ "m" ]]; then
                res=$res"b"
        fi
        # echo "$res"
        mkdir -p $dir/mustache_ICE_loops/$res
        # cd $dir/mustache_ICE_loops/$res
        chrs=($(seq 1 22))
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                echo "conda activate mustache
mustache -f $dir/cool/merge/$condition.$resolution.cool -ch $chr -r $res -pt 0.1 -o $dir/mustache_ICE_loops/$res/loops.$chr.bedpe
" > $dir/scripts/mustache.$res.$chr.sh
                cd $dir/scripts
                qsub -cwd -l h_vmem=16G mustache.$res.$chr.sh
        done
done


cd $dir/mustache_ICE_loops
resolutions=($(ls -d */))
for resol in ${resolutions[@]}; do
        resol=$(echo $resol | sed "s/\///")
        files=($(ls $dir/mustache_ICE_loops/$resol/*.bedpe))
        for file in ${files[@]}; do
                sed 1d $file
        done > mustache_ICE_loops.$resol.bedpe
        # rm -rf $resol/
done

# APA to verify loops
cd $dir/mustache_ICE_loops
echo "source ~/.bashrc
conda activate py36
cut -f 1-6 $dir/mustache_ICE_loops/mustache_ICE_loops.4000.bedpe > tmp.bedpe
coolpup.py $dir/cool/merge/$condition.4000.cool tmp.bedpe --nshifts 10 --pad 40 --mindist 100000 --outname $dir/mustache_ICE_loops/$condition\"_loop.4000.txt\"
plotpup.py $dir/mustache_ICE_loops/$condition\"_loop.4000.txt\" --row_names $condition --output $condition\"_loopsAPA.4000.pdf\"
rm tmp.bedpe
" > $dir/scripts/apa_loops.sh
sed -i '1i#!/bin/bash' $dir/scripts/apa_loops.sh
cd $dir/scripts
sbatch --mem 16G $dir/scripts/apa_loops.sh



# # summarize loops
# tss_file="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.nonpseudo.TSS.bed"
# gene_type_file="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt"
# ocr_file="/mnt/isilon/sfgi/pahlm/datasets_downloaded/GEO/GSE76268_Pancreas_Cells_ATAC/Acinar_extended_atac_peaks_hg19.bed"
# resols=("500" "1kb" "2000" "4000")
# bedpe_list=$(for resol in ${resols[@]}; do
#         echo "$dir/mustache_ICE_loops/mustache_ICE_loops.$resol.bedpe"
# done | tr "\n" "," | sed "s/,$//")
# 
# echo "module load R/3.3.2
# Rscript /home/suc1/hic/scripts/R/summarise_bedpe.R $tss_file $gene_type_file $ocr_file $bedpe_list
# " > $dir/scripts/summarise_bedpe.sh
# cd $dir/scripts/
# qsub -cwd -l h_vmem=8G $dir/scripts/summarise_bedpe.sh
# cd $dir/mustache_ICE_loops
# 
# mv bedpe_summary.csv $condition.bedpe_summary.csv

# ### merge loops for visualization
# # all loops
# resols=("500" "1kb" "2000" "4000")
# bedpe_list=$(for resol in ${resols[@]}; do
#         echo "$dir/mustache_ICE_loops/mustache_ICE_loops.$resol.bedpe"
# done | tr "\n" "," | sed "s/,$//")
# Rscript /home/suc1/hic/scripts/R/mergeBedpeResForViz.R --genome "hg19" -i $bedpe_list -o $dir/mustache_ICE_loops/$condition.mustache_loops.bedpe
# 
# perl ~/hic/scripts/perl/bedpe2ucscBigInteract.pl $dir/mustache_ICE_loops/$condition.mustache_loops.bedpe "no" "#7A67EE" | sort -k1,1 -k2,2n > tmp.bed
# 
# bedToBigBed -as=/mnt/isilon/sfgi/referenceSequences/interact.as -type=bed5+13 tmp.bed /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $dir/mustache_ICE_loops/$condition.mustache_loops.bb
# 
# rm tmp.bed
# 
# # gene_loops
# resols=("500" "1kb" "2000" "4000")
# bedpe_list=$(for resol in ${resols[@]}; do
#         echo "$dir/mustache_ICE_loops/mustache_ICE_loops.$resol.bedpe"
# done | tr "\n" "," | sed "s/,$//")
# Rscript /home/suc1/hic/scripts/R/mergeBedpeResForViz.R --genome "hg19" --tss T -i $bedpe_list -o $dir/mustache_ICE_loops/$condition.mustache_loops_anno.bedpe
# 
# perl ~/hic/scripts/perl/bedpe2ucscBigInteract.pl $dir/mustache_ICE_loops/$condition.mustache_loops_anno.bedpe "no" "0,64,255" | sort -k1,1 -k2,2n > tmp.bed
# 
# bedToBigBed -as=/mnt/isilon/sfgi/referenceSequences/interact.as -type=bed5+13 tmp.bed /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $dir/mustache_ICE_loops/$condition.mustache_loops_anno.bb
# 
# rm tmp.bed
# 
# # upload to webserver
# webserver_dir="/var/www/html/ucsc/sfgi/suc1/grant/hiC/pancreatic_cells_2reps_hic_loops"  # change
# scp -o PubkeyAuthentication=no -vvv *.bb suc1@reslncsfgweb01.research.chop.edu:$webserver_dir
# 


################################## fithic_loops #################################
mkdir -p $dir/fithic_loops
cd $dir/fithic_loops

resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        mkdir -p $dir/fithic_loops/$resol
        cool_file=$dir/cool/merge/$condition.$resol.cool
        chrs=(`seq 1 22` "X")
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                mkdir -p $dir/fithic_loops/$resol/$chr
                cd $dir/fithic_loops/$resol/$chr
                sbatch  -t 12:00:00 --mem 32G -J fithic ~/hic/scripts/bash/run_fithic.sh $chr $resol $cool_file $dir/fithic_loops/$resol/$chr
        done
        cd $dir/fithic_loops
done

resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        ls $dir/fithic_loops/$resol/*/FitHiC.spline_pass2.res$resol.significances.txt.gz | wc -l
        # 23
done

# merge all chr results, recal fdr and filter fdr < 1e-4
resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        # merge chr results
        echo "dir=$dir
resol=$resol
        " > $dir/fithic_loops/$resol/fithic_postprocess.sh
        echo 'chrs=($(seq 1 22) "X")
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                if [ $chr = "chr1" ]; then
                        zcat $dir/fithic_loops/$resol/$chr/FitHiC.spline_pass2.res$resol.significances.txt.gz
                else
                        zcat $dir/fithic_loops/$resol/$chr/FitHiC.spline_pass2.res$resol.significances.txt.gz | sed 1d
                fi
        done | gzip > $dir/fithic_loops/$resol/FitHiC.spline_pass2.res$resol.significances.txt.gz
        ' >> $dir/fithic_loops/$resol/fithic_postprocess.sh
        
        # recal fdr
        infile="$dir/fithic_loops/$resol/FitHiC.spline_pass2.res$resol.significances.txt.gz"
        outdir="$dir/fithic_loops/$resol"
        echo "module load R/4.0.2
Rscript ~/hic/scripts/R/fithic_padjust.R $infile
        " >> $dir/fithic_loops/$resol/fithic_postprocess.sh
        
        # filter fdr < 1e-4
        filter_col="fdr"
        cutoff="1e-4"
        infile="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/Acinar_2reps/fithic_loops/$resol/FitHiC.spline_pass2.res$resol.significances.fdr.txt.gz"
        outfile="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/Acinar_2reps/fithic_loops/$resol/FitHiC.spline_pass2.res$resol.sigInt.bedpe"
        echo "module load R/4.0.2
Rscript ~/hic/scripts/R/filter_fithic.R $filter_col $cutoff $resol $infile $outfile
        " >> $dir/fithic_loops/$resol/fithic_postprocess.sh

        sed -i '1i#!/bin/bash' $dir/fithic_loops/$resol/fithic_postprocess.sh
        cd $dir/fithic_loops/$resol
        sbatch -t 24:00:00 --mem 120G fithic_postprocess.sh
        # sbatch --mem 32G fithic_postprocess.sh
        # cd $dir/fithic_loops
done


### run hicACT (initial thred p < 1.8e-7) and filter ACT_pvalue < 1e-10 and write to bedpe
resols=("1000" "2000" "4000")
chrs=($(seq 1 22) "X")
for resol in ${resols[@]}; do
        for chr in ${chrs[@]}; do
                # run hicACT
                chr="chr"$chr
                resolKB=$(echo $resol | sed "s/000//")
                infile="$dir/fithic_loops/$resol/$chr/FitHiC.spline_pass2.res$resol.significances.txt.gz"
                outdir="$dir/fithic_loops/$resol/$chr"
                echo "module load R/4.0.2
Rscript ~/hic/scripts/R/hicACT.R $resolKB $infile $outdir
        " > $dir/fithic_loops/$resol/$chr/hicACT.sh
        
                # filter_fithic by ACT_adjusted < 1e-10
                filter_col="ACT_pvalue"
                cutoff="1e-10"
                infile="$dir/fithic_loops/$resol/$chr/ACT_adjusted.txt.gz"
                outfile="$dir/fithic_loops/$resol/$chr/ACT_adjusted.sigInt.bedpe"
                echo "Rscript ~/hic/scripts/R/filter_fithic.R $filter_col $cutoff $resol $infile $outfile
" >> $dir/fithic_loops/$resol/$chr/hicACT.sh
        sed -i '1i#!/bin/bash' $dir/fithic_loops/$resol/$chr/hicACT.sh
        cd $dir/fithic_loops/$resol/$chr/
        sbatch --mem 16G hicACT.sh
        cd $dir/fithic_loops
        done
done

resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        # ls $dir/fithic_loops/$resol/*/ACT_adjusted.txt.gz | wc -l
        ls $dir/fithic_loops/$resol/*/ACT_adjusted.sigInt.bedpe | wc -l
# 23
done

# module load R/4.0.2
# resols=("1000" "2000" "4000")
# chrs=($(seq 1 22) "X")
# for resol in ${resols[@]}; do
#         for chr in ${chrs[@]}; do
#                 chr="chr"$chr
#                 filter_col="ACT_pvalue"
#                 cutoff="1e-10"
#                 infile="$dir/fithic_loops/$resol/$chr/ACT_adjusted.txt.gz"
#                 outfile="$dir/fithic_loops/$resol/$chr/ACT_adjusted.sigInt.bedpe"
#                 cd $dir/fithic_loops/$resol/$chr/
#                 Rscript ~/hic/scripts/R/filter_fithic.R $filter_col $cutoff $resol $infile $outfile
#                 cd $dir/fithic_loops
#         done
# done

resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        cat $dir/fithic_loops/$resol/*/ACT_adjusted.sigInt.bedpe > $dir/fithic_loops/$resol/ACT_adjusted.res$resol.sigInt.bedpe
done

files=(ls */*.bedpe)
for file in ${files[@]}; do
        wc -l $file
done

# 763 1000/ACT_adjusted.res1000.sigInt.bedpe
# 551 1000/FitHiC.spline_pass2.res1000.sigInt.bedpe
# 25187 2000/ACT_adjusted.res2000.sigInt.bedpe
# 27216 2000/FitHiC.spline_pass2.res2000.sigInt.bedpe
# 227305 4000/ACT_adjusted.res4000.sigInt.bedpe
# 288181 4000/FitHiC.spline_pass2.res4000.sigInt.bedpe

############################# hic_loops ##############################
mkdir -p $dir/hic_loops
mkdir -p $dir/hic_loops/all
mkdir -p $dir/hic_loops/anno

### link loops
cd $dir/hic_loops/all
files=($(ls $dir/mustache_ICE_loops/*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 2)
        if [[ $resol == "1kb" ]]; then
                #echo "1000"
                ln -s $file $condition.mustache.pt01.res1000.bedpe
        elif [[ $resol == "2000" ]]; then
                #echo "2000"
                ln -s $file $condition.mustache.pt01.res2000.bedpe
        elif [[ $resol == "4000" ]]; then
                #echo "4000"
                ln -s $file $condition.mustache.pt01.res4000.bedpe
        fi
done

cd $dir/hic_loops/all
files=($(ls $dir/fithic_loops/*/ACT_adjusted.*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 2)
        ln -s $file $condition.fithic_ACT.pv1e10.$resol.bedpe
done

cd $dir/hic_loops/all
files=($(ls $dir/fithic_loops/*/FitHiC.*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 3)
        ln -s $file $condition.fithic_FDR.fdr1e4.$resol.bedpe 
done

### different cutoff for fithic_ACT and fithic_FDR
cd $dir/hic_loops/all
files=($(ls $dir/hic_loops/all/$condition.fithic_ACT.pv1e10.*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 4)
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($8 <= 1e-14) print $0}' $file > $dir/hic_loops/all/$condition.fithic_ACT.pv1e14.$resol.bedpe
done

cd $dir/hic_loops/all
files=($(ls $dir/hic_loops/all/$condition.fithic_FDR.fdr1e4.*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 4)
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($8 <= 1e-6) print $0}' $file > $dir/hic_loops/all/$condition.fithic_FDR.fdr1e6.$resol.bedpe
done

### generate annotated files - anno.bedpe, anno.gene2OCR.txt anno.summary.txt
mkdir -p $dir/hic_loops/scripts
# using concensus atac peak here
atac_peak_file="/mnt/isilon/sfgi/suc1/analyses/grant/scATAC/pancreaticCells/peaks/consensus_3cells.bed"

files=($(ls $dir/hic_loops/all/*.bedpe))
for file in ${files[@]}; do
        prefix=$(basename $file | sed "s/\.bedpe/\.anno/")
        outfile_prefix="$dir/hic_loops/anno/$prefix"
        echo "module load R/4.0.2
Rscript /mnt/isilon/sfgi/suc1/scripts/R/annotate_bedpe2geneOCR_hg19.R $file $atac_peak_file $outfile_prefix
" > $dir/hic_loops/scripts/$prefix.sh
        sed -i '1i#!/bin/bash' $dir/hic_loops/scripts/$prefix.sh
        cd $dir/hic_loops/scripts/
        sbatch --mem 8G $dir/hic_loops/scripts/$prefix.sh
done

### put summary together
cd $dir/hic_loops
files=($(ls $dir/hic_loops/anno/*.anno.summary.txt))
for i in `seq 1 ${#files[@]}`; do
        j=$((i-1))
        if [[ "$i" -eq "1" ]]; then
                cat ${files[$j]}
        else
                sed 1d ${files[$j]}
        fi
done > $dir/hic_loops/$condition.hic_loops.summary.txt

cd $dir/hic_loops
midfixes=($(ls all/*.bedpe | cut -d "/" -f 2 | cut -d "." -f 2-3 | sort -u))
for midfix in ${midfixes[@]}; do
        files=($(ls $dir/hic_loops/anno/*$midfix*.anno.gene2OCR.txt))
        geneOCR_pairN=$(cat ${files[@]} | sort -u | wc -l)
        gene_N=$(cat ${files[@]} | cut -f 3 | sort -u | wc -l)
        OCR_N=$(cat ${files[@]} | cut -f 1 | sort -u | wc -l)
        echo -e "$midfix\t$geneOCR_pairN\t$gene_N\t$OCR_N"
done > $dir/hic_loops/$condition.hic_loops_anno.noRes.summary.txt
sed -i "1i#prefix\tgeneOCR_pairN\tgene_N\tOCR_N" $dir/hic_loops/$condition.hic_loops_anno.noRes.summary.txt
