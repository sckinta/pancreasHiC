# v0.7.4 split
PROMOTEROME_PARENT_DIR="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/hicup"

REP="1" # change
PROMOTEROME_DIR="Acinar_"$REP # change
FASTQ_R1="Acinar_"$REP"_R1.fastq.gz" # change
FASTQ_R2="Acinar_"$REP"_R2.fastq.gz" # change
RAWDATA_DIR="/mnt/isilon/sfgi/rawData/grant/hiC/Acinar/Acinar_"$REP  # change

############## CREATE ANALYSIS DIR FOR EACH REP ################
mkdir $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/
mkdir $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

### CREATE FASTQ SYMBOLIC LINK FOR EACH REP
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR
ln -s $RAWDATA_DIR/$FASTQ_R1
ln -s $RAWDATA_DIR/$FASTQ_R2

### PREPARE FOR HICUP
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

#### CREATE FASTQ CHUNK FOR R1 AND R2
qsub -cwd -j y -o split_1_1.log ~/captureC/scripts/bash/splitFastqGz.sh ../$FASTQ_R1 200000000 1 R1
qsub -cwd -j y -o split_1_2.log ~/captureC/scripts/bash/splitFastqGz.sh ../$FASTQ_R2 200000000 1 R2

### AFTER CREAT CHUNK CHECK WHETHER TOTAL IS RIGHT
total_split=$(ls lane_1_R1* | wc -l)
total_split_0base=$((total_split-1))
files=($(ls -lhtr lane_1_R1* | awk '{print $9}'))
last_i=$((${#files[@]} - 1))
reminder=$(wc -l ${files[$last_i]} | awk '{num=$1/4; print num}')
total=$(cat $RAWDATA_DIR/readCount.txt)
cal_total=$((total_split_0base*50000000+reminder))
if [ $total = $cal_total ]; then
        echo "$cal_total is correct"
else
        echo "splitted read number not match"
fi

### CHANGE FILE NAME
files=($(ls -lhtr lane_1_R1* | awk '{print $9}'))
for i in `seq 1 ${#files[@]}`;do 
        j=$((i-1))
        newfile="lane_1_"$i"_R1.fq"
        if [[ $newfile != $files[$j] ]]; then
                # echo "$newfile ${files[$j]}"
                mv ${files[$j]} $newfile
        fi
done

files=($(ls -lhtr lane_1_R2* | awk '{print $9}'))
for i in `seq 1 ${#files[@]}`;do 
        j=$((i-1))
        newfile="lane_1_"$i"_R2.fq"
        if [[ $newfile != $files[$j] ]]; then
                # echo "$newfile ${files[$j]}"
                mv ${files[$j]} $newfile
        fi
done

### FOR EACH LANE, CREATE CHUNK DIRECTORIES
for i in `seq 1 $total_split`;do mkdir lane_1_${i};done

### GZIP FASTQ
for i in `seq 1 $total_split`;do qsub -cwd -j y -o logs -N "gzip_1_"${i}"_R1"  ~/captureC/scripts/bash/gzip_mv.sh "lane_1_"$i"_R1.fq";done
for i in `seq 1 $total_split`;do qsub -cwd -j y -o logs -N "gzip_1_"${i}"_R2"  ~/captureC/scripts/bash/gzip_mv.sh "lane_1_"$i"_R2.fq";done


### PREPARE JOB BASH AND CONFIG FILE FOR HICUP (hg19) - can run hg19|hg38 and Arima|DpnII
# more info `perl ~/captureC/scripts/perl/writeHicupScripts.v2.pl`

cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
total_split=$(ls -d lane_1_*/ | wc -l)
perl ~/captureC/scripts/perl/writeHicupScripts.v2.pl --genome hg19 --enzyme Arima $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup $total_split 1 8G

for i in `seq 1 $total_split`;do 
        cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/lane_1_${i}
        qsub -cwd -pe smp 2 -l h_vmem=8G -o hicup.log -N "hicup_"$i hicup.sh
        cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
done



### CHECK HICUP RUNNING (IN DIFFERENT CHUNK) ARE SUCCESSFULLY COMPLETED
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
grep "HiCUP processing complete" lane*/hicup.log | wc -l
# check whether it matches split chunk number
total_split=$(ls -d lane_1_*/ | wc -l)
echo $total_split

echo "files=(\$(ls $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/*/*.hicup.bam))
perl ~/captureC/scripts/perl/hicupMergeDedup.v0.7.4.pl --threads 4 \${files[@]}
" > hicupMergeDedup.sh
qsub -cwd -l h_vmem=4G -pe smp 4 -N hicupMergeDedup -o hicupMergeDedup.log hicupMergeDedup.sh
# qsub -cwd -l h_vmem=4G -pe smp 4 -N hicupMergeDedup -o hicupMergeDedup.log hicup_dedup.sh

# CLEAN-UP (kept filt.bam)
rm -f lane_*_*/*fastq.gz
rm -f lane_*/*fq.gz
rm -f lane_*/*pair.bam
# rm -f lane_*/*filt.bam
rm -f lane_*/*hicup.bam

### Virtual CAPUTRE
# CONVERT TO CHICAGO INPUT (1 FRAGMENT RESOLUTION)
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
DIR_1FRAG="virtualDpnII_"$REP"_1frag"

bash ~/captureC/scripts/bash/run_bam2chicago.sh \
-k DpnII -g hg19 -f 1frag \
-m 32G -s 1 \
$DIR_1FRAG $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/merged.dedup.bam

# CONVERT TO CHICAGO INPUT (4 FRAGMENT RESOLUTION)
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
DIR_1FRAG="virtualDpnII_"$REP"_4frag"

bash ~/captureC/scripts/bash/run_bam2chicago.sh \
-k DpnII -g hg19 -f 4frag \
-m 32G -s 1 \
$DIR_1FRAG $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/merged.dedup.bam


### COLLECT INFO AND OUTPUT AS HICUP SUMMARY
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
total_split=$(ls -d lane_1_*/ | wc -l)
deduplicator_FILE=$(ls hicup_deduplicator_summary*)
SUMMARY_FILE=$PROMOTEROME_DIR"_hicupSummary.txt"

DIR_1FRAG="virtualDpnII_"$REP"_1frag" # optional only used for capture C
chinput_FILE=$DIR_1FRAG".chinput" # optional only used for capture C

perl ~/captureC/scripts/perl/collectHicupResults.v2.pl 1 $total_split \
$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup \
$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/$deduplicator_FILE \
$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/$DIR_1FRAG/$chinput_FILE \
> $SUMMARY_FILE

DIR_4FRAG="virtualDpnII_"$REP"_4frag" # optional only used for capture C
chinput_FILE=$DIR_4FRAG".chinput" # optional only used for capture C
TAG_NUM=$(sed '6q;d' $SUMMARY_FILE | awk '{print $24}')
perl ~/captureC/scripts/perl/captureStatsFromChinput.pl $DIR_4FRAG/$chinput_FILE $TAG_NUM > $DIR_4FRAG"_captured.txt"


bash ~/captureC/scripts/bash/summarize_hicup.v5.sh -o $PROMOTEROME_DIR"_ShortSummary.txt" $PROMOTEROME_DIR"_hicupSummary.txt"

################# convert to matrix ####################
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR
mkdir -p $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix

# bam2pairs_pre.sh (24hr)
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix
ln -s $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/merged.dedup.bam $PROMOTEROME_DIR.hicup.bam

# cp /mnt/isilon/sfgi/suc1/analyses/grant/hiC/hicup/Acinar_1/hicup/merged.dedup.bam old.hicup.bam
# BAM_IN=$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR.hicup.bam
# dedup_time=$(samtools view -H $BAM_IN | cut -f 2 | grep -c "Deduplicator")
# if [[ $dedup_time > 1 ]]; then
#         # echo "dedup"
#         samtools view -H $BAM_IN  | awk 'BEGIN{FS="\t"}$2~/Deduplicator/&&c++>0 {next} 1' | samtools reheader - $BAM_IN > tmp.bam
#         mv tmp.bam $BAM_IN
# fi

BAM_IN=$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR.hicup.bam
OUT_DIR=$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix
OUT_PREFIX=$PROMOTEROME_DIR
bash ~/captureC/scripts/bash/run_bam2pairs_pre.sh \
-g hg19 -t 8 -s 1 \
$BAM_IN \
$OUT_DIR \
$OUT_PREFIX


# pairs2cool.sh
PAIR_IN="$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR.pairs.gz"
COOL_OUT_PREFIX="$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR"
bash ~/captureC/scripts/bash/run_pairs2cool.sh \
-g "hg19" \
-t 8 -m "4G" -s 1 \
$PAIR_IN \
$COOL_OUT_PREFIX

PAIR_IN="$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR.pairs.gz"
COOL_OUT_PREFIX="$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR"
bash ~/captureC/scripts/bash/run_pairs2cool.sh \
-g "hg19" \
-t 8 -m "4G" -s 1 \
-b 1500,2000,2500,4000,40K \
$PAIR_IN \
$COOL_OUT_PREFIX

PAIR_IN="$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR.pairs.gz"
COOL_OUT_PREFIX="$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR"
bash ~/captureC/scripts/bash/run_pairs2cool.sh \
-g "hg19" \
-t 8 -m "4G" -s 1 \
-b 8000 \
$PAIR_IN \
$COOL_OUT_PREFIX

# pre2hic.sh
PRE_IN="$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR.pre.txt.gz"
HIC_OUT="$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/matrix/$PROMOTEROME_DIR.hic"
bash ~/captureC/scripts/bash/run_pre2hic.sh \
-g "hg19" -t 4 -m "120G" -s 0 \
-b 2.5M,1M,500K,250K,100K,50K,25K,10K,5K,1K,500 \
$PRE_IN $HIC_OUT

rm pair2cool.*.sh*



