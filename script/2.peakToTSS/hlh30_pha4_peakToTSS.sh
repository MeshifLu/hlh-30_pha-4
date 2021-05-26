################ Common parameter
sraDownFolder=/media/meshif/sda/SRA/sra # Already finish the NGS process. Shift files to project folder.


# codeFolder
code_Folder_Path='/home/meshif/codeFolder/Project/HLH-30/20210414_reProcess_modEncode/peakToTSS'
myPyLibFolder='/home/meshif/codeFolder/myPyLib'

### bowtie2 index
# WS220
ws220Bt2Index=/media/meshif/sda/Genome/Worm/WS220/bowtie2Index/ws220_bt2Index
ws220Fasta=/media/meshif/sda/Genome/Worm/WS220/c_elegans.WS220.genomic_chrN.fa

#####################################################################################################################
### hlh30_pha4_chip
#####################################################################################################################

projectName='hlh30_pha4_chip' # used for creating project folders
# rawFastqFolder=/media/meshif/NAS/ngs_raw/PC/$projectName # Please check sraToolKits to know the path of download folder

# projectFolder=$sraDownFolder/$projectName  # Already finish the NGS process. Shift files to project folder.
projectFolder=/media/meshif/NAS/Project/HLH-30/$projectName
ngsFolder=$projectFolder/ngsData
sampleTableFile=$projectFolder/SraRunTable_pha4_hlh30.csv
chipTableFile=$projectFolder/SraRunTable_pha4_hlh30_chipTable.csv
sra_L=$(awk -F',' 'NR!=1{print $1}' $sampleTableFile)
# echo $sra_L

### conda env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ngs

### 1-3. Preprocess --> prepare folder, download SRA --> fastq-dump --> fastqc
preprocessPy=$myPyLibFolder/ngsPreprocess_fromSRA_pipe.py
python $preprocessPy $sra_L $projectName $sraDownFolder $ngsFolder $sampleTableFile

### 4. Trim
fastqFolder=$ngsFolder/2.fastq
trimFolder=$ngsFolder/4.trim_fastQC
mkdir -p $trimFolder

echo Start trim_galore trimming!
awk -F',' 'NR!=1{print $1}' $sampleTableFile | parallel -j 5 "echo {} && trim_galore -j 7 -o $trimFolder $fastqFolder/{}_R1.fastq.gz"
echo Finish trimming!

echo
echo Start fastQC
find $trimFolder -name "*.fq.gz" | parallel -j 20 fastqc {} -o $trimFolder
echo Finish FastQC!!!    


### Map reads by bowtie2
bt2Folder=$ngsFolder/5.bt2
samFolder=$bt2Folder/1.sam
bamFolder=$bt2Folder/2.bam
mkdir -p $bt2Folder $samFolder $bamFolder

# # build-genome
# bowtie2-build --threads 40 $ws220Fasta $ws220Bt2Index

# bowtie2 mapping
echo Start bowtie2 mapping
awk -F',' 'NR!=1{print $1}' $sampleTableFile | \
    parallel -j 1 "echo Map {} reads using bowtie2 && \
    bowtie2 -p 40 -x $ws220Bt2Index -U $trimFolder/{}_R1_trimmed.fq.gz -S $samFolder/{}.sam 2> $samFolder/{}.log"
echo Finish bowtie2


# sam to sorted bam and index
echo SAM to sorted BAM and Index files
awk -F',' 'NR!=1{print $1}' $sampleTableFile | \
    parallel -j 1 "echo Sort and Index bam files of {} && \
    samtools view -Shu $samFolder/{}.sam -@ 40 | samtools sort - -@ 40 -T ~/temp/{}.temp | \
    tee $bamFolder/{}_sorted.bam | samtools index - $bamFolder/{}_sorted.bam.bai"
echo Finish bam conversion
    

##################################################################
### Use peaks idenfied by the modENCODE project
##################################################################

### 5. Overlapping peaks, distance to TSS
modEncodePeakFolder=/media/meshif/sda/Resources/Worm/ChIP/modencode/C.elegans/Transcriptional-Factor/ChIP-seq/computed-peaks_bed
modifiedBedFolder=$projectFolder/modEncode/modifiedBed
overlapPeakFolder=$projectFolder/modEncode/overlapPeak
distanceToTssFolder=$projectFolder/modEncode/disToTSS
mkdir -p $overlapPeakFolder $distanceToTssFolder

echo Identify ChIP peaks co-occupied by PHA-4 and HLH-30
echo Min distance to TSS
echo $code_Folder_Path
python $code_Folder_Path/1.peakCenter_to_tss.py $modEncodePeakFolder $modifiedBedFolder $overlapPeakFolder $distanceToTssFolder
echo Finish check Peak-to-TSS distance





##################################################################
### Use peaks idenfied by the modENCODE project
### 
##################################################################

### Metagene
metageneFolder=$ngsFolder/6.metagene
mkdir -p $metageneFolder

### (1) deeptools: bamCoverage ==> files for visualization
echo bamCoverage bamFile files
awk -F',' 'NR!=1{print $1}' $sampleTableFile | \
    parallel -j 10 "echo bamCoverage {} && bamCoverage --bam $bamFolder/{}_sorted.bam --outFileName $bamFolder/{}.bigWig --outFileFormat bigwig --binSize 5 -p 4 \
    2> $bamFolder/{}_bamCoverge.stderr"
echo finish bamCoverage!


### (2) deeptools: bamCompare
echo bamCompare samples
awk -F ',' 'NR!=1{print $1" "$2}' $chipTableFile | \
    while read line; do set -- $line; # read the whole line 
    targetName="$1"
    controlName="$2"

    target_Bam=$bamFolder/$targetName"_sorted.bam"
    control_Bam=$bamFolder/$controlName"_sorted.bam"
    # subtract
    outSubtractBigwig=$metageneFolder/$targetName"_"$controlName"_subtract.bigWig"
    outRatioBigwig=$metageneFolder/$targetName"_"$controlName"_log2.bigWig"

    echo ============================
    echo target_Bam: $target_Bam
    echo control_Bam: $control_Bam
    echo outSubtractBigwig: $outSubtractBigwig
    echo outRatioBigwig: $outRatioBigwig
    echo ============================
    echo
    # bamCompare -p 40 --operation subtract --outFileFormat bigwig --binSize 1 -b1 $target_Bam -b2 $control_Bam -o $outSubtractBigwig
    bamCompare -p 40 --outFileFormat bigwig --binSize 1 -b1 $target_Bam -b2 $control_Bam -o $outRatioBigwig \
        2> $metageneFolder/$targetName"_"$controlName"_bamCompare_log2.stderr";done


### (2) deeptools: computeMatrix & plotProfile
tssBed='/media/meshif/sda/Resources/Worm/TSS/2013_GR_Chen/2013_GR_Chen_majorTSS.sorted.bed'


# touch $metageneFolder/all_bigWig.txt
find $metageneFolder -name "*_log2.bigWig" | sort - > $metageneFolder/all_bigWig.txt

allBigwig=$(cat $metageneFolder/all_bigWig.txt)
echo deeptools: computeMatrix reference-point
computeMatrix reference-point -p 40 -R $tssBed -S $allBigwig \
        -b 3000 -a 3000 -o $metageneFolder/pha4hlh30_log2.mtx 1> $metageneFolder/allChIP_computeMatrix_log2.stdout


plotProfile -m $metageneFolder/pha4hlh30_log2.mtx \
    -o $metageneFolder/pha4hlh30_log2_plot_name.png --dpi 1200 \
    --plotWidth 9 --legendLocation best \
    --perGroup --refPointLabel TSS --samplesLabel PHA-4_rep1 PHA-4_rep2 HLH-30_rep1 HLH-30_rep2 \
    1> $metageneFolder/allChIP__plotProfile_log2.stdout

echo Finish plotting
rm $metageneFolder/all_bigWig.txt

mkdir -p $projectFolder/modEncode/1.chipTss/metagene
cp $metageneFolder/pha4hlh30_log2_plot_name.png $projectFolder/modEncode/1.chipTss/metagene


