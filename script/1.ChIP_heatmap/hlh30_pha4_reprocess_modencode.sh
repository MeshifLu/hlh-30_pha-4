### conda env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ngs


################################################################################################
### Use data from modEncode
###
################################################################################################


### Parameters
hlh30pha4Folder=/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip
hlh30pha4ModEncode=$hlh30pha4Folder/modEncode
modEncodeChipFolder=/media/meshif/sda/Resources/Worm/ChIP/modencode/C.elegans/Transcriptional-Factor/ChIP-seq
modEncodeGff3Folder=$modEncodeChipFolder/computed-peaks_gff3
hlh30codeFolder=/home/meshif/codeFolder/Project/HLH-30/20210414_reProcess_modEncode/1.ChIP_heatmap
TFs="hlh-30 pha-4"


## 0. Download .gff3 and decompress

## Download data
## Download peak files from modEnocde
## ftp://data.modencode.org/C.elegans/Transcriptional-Factor/ChIP-seq/computed-peaks_gff3/

find $modEncodeGff3Folder -name "*.gz" | xargs -n1 -P3 -I{} gzip -d {}


## 1. Process .gff3 files
## Process .gff3 files to .bed files and sorted.bed files
bedFileFolder=$modEncodeChipFolder/computed-peaks_bed
sortedBedFileFolder=$modEncodeChipFolder/computed-peaks_bed_sorted
mkdir -p $bedFileFolder $sortedBedFileFolder

echo Process and sort bedFiles
find $modEncodeGff3Folder -name "*.gff3" | xargs -n1 -P20 -I{} python $hlh30codeFolder/1.parseGff3ToBed.py {} $bedFileFolder $sortedBedFileFolder
echo Finish processing


## 2. ChIP-TSS # all ChIP data from modEncode

################################################
### Use bedtools windows function to find TSS sites within 1kb of ChIP sites
################################################

chipTssFolder=$modEncodeChipFolder/chipTss
mkdir -p $chipTssFolder
hlh30pha4ChipTss=$hlh30pha4ModEncode/1.chipTss
mkdir -p $hlh30pha4ChipTss
# TSS file from: http://www.genome.org/cgi/doi/10.1101/gr.153668.112
majorTssBed=/media/meshif/sda/Resources/Worm/TSS/2013_GR_Chen/2013_GR_Chen_majorTSS.sorted.bed
winSize_L="1000 3000 5000"

for winSize in $winSize_L; do
    chipTssWinSizeFolder=$chipTssFolder/$winSize
    mkdir -p $chipTssWinSizeFolder

    echo winSize: $winSize
    python $hlh30codeFolder/2.hlh30_pha4_chipSite_to_tss.py $hlh30pha4ModEncode $bedFileFolder $chipTssWinSizeFolder $majorTssBed $winSize
done

## copy hlh-30 pha-4 ChIP-TSS.bed file 

for winSize in $winSize_L; do
    chipTssWinSizeFolder=$chipTssFolder/$winSize
    mkdir -p $hlh30pha4ChipTss/$winSize
    echo winSize: $winSize
    find $chipTssWinSizeFolder -name 'HLH-30*' | xargs -n100 -I{} cp {} $hlh30pha4ChipTss/$winSize
    find $chipTssWinSizeFolder -name 'PHA-4*' | xargs -n100 -I{} cp {} $hlh30pha4ChipTss/$winSize
done


## 3. Summarize gene lists for heatmap
matrixFile=$hlh30pha4ChipTss/chipSummary.csv
matrixNumFile=$hlh30pha4ChipTss/chipSummary_Num.csv

python $hlh30codeFolder/3.parse_tss-chip_forEachTf.py $hlh30pha4ChipTss $matrixFile $matrixNumFile

## 4. heatmap for hlh-30 & pha-4
conda activate base
heatmapFolder=$hlh30pha4ChipTss/heatmap
mkdir -p $heatmapFolder

echo Plot ChIP binding patterns using heatmap
Rscript $hlh30codeFolder/4.hlh30_pha4_heatmap.R $matrixFile $matrixNumFile $heatmapFolder
echo Finish plotting heatmap


#####################################################################################
### Unable to plot heatmaps containing all conditions ==> cause errors in the program
### ==> use binding patterns instead of clustering results


## 5. peak overlapped by HLH-30 and PHA-4 in the L4 stage
conda activate ngs
modifiedBedFolder=$hlh30pha4Folder/modEncode/2.modifiedBed
overlapPeakFolder=$hlh30pha4Folder/modEncode/3.overlapPeak
mkdir -p $modifiedBedFolder $overlapPeakFolder

echo 1. Modify peak Name
echo 2. Find overlapping peaks
python $hlh30codeFolder/5.peakOverlap.py $sortedBedFileFolder $modifiedBedFolder $overlapPeakFolder
echo Finish Modification
echo Finish idenfying overlapping peaks


### Summary of overlap peaks
# 	            HLH-30	PHA-4
# Total	        2791	3890
# Overlap	    1379    1475
# Non-overlap	1412	2415


# 5-2. vennDiagram
conda activate base
Rscript $hlh30codeFolder/5-2.hlh-30_pha-4_L4venn.R $overlapPeakFolder $overlapPeakFolder/Overlap_PHA-4_HLH-30.bed $overlapPeakFolder/pha4_noOverlap.list $overlapPeakFolder/hlh30_noOverlap.list




