### conda env
source ~/anaconda3/etc/profile.d/conda.sh
conda activate ngs

### code
hlh30pha4CodeFolder='/home/meshif/codeFolder/Project/HLH-30/20210414_reProcess_modEncode'
myPyLibFolder='/home/meshif/codeFolder/myPyLib'


################################################################################################
### Calculate distance between pha-4 and hlh-30 peaks
### Use data from modEncode
################################################################################################

### Parameters
hlh30pha4Folder=/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip
hlh30pha4ModEncode=$hlh30pha4Folder/modEncode
modEncodeChipFolder=/media/meshif/sda/Resources/Worm/ChIP/modencode/C.elegans/Transcriptional-Factor/ChIP-seq
hlh30codeRolder=/home/meshif/codeFolder/Project/HLH-30/20210414_reProcess_modEncode
bedFileFolder=$modEncodeChipFolder/computed-peaks_bed
TFs="HLH-30 PHA-4"


## 1. (1)Assign peak centers, (2) sort bed files, (3) calculate PHA-4-HLH-30 distances
peakCenterBedFolder=$modEncodeChipFolder/computed-peaks_bed_peakCenter
hlh30pha4DistanceCodeFoler=$hlh30pha4CodeFolder/3.hlh30_pha4_distance
sortedBedFolder=$modEncodeChipFolder/computed-peaks_bed_sorted
peakCenterSortedBedFolder=$modEncodeChipFolder/computed-peaks_bed_peakCenter_sorted
hlh30pha4DistanceFolder=$hlh30pha4ModEncode/4.hlh30_pha4_distance
mkdir -p $hlh30pha4DistanceCodeFoler $sortedBedFolder $peakCenterSortedBedFolder $hlh30pha4DistanceFolder

echo Assign peak centers
python $hlh30pha4DistanceCodeFoler/1.calculate_disToPeakCenter.py $bedFileFolder $sortedBedFolder $peakCenterBedFolder $peakCenterSortedBedFolder $hlh30pha4DistanceFolder
echo Finish assigning peak centers

## 2. Plot distance between PHA-4 and HLH-30 binding sites

conda activate base
echo Plot HLH-30 to PHA-4 distance
Rscript $hlh30pha4DistanceCodeFoler/2.hlh30_to_pha4_distance.R $hlh30pha4DistanceFolder
echo Finish plotting

