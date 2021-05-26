from os import system, mkdir, path
import sys
from multiprocessing import Pool
from functools import partial
import pandas as pd
from glob import glob
import csv
import numpy as np

###

### Use average value of coordinates as the peak center (column 2:3 represents the coordinates of peak centers)
def parse_peakCenter(bedFileFolder, peakCenterBedFolder):
    globBedFile=glob(f"{bedFileFolder}/*.bed")
    for bedFile in sorted(globBedFile):
        bedFileName=bedFile.split('/')[-1][:-4]
        peakCenterBedFile=f"{peakCenterBedFolder}/{bedFileName}_peakCenter.bed"

        with open(peakCenterBedFile, 'w') as OutFile1:
            writer1=csv.writer(OutFile1, delimiter='\t')
            with open(bedFile) as InFile1:
                reader1=csv.reader(InFile1, delimiter='\t')
                for row1 in reader1:
                    peakCenter1= int( (float(row1[1])+float(row1[2])) /2 )
                    peakCenterCoordR1=peakCenter1+1
                    rowOut1= [row1[0], str(peakCenter1), str(peakCenterCoordR1)]+\
                        row1[3:] + row1[1:3]
                    writer1.writerow(rowOut1)

### Sort bed files
def sortBed(bedFileFolder, sortedBedFolder, peakCenterBedFolder, peakCenterSortedBedFolder):
    # Sort bed
    globBed=glob(f"{bedFileFolder}/*.bed")
    for bedFile in sorted(globBed):
        bedFileName=bedFile.split('/')[-1][:-4]
        command=f"bedtools sort -i {bedFile} > {sortedBedFolder}/{bedFileName}_sorted.bed"
        system(command)

    # Sort peakCenterBed
    globPeakcenterBed=glob(f"{peakCenterBedFolder}/*.bed")
    for peakCenterBedFile in sorted(globPeakcenterBed):
        bedFileName=peakCenterBedFile.split('/')[-1][:-4]
        command=f"bedtools sort -i {peakCenterBedFile} > {peakCenterSortedBedFolder}/{bedFileName}_sorted.bed"
        system(command)


def nearnestTSS(distanceToTssFolder, modifiedBedFolder):

    for gene in ['HLH-30', 'PHA-4']:
        # print(gene)
        eachTfBedFile = '%s/%s_L4_combined.sorted.bed' % (modifiedBedFolder, gene)
        tssBedFile = '/media/meshif/sda/Resources/Worm/TSS/2013_GR_Chen/2013_GR_Chen_WS220_chrN.sorted.bed'

        closestCommand = 'bedtools closest -D ref -a %s -b %s > %s/%s_L4_toTSS.bed' % \
                         (eachTfBedFile, tssBedFile, distanceToTssFolder, gene)
        system(closestCommand)




if __name__ == '__main__':
    # print(sys.argv[1:])

    "Confirm the following parameters!!!"
    (modEncodePeakFolder, modifiedBedFolder, overlapPeakFolder, distanceToTssFolder) = sys.argv[1:]
    (bedFileFolder, sortedBedFolder, peakCenterBedFolder, peakCenterSortedBedFolder) = sys.argv[1:]

    # ### parse_peakCenter
    # print("Peak centers")
    # parse_peakCenter(bedFileFolder, peakCenterBedFolder)

    # ### sort bed files
    # print("Sort bed files")
    # sortBed(bedFileFolder, sortedBedFolder, peakCenterBedFolder, peakCenterSortedBedFolder)

    "Upper function has been executed and files are put in the following folder:" \
    "/media/meshif/sda/Resources/Worm/ChIP/modencode/C.elegans/Transcriptional-Factor/ChIP-seq/"



    "To-do list:" \
    "1. Plot the binding signal of PHA-4 and HLH-30 around the TSS (Need more thinking)"
    "2. The nearest distance between PHA-4 and HLH-30 peaks"

    "Finish"



    # # nearnestTSS
    # nearnestTSS(distanceToTssFolder, modifiedBedFolder)






