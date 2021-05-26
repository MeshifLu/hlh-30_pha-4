from os import system, mkdir, path
import sys
from multiprocessing import Pool
from functools import partial
import pandas as pd
from glob import glob
import csv
import numpy as np


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

### process distances between PHA-4 and HLH-30 peaks
def closestPeak(peakCenterSortedBedFolder, hlh30pha4DistanceFolder):
    for TF in ['HLH-30', 'PHA-4']:
        tfPeakcenterBed=glob(f'{peakCenterSortedBedFolder}/{TF}_*combined*')[0]
        # print(tfPeakcenterBed)

        if TF == 'HLH-30': hlh30peakBed=tfPeakcenterBed
        else: pha4peakBed=tfPeakcenterBed

    command=f'bedtools closest -a {pha4peakBed} -b {hlh30peakBed} > ' \
            f'{hlh30pha4DistanceFolder}/closeset_PHA-4_HLH-30.bed'
    system(command)


### calculate distances between PHA-4 and HLH-30 peaks
### Use peak center to calculate peak distances
def calDistance(hlh30pha4DistanceFolder):
    peakDis_D={}

    with open(f'{hlh30pha4DistanceFolder}/distance_PHA-4_HLH-30.csv', 'w') as OutFile1:
        writer1=csv.writer(OutFile1)
        writer1.writerow( ['peakName', 'Distance'] )

        with open(f'{hlh30pha4DistanceFolder}/closeset_PHA-4_HLH-30.bed') as InFile1:
            reader1=csv.reader(InFile1, delimiter='\t')
            for row1 in reader1:
                pha4PeakName='_'.join(row1[:4])
                # print(row1)
                distance=int(row1[-2]) - int(row1[1]) +1

                if pha4PeakName not in peakDis_D.keys():
                    peakDis_D[pha4PeakName] = distance
                else:
                    if peakDis_D[pha4PeakName] > distance:
                        peakDis_D[pha4PeakName] = distance

            for pha4PeakName in sorted(peakDis_D.keys()):
                writer1.writerow([ pha4PeakName, str(peakDis_D[pha4PeakName]) ])




if __name__ == '__main__':
    # print(sys.argv[1:])
    (bedFileFolder, sortedBedFolder, peakCenterBedFolder, peakCenterSortedBedFolder,
     hlh30pha4DistanceFolder) = sys.argv[1:]

    ### parse_peakCenter
    print("Peak centers")
    parse_peakCenter(bedFileFolder, peakCenterBedFolder)

    ### sort bed files
    print("Sort bed files")
    sortBed(bedFileFolder, sortedBedFolder, peakCenterBedFolder, peakCenterSortedBedFolder)

    ### closestPeak
    print('Find HLH-30 peaks which are most close to the PHA-4 peaks')
    closestPeak(peakCenterSortedBedFolder, hlh30pha4DistanceFolder)

    ### calDistance
    print('Calculate min distances between PHA-4 and HLH-30 peaks')
    calDistance(hlh30pha4DistanceFolder)




