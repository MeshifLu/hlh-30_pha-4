from os import system, mkdir, path
import sys
from multiprocessing import Pool
from functools import partial
import pandas as pd
from glob import glob
import csv
import numpy as np

###

def modifyPeakName( sortedBedFileFolder, modifiedBedFolder):
    for gene in ['HLH-30', 'PHA-4']:
        print(gene)

        globBedFile = glob(f'{sortedBedFileFolder}/{gene}*L4*_combined*bed')
        for bedFile in sorted(globBedFile):
            print(bedFile)

            # if "HLH-30" in bedFile:
            #     hlh30BedPath = bedFile
            # else:
            #     pha4BedPath = bedFile

            with open(f"{modifiedBedFolder}/{gene}_L4_combined.sorted.bed", 'w') as OutFile1:
                writer1=csv.writer(OutFile1, delimiter='\t')

                with open(bedFile) as InFile1:
                    reader1=csv.reader(InFile1, delimiter='\t')
                    for row1 in reader1:
                        # print(row1)
                        row1[3] = '_'.join([row1[3]]+row1[:3])
                        # print(row1[3])
                        writer1.writerow(row1)



def overlappingPeak( modifiedBedFolder, overlapPeakFolder ):
    for gene in ['HLH-30', 'PHA-4']:
        print(gene)

        bedFile=f'{modifiedBedFolder}/{gene}_L4_combined.sorted.bed'
        # print(bedFile)
        if "HLH-30_L4" in bedFile: hlh30BedPath=bedFile
        else: pha4BedPath=bedFile

    ## overlapPha4Hlh30
    hlh30pha4OverlapBed = '%s/Overlap_PHA-4_HLH-30.bed' % (overlapPeakFolder)
    print(hlh30BedPath)
    print(pha4BedPath)

    print('Intersect Two factors')
    bedIntersectCommnad=f'bedtools intersect -wa -wb -a "{pha4BedPath}" -b "{hlh30BedPath}" > "{hlh30pha4OverlapBed}"'
    # f'bedtools intersect -wa -wb -a {pha4BedPath} -b {hlh30BedPath} > {hlh30pha4OverlapBed}'
    print(bedIntersectCommnad)
    system(bedIntersectCommnad)



def summerizeOverlap(modifiedBedFolder, overlapPeakFolder):
    print('Summarize results')

    # overlapPha4Hlh30
    pha4BedPath = '%s/PHA-4_L4_combined.sorted.bed' % (modifiedBedFolder)
    hlh30BedPath='%s/HLH-30_L4_combined.sorted.bed' % (modifiedBedFolder)
    hlh30pha4OverlapBed = '%s/Overlap_PHA-4_HLH-30.bed' % (overlapPeakFolder)

    (pha4_overlapped_L, pha4_noOverlap_L, hlh30_overlapped_L, hlh30_noOverlap_L)= ([], [], [], [])

    ### Use bedtools-intersect positive peaks as overlapping peaks
    with open(hlh30pha4OverlapBed) as inOverlapPeak:
        for lineOverlap in inOverlapPeak:
            pha4PeakName = lineOverlap.rstrip().split('\t')[3]
            hlh30PeakName = lineOverlap.rstrip().split('\t')[-4]
            # print(pha4PeakName)
            # print(hlh30PeakName)
            # print()

            if pha4PeakName not in pha4_overlapped_L: pha4_overlapped_L.append(pha4PeakName)
            if hlh30PeakName not in hlh30_overlapped_L: hlh30_overlapped_L.append(hlh30PeakName)


    # # Check peak numbers in overlapping list
    # for gene in ['HLH-30', 'PHA-4']:
    #     # print(gene)
    #     eachTfBedFile = '%s/%s_L4_combined.sorted.bed' % (modifiedBedFolder, gene)


    ### Non-overlapping peaks
    for bedFile in [pha4BedPath, hlh30BedPath]:
        with open(bedFile) as tfPeakFile:
            for lineMergedPeak in tfPeakFile:
                peakName = lineMergedPeak.rstrip().split('\t')[3]
                # print(peakName)
                if 'HLH-30' in peakName:
                    if peakName not in hlh30_overlapped_L:
                        hlh30_noOverlap_L.append(peakName)
                else:
                    if peakName not in pha4_overlapped_L:
                        pha4_noOverlap_L.append(peakName)

    #
    outFile_D = {}
    outFile_D['pha4_overlapped_L'] = pha4_overlapped_L
    outFile_D['pha4_noOverlap_L'] = pha4_noOverlap_L
    outFile_D['hlh30_overlapped_L'] = hlh30_overlapped_L
    outFile_D['hlh30_noOverlap_L'] = hlh30_noOverlap_L

    for outFileName in sorted(outFile_D.keys()):
        outFilePath = '%s/%s.list' % (overlapPeakFolder, outFileName[:-2])
        outFile = open(outFilePath, 'w')
        for chipRecord in sorted(outFile_D[outFileName]): outFile.write(chipRecord + '\n')



if __name__ == '__main__':
    # print(sys.argv[1:])
    (sortedBedFileFolder, modifiedBedFolder, overlapPeakFolder) = sys.argv[1:]

    # modifyPeakName
    modifyPeakName(sortedBedFileFolder, modifiedBedFolder)

    # overlappingPeak
    overlappingPeak(modifiedBedFolder, overlapPeakFolder)

    # summerizeOverlap
    summerizeOverlap(modifiedBedFolder, overlapPeakFolder)












