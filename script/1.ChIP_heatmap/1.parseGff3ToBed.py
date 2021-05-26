from os import system, chdir, getcwd, mkdir, listdir, path, stat
import csv
from glob import glob
import sys
from functools import partial
from multiprocessing import Pool

def gffToBed(gff3File, bedFileFolder):

    fileName = gff3File.split('/')[-1]
    geneName = fileName.split('_')[0]
    if '#' in fileName:
        stage, strain, condition= fileName.split('_')[1].split('#')
    if '%23' in fileName:
        # File:
        # ZTF-7_Developmental-Stage=Larvae-L4-stage%23Strain=OP332%23temperature=20-degree-celsius_ChIP-seq_Rep-1__Cele_WS220_modENCODE_3219_combined.bed
        stage, strain, condition = fileName.split('_')[1].split('%23')
    if '%25252525252525252523' in fileName:
        stage, strain, condition = fileName.split('_')[1].split('%25252525252525252523')

    fileName='_'.join([geneName, '#'.join([stage, strain, condition])] +
              gff3File.split('/')[-1].split('_')[2:] )

    # print(geneName)
    # print(stage,strain, condition)
    # print()

    bedFilePath=f'{bedFileFolder}/{fileName}'[:-5] + '.bed'
    # print(bedFilePath)
    with open(bedFilePath, 'w') as OutFile1:
        with open(gff3File) as InFile1:
            for line1 in InFile1:
                if not line1.startswith('#'):
                    (chrN1, annot1, annot2, coordL1, coordR1, pvalue1, other1, other2, annot3) = line1.strip().split('\t')
                    bedOut='\t'.join([chrN1, coordL1, coordR1, geneName, pvalue1, '.',
                                      '#'.join([stage, strain, condition])])
                    OutFile1.write(bedOut + '\n')

    return fileName,bedFilePath

def sortBed(fileName, bedFilePath, sortedBedFileFolder):
    sortedBedFile=f"{sortedBedFileFolder}/{fileName[:-5]}_sorted.bed"

    command=f"bedtools sort -i {bedFilePath} > {sortedBedFile}"
    system(command)

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

if __name__ == '__main__':
    # print(sys.argv[1:])
    (gff3File, bedFileFolder, sortedBedFileFolder) = sys.argv[1:]

    ### gffToBed
    (fileName,bedFilePath) = gffToBed(gff3File, bedFileFolder)

    ### sortBed
    sortBed(fileName, bedFilePath, sortedBedFileFolder)


