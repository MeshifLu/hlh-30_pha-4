from os import system, chdir, getcwd, mkdir, listdir, path, stat
import csv
from glob import glob
import sys
from functools import partial
from multiprocessing import Pool

def bedtools_window_fun(chipBedFile, winSize, chipTssFolder, majorTssBed):
    # print('bedtools_window_fun, chipTssFolder: ', chipTssFolder)

    fileName = chipBedFile.split('/')[-1]
    TF = fileName.split('_')[0]
    # print(fileName)
    # print(fileName.split('_')[1].split('#'))
    stage, strain, condition = fileName.split('_')[1].split('#')

    ## Use bedtools window function to find TSS-proximal ChIP records
    windowBed=f'{chipTssFolder}/{fileName}'[:-4] + f'_TSS_{winSize}.bed'
    command=f'bedtools window -a {majorTssBed} -b "{chipBedFile}" -w {winSize} > "{windowBed}"'
    system(command)

def bedtools_window(bedFileFolder, chipTssFolder, majorTssBed, winSize):
    # print('bedtools_window, chipTssFolder: ', chipTssFolder)
    # print('winSize_L: ', winSize_L)

    globChipBedFile=glob(f'{bedFileFolder}/*.bed')

    p = Pool(processes=10)
    func = partial(bedtools_window_fun, chipTssFolder=chipTssFolder, majorTssBed=majorTssBed, winSize=winSize)
    result = p.map(func, sorted(globChipBedFile))
    p.close()
    p.join()

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




# def bedFilesIn(modEncodeChipFolder, bedFileFolder, chipTssFolder, majorTssBed, winSize):
#     # tfStageFolders=glob(f"{modEncodeChipFolder}/ce10/{TF}/*/")
#     # for stageFolder in sorted(tfStageFolders):
#     #     # print(stageFolder)
#     #     stageName=stageFolder.split('/')[-2]
#     #     # print(stageName)
#     #     chipBedGlobs=glob(f'{stageFolder}optimal IDR thresholded peaks/{TF}*.bed')
#     #
#     #     for chipBed in sorted(chipBedGlobs):
#     #         print(chipBed)
#     #         if 'chrName' not in chipBed:
#     #             stageName2 = chipBed.split('/')[-1].split('_')[1].split('.')[0]
#     #             bedtools_window(chipTssFolder, TF, stageName2, majorTssBed, chipBed, winSize_L)
#
#     globChipBedFile=glob(f'{bedFileFolder}/*.bed')
#     for chipBedFile in sorted(globChipBedFile):
#         fileName = chipBedFile.split('/')[-1]
#         TF = fileName.split('_')[0]
#         stage, strain, condition = fileName.split('_')[1].split('#')
#
#         bedtools_window(chipBedFile, chipTssFolder, TF, fileName, majorTssBed, winSize)
#     print()

if __name__ == '__main__':
    # print(sys.argv[1:])
    (modEncodeChipFolder, bedFileFolder, chipTssFolder, majorTssBed, winSize) = sys.argv[1:]
    # print(modEncodeChipFolder, bedFileFolder, chipTssFolder, majorTssBed)

    # bedFilesIn(modEncodeChipFolder, bedFileFolder, chipTssFolder, majorTssBed, winSize)
    bedtools_window(bedFileFolder, chipTssFolder, majorTssBed, winSize)



