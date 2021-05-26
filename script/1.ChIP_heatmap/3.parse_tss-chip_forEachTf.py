from os import system, chdir, getcwd, mkdir, listdir, path, stat
import csv
from glob import glob
import sys
from functools import partial
from multiprocessing import Pool


if __name__ == '__main__':
    # print(sys.argv[1:])
    (hlh30pha4ChipTss, matrixFile, matrixNumFile) = sys.argv[1:]
    # print(modErnFolder, TF, hlh30pha4ChipTss, majorTssBed)

    ###
    ### Use short name to represent each ChIP file
    ###
    gene_L=[]
    chipFile_D={}
    globchipFile=glob(f'{hlh30pha4ChipTss}/1000/*_1000.bed') # Use 1kb file
    for chipFile in sorted(globchipFile):
        # print(chipFile)
        TF=chipFile.split('/')[-1].split('_')[0]
        stageNameLong, strain, condition = chipFile.split('/')[-1].split('_')[1].split('#')
        # print(TF, stageNameLong)

        stageNameLong=stageNameLong.split('=')[-1]

        stageNameShort=''
        while stageNameShort == '':
            if 'Embryo' in stageNameLong:
                if 'Late' in stageNameLong:
                    stageNameShort='LE'
                else:
                    stageNameShort = 'Em'
            else:
                if 'adult' in stageNameLong:
                    stageNameShort='YA'
                else:
                    if 'Starved' in stageNameLong or 'Fed' in stageNameLong:
                        if 'Fed' in stageNameLong: stageNameShort='L1-Fed'
                        else: stageNameShort='L1-Starved'
                    else:
                        stageNameShort=stageNameLong.split('-')[1]

        chipDicName=f'{TF}_{stageNameShort}'
        if chipDicName not in chipFile_D.keys():
            chipFile_D[chipDicName]=[chipFile]
        else:
            chipFile_D[chipDicName].append(chipFile)

        # print('stageNameLong: ', stageNameLong)
        # print('stageNameShort: ', stageNameShort)
        # print('chipDicName: ', chipDicName)
        # print(chipFile.split('/')[-1])

    ###
    ### Use combined peak to represent the overall signal
    ### If there's no combined file, choose one file to represent the whole data set
    ### PHA-4 L3: use Rep-2
    for chipDicName in sorted(chipFile_D.keys()):
        checkMergeFile='N'
        # while checkMergeFile =='N':

        if len(chipFile_D[chipDicName]) != 1:
            # print(chipFile_D[chipDicName])
            # for i in sorted(chipFile_D[chipDicName]):
            #     print(i)

            if 'combined' not in str(chipFile_D[chipDicName]):
                # print(chipFile_D[chipDicName])
                chipFile_D[chipDicName]=chipFile_D[chipDicName][1]
                # print(chipFile_D[chipDicName])
            else:
                for chipFile in chipFile_D[chipDicName]:
                    if 'combined' in chipFile:
                        chipFile_D[chipDicName] = chipFile
        else:
            chipFile_D[chipDicName]=chipFile_D[chipDicName][0]

        # print(chipDicName)
        # print(chipFile_D[chipDicName])
        # print()


    ###
    ### Summarise all ChIP records into one single table
    ###
    chipRecordGene_L=[]
    chipRecordEachFile_D={} # To trace genes reported in each stage(condition)

    ### Report all genes found in ChIP exp
    for chipDicName in sorted(chipFile_D.keys()):
        chipFile = chipFile_D[chipDicName]
        chipRecordEachFile_D[chipDicName] = []

        with open(chipFile) as InFile1:
            for line1 in InFile1:
                line1_L = line1.strip().split('\t')
                geneID1 = line1_L[9]
                if geneID1 not in chipRecordGene_L:
                    chipRecordGene_L.append(geneID1)
                if geneID1 not in chipRecordEachFile_D[chipDicName]:
                    chipRecordEachFile_D[chipDicName].append(geneID1)

    ### output result into one .csv file
    OutFile1=open(f'{matrixFile}', 'w')
    OutFile2 = open(f'{matrixNumFile}', 'w')
    writer1 = csv.writer(OutFile1)
    writer2 = csv.writer(OutFile2)

    writer1.writerow(['seqName'] + sorted(chipRecordEachFile_D.keys()))
    writer2.writerow(['seqName'] + sorted(chipRecordEachFile_D.keys()))

    for geneID1 in sorted(chipRecordGene_L):
        out1_L = [geneID1]
        out2_L = [geneID1]
        for chipRecordEachFile in sorted(chipRecordEachFile_D.keys()):
            # print(geneID1)
            # print(chipRecordEachFile_D[chipRecordEachFile])
            # print()
            if geneID1 in chipRecordEachFile_D[chipRecordEachFile]:
                out1_L.append('Bound')
                out2_L.append('1')
            else:
                out1_L.append('Unbound')
                out2_L.append('0')

        writer1.writerow(out1_L)
        writer2.writerow(out2_L)

