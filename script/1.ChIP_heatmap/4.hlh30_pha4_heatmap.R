# load library
suppressMessages({
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(ComplexHeatmap)
  library(pheatmap)
}) 

filenames <- commandArgs(trailingOnly=TRUE)
# print(filenames)

allChipFile <- filenames[1]
allChipNumFile <- filenames[2]
heatmapFolder <- filenames[3]

# import file
allChip.df = read.csv(allChipFile, header = TRUE)
allChipNum.df <- read.csv(allChipNumFile, header = TRUE, row.names = 1)
# (temp)
allChip.df = read.csv('/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip/modEncode/1.chipTss/chipSummary.csv', header = TRUE)
allChipNum.df <- read.csv('/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip/modEncode/1.chipTss/chipSummary_Num.csv', header = TRUE, row.names = 1)

# head(allChip.df)
rownames(allChip.df) <- allChip.df$seqName
# head(allChip.df)

allChip.df.mtx = as.matrix(allChip.df)
allChipNum.df.mtx <- as.matrix(allChipNum.df)
# head(allChipNum.df.mtx)
col = c(Bound = "Blue", Unbound = "grey")

# ## Normal complexheatmap
# Heatmap(allChipNum.df.mtx, show_row_names = FALSE, col = col ) #Error: data length exceeds size of matrix


## select L3, L4, YA stages for plotting heatmap
### Use Characters
print("Select L3, L4, YA stages for plotting heatmap")
afterL2Chip.df = allChip.df[, c(1,2,8,9,11)]
rownames(afterL2Chip.df) <- afterL2Chip.df$seqName
afterL2Chip.df$seqName <- NULL
# head(afterL2Chip.df %>% filter(HLH.30_L4 == 'Bound' | PHA.4_L3 == 'Bound' | PHA.4_L4 == 'Bound' | PHA.4_YA == 'Bound'))
afterL2Chip.df = afterL2Chip.df %>% filter(HLH.30_L4 == 'Bound' | PHA.4_L3 == 'Bound' | PHA.4_L4 == 'Bound' | PHA.4_YA == 'Bound')

png(filename = paste0(heatmapFolder, '/hlh-30_pha-4_binding_afterL2.png'),
    width = 2000, height = 3000, res = 300)
Heatmap(as.matrix(afterL2Chip.df), show_row_names = FALSE, col = col, name = 'TF_binding')
dev.off()

### Select L4 stage for plotting heatmap
print("Select L4 stage for plotting heatmap")
L4Chip.df = allChip.df[, c(1,2,9)]
rownames(L4Chip.df) <- L4Chip.df$seqName


#################################################################
### Sort L4 HLH-30 and PHA-4 signals for plotting simple heatmap
#################################################################
L4Chip.hlh30.pha4.df = L4Chip.df %>% filter(HLH.30_L4 == 'Bound' & PHA.4_L4 == 'Bound') %>% arrange(seqName)
L4Chip.hlh30.only.df = L4Chip.df %>% filter(HLH.30_L4 == 'Bound' & PHA.4_L4 == 'Unbound') %>% arrange(seqName)
L4Chip.pha4.only.df  = L4Chip.df %>% filter(HLH.30_L4 == 'Unbound' & PHA.4_L4 == 'Bound') %>% arrange(seqName)
allL4Chip.df = L4Chip.hlh30.pha4.df %>% bind_rows(L4Chip.hlh30.only.df) %>% bind_rows(L4Chip.pha4.only.df)
allL4Chip.df$seqName <- NULL

# heatmapFolder='/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip/modEncode/1.chipTss/heatmap'
png(filename = paste0(heatmapFolder, '/hlh-30_pha-4_L4_sorted.png'),
    width = 2000, height = 3000, res = 300)
Heatmap(as.matrix(allL4Chip.df), col = col, show_row_names = FALSE, name = 'TF_binding')
dev.off()



## Unsorted L4 heatmap
L4Chip.df$seqName <- NULL
# head(L4Chip.df)
# L4Chip.df = L4Chip.df %>% filter(HLH.30_L4 == 'Bound' | PHA.4_L3 == 'Bound' | PHA.4_L4 == 'Bound' | PHA.4_YA == 'Bound')

png(filename = paste0(heatmapFolder, '/hlh-30_pha-4_binding_L4.png'),
    width = 2000, height = 3000, res = 300)
Heatmap(as.matrix(L4Chip.df), show_row_names = FALSE, col = col, name = 'TF_binding')
dev.off()

# Use 0,1 # cause error
L4ChipNum.df = allChipNum.df[,c(1,8)]
# head(L4ChipNum.df %>% filter(HLH.30_L4 > 0 | PHA.4_L3 > 0 | PHA.4_L4>0 | PHA.4_YA>0 ))
L4ChipNum.df = L4ChipNum.df %>% filter(HLH.30_L4 > 0 | PHA.4_L4>0)
# Heatmap(as.matrix(L4ChipNum.df), show_row_names = FALSE, col = col) # Errors



### All ChIP records
print('All ChIP records')
allChip.df$seqName <- NULL
png(filename = paste0(heatmapFolder, '/hlh-30_pha-4_binding.png'),
    width = 2000, height = 3000, res = 300)
Heatmap(as.matrix(allChip.df), show_row_names = FALSE, col = col, name = 'TF_binding')
dev.off()

