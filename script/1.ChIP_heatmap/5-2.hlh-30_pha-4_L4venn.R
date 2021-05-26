suppressMessages({
  library(eulerr)
  library(tidyverse)
  library(cowplot)
})

###
filenames <- commandArgs(trailingOnly = TRUE)
# print(filenames)
overlapPeakFolder <- filenames[1]
hlh30pha4OverlapBed <- filenames[2]
pha4_noOverlapFile  <- filenames[3]
hlh30_noOverlapFile <- filenames[4]

print(hlh30pha4OverlapBed)

### overlap
print('Process overlap')
overlap.df = read.delim(hlh30pha4OverlapBed, header = FALSE)
# # (temp)
# overlap.df = read.delim('/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip/modEncode/3.overlapPeak/Overlap_PHA-4_HLH-30.bed', header = FALSE)
overlap.df = overlap.df %>% distinct(V4, .keep_all = TRUE) # 1491 to 1475
overlap.df = overlap.df %>% select(V4, V11)
colnames(overlap.df) = c("PHA-4_peaks", "HLH-30_peaks")
overlap.df$`PHA-4_peaks` = TRUE
overlap.df$`HLH-30_peaks` = TRUE

### pha4-only
print('Process PHA-4')
pha4Only.df = read.csv(pha4_noOverlapFile, header = FALSE)
# # (temp)
# pha4Only.df = read.csv('/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip/modEncode/3.overlapPeak/pha4_noOverlap.list', header = FALSE)
colnames(pha4Only.df) = "PHA-4_peaks"
pha4Only.df$`PHA-4_peaks` = TRUE
pha4Only.df["HLH-30_peaks"] = FALSE

### hlh30-only
print('Process HLH-30')
hlh30Only.df = read.csv(hlh30_noOverlapFile, header = FALSE)
# # (temp)
# hlh30Only.df = read.csv('/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip/modEncode/3.overlapPeak/hlh30_noOverlap.list', header = FALSE)
colnames(hlh30Only.df) = "HLH-30_peaks"
hlh30Only.df$`HLH-30_peaks`=TRUE
hlh30Only.df["PHA-4_peaks"] = FALSE
hlh30Only.df = hlh30Only.df[,c(2,1)]

allPeaks.df = overlap.df %>% bind_rows(pha4Only.df) %>% bind_rows(hlh30Only.df)

print('Output')
png(filename = paste0( overlapPeakFolder, '/vennDiagram.png' ), 
    res = 300, width = 1800, height = 1000)
plot(euler(allPeaks.df), quantities = TRUE)
dev.off()



