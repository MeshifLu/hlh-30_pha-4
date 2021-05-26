suppressMessages({
  library(tidyverse)
  library(cowplot)
})

filenames <- commandArgs(trailingOnly=TRUE)
hlh30pha4DistanceFolder <- filenames[1]
# print(filenames)

theme_set(theme_cowplot())

pha4_hlh30_dis.df<- read.csv(file = paste0(hlh30pha4DistanceFolder, '/distance_PHA-4_HLH-30.csv'))
# pha4_hlh30_dis.df=read.csv(file = '/media/meshif/NAS/Project/HLH-30/hlh30_pha4_chip/modEncode/4.hlh30_pha4_distance/distance_PHA-4_HLH-30.csv')
pha4_hlh30_dis.df$Distance = pha4_hlh30_dis.df$Distance/1000 # Use kb to represent the distance


png(filename = paste0(hlh30pha4DistanceFolder, '/distance_PHA-4_HLH-30_L4.png'),
    width = 1000, height = 1000, res = 300)
print(
  pha4_hlh30_dis.df %>%
  ggplot(aes(Distance))+
  geom_histogram(binwidth = 0.1)+
  coord_cartesian(xlim = c(-3, 3))+
  labs(x='Distance (kb)')+
  scale_x_continuous(n.breaks = 7)
)

dev.off()

pha4_hlh30_dis.df$within250=abs(pha4_hlh30_dis.df$Distance) < 0.250
sum250.df <- pha4_hlh30_dis.df %>% group_by(within250) %>% summarise(n250=n())
pha4_hlh30_dis.df$within500=abs(pha4_hlh30_dis.df$Distance) < 0.500
sum500.df <- pha4_hlh30_dis.df %>% group_by(within500) %>% summarise(n500=n())
pha4_hlh30_dis.df$within1000=abs(pha4_hlh30_dis.df$Distance) < 1
sum1000.df <- pha4_hlh30_dis.df %>% group_by(within1000) %>% summarise(n1000=n())

allSummary.df = bind_cols(sum250.df, sum500.df, sum1000.df)
write.csv(allSummary.df, file = paste0(hlh30pha4DistanceFolder, '/distance_PHA-4_HLH-30_summary.csv'), row.names = FALSE)



