library(tidyverse)

dwpc.data <- read.csv('out2_11Mar19.csv')

dwpc.data %>%
  ggplot(aes(x=metapath, y=DWPC, color=source)) +
  geom_jitter(width=0.2, size=2.5, alpha=0.7) +
  theme_bw() +
  theme(legend.pos='none') +
  labs(x=NULL, y='DWPC', title='Query 1 Gene Metapaths')
