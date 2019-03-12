library(tidyverse)

dwpc.data <- read.csv('query1_DWPC_11Mar19.csv')

dwpc.data %>%
  mutate(log2_norm_dwpc = log2(DWPC/median(DWPC))) %>%
  ggplot(aes(x=metapath, y=DWPC, color=source)) +
  geom_jitter(width=0.1, size=2.5, alpha=0.7) +
  theme_bw() +
  theme(legend.pos='none', 
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x='Metapath', y='log2 Normalized DWPC', title='Query 1 Gene Metapaths')


dwpc.data2 <- read.csv('PC1_DWPC_11Mar19.csv')

dwpc.data2 %>%
  mutate(log2_norm_dwpc = log2(DWPC/median(DWPC))) %>%
  ggplot(aes(x=metapath, y=DWPC, color=source)) +
  geom_jitter(width=0.1, size=2.5, alpha=0.7) +
  theme_bw() +
  theme(legend.pos='none', 
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x='Metapath', y='log2 Normalized DWPC', title='PC 1 Gene Metapaths')

dwpc.data3 <- read.csv('random_DWPC_11Mar19.csv')

dwpc.data3 %>%
  mutate(log2_norm_dwpc = log2(DWPC/median(DWPC))) %>%
  ggplot(aes(x=metapath, y=DWPC, color=source)) +
  geom_jitter(width=0.1, size=2.5, alpha=0.7) +
  theme_bw() +
  theme(legend.pos='none', 
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x='Metapath', y='log2 Normalized DWPC', title='Random Gene Metapaths')

dwpc.total <- rbind(dwpc.data, dwpc.data2) %>%
  rbind(dwpc.data3) %>%
  left_join(.,gene.dict, by = c('source'='ID')) %>%
  rename('source.name'='name') %>%
  select(-query) %>%
  left_join(.,gene.dict, by = c('target'='ID')) %>%
  rename('target.name'='name') %>%
  select(-X.x, -X.y, -X, -pair)

write.csv(dwpc.total, 'DWPC_12Mar19.csv')

dwpc.total %>%
  mutate(log2_norm_dwpc = log2(DWPC/median(DWPC))) %>%
  ggplot(aes(x=metapath, y=DWPC, color=query)) +
  geom_jitter(width=0.1, size=2.5, alpha=0.4) +
  theme_bw() +
  theme(legend.pos='right', 
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x='Metapath', y='log2 Normalized DWPC', title='Gene Metapaths', color='Query')
