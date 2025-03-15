library(tidyverse)
library(ggplot2)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

sol1 <- read.csv('../output_table/scenario1.txt') %>%
  pivot_longer(cols = -t, names_to = 'type', values_to = 'value') %>%
  mutate(category = str_sub(type, 1, 1),
         community = str_sub(type, 2, -1) %>% as.numeric(),
         society = floor(community/4) + 1,
         category = case_when(
           category == 'S' ~ 'susceptible',
           category == 'E' ~ 'exposed',
           category == 'I' ~ 'infectious',
           category == 'R' ~ 'recovered'
         ) %>%
           factor(levels = c('susceptible', 'exposed', 'infectious', 'recovered')),
         scenario = 1)

sol2 <- read.csv('../output_table/scenario2.txt') %>%
  pivot_longer(cols = -t, names_to = 'type', values_to = 'value') %>%
  mutate(category = str_sub(type, 1, 1),
         community = str_sub(type, 2, -1) %>% as.numeric(),
         society = floor(community/4) + 1,
         category = case_when(
           category == 'S' ~ 'susceptible',
           category == 'E' ~ 'exposed',
           category == 'I' ~ 'infectious',
           category == 'R' ~ 'recovered'
         ) %>%
           factor(levels = c('susceptible', 'exposed', 'infectious', 'recovered')), 
         scenario = 2)

sol3 <- read.csv('../output_table/scenario3.txt') %>%
  pivot_longer(cols = -t, names_to = 'type', values_to = 'value') %>%
  mutate(category = str_sub(type, 1, 1),
         community = str_sub(type, 2, -1) %>% as.numeric(),
         society = floor(community/4) + 1,
         category = case_when(
           category == 'S' ~ 'susceptible',
           category == 'E' ~ 'exposed',
           category == 'I' ~ 'infectious',
           category == 'R' ~ 'recovered'
         ) %>%
           factor(levels = c('susceptible', 'exposed', 'infectious', 'recovered')), 
         scenario = 3)

sol4 <- read.csv('../output_table/scenario4.txt') %>%
  pivot_longer(cols = -t, names_to = 'type', values_to = 'value') %>%
  mutate(category = str_sub(type, 1, 1),
         community = str_sub(type, 2, -1) %>% as.numeric(),
         society = floor(community/4) + 1,
         category = case_when(
           category == 'S' ~ 'susceptible',
           category == 'E' ~ 'exposed',
           category == 'I' ~ 'infectious',
           category == 'R' ~ 'recovered'
         ) %>%
           factor(levels = c('susceptible', 'exposed', 'infectious', 'recovered')), 
         scenario = 4)

sol5 <- read.csv('../output_table/scenario5.txt') %>%
  pivot_longer(cols = -t, names_to = 'type', values_to = 'value') %>%
  mutate(category = str_sub(type, 1, 1),
         community = str_sub(type, 2, -1) %>% as.numeric(),
         society = floor(community/4) + 1,
         category = case_when(
           category == 'S' ~ 'susceptible',
           category == 'E' ~ 'exposed',
           category == 'I' ~ 'infectious',
           category == 'R' ~ 'recovered'
         ) %>%
           factor(levels = c('susceptible', 'exposed', 'infectious', 'recovered')), 
         scenario = 5)

sol6 <- read.csv('../output_table/scenario6.txt') %>%
  pivot_longer(cols = -t, names_to = 'type', values_to = 'value') %>%
  mutate(category = str_sub(type, 1, 1),
         community = str_sub(type, 2, -1) %>% as.numeric(),
         society = floor(community/4) + 1,
         category = case_when(
           category == 'S' ~ 'susceptible',
           category == 'E' ~ 'exposed',
           category == 'I' ~ 'infectious',
           category == 'R' ~ 'recovered'
         ) %>%
           factor(levels = c('susceptible', 'exposed', 'infectious', 'recovered')), 
         scenario = 6)


p1 <- ggplot(sol1, aes(x=t, y=value, group=category, color=category)) +
  facet_wrap(society ~ community, labeller = labeller(society = label_both, community = label_both), ncol = 4) +
  geom_line() +
  scale_color_manual(values=c('#3c7a6e', '#ffa500', '#990000', 'navy'))

p2 <- ggplot(sol2, aes(x=t, y=value, group=category, color=category)) +
  facet_wrap(society ~ community, labeller = labeller(society = label_both, community = label_both), ncol = 4) +
  geom_line() +
  scale_color_manual(values=c('#3c7a6e', '#ffa500', '#990000', 'navy'))


pdf('../figs/scenarios_12.pdf', width=8, height=12)
ggarrange(p1, p2, ncol = 1, labels = c('A', 'B'))
dev.off()

p3 <- ggplot(sol3, aes(x=t, y=value, group=category, color=category)) +
  facet_wrap(society ~ community, labeller = labeller(society = label_both, community = label_both), ncol = 4) +
  geom_line() +
  scale_color_manual(values=c('#3c7a6e', '#ffa500', '#990000', 'navy'))

p4 <- ggplot(sol4, aes(x=t, y=value, group=category, color=category)) +
  facet_wrap(society ~ community, labeller = labeller(society = label_both, community = label_both), ncol = 4) +
  geom_line() +
  scale_color_manual(values=c('#3c7a6e', '#ffa500', '#990000', 'navy'))

pdf('../figs/scenarios_34.pdf', width=8, height=12)
ggarrange(p3, p4, ncol = 1, labels = c('A', 'B'))
dev.off()

p5 <- ggplot(sol5, aes(x=t, y=value, group=category, color=category)) +
  facet_wrap(society ~ community, labeller = labeller(society = label_both, community = label_both), ncol = 4) +
  geom_line() +
  scale_color_manual(values=c('#3c7a6e', '#ffa500', '#990000', 'navy'))

p6 <- ggplot(sol6, aes(x=t, y=value, group=category, color=category)) +
  facet_wrap(society ~ community, labeller = labeller(society = label_both, community = label_both), ncol = 4) +
  geom_line() +
  scale_color_manual(values=c('#3c7a6e', '#ffa500', '#990000', 'navy'))

pdf('../figs/scenarios_56.pdf', width=8, height=12)
ggarrange(p5, p6, ncol = 1, labels = c('A', 'B'))
dev.off()

halfpoints <- read.csv('../output_table/half_points_scenario1-6.txt') %>%
   mutate(community = row_number(),
          society = ceiling(community/4) %>% as.factor())  %>%
   pivot_longer(cols = starts_with('Scenario'), names_to = 'Scenario', values_to = 'halfpoint')


pdf('../figs/scenarios_1-6_hp.pdf', width=8, height=6)
ggplot(halfpoints, aes(x = halfpoint, y = Scenario, color = society)) +
  geom_jitter(width=0, height=0.3, size=3) +
  scale_color_brewer(palette = 'Dark2') +
  xlim(0,100) +
  xlab('Halfpoint (days)') +
  ylab('') +
  theme_bw(18)
dev.off()

# df1 <- read.csv('../output_table/trudgill.txt') %>%
#   mutate(group = row_number(),
#          population = ceiling(group/4) %>% as.factor()) %>%
#   pivot_longer(cols = starts_with('Type'), names_to = 'type', values_to = 'halfpoint')
# 
# df2 <- read.csv('../output_table/simple.txt') %>%
#   mutate(population = row_number() %>% as.factor(),
#          group = 1) %>%
#   pivot_longer(cols = starts_with('Type'), names_to = 'type', values_to = 'halfpoint')
# 
# df3 <- read.csv('../output_table/XS.txt') %>%
#   mutate(group = row_number(),
#          population = ceiling(group/4) %>% as.factor()) %>%
#   pivot_longer(cols = starts_with('Scenario'), names_to = 'type', values_to = 'halfpoint')
# 
# df <- rbind(df1, df2, df3) %>%
#   mutate(type = factor(type, 
#                        levels = c(paste0('Type', 1:6), paste0('Type', 1:3, '_simple'), paste0('Scenario', 1:4)))) %>%
#   filter(!(type == 'Type1_simple' & population != 1))
# 
# ggplot(df1, aes(x = halfpoint, y = type, color = population)) +
#   geom_jitter(width=0, height=0.3, size=3) +
#   scale_color_brewer(palette = 'Dark2') +
#   xlim(0,100) +
#   xlab('Halfpoint (days)') +
#   ylab('Type') +
#   theme_bw(18)

