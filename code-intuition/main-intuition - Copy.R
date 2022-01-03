# 2021 11 Andrew
# Just making a figure for the intro

# ENVIRONMENT ====

setwd('D:/OneDrive/t-hurdles/2021-07-PostOpenAP/code-intuition')

rm(list=ls())
library(data.table)
library(tidyverse)
library(ggplot2)
#library(ggthemes)
library(gridExtra)

# parameters
pnull = 0.5
N = 200 #1e4
alt_mu  = 2
alt_sig = 2

n_simulations <- 10

# function for creating histogram data hdat
sim = function(pnull, alt_mu, alt_sig){
  dat = data.frame(
    t = abs(rnorm(pnull*N*n_simulations)), group = F
  ) %>% rbind(
    data.frame(
      t = abs(rnorm((1-pnull)*N*n_simulations) + rexp((1-pnull)*N*n_simulations, rate = 1 / alt_sig) + alt_mu)
      , group = T
    )
  ) %>% 
    mutate(
      group = factor(
        group
        , levels = c(T,F)
      )
    )
  
} # end function sim

makehist = function(dat){
  edge = seq(0,10, 0.5)
  tmid = edge[1:(length(edge)-1)] + 0.25
  hdat2 = data.frame(
    bin = 1:length(tmid), tmid
  )
  
  hdat = dat %>% 
    filter(t>min(edge), t<max(edge)) %>% 
    group_by(group) %>% 
    mutate(
      bin = findInterval(t,edge)
    ) %>% 
    group_by(group,bin) %>% 
    summarize(
      n = n()/n_simulations
    ) %>% 
    left_join(hdat2)  
} # end function makehist

# HIGH FDR ====

dat = sim(pnull = 0.66, alt_mu, alt_sig)
hdat = makehist(dat)

fdr = dat %>% filter(t>1.96) %>% 
  summarize(fdr = round(mean(group == F)*100))


plot = ggplot(hdat, aes(x=tmid, y=n, fill=group)) + 
  geom_bar(stat='identity', position='stack') +
  geom_vline(xintercept = 1.96, size = 1) + 
  labs(title = "", x = "t-statistic", y = "Factors") +
  
  scale_fill_manual(
    labels = c("True Factor", "False Factor"), 
    values = c(rgb(0,0.4470,0.7410), rgb(0.8500, 0.3250, 0.0980)),
    name = ""
  ) +
  theme_light() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.title = NULL
  )  +
  annotate(
    geom = "text",
    label = "t = 1.96",
    x = 2.5,
    y = 40
  ) +
  geom_curve(aes(x = 3, y = 16, xend = 2.3, yend = 6),
             arrow = arrow(length = unit(0.03, "npc")),
             colour = "black", size = 0.3
  ) +
  annotate(
    geom = "text",
    label = paste0("FDR = ", fdr$fdr, '%'),
    x = 3.5,
    y = 16
  )

plot

ggsave(plot, filename = 'output/intro-high-FDR.pdf', width = 5, height = 4)



# plot all false ====

plot <- ggplot(
  dat %>% 
    filter(group | (!group & runif(1) < 0.5))
  , aes(x = t, fill = group)) + 
  
  geom_histogram(breaks = seq(0,4, 0.5)) +
  geom_vline(xintercept = 1.96, size = 1) + 
  labs(title = "", x = "t-statistic", y = "Factors") +
  
  scale_fill_manual(
    labels = c("False Factor"), 
    values = rgb(0.8500, 0.3250, 0.0980),
    name = ""
  ) +
  
  theme_light() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.title = NULL
  ) +
annotate(
  geom = "text",
  label = "t = 1.96",
  x = 2.5,
  y = 40
) + 
  geom_curve(aes(x = 3, y = 16, xend = 2.3, yend = 6),
             arrow = arrow(length = unit(0.03, "npc")), 
             colour = "black", size = 0.3
  ) + 
  annotate(
    geom = "text",
    label = "FDR = 100%",
    x = 3.5,
    y = 16
  )


plot


ggsave(plot, filename = 'output/intro-fdr-100.pdf', width = 5, height = 4)
