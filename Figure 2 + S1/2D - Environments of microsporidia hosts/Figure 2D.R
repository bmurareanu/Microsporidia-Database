
#
# Creating venn diagram of microsporidia host environments
#
# Jason Jiang - Created: 2020/07/07
#                     Last Edit: 2021/05/03
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code to clean up host environment data and create a venn diagram
# showing the overlap and frequency of microsporidia host habitats.
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------


# Load useful packages
library(tidyverse)
library(readxl)
library(VennDiagram)


# Set the working directory for R. This is what I use
setwd("P:/Shared/Microsporidia database/Figure 2/2D - Environments of microsporidia hosts")


# Load in host masterlist spreadsheet, which contains host environment data
host_envs <- read_xlsx("Table S2 Host Masterlist.xlsx")

host_envs <- host_envs %>%
  select(Environment) %>% # Only consider host environment data
  filter(!is.na(Environment))


# Create table of all unique environment combinations, and their frequencies
host_env_frequencies <- as.data.frame(table(host_envs))

# arbitrary areas were assigned to venn diagram, and area labels were hidden
env_venn_diagram <- draw.quad.venn(area1 = 100, # Terrestrial
                                  area2 = 200,
                                  area3 = 300,
                                  area4 = 400,
                                  n12 = 2,
                                  n13 = 2,
                                  n14 = 10,
                                  n23 = 2,
                                  n24 = 2,
                                  n34 = 2,
                                  n123 = 2,
                                  n124 = 2,
                                  n134 = 2,
                                  n234 = 2,
                                  n1234 = 2,
                                  col = 'black',
                                  fill = c('green3', 'blue3', 'gray92', 'cyan1'),
                                  cex = 0,
                                  category = c('Terrestrial', 'Marine', 'Brackish', 'Freshwater'),
                                  cat.cex = c(1.5, 1.5, 1.5, 1.5),
                                  cat.fontfamily = c('sans', 'sans', 'sans', 'sans'),
                                  cat.fontface = c('bold', 'bold', 'bold', 'bold'))


# The frequency of each environment combination was counted without any overlap
# between the groups.

# As such, I cannot add the labels of the frequency of each environment
# combination through R.

# Instead, I will take a screenshot of the blank venn diagram, and I will
# manually add the data labels to the venn diagram using image editing software.

test <- host_envs %>%
  filter(grepl(",", Environment, fixed = TRUE))

brackish <- host_envs %>%
  filter(grepl("brackish", Environment, fixed = TRUE))

fresh <- host_envs %>%
  filter(grepl("fresh", Environment, fixed = TRUE))

terrestrial <- host_envs %>%
  filter(grepl("terres", Environment, fixed = TRUE))

marine <- host_envs %>%
  filter(grepl("marine", Environment, fixed = TRUE))