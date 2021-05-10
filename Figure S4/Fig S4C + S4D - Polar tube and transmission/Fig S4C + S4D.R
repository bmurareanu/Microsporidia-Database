# -----------------------------------------------------------------------------
#
# Compare tube lengths to infection
#
# Jason Jiang - Created: 2020/08/14
#                     Last Edit: 2021/05/09
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code to compare tube lengths across infection strategies (i.e: 
# somatic, germline or both infection) for microsporidia species.
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------

# Load in useful packages
library(tidyverse)
library(readxl)


# Set working directory
setwd("P:/Shared/Microsporidia database/Figure S4/Fig S4C + S4D - Polar tube and transmission")


# Load in spreadsheets containing relevant data
transmission_tube_db <- read_xlsx('polar tube and transmission.xlsx')


# Filter out microsporidia species without polar tube data, or without
# transmission data
transmission_tube_db <- transmission_tube_db %>%
  filter(!(is.na(Tube_act) & is.na(Tube_calc))) %>%
  filter(!is.na(Transmission))


# To compare tube length to infection type, we want to select which tube values
# we will use for comparison

# First, write functions to strip parentheses from the tube data.
# This will be used later to format the tube data for analysis.
strip_parentheses <- function(tube_value) {
  strsplit(trimws(strsplit(tube_value, "\\(")[[1]][1]), " ")[[1]][1]
}

strip_parentheses <- Vectorize(strip_parentheses)


# Write functions for formatting the tube data, so all parenthesized text is
# stripped from the tube data.
clean_tube_data <- function(tube_data) {
  raw_values <- strsplit(tube_data, "; ")[[1]]
  paste(strip_parentheses(raw_values), collapse = "; ")
}

clean_tube_data <- Vectorize(clean_tube_data)

# Separate transmission_tube_db into 2 dataframes, one for calculated tube values
# only and the other for experimental tube values
exp <- transmission_tube_db %>%
  filter(!is.na(Tube_act)) %>%
  mutate(Tube_act = clean_tube_data(Tube_act)) %>%
  separate_rows(Tube_act, sep = "; ") %>%
  mutate(Tube_act = as.numeric(Tube_act)) %>%
  select(-Tube_calc)

calc <- transmission_tube_db %>%
  filter(!is.na(Tube_calc)) %>%
  mutate(Tube_calc = clean_tube_data(Tube_calc)) %>%
  separate_rows(Tube_calc, sep = "; ") %>%
  mutate(Tube_calc = as.numeric(Tube_calc)) %>%
  select(-Tube_act)


# Statistical testing for experimental tube values
kruskal.test(data = exp, Tube_act ~ Transmission) # NS
kruskal.test(data = calc, Tube_calc ~ Transmission) # NS


# Make boxplot of all infection groups, for experimental and calculated
ggplot(exp, aes(x = Transmission, y = Tube_act, fill = Transmission)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Transmission", y = "Polar Tube Length (\u03bcm)", title = "Experimental Polar Tubes\np = 0.12") +
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),???
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(calc, aes(x = Transmission, y = Tube_calc, fill = Transmission)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Transmission", y = "Polar Tube Length (\u03bcm)", title = "Calculated Polar Tubes\np = 0.30") +
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())