# -----------------------------------------------------------------------------
#
# Making tissues frequency histogram
#
# Jason Jiang - Created: 2020/07/10
#                     Last Edit: 2021/02/12
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code for making a frequency histogram of the number of tissues
# infected by each microsporidia species.
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------


# Load useful packages
library(tidyverse)
library(readxl)


# Set working directory
setwd("P:/Shared/Microsporidia database/Figure 3 (Jason)/3B - No. tissues infected histogram")


# Open spreadsheet containing tissue data for analysis
tissue_data <- read_xlsx('Infected Hosts and Tissues.xlsx')


# Create a vector containing terms for systemic and ambiguous infections, from
# supplementary table X.
# Species with systemic or ambiguous infections listed infect an unclear number
# of tissues.
infections_to_exclude <-
  tissue_data$"Systemic and Ambiguous Infections"[!is.na(tissue_data$"Systemic and Ambiguous Infections")]


# Write code to clean up tissue_data, removing parentheses from tissue data entries
clean_str <- function(str) {
  paste(sapply(strsplit(str, "; ")[[1]], function(x) {trimws(strsplit(x, "\\(")[[1]][1])}),
        collapse = "; ")
}

check_if_ambiguous_or_systemic_infection <- function(Tissues) {
  any(infections_to_exclude %in% strsplit(Tissues, "; ")[[1]])
}

cleaned_tissues <- tissue_data %>%
  rowwise() %>%
  filter(!is.na(Tissues)) %>% # Exclude species w/out tissue data
  mutate(Tissues = clean_str(Tissues)) %>% # Strip parentheses from tissue data
  filter(!check_if_ambiguous_or_systemic_infection(Tissues)) # Exclude species w/ general infections


# Add new column with number of tissues infected per species, and a new column
# sorting each species into bins of no. tissues infected.
determine_tissue_bin <- function(num_tissues) {
  if (num_tissues == 1) {
    return('1 Tissue')
  }
  else if (num_tissues >= 2 & num_tissues <= 3) {
    return('2 - 3 Tissues')
  }
  else if (num_tissues >= 4 & num_tissues <= 7) {
    return('4 - 7 Tissues')
  }
  else {
    return('>7 Tissues')
  }
}

counted_tissues <- cleaned_tissues %>%
  mutate(num_tissues = lengths(strsplit(Tissues, '; '))) %>%
  rowwise() %>%
  mutate(Infection_bin = determine_tissue_bin(num_tissues))


# Create a table of the tissue bin frequences, and manually move this data
# to Excel to create a bar chart of tissue infection frequencies (I think bar
# graphs look better in Excel)
tissue_bin_freqs <-
  as.data.frame(table(counted_tissues$Infection_bin)) # Manually write as .xlsx file