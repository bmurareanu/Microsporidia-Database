# -----------------------------------------------------------------------------
#
# Count animal tissues infected by microsporidia
#
# Jason Jiang - Created: 2020/09/02
#                     Last Edit: 2021/02/19
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code for counting the number of times distinct animal tissues are
# infected by unique microsporidia species. This is a reformatted version of an
# older Python script for doing the same thing.
#
# Thanks to Brandon for this script layout
#
# -----------------------------------------------------------------------------


# First, load in useful packages
library(tidyverse)
library(hash) # For creating Python-style dictionaries in R
library(readxl)


# Set working directory
setwd("P:/Shared/Microsporidia database/Figure 3 (Jason)/3A - Frequency of tissues infected by microsporidia")


# Load in spreadsheet containing sorted animal tissue data.
# Each animal tissue infected by microsporidia is sorted into different columns,
# in which each column represents a different tissue type.
sorted_tissues <- read_xlsx("Table S3 Sorted Infected Tissues.xlsx")

systemic <- sorted_tissues$"systemic infection"[!is.na(sorted_tissues$"systemic infection")]


# Create a "dictionary" from sorted_tissues, in which each tissue type is a
# key, and vectors of all the tissues associated with each tissue type are the
# values.
sorted_tissues <- hash(sorted_tissues)

for (key in keys(sorted_tissues)) { # Remove all NA values from tissue vectors
  sorted_tissues[[key]] <- sorted_tissues[[key]][!is.na(sorted_tissues[[key]])]
}

del(c('ambiguous', 'protists'), sorted_tissues) # Exclude ambiguous tissues and protist tissues


# Load in a spreadsheet matching microsporidia species to the tissues they infect,
# as a dataframe
tissues_db <- read_xlsx("Infected Hosts and Tissues.xlsx") %>%
  select(Species, Tissues) %>%
  filter(!is.na(Tissues)) # Exclude species without infected tissue data


# Strip parenthesized text from tissue data entries, as the sorted tissue names
# in sorted_tissues have parenthesized text.
strip_parentheses <- function(Tissues) {
  split_str <- strsplit(Tissues, "; ")[[1]]
  unname(vapply(split_str, function(x) {strsplit(x, " \\(")[[1]][1]}, character(1)))
}

vstrip_parentheses <- Vectorize(strip_parentheses)

tissues_db <- tissues_db %>% mutate(Tissues = vstrip_parentheses(Tissues))


# Create a new dictionary with all the tissue types from sorted_tissues as the
# keys, and 0 as all of the values.
# We will use this dictionary to keep track of how many times each tissue type
# is infected by individual microsporidia species.
counted_tissues <- copy(sorted_tissues)

for (key in keys(counted_tissues)) { # Initialize all values to zero
  counted_tissues[[key]] = 0
}


# Count how many times each tissue type is infected by individual microsporidia
# species.
# Each species may only infect a single tissue type once.
# This is because authors can be inconsistent with how they specific they are
# with microsporidia infection descriptions (ex: just listing gut vs. listing
# intestinal enterocytes, duodenum, ..., for the digestive tract).
for (i in 1:length(tissues_db$Tissues)) {
  species_tissues <- tissues_db$Tissues[[i]]
  local_dict <- hash()
  for (tissue in species_tissues) {
    for (key in keys(sorted_tissues)) {
      if (tissue %in% sorted_tissues[[key]] & !(key %in% keys(local_dict))) {
        counted_tissues[[key]] = counted_tissues[[key]] + 1
        local_dict[[key]] = 1
      }
    }
  }
}


# Manually transfer data of counted_tissues into an Excel spreadsheet, too lazy
# to figure out how to convert a hash object to a dataframe.
counted_tissues