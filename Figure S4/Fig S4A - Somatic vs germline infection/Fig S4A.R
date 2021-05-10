
#
# Making venn diagram of species infecting soma vs germline tissue
#
# Jason Jiang - Created: 2020/07/17
#                     Last Edit: 2021/05/03

#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code to construct a venn diagram showing the frequency of species
# infecting somatic tissue, germline tissue or both.
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------


# Load useful packages
library(tidyverse)
library(VennDiagram)
library(readxl)


# Set working directory, this is what I use
setwd('P:/Shared/Microsporidia database/Figure S4/Fig S4A - Somatic vs germline infection')


# Create dataframe containing our data of interest
tissue_transmission_data <-
  read_xlsx('Tissues and Transmission Masterlist (NEW).xlsx')


# Create a vector of all the germline tissue names
# We will use this vector to determine if a species infects germline tissue,
# by checking if any of the tissues it infects is in this vector
germline_tissue_names <-
  tissue_transmission_data$Germline_tissues[!is.na(tissue_transmission_data$Germline_tissues)]

systemic_infections <-
  tissue_transmission_data$Systemic_infection[!is.na(tissue_transmission_data$Systemic_infection)]


# Clean up our transmission and tissue columns, by removing parentheses from
# data entries
remove_parenth_str <- function(string) {
  if (!is.na(string)) { # Remove parentheses from all entries in semicolon separated string
    split_string = strsplit(string, '; ')[[1]]
    for (i in 1:length(split_string)) {
      split_string[i] = strsplit(split_string[i], ' \\(')[[1]][1]
    }
    return(paste(split_string, collapse = '; '))
  }
  else {
    return(NA)
  }
}

remove_parenth_str <- Vectorize(remove_parenth_str)


# Filter out rows without tissue AND transmission data from our dataframe, and
# exclude the columns we don't need anymore.
# Format our Tissues and Transmission columns for analysis, using the functions
# above to strip away parenthesized text
tissue_transmission_data <- tissue_transmission_data %>%
  select(-Germline_tissues) %>%
  filter(!(is.na(Tissues) & is.na(Transmission))) %>%
  mutate(Tissues_formatted = remove_parenth_str(Tissues),
         Transmission_formatted = remove_parenth_str(Transmission))


# Also, exclude microsporidia with systemic infections described, as it's unclear
# what host tissues they infect (see methods).
# However, include species with systemic infections if they have transmission
# data described from literature, as that resolves the ambiguity of whether they
# are horizontally or vertically transmitted.
check_for_systemic_infection <- function(Tissues_formatted) {
  any(str_split(Tissues_formatted, "; ")[[1]] %in% systemic_infections)
}

check_for_systemic_infection <- Vectorize(check_for_systemic_infection)

tissue_transmission_data <- tissue_transmission_data %>%
  filter(!(check_for_systemic_infection(Tissues_formatted) & is.na(Transmission_formatted)))


# Write code to determine whether a species infects somatic, germline or both
# I will consider two cases: rows with tissue data, and rows without tissue data
# I will first write helper functions that determines what tissues a species
# infects, one for each case
#
# See methods on "Determining microsporidia transmission mode" for more details
#
# UPDATE, MAY 3 2021: Replaced NA Tissue_formatted values with 'NA' string instead,
# as the NAs were throwing errors in str_detect
get_infections_if_tissues_NOT_empty <- function(Tissues_formatted, Transmission_formatted) {
  if (!(any(strsplit(Tissues_formatted, '; ')[[1]] %in% germline_tissue_names))) { # No germline tissues listed in Tissues
    if (str_detect(Transmission_formatted, 'vertical')) { # Vertical infection listed in Transmission
      return('both')
    }
    else {
      return('horizontal')
    }
  }
  else { # At least one germline tissue listed in Tissues
    if (all(strsplit(Tissues_formatted, '; ')[[1]] %in% germline_tissue_names) & # ONLY germline tissues listed in Tissues
        !(str_detect(Transmission_formatted, 'horizontal'))) { # No horizontal infection listed in Transmission
      return('vertical')
    }
    else {
      return('both')
    }
  }
}

get_infections_if_tissues_empty <- function(Transmission_formatted) {
  if (str_detect(Transmission_formatted, 'vertical') &
      str_detect(Transmission_formatted, 'horizontal')) {
    return('both')
  }
  else if (str_detect(Transmission_formatted, 'vertical')) {
    return('vertical')
  }
  else if (str_detect(Transmission_formatted, 'autoinfection') & !str_detect(Transmission_formatted, 'horizontal')) {
    return('NA') # Autoinfection data w/out horizontal transmission, infection status is uncelar
  }
  else {
    return('horizontal')
  }
}

get_infections_if_tissues_NOT_empty <- Vectorize(get_infections_if_tissues_NOT_empty,
                                                 vectorize.args = c('Tissues_formatted',
                                                                    'Transmission_formatted'))
get_infections_if_tissues_empty <- Vectorize(get_infections_if_tissues_empty)


# Mutate a new column into the dataframe, 'Infection_type', which will tell us
# whether a species infects soma, germline or both.
tissue_transmission_data <- tissue_transmission_data %>%
  mutate(Transmission_formatted = ifelse(is.na(Transmission_formatted),
                                         'NA',
                                         Transmission_formatted)) %>%
  mutate(Infection_type =
           ifelse(!is.na(Tissues_formatted),
                  get_infections_if_tissues_NOT_empty(Tissues_formatted, Transmission_formatted),
                  get_infections_if_tissues_empty(Transmission_formatted))) %>%
  # Exclude species with only autoinfection data or ambiguous literature transmission data
  filter(!(Infection_type == 'NA' | str_detect(Transmission_formatted, 'ambiguous')))


# Create a new dataframe containing just the frequencies of each infection type
# We will use this dataframe to create the venn diagram
infection_freqs <-
  as.data.frame(table(select(tissue_transmission_data, Infection_type)))


# Create constants for venn diagram area sizes, using the infection_freqs
# venn diagram
both_freq = as.numeric(infection_freqs[1, 'Freq'])
somatic_freq = as.numeric(infection_freqs[3, 'Freq']) + as.numeric(infection_freqs[1, 'Freq'])
germline_freq = as.numeric(infection_freqs[2, 'Freq']) + as.numeric(infection_freqs[1, 'Freq'])
# We need to add the frequency of 'both' to 'somatic' and 'germline', to account
# for the overlap between the sets


# Let's now create the venn diagram of infection frequencies!
infection_venn <-
  draw.pairwise.venn(area1 = somatic_freq,
                     area2 = germline_freq,
                     cross.area = both_freq,
                     # fill = c('blue', 'red'),
                     # alpha = 0.8,
                     cex = 2,
                     fontfamily = 'sans',
                     fontface = 'bold',
                     cat.cex = 2,
                     cat.fontfamily = c('sans', 'sans'),
                     cat.fontface = c('bold', 'bold'),
                     cat.just = list(c(1.7, -9), c(-2.5, -20)))


# Perfect, now finish off by saving our tissue_transmission_data dataframe as an
# Excel spreadsheet, for use in later figures requiring transmission data

# openxlsx::write.xlsx(tissue_transmission_data, 'Transmission Database.xlsx')


# As a sidenote, I also replaced all the 'NA' text in the transmission columns
# as actual NAs in Excel, to simplify downstream analysis.