
#
# Calculate Volume and Polar Filament
#
# Jason Jiang - Created: 2021/01/20
#                     Last Edit: 2021/02/26
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code to estimate volume and polar filament lengths of spores
#
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------

# First. load in useful packages.
library(tidyverse)
library(readxl)


# Set working directory, here's what I use
setwd("C:/Users/Jason Jiang/Desktop/temp work")


# Now, load in the dataframe containing spore measurements we can use for
# estimating volume and polar tube length (see methods).
# This data has been cleaned beforehand, so that spore dimensions properly match
# up to their spore classes.
new_db <- read_xlsx("Volume + PF data.xlsx")


# OK, let's write functions for calculating spore volume
clean <- function(x) {
  # Format spore dimensions for calculation
  if (str_detect(x, "^\\d+?\\.")) {
    return(as.numeric(str_extract(x, "^\\d+\\.\\d+")))
  }
  return(as.numeric(str_extract(x, "^\\d+")))
}

calc_vol <- function(L, W) {
  # Estimate spore volume from length and width
  if(is.na(L) | is.na(W)) {
    return(NA)
  }
  spore_class <- str_extract(L, "\\(.+?\\)") # Retrieve spore class in parentheses
  volume <- (4/3)*pi*((clean(W)/2)^2)*(clean(L)/2)
  return(paste(str_c(volume, if (!is.na(spore_class)) spore_class, sep = " "),
               collapse = "; ")) # Return calculated volume with spore class
}

vcalc_vol <- Vectorize(calc_vol, vectorize.args = c("L", "W"))


# Mutate in a new column for volume
new_db_vols <- new_db %>%
  rowwise() %>%
  mutate(Volume = str_c(vcalc_vol(str_split(Length, "; ")[[1]],
                                  str_split(Width, "; ")[[1]]),
                        collapse = "; "))


# Great, now let's calculate polar filament lengths!
calc_pf <- function(W, Coil) { # Basically the same as calc_vol
  if(!(!is.na(W) & !is.na(Coil))) {
    return(NA)
  }
  spore_class <- str_extract(W, "\\(.+?\\)")
  calc_pf <- clean(W)*clean(Coil)*pi
  return(paste(str_c(calc_pf, if (!is.na(spore_class)) spore_class, sep = " "),
               collapse = "; "))
}

vcalc_pf <- Vectorize(calc_pf, vectorize.args = c("W", "Coil"))


# Awesome, now let's mutate in a new column for calculated polar filament length
new_db_pf_calc <- new_db_vols %>%
  rowwise() %>%
  # Entries with "unclear" written in coil data means that the coils don't
  # clearly match up to a spore class, so we don't calculate PF in these cases.
  mutate(PF_calc = ifelse(str_detect(Coils_for_PF_calc, "unclear"), NA,
    str_c(vcalc_pf(str_split(Width_for_PF_calc, "; ")[[1]],
                                  str_split(Coils_for_PF_calc, "; ")[[1]]),
                         collapse = "; ")))


# Alrighty, we can now save the new_db_pf_calc as a spreadsheet (.xlsx, .csv,
# whatever) and add these volumes + calculated PFs to Table S1!