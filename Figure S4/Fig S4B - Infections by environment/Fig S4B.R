# -----------------------------------------------------------------------------
#
# Analyze infections across environments
#
# Jason Jiang - Created: 2020/07/25
#                     Last Edit: 2021/05/03
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code to analyze patterns of soma and germline infections
# across environments by microsporidia, eventually creating bar graphs
# to represent this data.
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------


# First, I need to modify our microsporidia database so that each row
# corresponds to a single host.
# Afterwards, I'll manually adjust the transmission data for each row,
# so it is reflective of each individual host species.
# Finally, we want to update the host species names to the accepted scientific
# name and add in host environment data, using our host masterlist (Table SX).

# First, load in useful packages.
library(tidyverse)
library(readxl)


# Set working directory, this is what I use.
setwd("P:/Shared/Microsporidia database/Figure S4/Fig S4B - Infections by environment")


# Load in the dataframe created from making Figure S4A
# This dataframe contains data on mode of transmission for all microsporidia
infection_data <- read_xlsx('Transmission Database.xlsx')


# Load in host species list, as this contains up-to-date scientific names of all
# host species, as well as taxonomic + environment data for each host.
host_masterlist <- read_xlsx("Table S2 Host Masterlist.xlsx")


# Write a helper function to help remove parentheses from host entries
remove_parenth_str <- function(x) {
  ifelse(str_sub(x, -1, -1) == ")", trimws(strsplit(x, "\\(")[[1]][1]), x)
}

vremove_parenth_str <- Vectorize(remove_parenth_str)


# Separate each host species into its own row, so each row represents
# a single host species
# Then, mutate new column w/ formatted host data, in which parentheses are stripped
infection_data <- infection_data %>%
  mutate(to_separate = Hosts) %>%
  separate_rows(to_separate, sep = '; ', convert = FALSE) %>%
  mutate(Hosts_formatted = vremove_parenth_str(to_separate)) %>%
  select(Species, Hosts, Hosts_formatted, Tissues, Tissues_formatted, Transmission,
         Transmission_formatted, Infection_type) # Reorganize column ordering


# For microsporidia species infecting multiple hosts, we want to make sure the
# infection status data (somatic, germline, both) is correct for each individual
# host.
# This is because infection status was originally determined for each microsporidia
# species as a whole, using the tissue/transmission data from ALL their hosts.

# First, find all microsporidia species infecting multiple hosts
multiple_hosts <- infection_data %>%
  filter(str_detect(Hosts, ";"))


# Then, save the multiple_hosts dataframe as an Excel file.
# I will manually go over each host entry, and make sure the tissue/transmission
# data is accurate for each individual host.
# Then, I will manually adjust the infection status data of each host, if
# necessary.


# Please refer to the methods for more details about this ("Analysis of Host Tissues").
# openxlsx::write.xlsx(multiple_hosts, "Manually Corrected Species (NEW).xlsx")


# Great, we've manually corrected all the individual host entries, ensuring that
# their infection status data matches up with their tissue infection/transmission
# data,
# Let's first load in our manually corrected host entries.

# 252 host infections excluded
multiple_hosts_fixed <- read_xlsx("Manually Corrected Species (NEW).xlsx")


# Some host species are being excluded from analysis, as we cannot match clear
# tissue infection/transmission data to them, and thus can't reliably assign an
# infection status.
# Host species to be excluded had their Host_formatted cell turned blank, so
# we can use that to find the corresponding host entries in the multiple_hosts
# dataframe.
entries_to_exclude <- multiple_hosts[which(is.na(multiple_hosts_fixed$Hosts_formatted)),]


# Now remove the hosts to exclude from our fixed multiple_hosts dataframe.
multiple_hosts_fixed <- multiple_hosts_fixed %>% filter(!is.na(Hosts_formatted))


# Then, remove the hosts to exclude from our overall infection_data dataframe.
infection_data_good <- anti_join(infection_data, entries_to_exclude)


# Now, add our manually fixed host data (multiple_hosts_fixed) to the
# infection_data_good dataframe.
# Then, remove the original unfixed host entries from this dataframe.
infection_data_good <- rbind(multiple_hosts_fixed, infection_data_good)

infection_data_good <-
  infection_data_good[!duplicated(infection_data_good[c("Species", "Hosts_formatted")]),]


# Cool, now let's create column of host taxonomic data.
# We'll also change all host names to its most up-to-date-name from the host masterlist
# Also, let's add in a new column of the environment of each host, using the
# host masterlist.
infection_data_good$Host_updated <- "" # Create new column for updated host names
infection_data_good$Host_taxonomy <- "" # Create new column for host taxo. data
infection_data_good$Host_env <- "" # Create new column for host env data

for(i in 1:nrow(infection_data_good)) {
  host <- infection_data_good$Hosts_formatted[i]
  host_row <- filter(host_masterlist, grepl(host, Host_original, fixed = TRUE))
  infection_data_good[i, "Host_updated"] <- host_row$Host_formatted[1]
  infection_data_good[i, "Host_taxonomy"] <- host_row$Taxonomy[1]
  infection_data_good[i, "Host_env"] <- host_row$Environment[1]
}


################################################################################


# With all our data cleanup/prep done, let's get to making the stacked bar
# charts of transmission mode frequencies across environments.


# Exclude any host species w/out environment data, not useful for our analysis.
# We are also only interested in studying purely terrestrial, freshwater
# and marine species (no pure brackish species were available after data cleanup).
envs_of_interest <- c('terrestrial', 'marine', 'freshwater')


# Filter down to hosts in environments of interest, and not a non-species entry
infection_data_good <- infection_data_good %>%
  filter(!is.na(Host_env) & Host_env %in% envs_of_interest)


# Create bar graphs of soma and germline infection frequencies across
# terrestrial, freshwater and marine environments
ggplot(data = infection_data_good) +
  geom_bar(mapping = aes(x = Host_env, fill = Infection_type),
           position = 'fill') +
  ylab('Fraction of Hosts') +
  guides(fill = guide_legend(title = 'Infection type')) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        legend.title=element_text(size = 18),
        legend.text=element_text(size = 14),
        panel.background = element_rect(fill = 'white', colour = 'black'))