# -----------------------------------------------------------------------------
#
# Making species frequency histogram
#
# Jason Jiang - Created: 2020/07/10
#                     Last Edit: 2020/12/25
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code for making a frequency histogram of the number of hosts
# infected by each microsporidia species.
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------


# Load useful packages
library(tidyverse)
library(readxl)


# Set working directory
setwd('P:/Shared/Microsporidia database/Figure 3 (Jason)/3C - No. hosts infected histogram')


# Open spreadsheet containing host data for analysis
host_data = read_xlsx("Infected Hosts and Tissues.xlsx")


# Create new column in spreadsheet w/ no. hosts infected by each microsporidium,
# then create another column grouping species by no. hosts infected
determine_host_bin <- function(num_hosts) {
  if (num_hosts == 1) {
    return('1 Host')
  }
  else if (num_hosts >= 2 & num_hosts <= 4) {
    return('2 - 4 Hosts')
  }
  else if (num_hosts >= 5 & num_hosts <= 10)
    return('5 - 10 Hosts')
  else {
    return('>10 Hosts')
  }
}

hosts_counted <- host_data %>%
  filter(!is.na(Hosts)) %>% # exclude rows w/out host data
  rowwise() %>%
  mutate(num_hosts =
           lengths(strsplit(Hosts, '; '))) %>%
  mutate(host_bin =
           determine_host_bin(num_hosts))


# Manually transfer bin frequency data to Excel, to create a frequency bar
# chart (I think the graph looks better in Excel)
view(as.data.frame(table(hosts_counted$host_bin)))
