
#
# Comparing hosts infected to tissues infected
#
# Jason Jiang - Created: 2020/07/11
#                     Last Edit: 2021/02/12
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write code to demonstrate the relationship between number of hosts
# infected by microsporidia and number of tissues infected by microsporidia.
#
#
# Thanks to Brandon for this RScript layout
#
# -----------------------------------------------------------------------------

# Load useful packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("readxl")) install.packages("readxl")
if (!require("ggsignif")) install.packages("ggsignif")
library(tidyverse)
library(readxl)
library(ggsignif)


# Set working directory, this is the directory I use
setwd("P:/Shared/Microsporidia database/Figure 3 (Jason)/3D - No. tissues vs no. hosts")


# Load spreadsheet containing tissue + host data into R
host_tissue_data <- read_xlsx("Infected Hosts and Tissues.xlsx")


# Create a vector of terms describing generalized/systemic infections.
# Create a vector of terms describing ambiguous tissue infections.
general_terms <-
  host_tissue_data$General_terms[!is.na(host_tissue_data$General_terms)]

ambiguous_terms <-
  host_tissue_data$Ambiguous_terms[!is.na(host_tissue_data$Ambiguous_terms)]


# Clean up tissue data by separating parenthesized data from entries,
# so we can recognize and remove rows w/ general infections.
# Also, determine whether or not each species infects >1 tissues.
clean_str <- function(string) {
  paste(vapply(strsplit(string, "; ")[[1]],
               function(x) {trimws(strsplit(x, "\\(")[[1]][1])}, character(1)),
        collapse = "; ")
}

vclean_str <- Vectorize(clean_str)


# Species with more than one tissue listed, and/or have systemic infections are
# considered to infect more than one tissue.
check_more_than_one_tissue <- function(Tissues) {
  any(general_terms %in% strsplit(Tissues, "; ")[[1]]) |
    lengths(strsplit(Tissues, '; ')) > 1
}

vcheck_more_than_one_tissue <- Vectorize(check_more_than_one_tissue)


# Also, write a function to check if a species has ambiguous infections listed.
# We will exclude these species, as its unclear what and how many tissues they're
# infecting.
check_ambiguous <- function(Tissues) {
  any(strsplit(Tissues, "; ")[[1]] %in% ambiguous_terms)
}

vcheck_ambiguous <- Vectorize(check_ambiguous)


# Take our dataframe and filter out microsporidia species without both host and
# tissue data.
# Mutate in a new column with the number of infected hosts for each microsporidia.
# Mutate in a new column determining whether or not a species infects more than
# one tissue.
host_tissue_data <- host_tissue_data %>%
  select(-General_terms, -Ambiguous_terms) %>%
  filter(!is.na(Hosts) & !is.na(Tissues)) %>%
  mutate(Tissues = vclean_str(Tissues)) %>%
  filter(!vcheck_ambiguous(Tissues)) %>%
  mutate(num_hosts = lengths(strsplit(Hosts, '; ')),
         More_1_tissue = ifelse(vcheck_more_than_one_tissue(Tissues), ">1 tissues", "1 tissue"))


# Group species into bins of how many hosts they infect.
# Bins: 1 host, 2 hosts, 3 hosts, 4+ hosts
determine_host_bin <- function(num_hosts) {
  if (num_hosts %in% 1:3) { # EDITED
    return(stringr::str_c(num_hosts, " Host"))
  }
  else {
    return("4+ Host") # EDITED
  }
}

vdetermine_host_bin <- Vectorize(determine_host_bin)

host_tissue_data <- host_tissue_data %>%
  mutate(host_bin = vdetermine_host_bin(num_hosts)) %>%
  group_by(host_bin) %>%
  # Fixed error with "`n()` must only be used inside dplyr verbs."
  dplyr::mutate(host_bin_freq = n())

host_bins_counted <- as.matrix(table(host_tissue_data$More_1_tissue,
                                     host_tissue_data$host_bin))


# Chi-square test of independence
chisq.test(host_bins_counted)


# We got a significant p-value from the chi square test (p-value = 0.0003977), so now
# let's run a bunch of pairwise chi square tests for post-hoc testing.
# Set correct = FALSE, we will be applying a Bonferroni correction later.

# I apologize in advance for this tedious monstrosity :/

# Bonferroni corrected p-value: 0.05/6 comparisons
p_bonferroni <- 0.05/6 # 0.008333333

# 1 host to 2 host
p12 <- chisq.test(host_bins_counted[, c(1, 2)], correct = FALSE)$p.value
p12 < p_bonferroni # TRUE
# p = 0.006714288

# 1 host to 3 host
p13 <- chisq.test(host_bins_counted[, c(1, 3)], correct = FALSE)$p.value
p13 < p_bonferroni # FALSE, lack of statistical power?
# p = 0.02979858


# 1 host to 4+ host
p14 <- chisq.test(host_bins_counted[, c(1, 4)], correct = FALSE)$p.value
p14 < p_bonferroni # TRUE
# p = 0.001005971

# 2 host to 3 host
p23 <- chisq.test(host_bins_counted[, c(2, 3)], correct = FALSE)$p.value
p23 < p_bonferroni # FALSE
# p = 0.6659246

# 2 host to 4+ host
p24 <- chisq.test(host_bins_counted[, c(2, 4)], correct = FALSE)$p.value
p24 < p_bonferroni # FALSE
# p = 0.2041254

# 3 host to 4+ host
p34 <- chisq.test(host_bins_counted[, c(3, 4)], correct = FALSE)$p.value
p34 < p_bonferroni # FALSE
# p = 0.5111377


# Turn the table of frequencies of species infecting 1 or more tissues into
# a dataframe, so we can make a bar graph of the proportions
frequency_df <- as.data.frame(host_bins_counted) %>%
  group_by(Var2) %>%
  mutate(host_bin_count = sum(Freq)) %>%
  group_by(host_bin_count) %>%
  mutate(bin_proportion = round(Freq / host_bin_count, 2)) %>%
  filter(Var1 == "1 tissue") %>%
  select(Var2, Freq, host_bin_count, bin_proportion)


# Perfect, let's make that bar graph of proportions!
# NOTE: sometimes, extremely low frequencies get plotted instead of what is
# seen in frequency_df
# I don't know why this happens, but just a heads-up
ggplot(data = frequency_df, aes(x = Var2, y = bin_proportion)) +
  geom_bar(stat = "identity", fill = "black") +
  labs(y = "Fraction Infecting 1 Tissue") +
  theme_bw() +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18),
        axis.title.x = element_blank()) +
  geom_text(aes(label = bin_proportion), nudge_y = 0.02, size = 5) +
  geom_signif(comparisons = list(c("1 Host", "2 Host"), c("1 Host", "4+ Host")),
              annotations = c("p = 5.6e-3", p = "    p = 8.0e-4"), size = 1, tip_length = 0.1,
              y_position = c(0.72, 0.78), textsize = 5)