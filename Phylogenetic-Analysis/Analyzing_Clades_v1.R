# -----------------------------------------------------------------------------
#
# Analyzing Clades - v1
#
# Brandon Murareanu - Created: 2020/09/23 
#                     Last Edit: 2021/02/19 
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write a script that can organize the four major microsporidia clades
# from the phylogenetic trees, extract the species' characteristics from the
# master list, and form graphs and charts based on that info for analysis.
#
# -----------------------------------------------------------------------------
#
# As usual, load helpful packages first.

if (! requireNamespace("ape", quietly=TRUE)) {
  install.packages("ape")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("ggtree", quietly = TRUE)) {
  install.packages("ggtree")
}
if (!requireNamespace("treeio", quietly = TRUE)) {
  install.packages("treeio")
}
if (!requireNamespace("tidytree", quietly = TRUE)) {
  install.packages("tidytree")
}
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("DECIPHER", quietly = TRUE)) {
  install.packages("DECIPHER")
}
if (!requireNamespace("seqinr", quietly = TRUE)) {
  install.packages("seqinr")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
}

library(BiocManager)
library(ape)
library(ggtree)
library(treeio)
library(tidytree)
library(ggnewscale)
library(ggplot2)
library(DECIPHER)
library(seqinr)
library(dplyr)
library(ggpubr)

# Set working directory. Modify as need.
# For 270 Species Tree, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_Tree/Source"
anadir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_Tree/Analysis"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_tree"

setwd(anadir)

# Now, we need to read all of the species from the clades into R objects:

cla <- read.csv("./Clade_Species.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
colnames(cla) <- cla[1,]
cla <- cla[-1,]

clade1 <- cla[1]
clade1 <- data.frame(clade1[!apply(is.na(clade1) | clade1 == "", 1, all),])
clade1[2:6] <- ""
colnames(clade1) <- c("species", "hosts", "environments", "spore volume", "polar tube length", "multiple host environments tag")

clade3 <- cla[2]
clade3 <- data.frame(clade3[!apply(is.na(clade3) | clade3 == "", 1, all),])
clade3[2:6] <- ""
colnames(clade3) <- c("species", "hosts", "environments", "spore volume", "polar tube length", "multiple host environments tag")

clade4 <- cla[3]
clade4 <- data.frame(clade4[!apply(is.na(clade4) | clade4 == "", 1, all),])
clade4[2:6] <- ""
colnames(clade4) <- c("species", "hosts", "environments", "spore volume", "polar tube length", "multiple host environments tag")

clade5 <- cla[4]
clade5 <- data.frame(clade5[!apply(is.na(clade5) | clade5 == "", 1, all),])
clade5[2:6] <- ""
colnames(clade5) <- c("species", "hosts", "environments", "spore volume", "polar tube length", "multiple host environments tag")

noclade <- cla[5]
noclade <- data.frame(noclade[!apply(is.na(noclade) | noclade == "", 1, all),])
noclade[2:6] <- ""
colnames(noclade) <- c("species", "hosts", "environments", "spore volume", "polar tube length", "multiple host environments tag")

# --- ASIDE: Comparison of Clades between Bayesian and ML Phylogenies ---

# Now, do this again, but for the ML clades, for comparison purposes:

claml <- read.csv("./Clade_Species_ML.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
colnames(claml) <- claml[1,]
claml <- claml[-1,]

clade1ml <- claml[1]
clade1ml <- data.frame(clade1ml[!apply(is.na(clade1ml) | clade1ml == "", 1, all),])
clade1ml[2] <- "clade 1"
colnames(clade1ml) <- c("species", "clade")

clade3ml <- claml[2]
clade3ml <- data.frame(clade3ml[!apply(is.na(clade3ml) | clade3ml == "", 1, all),])
clade3ml[2] <- "clade 3"
colnames(clade3ml) <- c("species", "clade")

clade4ml <- claml[3]
clade4ml <- data.frame(clade4ml[!apply(is.na(clade4ml) | clade4ml == "", 1, all),])
clade4ml[2] <- "clade 4"
colnames(clade4ml) <- c("species", "clade")

clade5ml <- claml[4]
clade5ml <- data.frame(clade5ml[!apply(is.na(clade5ml) | clade5ml == "", 1, all),])
clade5ml[2] <- "clade 5"
colnames(clade5ml) <- c("species", "clade")

noclademl <- claml[5]
noclademl <- data.frame(noclademl[!apply(is.na(noclademl) | noclademl == "", 1, all),])
noclademl[2] <- "no clade"
noclademl[7,2] <- "outgroup"
colnames(noclademl) <- c("species", "clade")

cladestotalml <- rbind(clade1ml, clade3ml, clade4ml, clade5ml, noclademl)

# Now, compare clades between both cladestotalml and cladestotal. Note that cladestotal is formed
# at the end of this script. Go there first before running this next block of code.

cbacla <- as.vector(cladestotal[,7])
cbaspe <- as.vector(cladestotal[,1])
cmlcla <- as.vector(cladestotalml[,2])
cmlspe <- as.vector(cladestotalml[,1])

comp <- character(length(cbaspe))
for (i in 1:length(cbaspe)) {
  spenum <- which(cmlspe==cbaspe[i])
  bacla <- cbacla[i]
  mlcla <- cmlcla[spenum]
  comp[i] <- as.character(mlcla==bacla)
}

prop.table(table(comp))*100

compf <- cbaspe[which(comp=="FALSE")]

# --- END OF ASIDE: Comparison of Clades between Bayesian and ML Phylogenies ---

# Now, we need to download the Master List, and assign hosts to every species.

MLhosts1 <- read.csv("./Master_List_hosts.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
MLhosts1 <- data.frame(MLhosts1[!apply(is.na(MLhosts1) | MLhosts1 == "", 1, all),])
colnames(MLhosts1) <- MLhosts1[1,]
MLhosts1 <- MLhosts1[-1,]

MLspe <- as.vector(MLhosts1$species)
MLhos <- as.vector(MLhosts1$hosts)

MLspe <- gsub("unnamed ", "", MLspe)
MLspe <- gsub("Unnamed ", "", MLspe)
MLspe <- gsub("\\s*\\([^\\)]+\\)","", MLspe)
MLspe <- gsub("\\.", "", MLspe)       # "." normally refers to "any character", so, we escape it with "//" to
MLspe <- gsub(" ", "_", MLspe)        # search for an actual "."
MLspe <- gsub(";_", ";", MLspe)

MLhos <- gsub("unnamed ", "", MLhos)
MLhos <- gsub("Unnamed ", "", MLhos)
MLhos <- gsub("\\s*\\([^\\)]+\\)","", MLhos)
MLhos <- gsub("\\.", "", MLhos)       # "." normally refers to "any character", so, we escape it with "//" to
MLhos <- gsub(" ", "_", MLhos)        # search for an actual "."
MLhos <- gsub(";_", ";", MLhos)

MLhosts2 <- data.frame(MLspe, MLhos)
colnames(MLhosts2) <- colnames(MLhosts1)

# Now, assign hosts to each species by reffering to the Master List.

# --- Clade 1 ---

for (i in 1:length(as.vector(clade1$species))) {
  xtemp <- as.vector(clade1[i,1])
  num <- which(xtemp == MLspe)
  ytemp <- MLhos[num]
  ztemp <- ytemp[1]
  clade1[i,2] <- ztemp
}

# --- Clade 3 ---

for (i in 1:length(as.vector(clade3$species))) {
  xtemp <- as.vector(clade3[i,1])
  num <- which(xtemp == MLspe)
  ytemp <- MLhos[num]
  ztemp <- ytemp[1]
  clade3[i,2] <- ztemp
}

# --- Clade 4 ---

for (i in 1:length(as.vector(clade4$species))) {
  xtemp <- as.vector(clade4[i,1])
  num <- which(xtemp == MLspe)
  ytemp <- MLhos[num]
  ztemp <- ytemp[1]
  clade4[i,2] <- ztemp
}

# --- Clade 5 ---

for (i in 1:length(as.vector(clade5$species))) {
  xtemp <- as.vector(clade5[i,1])
  num <- which(xtemp == MLspe)
  ytemp <- MLhos[num]
  ztemp <- ytemp[1]
  clade5[i,2] <- ztemp
}

# ---- No Clade ---

for (i in 1:length(as.vector(noclade$species))) {
  xtemp <- as.vector(noclade[i,1])
  num <- which(xtemp == MLspe)
  ytemp <- MLhos[num]
  ztemp <- ytemp[1]
  noclade[i,2] <- ztemp
}

# ------------------------------------------------------------------------------
#
# Part 1a: Environments by Clade 1.0
#
# ------------------------------------------------------------------------------

# Note, make sure to go back through and check for any species that may have different genus names,
# or incorrect spelling. Such species will likely yield "NA" in the hosts column. Correct and annotate 
# these discrepancies in Master_List_hosts.csv, and rerun this section of code.Awesome!

# Now, we need to read the Host Data containing phylum, order and environment into R, and format it
# so that we can easily extract the information we need. 

Hosdat1 <- read.csv("./Host_Data.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
Hosdat1 <- data.frame(Hosdat1[!apply(is.na(Hosdat1) | Hosdat1 == "", 1, all),])
colnames(Hosdat1) <- c("host", "taxonomy", "environment")
Hosdat1 <- Hosdat1[-1,]

HDhos <- as.vector(Hosdat1$host)
HDtax <- as.vector(Hosdat1$taxonomy)
HDenv <- as.vector(Hosdat1$environment)

HDhos <- gsub("unnamed ", "", HDhos)
HDhos <- gsub("Unnamed ", "", HDhos)
HDhos <- gsub("\\s*\\([^\\)]+\\)","", HDhos)
HDhos <- gsub("\\.", "", HDhos)
HDhos <- gsub(" ", "_", HDhos) 
HDhos <- gsub(";_", ";", HDhos)

HDenv <- gsub(",", ";", HDenv)

Hosdat2 <- data.frame(HDhos, HDtax, HDenv)
colnames(Hosdat2) <- colnames(Hosdat1)

# Awesome! Now, we can assign environments to every species based on their host data.

# --- Clade 1 ---

for (i in 1:length(as.vector(clade1$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(clade1[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  clade1[i,3] <- qtemp
  clade1[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

C1env <- as.vector(clade1$environments)
C1env <- unlist(strsplit(C1env, ";"))


  

# --- Clade 3 ---

for (i in 1:length(as.vector(clade3$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(clade3[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  clade3[i,3] <- qtemp
  clade3[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

C3env <- as.vector(clade3$environments)
C3env <- unlist(strsplit(C3env, ";"))

# --- Clade 4 ---

for (i in 1:length(as.vector(clade4$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(clade4[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  clade4[i,3] <- qtemp
  clade4[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

C4env <- as.vector(clade4$environments)
C4env <- unlist(strsplit(C4env, ";"))

# --- Clade 5 ---

for (i in 1:length(as.vector(clade5$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(clade5[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  clade5[i,3] <- qtemp
  clade5[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

C5env <- as.vector(clade5$environments)
C5env <- unlist(strsplit(C5env, ";"))

# --- No Clade ---

for (i in 1:length(as.vector(noclade$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(noclade[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  noclade[i,3] <- qtemp
  noclade[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

nCenv <- as.vector(clade5$environments)
nCenv <- unlist(strsplit(C5env, ";"))

# Now to table the results, and compact them into a single data frame:

# --- Data frame formatted for visualization ---

prop.table(table(C1env))
prop.table(table(C3env))
prop.table(table(C4env))
prop.table(table(C5env))

cladeData <- data.frame(prop.table(table(C1env)),
                        prop.table(table(C3env)),
                        prop.table(table(C4env)),
                        prop.table(table(C5env)))
cladeData <- cladeData[c(2,4,6,8)]
colnames(cladeData) <- c("Clade 1", "Clade 3", "Clade 4", "Clade 5")
rownames(cladeData) <- c("Brackish", "Freshwater", "Marine", "Terrestrial")
cladeData <- (cladeData*100)

# --- Data frame formatted for graphing ---

cladeData2 <- data.frame(1:16,1:16,1:16)
colnames(cladeData2) <- c("Clade", "Environment", "Frequency")
cladeData2[1:4,1] <- colnames(cladeData[1])
cladeData2[1:4,2] <- rownames(cladeData)
cladeData2[1:4,3] <- cladeData[,1]
cladeData2[5:8,1] <- colnames(cladeData[2])
cladeData2[5:8,2] <- rownames(cladeData)
cladeData2[5:8,3] <- cladeData[,2]
cladeData2[9:12,1] <- colnames(cladeData[3])
cladeData2[9:12,2] <- rownames(cladeData)
cladeData2[9:12,3] <- cladeData[,3]
cladeData2[13:16,1] <- colnames(cladeData[4])
cladeData2[13:16,2] <- rownames(cladeData)
cladeData2[13:16,3] <- cladeData[,4]

# --- Graph (Stacked Bar Chart) ---

# --- Barplot 1 ---

ggplot(cladeData2, aes(fill=Environment, y=Frequency, x=Clade)) + 
  geom_bar(position="stack", stat="identity", col="black") +
  scale_fill_manual(breaks=c("Brackish", "Freshwater", "Marine", "Terrestrial"),
                    values=c("Brackish"="#76D7C4",
                             "Freshwater"="#AED6F1",
                             "Marine"="#5DADE2",
                             "Terrestrial"="#1E8449"),
                    name="") +
  theme_classic()

# --- Barplot 2 ---

ggplot(cladeData2, aes(fill=Environment, y=Frequency, x=Clade)) + 
  geom_col(position=position_dodge(0.7), width=0.7, col="black") +
  scale_fill_manual(breaks=c("Brackish", "Freshwater", "Marine", "Terrestrial"),
                    values=c("Brackish"="#76D7C4",
                             "Freshwater"="#AED6F1",
                             "Marine"="#5DADE2",
                             "Terrestrial"="#1E8449"),
                    name="") +
  theme_classic()

# ------------------------------------------------------------------------------
#
# Part 1b: Environments by Clade 2.0
#
# ------------------------------------------------------------------------------

# So, we are taking Aaron's suggestions, and calculating this in a different way.
# Calculate the number of species in each environment / the total number of species
# in the clade. not, as we did before, the frequency of each environment in the clade.

# First part where we collect host data is exactly the same as above:

Hosdat1 <- read.csv("./Host_Data.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
Hosdat1 <- data.frame(Hosdat1[!apply(is.na(Hosdat1) | Hosdat1 == "", 1, all),])
colnames(Hosdat1) <- c("host", "taxonomy", "environment")
Hosdat1 <- Hosdat1[-1,]

HDhos <- as.vector(Hosdat1$host)
HDtax <- as.vector(Hosdat1$taxonomy)
HDenv <- as.vector(Hosdat1$environment)

HDhos <- gsub("unnamed ", "", HDhos)
HDhos <- gsub("Unnamed ", "", HDhos)
HDhos <- gsub("\\s*\\([^\\)]+\\)","", HDhos)
HDhos <- gsub("\\.", "", HDhos)
HDhos <- gsub(" ", "_", HDhos) 
HDhos <- gsub(";_", ";", HDhos)

HDenv <- gsub(",", ";", HDenv)

Hosdat2 <- data.frame(HDhos, HDtax, HDenv)
colnames(Hosdat2) <- colnames(Hosdat1)

# Awesome! Now, we can assign environments to every species based on their host data.
# Except this time, we are not collapsing the species-specific environmental data
# at the end of each clade analysis. Instead, we are sticking with it, and then
# calculating what percentage of species in the clade are found in each of the four
# environments. 

gtemp <- character(0)

# --- Clade 1 ---

for (i in 1:length(as.vector(clade1$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(clade1[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  clade1[i,3] <- qtemp
  clade1[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

C1env <- as.vector(clade1$environments)

C1envbrac <- character(length(C1env))
for (i in 1:length(C1env)) {
  gtemp <- C1env[i]
  envbrac <- "brackish"
  C1envbrac[i] <- grepl(envbrac, gtemp, fixed=TRUE)
}
C1envfres <- character(length(C1env))
for (i in 1:length(C1env)) {
  gtemp <- C1env[i]
  envfres <- "freshwater"
  C1envfres[i] <- grepl(envfres, gtemp, fixed=TRUE)
}
C1envmar <- character(length(C1env))
for (i in 1:length(C1env)) {
  gtemp <- C1env[i]
  envmar <- "marine"
  C1envmar[i] <- grepl(envmar, gtemp, fixed=TRUE)
}
C1envterr <- character(length(C1env))
for (i in 1:length(C1env)) {
  gtemp <- C1env[i]
  envterr <- "terrestrial"
  C1envterr[i] <- grepl(envterr, gtemp, fixed=TRUE)
}

table(C1envbrac)  # visualization.
table(C1envfres)
table(C1envmar)
table(C1envterr)

# --- Clade 3 ---

for (i in 1:length(as.vector(clade3$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(clade3[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  clade3[i,3] <- qtemp
  clade3[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

C3env <- as.vector(clade3$environments)

C3envbrac <- character(length(C3env))
for (i in 1:length(C3env)) {
  gtemp <- C3env[i]
  envbrac <- "brackish"
  C3envbrac[i] <- grepl(envbrac, gtemp, fixed=TRUE)
}
C3envfres <- character(length(C3env))
for (i in 1:length(C3env)) {
  gtemp <- C3env[i]
  envfres <- "freshwater"
  C3envfres[i] <- grepl(envfres, gtemp, fixed=TRUE)
}
C3envmar <- character(length(C3env))
for (i in 1:length(C3env)) {
  gtemp <- C3env[i]
  envmar <- "marine"
  C3envmar[i] <- grepl(envmar, gtemp, fixed=TRUE)
}
C3envterr <- character(length(C3env))
for (i in 1:length(C3env)) {
  gtemp <- C3env[i]
  envterr <- "terrestrial"
  C3envterr[i] <- grepl(envterr, gtemp, fixed=TRUE)
}

table(C3envbrac)  # visualization.
table(C3envfres)
table(C3envmar)
table(C3envterr)

# --- Clade 4 ---

for (i in 1:length(as.vector(clade4$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(clade4[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  clade4[i,3] <- qtemp
  clade4[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

C4env <- as.vector(clade4$environments)

C4envbrac <- character(length(C4env))
for (i in 1:length(C4env)) {
  gtemp <- C4env[i]
  envbrac <- "brackish"
  C4envbrac[i] <- grepl(envbrac, gtemp, fixed=TRUE)
}
C4envfres <- character(length(C4env))
for (i in 1:length(C4env)) {
  gtemp <- C4env[i]
  envfres <- "freshwater"
  C4envfres[i] <- grepl(envfres, gtemp, fixed=TRUE)
}
C4envmar <- character(length(C4env))
for (i in 1:length(C4env)) {
  gtemp <- C4env[i]
  envmar <- "marine"
  C4envmar[i] <- grepl(envmar, gtemp, fixed=TRUE)
}
C4envterr <- character(length(C4env))
for (i in 1:length(C4env)) {
  gtemp <- C4env[i]
  envterr <- "terrestrial"
  C4envterr[i] <- grepl(envterr, gtemp, fixed=TRUE)
}

table(C4envbrac)  # visualization.
table(C4envfres)
table(C4envmar)
table(C4envterr)

# --- Clade 5 ---

for (i in 1:length(as.vector(clade5$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(clade5[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  clade5[i,3] <- qtemp
  clade5[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

C5env <- as.vector(clade5$environments)

C5envbrac <- character(length(C5env))
for (i in 1:length(C5env)) {
  gtemp <- C5env[i]
  envbrac <- "brackish"
  C5envbrac[i] <- grepl(envbrac, gtemp, fixed=TRUE)
}
C5envfres <- character(length(C5env))
for (i in 1:length(C5env)) {
  gtemp <- C5env[i]
  envfres <- "freshwater"
  C5envfres[i] <- grepl(envfres, gtemp, fixed=TRUE)
}
C5envmar <- character(length(C5env))
for (i in 1:length(C5env)) {
  gtemp <- C5env[i]
  envmar <- "marine"
  C5envmar[i] <- grepl(envmar, gtemp, fixed=TRUE)
}
C5envterr <- character(length(C5env))
for (i in 1:length(C5env)) {
  gtemp <- C5env[i]
  envterr <- "terrestrial"
  C5envterr[i] <- grepl(envterr, gtemp, fixed=TRUE)
}

table(C5envbrac)  # visualization.
table(C5envfres)
table(C5envmar)
table(C5envterr)

# --- No Clade ---

for (i in 1:length(as.vector(noclade$species))) {
  
  # 1. Extract environment data based on hosts.
  xtemp <- character(0)
  xtemp <- as.vector(noclade[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDenv[num]
    }
  }
  
  # 2. Tag species that have multiple hosts in different environments, 
  # or multiple hosts in the same environment(s).
  ftemp <- character(0)
  ftemp <- ztemp[ztemp != ""]    # If host has no environmental data, don't include it.
  mhtag <- character(0)
  if ((length(ftemp) > 1)) {
    ftemp <- unique(ftemp)
    if ((length(ftemp) > 1)) {
      mhtag <- "Multiple_Hosts-Different_Environments" 
    } else {
      mhtag <- "Multiple_Hosts-Same_Environments"
    }
  } else {
    mhtag <- ""
  }
  
  # 3. Add data to clade data frame.
  qtemp <- character(0)
  qtemp <- unlist(strsplit(ztemp, ";"))
  qtemp <- unique(qtemp)
  qtemp <- paste(qtemp, collapse=";")
  noclade[i,3] <- qtemp
  noclade[i,6] <- mhtag
  rm(ztemp, qtemp, ytemp, ftemp, mhtag)
}

nCenv <- as.vector(clade5$environments)
nCenv <- unlist(strsplit(C5env, ";"))

# ---

rm(gtemp)

# Now to table the results, and compact them into a single data frame:

(length(which(C1envbrac==TRUE))/length(C1envbrac))*100

cladeData3 <- data.frame(1:16,1:16,1:16)
colnames(cladeData3) <- c("Clade", "Environment", "Frequency")

cladeData3[1:4,1] <- colnames(cladeData[1])
cladeData3[1:4,2] <- rownames(cladeData)
cladeData3[1,3] <- (length(which(C1envbrac==TRUE))/length(C1envbrac))*100
cladeData3[2,3] <- (length(which(C1envfres==TRUE))/length(C1envfres))*100
cladeData3[3,3] <- (length(which(C1envmar==TRUE))/length(C1envmar))*100
cladeData3[4,3] <- (length(which(C1envterr==TRUE))/length(C1envterr))*100

cladeData3[5:8,1] <- colnames(cladeData[2])
cladeData3[5:8,2] <- rownames(cladeData)
cladeData3[5,3] <- (length(which(C3envbrac==TRUE))/length(C3envbrac))*100
cladeData3[6,3] <- (length(which(C3envfres==TRUE))/length(C3envfres))*100
cladeData3[7,3] <- (length(which(C3envmar==TRUE))/length(C3envmar))*100
cladeData3[8,3] <- (length(which(C3envterr==TRUE))/length(C3envterr))*100

cladeData3[9:12,1] <- colnames(cladeData[3])
cladeData3[9:12,2] <- rownames(cladeData)
cladeData3[9,3] <- (length(which(C4envbrac==TRUE))/length(C4envbrac))*100
cladeData3[10,3] <- (length(which(C4envfres==TRUE))/length(C4envfres))*100
cladeData3[11,3] <- (length(which(C4envmar==TRUE))/length(C4envmar))*100
cladeData3[12,3] <- (length(which(C4envterr==TRUE))/length(C4envterr))*100

cladeData3[13:16,1] <- colnames(cladeData[4])
cladeData3[13:16,2] <- rownames(cladeData)
cladeData3[13,3] <- (length(which(C5envbrac==TRUE))/length(C5envbrac))*100
cladeData3[14,3] <- (length(which(C5envfres==TRUE))/length(C5envfres))*100
cladeData3[15,3] <- (length(which(C5envmar==TRUE))/length(C5envmar))*100
cladeData3[16,3] <- (length(which(C5envterr==TRUE))/length(C5envterr))*100

# Now graph!

# --- Barplot ---

ggplot(cladeData3, aes(fill=Environment, y=Frequency, x=Clade)) + 
  geom_col(position=position_dodge(0.7), width=0.7, col="black") +
  ylim(0,80) +
  scale_fill_manual(breaks=c("Brackish", "Freshwater", "Marine", "Terrestrial"),
                    values=c("Brackish"="#76D7C4",
                             "Freshwater"="#AED6F1",
                             "Marine"="#5DADE2",
                             "Terrestrial"="#1E8449"),
                    name="") +
  theme_classic()

# ------------------------------------------------------------------------------
#
# Part 2a: Spore Volume by Clade
#
# ------------------------------------------------------------------------------

# First, read the spore volume spreadsheet into an R object, and format it.

Voldat1 <- read.csv("./Spore_Volume.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
Voldat1 <- data.frame(Voldat1[!apply(is.na(Voldat1) | Voldat1 == "", 1, all),])
colnames(Voldat1) <- c("species", "spore volume")
Voldat1 <- Voldat1[-1,]

VDspe <- as.vector(Voldat1$species)
VDvol <- as.vector(Voldat1$`spore volume`)

VDspe <- gsub("unnamed ", "", VDspe)
VDspe <- gsub("Unnamed ", "", VDspe)
VDspe <- gsub("\\s*\\([^\\)]+\\)","", VDspe)
VDspe <- gsub("\\.", "", VDspe)
VDspe <- gsub(" ", "_", VDspe) 
VDspe <- gsub(";_", ";", VDspe)

Voldat2 <- data.frame(VDspe, VDvol)
colnames(Voldat2) <- colnames(Voldat1)

# Awesome! Now, we can assign spore volumes to every species.

# --- Clade 1 ---

for (i in 1:length(as.vector(clade1$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade1[i,1])
  num <- which(xtemp == VDspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  clade1[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# --- Clade 3 ---

for (i in 1:length(as.vector(clade3$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade3[i,1])
  num <- which(xtemp == VDspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  clade3[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# --- Clade 4 ---

for (i in 1:length(as.vector(clade4$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade4[i,1])
  num <- which(xtemp == VDspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  clade4[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# --- Clade 5 ---

for (i in 1:length(as.vector(clade5$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade5[i,1])
  num <- which(xtemp == VDspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  clade5[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# --- No Clade ---

for (i in 1:length(as.vector(noclade$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(noclade[i,1])
  num <- which(xtemp == VDspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  noclade[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# Okay cool, so now we need to format this data into a data frame suitable for
# making a violin plot with ggplot:

C1vol <- clade1[4]
C1vol <- as.vector(C1vol[,1])
NAtemp <- which(C1vol == "NA")
C1vol <- C1vol[-NAtemp]
C1vol <- unlist(strsplit(C1vol, ";"))
C1 <- character(length(C1vol))
C1[1:length(C1)] <- "Clade 1"

C3vol <- clade3[4]
C3vol <- as.vector(C3vol[,1])
NAtemp <- which(C3vol == "NA")
C3vol <- C3vol[-NAtemp]
C3vol <- unlist(strsplit(C3vol, ";"))
C3 <- character(length(C3vol))
C3[1:length(C3)] <- "Clade 3"

C4vol <- clade4[4]
C4vol <- as.vector(C4vol[,1])
NAtemp <- which(C4vol == "NA")
C4vol <- C4vol[-NAtemp]
C4vol <- unlist(strsplit(C4vol, ";"))
C4 <- character(length(C4vol))
C4[1:length(C4)] <- "Clade 4"

C5vol <- clade5[4]
C5vol <- as.vector(C5vol[,1])
NAtemp <- which(C5vol == "NA")
C5vol <- C5vol[-NAtemp]
C5vol <- unlist(strsplit(C5vol, ";"))
C5 <- character(length(C5vol))
C5[1:length(C5)] <- "Clade 5"

clade <- c(C1, C3, C4, C5)
sporevol <- as.numeric(c(C1vol, C3vol, C4vol, C5vol))

sporevolData <- data.frame(clade, sporevol)
colnames(sporevolData) <- c("Clade", "Spore Volume")

# Now the graphs.

# --- Violin Plot ---

ggplot(sporevolData, aes(fill=Clade, x=Clade, y=`Spore Volume`)) + 
  geom_violin(col=FALSE) +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# --- Boxplot ---

ggplot(sporevolData, aes(fill=Clade, x=Clade, y=`Spore Volume`)) + 
  geom_boxplot() +
  ylim(0,800) +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()
  
# ------------------------------------------------------------------------------
#
# Part 2b: Spore Volume by Clade Updated
#
# ------------------------------------------------------------------------------

# We need to redo these graphs because of the new method of calculating spore volume,
# and because of some corrected errors within the spore volume and polar tube masterlist data.

# First, read the new spore volume spreadsheet into an R object, and format it.

VoldatU1 <- read.csv("./Spore_Volume_v3.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
VoldatU1 <- data.frame(VoldatU1[!apply(is.na(VoldatU1) | VoldatU1 == "", 1, all),])
colnames(VoldatU1) <- c("species", "spore volume")
VoldatU1 <- VoldatU1[-1,]

VDUspe <- as.vector(VoldatU1$species)
VDUvol <- as.vector(VoldatU1$`spore volume`)

VDUspe <- gsub("unnamed ", "", VDUspe)
VDUspe <- gsub("Unnamed ", "", VDUspe)
VDUspe <- gsub("\\s*\\([^\\)]+\\)","", VDUspe)
VDUspe <- gsub("\\.", "", VDUspe)
VDUspe <- gsub(" ", "_", VDUspe) 
VDUspe <- gsub(";_", ";", VDUspe)

VDUvol <- gsub("\\s*\\([^\\)]+\\)","", VDUvol)
VDUvol <- gsub(" ", "_", VDUvol) 
VDUvol <- gsub(";_", ";", VDUvol)

VoldatU2 <- data.frame(VDUspe, VDUvol)
colnames(VoldatU2) <- colnames(VoldatU1)

# Awesome! Now, we can assign spore volumes to every species.

# --- Clade 1 ---

for (i in 1:length(as.vector(clade1$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade1[i,1])
  num <- which(xtemp == VDUspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDUvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  clade1[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# --- Clade 3 ---

for (i in 1:length(as.vector(clade3$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade3[i,1])
  num <- which(xtemp == VDUspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDUvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  clade3[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# --- Clade 4 ---

for (i in 1:length(as.vector(clade4$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade4[i,1])
  num <- which(xtemp == VDUspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDUvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  clade4[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# --- Clade 5 ---

for (i in 1:length(as.vector(clade5$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade5[i,1])
  num <- which(xtemp == VDUspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDUvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  clade5[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# --- No Clade ---

for (i in 1:length(as.vector(noclade$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(noclade[i,1])
  num <- which(xtemp == VDUspe)
  ytemp <- character(0)
  for (j in 1:length(num)) {
    ytemp[j] <- VDUvol[num[j]]
  }
  qtemp <- character(0)
  qtemp <- paste(ytemp, collapse=";")
  noclade[i,4] <- qtemp
  rm(xtemp, num, qtemp, ytemp)
}

# Okay cool, so now we need to format this data into a data frame suitable for
# making a violin plot with ggplot:

C1volU <- clade1[4]
C1volU <- as.vector(C1volU[,1])
emtemp <- which(C1volU == "")
C1volU <- C1volU[-emtemp]
C1volU <- unlist(strsplit(C1volU, ";"))
C1U <- character(length(C1volU))
C1U[1:length(C1U)] <- "Clade 1"

C3volU <- clade3[4]
C3volU <- as.vector(C3volU[,1])
emtemp <- which(C3volU == "")
C3volU <- C3volU[-emtemp]
NAtemp <- which(C3volU == "NA")
C3volU <- C3volU[-NAtemp]
C3volU <- unlist(strsplit(C3volU, ";"))
C3U <- character(length(C3volU))
C3U[1:length(C3U)] <- "Clade 3"

C4volU <- clade4[4]
C4volU <- as.vector(C4volU[,1])
emtemp <- which(C4volU == "")
C4volU <- C4volU[-emtemp]
NAtemp <- which(C4volU == "NA")
C4volU <- C4volU[-NAtemp]
C4volU <- unlist(strsplit(C4volU, ";"))
C4U <- character(length(C4volU))
C4U[1:length(C4U)] <- "Clade 4"

C5volU <- clade5[4]
C5volU <- as.vector(C5volU[,1])
emtemp <- which(C5volU == "")
C5volU <- C5volU[-emtemp]
NAtemp <- which(C5volU == "NA")
C5volU <- C5volU[-NAtemp]
C5volU <- unlist(strsplit(C5volU, ";"))
C5U <- character(length(C5volU))
C5U[1:length(C5U)] <- "Clade 5"

cladeU <- c(C1U, C3U, C4U, C5U)
sporevolU <- as.numeric(c(C1volU, C3volU, C4volU, C5volU))

sporevolDataU <- data.frame(cladeU, sporevolU)
colnames(sporevolDataU) <- c("Clade", "Spore Volume")

# We want to do some statistical test to determine the difference. Thinking either
# unpaired two-sample T-tests, or two-sample U-tests. 

# First, make each dataset numeric:

C1voln <- as.numeric(C1volU)
C3voln <- as.numeric(C3volU)
C4voln <- as.numeric(C4volU)
C5voln <- as.numeric(C5volU)

# Now apply visual inspection and Shaprio-Wilk's test for normality, to see if 
# T-test can be applied.

ggdensity(C1voln)
ggqqplot(C1voln)
ggdensity(C3voln)
ggqqplot(C3voln)
ggdensity(C4voln)
ggqqplot(C4voln)
ggdensity(C5voln)
ggqqplot(C5voln)

shapiro.test(C1voln)
shapiro.test(C3voln)
shapiro.test(C4voln)
shapiro.test(C5voln)

# p-values indicate these values are significantly different from normal. So,
# probably can't use t-test. Use U-test (Mann-Whitney-Wilcox test) instead.

sporvolstats <- data.frame("clade comparison"=c("Clade 1; Clade 3",
                                                "Clade 1; Clade 4",
                                                "Clade 1; Clade 5",
                                                "Clade 3; Clade 4",
                                                "clade 3; Clade 5",
                                                "Clade 4; Clade 5"), 
                           "U-test p-value"=NA, "significant"=NA, stringsAsFactors=FALSE)
colnames(sporvolstats) <- c("clade comparison", "U-test p-value", "significant?")


sporvolstats[1,2] <- wilcox.test(C1voln, C3voln)$p.value*6
sporvolstats[2,2] <- wilcox.test(C1voln, C4voln)$p.value*6
sporvolstats[3,2] <- wilcox.test(C1voln, C5voln)$p.value*6
sporvolstats[4,2] <- wilcox.test(C3voln, C4voln)$p.value*6
sporvolstats[5,2] <- wilcox.test(C3voln, C5voln)$p.value*6
sporvolstats[6,2] <- wilcox.test(C4voln, C5voln)$p.value*6

for (i in 1:length(sporvolstats[,2])) {
  ptemp <- sporvolstats[i,2]
  if (ptemp > (0.05)) {
    sporvolstats[i,3] <- FALSE
  } else {
    sporvolstats[i,3] <- TRUE
  }
}

# Because we are doing multiple comparisons simultaneously, we might need to apply the Bonferroni correction, 
# in which case the p-value cutoff goes from 0.05 to 0.05/6 = 0.0083 = 8.3e-3. Or, we can multiply all the
# p-values by 6 instead, and keep the alpha value at 0.05.

# We can also use the pairwise U-test function to compare all simultaneously with a correction:

x <- pairwise.wilcox.test(sporevolDataU$`Spore Volume`, sporevolDataU$Clade, p.adjust.method="bonferroni")

# Now the graphs.

# --- Violin Plot ---

ggplot(sporevolDataU, aes(fill=Clade, x=Clade, y=`Spore Volume`)) + 
  geom_violin(col=FALSE) +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# --- Boxplot 1 ---

ggplot(sporevolDataU, aes(fill=Clade, x=Clade, y=`Spore Volume`)) + 
  geom_boxplot() +
  ylim(0,250) +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# --- Boxplot 2 ---

ggplot(sporevolDataU, aes(fill=Clade, x=Clade, y=`Spore Volume`)) + 
  geom_boxplot() +
  ylim(0,300) + 
  scale_y_continuous(breaks=round(seq(0,250,by=50))) +
  geom_signif(comparisons = list(c("Clade 1", "Clade 3"),
                                 c("Clade 1", "Clade 4"),
                                 c("Clade 3", "Clade 4"),
                                 c("Clade 3", "Clade 5"),
                                 c("Clade 4", "Clade 5")),
              annotations = c("p = 2.8e-9",
                              "p = 6.3e-18",
                              "p = 0.00023",
                              "p = 0.00312",
                              "p = 1.0e-8"),
              y_position = c(280, 260, 240, 300, 280)) +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# --- Boxplot 3 ---

ggplot(sporevolDataU, aes(fill=Clade, x=Clade, y=`Spore Volume`)) + 
  geom_boxplot() +
  ylim(0,300) + 
  scale_y_continuous(breaks=round(seq(0,250,by=50))) +
  geom_signif(comparisons = list(c("Clade 1", "Clade 3"),
                                 c("Clade 1", "Clade 4"),
                                 c("Clade 3", "Clade 4"),
                                 c("Clade 3", "Clade 5"),
                                 c("Clade 4", "Clade 5")),
              annotation=c("", "*", "*", "*", ""),
              y_position = c(275, 275, 250, 300, 300)) +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# --- Dotplot ---

ggplot(sporevolDataU, aes(fill=Clade, x=Clade, y=`Spore Volume`)) + 
  geom_boxplot() +
  geom_jitter(color="black", size=1) +
  ylim(0,60) +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# ------------------------------------------------------------------------------
#
# Part 3: Polar Tube Length by Clade
#
# ------------------------------------------------------------------------------

# First, read the polar tube length spreadsheet into an R object, and format it.

PTdat1 <- read.csv("./Polar_Tube_Length_v3.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
PTdat1 <- data.frame(PTdat1[!apply(is.na(PTdat1) | PTdat1 == "", 1, all),])
colnames(PTdat1) <- c("species", "all", "calculated", "actual")
PTdat1 <- PTdat1[-1,]

PTspe <- as.vector(PTdat1$species)
PTcal <- as.vector(PTdat1$calculated)
PTact <- as.vector(PTdat1$actual)

PTspe <- gsub("unnamed ", "", PTspe)
PTspe <- gsub("Unnamed ", "", PTspe)
PTspe <- gsub("\\s*\\([^\\)]+\\)","", PTspe)
PTspe <- gsub("\\.", "", PTspe)
PTspe <- gsub(" ", "_", PTspe) 
PTspe <- gsub(";_", ";", PTspe)

PTcal <- gsub("unnamed ", "", PTcal)
PTcal <- gsub("Unnamed ", "", PTcal)
PTcal <- gsub("\\s*\\([^\\)]+\\)","", PTcal)
PTcal <- gsub(" ", "_", PTcal) 
PTcal <- gsub(";_", ";", PTcal)

PTact <- gsub("unnamed ", "", PTact)
PTact <- gsub("Unnamed ", "", PTact)
PTact <- gsub("\\s*\\([^\\)]+\\)","", PTact)
PTact <- gsub(" ", "_", PTact) 
PTact <- gsub(";_", ";", PTact)

PTdat2 <- data.frame(PTspe, PTcal, PTact)
colnames(PTdat2) <- c("species", "calculated", "actual")

# Awesome! Now, we can assign spore volumes to every species.

# --- Clade 1 ---

for (i in 1:length(as.vector(clade1$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade1[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  if (all(acttemp=="")) {
    numx <- which(caltemp!="")
    qtemp <- caltemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  } else {
    numx <- which(acttemp!="")
    qtemp <- acttemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  }
  clade1[i,5] <- qtemp
  rm(xtemp, num, qtemp, numx, acttemp, caltemp)
}

# --- Clade 3 ---

for (i in 1:length(as.vector(clade3$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade3[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  if (all(acttemp=="")) {
    numx <- which(caltemp!="")
    qtemp <- caltemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  } else {
    numx <- which(acttemp!="")
    qtemp <- acttemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  }
  clade3[i,5] <- qtemp
  rm(xtemp, num, qtemp, numx, acttemp, caltemp)
}

# --- Clade 4 ---

for (i in 1:length(as.vector(clade4$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade4[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  if (all(acttemp=="")) {
    numx <- which(caltemp!="")
    qtemp <- caltemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  } else {
    numx <- which(acttemp!="")
    qtemp <- acttemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  }
  clade4[i,5] <- qtemp
  rm(xtemp, num, qtemp, numx, acttemp, caltemp)
}

# --- Clade 5 ---

for (i in 1:length(as.vector(clade5$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade5[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  if (all(acttemp=="")) {
    numx <- which(caltemp!="")
    qtemp <- caltemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  } else {
    numx <- which(acttemp!="")
    qtemp <- acttemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  }
  clade5[i,5] <- qtemp
  rm(xtemp, num, qtemp, numx, acttemp, caltemp)
}

# --- No Clade ---

for (i in 1:length(as.vector(noclade$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(noclade[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  if (all(acttemp=="")) {
    numx <- which(caltemp!="")
    qtemp <- caltemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  } else {
    numx <- which(acttemp!="")
    qtemp <- acttemp[numx]
    qtemp <- paste(qtemp, collapse=";")
  }
  noclade[i,5] <- qtemp
  rm(xtemp, num, qtemp, numx, acttemp, caltemp)
}

# Sweet, now to format the data for graphing:

C1pt <- clade1[5]
C1pt <- as.vector(C1pt[,1])
C1pt <- unlist(strsplit(C1pt, ";"))
C1 <- character(length(C1pt))
C1[1:length(C1)] <- "Clade 1"

C3pt <- clade3[5]
C3pt <- as.vector(C3pt[,1])
C3pt <- unlist(strsplit(C3pt, ";"))
C3 <- character(length(C3pt))
C3[1:length(C3)] <- "Clade 3"

C4pt <- clade4[5]
C4pt <- as.vector(C4pt[,1])
C4pt <- unlist(strsplit(C4pt, ";"))
C4 <- character(length(C4pt))
C4[1:length(C4)] <- "Clade 4"

C5pt <- clade5[5]
C5pt <- as.vector(C5pt[,1])
C5pt <- unlist(strsplit(C5pt, ";"))
C5 <- character(length(C5pt))
C5[1:length(C5)] <- "Clade 5"

clade <- c(C1, C3, C4, C5)
poltub <- as.numeric(c(C1pt, C3pt, C4pt, C5pt))

poltubData <- data.frame(clade, poltub)
colnames(poltubData) <- c("Clade", "Polar Tube Length")

# Now the graphs:

# --- Boxplot ---

ggplot(poltubData, aes(fill=Clade, x=Clade, y=`Polar Tube Length`)) + 
  geom_boxplot() +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# After meeting, what we'll want to do is the same thing, but make separate graphs for both
# "actual" or "calculated" polar filament lengths so we can compare how they look, and decide
# which graph we are going to show. Same code, but slightly modified:

# --- Clade 1 ---

clade1a <- clade1
clade1a[5] <- ""
clade1c <- clade1
clade1c[5] <- ""

for (i in 1:length(as.vector(clade1$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade1[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  rtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  numx <- which(caltemp!="")
  qtemp <- caltemp[numx]
  qtemp <- paste(qtemp, collapse=";")
  
  numy <- which(acttemp!="")
  rtemp <- acttemp[numy]
  rtemp <- paste(rtemp, collapse=";")

  clade1c[i,5] <- qtemp
  clade1a[i,5] <- rtemp
  
  rm(xtemp, num, qtemp, rtemp, numx, numy, acttemp, caltemp)
}

# --- Clade 3 ---

clade3a <- clade3
clade3a[5] <- ""
clade3c <- clade3
clade3c[5] <- ""

for (i in 1:length(as.vector(clade3$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade3[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  rtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  numx <- which(caltemp!="")
  qtemp <- caltemp[numx]
  qtemp <- paste(qtemp, collapse=";")
  
  numy <- which(acttemp!="")
  rtemp <- acttemp[numy]
  rtemp <- paste(rtemp, collapse=";")
  
  clade3c[i,5] <- qtemp
  clade3a[i,5] <- rtemp
  
  rm(xtemp, num, qtemp, rtemp, numx, numy, acttemp, caltemp)
}

# --- Clade 4 ---

clade4a <- clade4
clade4a[5] <- ""
clade4c <- clade4
clade4c[5] <- ""

for (i in 1:length(as.vector(clade4$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade4[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  rtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  numx <- which(caltemp!="")
  qtemp <- caltemp[numx]
  qtemp <- paste(qtemp, collapse=";")
  
  numy <- which(acttemp!="")
  rtemp <- acttemp[numy]
  rtemp <- paste(rtemp, collapse=";")
  
  clade4c[i,5] <- qtemp
  clade4a[i,5] <- rtemp
  
  rm(xtemp, num, qtemp, rtemp, numx, numy, acttemp, caltemp)
}

# --- Clade 5 ---

clade5a <- clade5
clade5a[5] <- ""
clade5c <- clade5
clade5c[5] <- ""

for (i in 1:length(as.vector(clade5$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(clade5[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  rtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  numx <- which(caltemp!="")
  qtemp <- caltemp[numx]
  qtemp <- paste(qtemp, collapse=";")
  
  numy <- which(acttemp!="")
  rtemp <- acttemp[numy]
  rtemp <- paste(rtemp, collapse=";")
  
  clade5c[i,5] <- qtemp
  clade5a[i,5] <- rtemp
  
  rm(xtemp, num, qtemp, rtemp, numx, numy, acttemp, caltemp)
}

# --- No Clade ---

nocladea <- noclade
nocladea[5] <- ""
nocladec <- noclade
nocladec[5] <- ""

for (i in 1:length(as.vector(noclade$species))) {
  xtemp <- character(0)
  xtemp <- as.vector(noclade[i,1])
  num <- which(xtemp==PTspe)
  acttemp <- character(0)
  caltemp <- character(0)
  qtemp <- character(0)
  rtemp <- character(0)
  if (isEmpty(num)) {
    acttemp <- ""
    caltemp <- ""
  } else {
    for (j in 1:length(num)) {
      acttemp[j] <- PTact[num[j]]
      caltemp[j] <- PTcal[num[j]]
    }
  }
  numx <- which(caltemp!="")
  qtemp <- caltemp[numx]
  qtemp <- paste(qtemp, collapse=";")
  
  numy <- which(acttemp!="")
  rtemp <- acttemp[numy]
  rtemp <- paste(rtemp, collapse=";")
  
  nocladec[i,5] <- qtemp
  nocladea[i,5] <- rtemp
  
  rm(xtemp, num, qtemp, rtemp, numx, numy, acttemp, caltemp)
}

# Sweet, now to format the data for graphing:

C1pta <- clade1a[5]
C1pta <- as.vector(C1pta[,1])
C1pta <- unlist(strsplit(C1pta, ";"))
C1a <- character(length(C1pta))
C1a[1:length(C1a)] <- "Clade 1"
C1ptc <- clade1c[5]
C1ptc <- as.vector(C1ptc[,1])
C1ptc <- unlist(strsplit(C1ptc, ";"))
C1c <- character(length(C1ptc))
C1c[1:length(C1c)] <- "Clade 1"

C3pta <- clade3a[5]
C3pta <- as.vector(C3pta[,1])
C3pta <- unlist(strsplit(C3pta, ";"))
C3a <- character(length(C3pta))
C3a[1:length(C3a)] <- "Clade 3"
C3ptc <- clade3c[5]
C3ptc <- as.vector(C3ptc[,1])
C3ptc <- unlist(strsplit(C3ptc, ";"))
C3c <- character(length(C3ptc))
C3c[1:length(C3c)] <- "Clade 3"

C4pta <- clade4a[5]
C4pta <- as.vector(C4pta[,1])
C4pta <- unlist(strsplit(C4pta, ";"))
C4a <- character(length(C4pta))
C4a[1:length(C4a)] <- "Clade 4"
C4ptc <- clade4c[5]
C4ptc <- as.vector(C4ptc[,1])
C4ptc <- unlist(strsplit(C4ptc, ";"))
C4c <- character(length(C4ptc))
C4c[1:length(C4c)] <- "Clade 4"

C5pta <- clade5a[5]
C5pta <- as.vector(C5pta[,1])
C5pta <- unlist(strsplit(C5pta, ";"))
C5a <- character(length(C5pta))
C5a[1:length(C5a)] <- "Clade 5"
C5ptc <- clade5c[5]
C5ptc <- as.vector(C5ptc[,1])
C5ptc <- unlist(strsplit(C5ptc, ";"))
C5c <- character(length(C5ptc))
C5c[1:length(C5c)] <- "Clade 5"

cladea <- c(C1a, C3a, C4a, C5a)
poltuba <- as.numeric(c(C1pta, C3pta, C4pta, C5pta))

cladec <- c(C1c, C3c, C4c, C5c)
poltubc <- as.numeric(c(C1ptc, C3ptc, C4ptc, C5ptc))

poltubDataa <- data.frame(cladea, poltuba)
colnames(poltubDataa) <- c("Clade", "Polar Tube Length")

poltubDatac <- data.frame(cladec, poltubc)
colnames(poltubDatac) <- c("Clade", "Polar Tube Length")

# We want to do some statistical test to determine the difference. Thinking either
# unpaired two-sample T-tests, or two-sample U-tests. 

# First, make each dataset numeric:

C1ptcn <- as.numeric(C1ptc)
C3ptcn <- as.numeric(C3ptc)
C4ptcn <- as.numeric(C4ptc)
C5ptcn <- as.numeric(C5ptc)

# Now apply visual inspection and Shaprio-Wilk's test for normality, to see if 
# T-test can be applied.

ggdensity(C1ptcn)
ggqqplot(C1ptcn)
ggdensity(C3ptcn)
ggqqplot(C3ptcn)
ggdensity(C4ptcn)
ggqqplot(C4ptcn)
ggdensity(C5ptcn)
ggqqplot(C5ptcn)

shapiro.test(C1ptcn)
shapiro.test(C3ptcn)
shapiro.test(C4ptcn)
shapiro.test(C5ptcn)

# p-values indicate these values are significantly different from normal. So,
# probably can't use t-test. Use U-test (Mann-Whitney-Wilcox test) instead.

ptstats <- data.frame("clade comparison"=c("Clade 1; Clade 3",
                                                "Clade 1; Clade 4",
                                                "Clade 1; Clade 5",
                                                "Clade 3; Clade 4",
                                                "Clade 3; Clade 5",
                                                "Clade 4; Clade 5"), 
                           "U-test p-value"=NA, "significant"=NA, stringsAsFactors=FALSE)
colnames(ptstats) <- c("clade comparison", "U-test p-value", "significant?")

ptstats[1,2] <- wilcox.test(C1ptcn, C3ptcn)$p.value*6
ptstats[2,2] <- wilcox.test(C1ptcn, C4ptcn)$p.value*6
ptstats[3,2] <- wilcox.test(C1ptcn, C5ptcn)$p.value*6
ptstats[4,2] <- wilcox.test(C3ptcn, C4ptcn)$p.value*6
ptstats[5,2] <- wilcox.test(C3ptcn, C5ptcn)$p.value*6
ptstats[6,2] <- wilcox.test(C4ptcn, C5ptcn)$p.value*6

for (i in 1:length(ptstats[,2])) {
  ptemp <- ptstats[i,2]
  if (ptemp > (0.05)) {
    ptstats[i,3] <- FALSE
  } else {
    ptstats[i,3] <- TRUE
  }
}

# Because we are doing multiple comparisons simultaneously, we need to apply the Bonferroni correction, 
# in which case the p-value cutoff goes from 0.05 to 0.05/6 = 0.0083 = 8.3e-3. Or, we can multiply all the
# p-values by 6 instead, and keep the alpha value at 0.05.

# --- Boxplot "actual" ---

ggplot(poltubDataa, aes(fill=Clade, x=Clade, y=`Polar Tube Length`)) + 
  geom_boxplot() +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# --- Boxplot "calculated" 1 ---

ggplot(poltubDatac, aes(fill=Clade, x=Clade, y=`Polar Tube Length`)) + 
  geom_boxplot() +
  ylim(0,900) +
  scale_y_continuous(breaks = round(seq(0,800,by=200))) +
  geom_signif(comparisons = list(c("Clade 1", "Clade 4"),
                                 c("Clade 3", "Clade 4"),
                                 c("Clade 3", "Clade 5"),
                                 c("Clade 4", "Clade 5")),
              annotation=c("8.1e-9", 
                           "5.3e-6", 
                           "0.01326", 
                           "8.0e-8"),
              y_position = c(780, 730, 840, 900)) +
  scale_fill_manual(breaks=c("Clade 1", "Clade 3", "Clade 4", "Clade 5"),
                    values=c("Clade 1"="#F73CFA",
                             "Clade 3"="#37D132",
                             "Clade 4"="#35B9EE",
                             "Clade 5"="#F89714"),
                    name="") +
  theme_classic()

# ------------------------------------------------------------------------------
#
# Part 4: Hosts Infected by Species in Multiple Clades
#
# ------------------------------------------------------------------------------

# To do this analysis, it would be easier to bind all the clades data frames together,
# so we can get a complete list of all unique hosts:

clades <- rbind(clade1, clade3, clade4, clade5)
hosts <- clades$hosts
hosts <- unlist(strsplit(hosts, ';'))
hosts <- unique(hosts)
hosts <- hosts[-176]     # Remove "NA"

hosts1 <- clade1$hosts
hosts1 <- unlist(strsplit(hosts1, ';'))
hosts1 <- unique(hosts1)

hosts3 <- clade3$hosts
hosts3 <- unlist(strsplit(hosts3, ';'))
hosts3 <- unique(hosts3)

hosts4 <- clade4$hosts
hosts4 <- unlist(strsplit(hosts4, ';'))
hosts4 <- unique(hosts4)

hosts5 <- clade5$hosts
hosts5 <- unlist(strsplit(hosts5, ';'))
hosts5 <- unique(hosts5)

# Now, we take each unique host, and find out which clades it appears in:

hostnum <- numeric(length(hosts))
for (i in 1:length(hosts)) {
  w <- which(hosts[i]==hosts1)
  x <- which(hosts[i]==hosts3)
  y <- which(hosts[i]==hosts4)
  z <- which(hosts[i]==hosts5)
  num <- c(w,x,y,z)
  hostnum[i] <- length(num)
  rm(num,w,x,y,z)
}

HDcla1 <- data.frame(hosts, hostnum)
colnames(HDcla1) <- c("hosts", "cladefreq")

HDcla2 <- data.frame(table(HDcla1[2]))

# Now the graphs.

# --- Barplot 1 ---

ggplot(data=HDcla2, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity",
           col="black",
           fill="black") +
  labs(subtitle="n = 408",
       x="Number of Clades",
       y="Number of Hosts") +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_classic()

# ------------------------------------------------------------------------------
#
# Part 5: Combine Data
#
# ------------------------------------------------------------------------------

# Combine all clades, even the noclade species:

cladestotal <- rbind(clade1c, clade3c, clade4c, clade5c, nocladec)
NAtemp <- which(is.na(cladestotal[,2]))
cladestotal[NAtemp,2] <- ""


cladestotal[,7] <- ""
cladestotal[,8] <- ""
colnames(cladestotal) <- c("species", 
                           "hosts", 
                           "environments", 
                           "spore volume", 
                           "polar tube length", 
                           "multiple host environments tag",
                           "clade",
                           "accession number")

cladestotal[1:63,7] <- "clade 1"
cladestotal[64:145,7] <- "clade 3"
cladestotal[146:227,7] <- "clade 4"
cladestotal[228:259,7] <- "clade 5"
cladestotal[260:273,7] <- "no clade"
cladestotal[274,7] <- "outgroup"


for (i in 1:length(files1)) {
  name <- as.character(cladestotal[i,1])
  num <- which(files1==name)
  cladestotal[i,8] <- accs[num]
  rm(name, num)
}
