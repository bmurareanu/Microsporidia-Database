# -----------------------------------------------------------------------------
#
# Creating Pairwise Sequence Identity Matrices - v1
#
# Brandon Murareanu - Created: 2020/08/04
#                     Last Edit: 2021/02/19 
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write some code that can create a pairwise sequence identity matrix from
# an alignment object.
#
# -----------------------------------------------------------------------------
#
# This analysis builds off data compiled in the Analyzing_Clades_v1.R script. Run through
# that script before running through this one. 
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
if (!requireNamespace("bio3d", quietly = TRUE)) {
  install.packages("bio3d")
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
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

library(BiocManager)
library(ape)
library(ggtree)
library(treeio)
library(tidytree)
library(ggnewscale)
library(ggplot2)
library(DECIPHER)
library(bio3d)
library(seqinr)
library(dplyr)
library(ggpubr)
library(writexl)

# Set working directory. Modify as need. For me, it's this (Big Tree 1 directory):

sourcedir <- "C:/Users/bmura/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_1/Source"
maindir <- "C:/Users/bmura/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_1"

# For 270 Species Tree, it's this:

sourcedir <- "C:/Users/bmura/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_Tree/Source"
maindir <- "C:/Users/bmura/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_tree"

setwd(maindir)

# ------------------------------------------------------------------------------
#
# Part 1: Making Pairwise Identity Matrices
#
# ------------------------------------------------------------------------------

# Read our alignment files in FASTA format into R (we will do both masked and unmaksed forms to compare):

MSA13 <- read.fasta("./MSA13.0_MAFFT_(FASTA).mafft")
MSA13masked <- read.fasta("./MSA13.0_MAFFT_(FASTA)_masked.mafft")

MSA13 <- read.alignment("./MSA13.0_MAFFT_(FASTA).mafft", format="fasta")
MSA13masked <- read.alignment("./MSA13.0_MAFFT_(FASTA)_masked.mafft", format="fasta")

# For ease-of-viewing purposes, we want to change the accession numbers into the species names. Luckily, all 
# species names are stored as the file names in the "Source" folder. To retirve those names into a vector:

setwd(sourcedir)
files1 <- list.files()
setwd(maindir)

# Awesome. Now we substitute the accession numbers for the species names in our alignment object:

MSA13[["id"]] <- files1
MSA13masked[["id"]] <- files1

# Now we compute a pairwise identity matrix of sequence similarity using the "seqidentity" function from the 
# bio3d package. Edit: NVM, use the dist.alignment function instead.

#pim13 <- seqidentity(MSA13)
#pim13masked <- seqidentity(MSA13masked)

pim13 <- dist.alignment(MSA13, matrix="identity", gap=1)
pim13 <- as.matrix(pim13)
pim13 <- (1-pim13^2)
colnames(pim13) <- files1
rownames(pim13) <- files1

pim13masked <- dist.alignment(MSA13masked, matrix="identity", gap=1)
pim13masked <- as.matrix(pim13masked)
pim13masked <- (1-pim13masked^2)
colnames(pim13masked) <- files1
rownames(pim13masked) <- files1

# Cool. Now save as an R object.

save(pim13, file="MSA13_pim.RData")
save(pim13masked, file="MSA13_masked_pim.RData")

# We will conduct the rest of the analysis using only the unmasked PIM. Can
# replicate the code for the masked PIM if need be later.

# ------------------------------------------------------------------------------
#
# Part 2: Extracting Pairwise Sequence Identity Data
#
# ------------------------------------------------------------------------------

# First, create data frame to store pairwise data:

pim13data <- data.frame("Species 1", "Species 2", "Similarity", stringsAsFactors = FALSE)
colnames(pim13data) <- pim13data[1,]
pim13data[1,] <- NA

# Awesome. Now, write a loop that parses through the PIM to extract sequence similarity for
# every possible unique species pair.

pim13rows <- as.vector(rownames(pim13))
pim13cols <- as.vector(colnames(pim13))

counter1 <- 1
counter2 <- 1
for (i in 1:length(pim13rows)){
  counter1 <- counter1+1
  for (j in counter1:length(pim13cols)) {
    pim13data[counter2,3] <- pim13[i,j]
    pim13data[counter2,1] <- pim13rows[i]
    pim13data[counter2,2] <- pim13cols[j]
    counter2 <- counter2+1
  }
}
rm(counter1, counter2)

# Amazing! Works just as I intended.

# Now, modify the PIM data frame to add additional columns for analysis, and remove all rows containing
# the outgroup Rozella allomycis:

pim13data <- data.frame(pim13data, "Shared Phylum?", "Shared Class?", "Shared Order?","Shared Family?", "Shared Environment?", 
                        "% Spore Volume Similarity", "% Polar Tube Length Similarity",
                        stringsAsFactors = FALSE)
pim13data[,4:10] <- NA
colnames(pim13data) <- c("Species 1", "Species 2", "Similarity",
                         "Shared Phylum?", "Shared Class?", "Shared Order?", "Shared Family?", "Shared Environment?",
                         "% Spore Volume Similarity", "% Polar Tube Length Similarity")

rozrows1 <- which(pim13data[,1]=="Rozella_allomycis")
rozrows2 <- which(pim13data[,2]=="Rozella_allomycis")
rozrows <- c(rozrows1, rozrows2)
pim13data <- pim13data[-rozrows,]
rm(rozrows, rozrows1, rozrows2)

# Remove all rows containing Pleistophora sp. 1, 2 or 3, according to Aaron's recommendations.

pleirows1 <- which(pim13data[,1]=="Pleistophora_sp_1" | pim13data[,1]=="Pleistophora_sp_2" | pim13data[,1]=="Pleistophora_sp_3")
pleirows2 <- which(pim13data[,2]=="Pleistophora_sp_1" | pim13data[,2]=="Pleistophora_sp_2" | pim13data[,2]=="Pleistophora_sp_3")
pleirows <- c(pleirows1, pleirows2)
pim13data <- pim13data[-pleirows,]
rm(pleirows1, pleirows, pleirows2)

# ------------------------------------------------------------------------------
#
# Part 3: Adding Host Taxonomy Data for Analysis
#
# ------------------------------------------------------------------------------

# Building off the "cladestotal" data frame from the "Analyzing_Clades_v1.R" script,
# add host phyla and order to two new columns.

cladestotal2 <- data.frame(cladestotal, "host phyla", "host classes", "host orders", "host families", stringsAsFactors = FALSE)
cladestotal2[,9:12] <- NA
colnames(cladestotal2) <- c("species", "hosts",
                            "environments", "spore volume",
                            "polar tube length", "multiple host environments tag",
                            "clade", "accession number", "host phyla", "host classes", "host orders", "host families")

# Using the "Hosdat2" data frame, extract taxonomy data for a desired host, and parse out
# the phylum, class, order and family for that host. Repeat for all hosts of a given species, then
# repeat for all species:

for (i in 1:length(files1)){
  # Select a microsporidia species:
  xtemp <- character(0)
  xtemp <- as.vector(cladestotal2[i,2])
  ytemp <- character(0)
  ytemp <- unlist(strsplit(xtemp, ";"))
  ztemp <- character(0)
  
  # Extract taxonomy data for each host for that species:
  for (j in 1:length(ytemp)) {
    num <- which(ytemp[j] == HDhos)
    if (isEmpty(num)) {
      ztemp[j] <- ""
    } else {
      ztemp[j] <- HDtax[num]
    }
  }
  
  # Parse phyla, classes, orders and families from host taxonomy data:
  phytemp <- character(0)
  clstemp <- character(0)
  ordtemp <- character(0)
  famtemp <- character(0)
  for (k in 1:length(ztemp)) {
    qtemp <- character(0)
    qtemp <- unlist(strsplit(ztemp[k], ";"))
    phytemp[k] <- qtemp[2]
    clstemp[k] <- qtemp[3]
    ordtemp[k] <- qtemp[4]
    famtemp[k] <- qtemp[5]
  }
  
  phytemp <- unique(phytemp)
  phytemp <- phytemp[!is.na(phytemp)]
  phytemp <- paste(phytemp, collapse=";")
  
  clstemp <- unique(clstemp)
  clstemp <- clstemp[!is.na(clstemp)]
  clstemp <- paste(clstemp, collapse=";")
  
  ordtemp <- unique(ordtemp)
  ordtemp <- ordtemp[!is.na(ordtemp)]
  ordtemp <- paste(ordtemp, collapse=";")
  
  famtemp <- unique(famtemp)
  famtemp <- famtemp[!is.na(famtemp)]
  famtemp <- paste(famtemp, collapse=";")
  
  # Add data to "cladestotal2" dataframe:
  cladestotal2[i,9] <- phytemp
  cladestotal2[i,10] <- clstemp
  cladestotal2[i,11] <- ordtemp
  cladestotal2[i,12] <- famtemp
}
rm(ztemp, ytemp, xtemp, phytemp, ordtemp, clstemp, famtemp, qtemp, num)

# Remove species with missing taxonomy data!!!

cladestotal3 <- cladestotal2
empt <- which(cladestotal3[,9]=="")
emptspe <- as.vector(cladestotal3[empt,1])
emptspe <- emptspe[-c(2,3,4,9)]

pim13data2 <- pim13data

for (i in 1:length(emptspe)) {
  emptrows1 <- which(pim13data2[,1]==emptspe[i])
  emptrows2 <- which(pim13data2[,2]==emptspe[i])
  emptrows <- c(emptrows1, emptrows2)
  pim13data2 <- pim13data2[-emptrows,]
  rm(emptrows, emptrows1, emptrows2)
}

# ------------------------------------------------------------------------------
#
# Part 4: Comparing Host Taxonomy within Species Pairs (Figure 6B)
#
# ------------------------------------------------------------------------------

# So now that we have all the host taxonomy data, we need some code that
# can determine whether or not the two species in each pair share at least one 
# unique host phylum, class, order or family.

# --- Phyla ---

#pim13data2
for(i in 1:length(pim13data2[,1])) {
  spe1 <- pim13data2[i,1]
  spe2 <- pim13data2[i,2]
  num1 <- which(cladestotal3[,1]==spe1)
  num2 <- which(cladestotal3[,1]==spe2)
  
  phy1 <- cladestotal3[num1,9]
  phy1 <- unlist(strsplit(phy1, ";"))
  phy2 <- cladestotal3[num2,9]
  phy2 <- unlist(strsplit(phy2, ";"))
  
  xtemp <- character(0)
  if (isEmpty(phy1)==TRUE) {
    phy1 <- ""
  }
  for (j in 1:length(phy1)) {
    xtemp[j] <- any(phy2==phy1[j])
  }
  if (any(xtemp==TRUE)==TRUE) {
    pim13data2[i,4] <- TRUE
  } else {
    pim13data2[i,4] <- FALSE
  }
 
  rm(spe1, spe2, num1, num2, phy1, phy2, xtemp)
}

# Little calculation for paper... calculating proportion of species pairs with more or less than 80%
# similarity that infect an arthropod or a chordate.

o80 <- which(pim13data2[,3]>0.80)
pim13data2o80 <- pim13data2[o80,]

for(i in 1:length(pim13data2o80[,1])) {
  spe1 <- pim13data2o80[i,1]
  spe2 <- pim13data2o80[i,2]
  num1 <- which(cladestotal3[,1]==spe1)
  num2 <- which(cladestotal3[,1]==spe2)
  
  phy1 <- cladestotal3[num1,9]
  phy1 <- unlist(strsplit(phy1, ";"))
  phy2 <- cladestotal3[num2,9]
  phy2 <- unlist(strsplit(phy2, ";"))
  
  if (isEmpty(phy1)==TRUE) {
    phy1 <- ""
  }
  arth <- any(c(phy1, phy2)=="Arthropoda")
  chor <- any(c(phy1, phy2)=="Chordata")
  if (arth==TRUE) {
    pim13data2o80[i,4] <- "Arthropoda"
  } else {
    pim13data2o80[i,4] <- "NO"
  }
  if (chor==TRUE) {
    pim13data2o80[i,5] <- "Chordata"
  } else {
    pim13data2o80[i,5] <- "NO"
  }
  if (chor==TRUE & arth==TRUE) {
    pim13data2o80[i,6] <- "Both"
  } else {
    pim13data2o80[i,6] <- "NO"
  }
 
  rm(spe1, spe2, num1, num2, phy1, phy2, arth, chor)
}

o80arth <- as.vector(pim13data2o80[,4])
o80chor <- as.vector(pim13data2o80[,5])
o80both <- as.vector(pim13data2o80[,6])

u80 <- which(pim13data2[,3]<0.80)
pim13data2u80 <- pim13data2[u80,]

for(i in 1:length(pim13data2u80[,1])) {
  spe1 <- pim13data2u80[i,1]
  spe2 <- pim13data2u80[i,2]
  num1 <- which(cladestotal3[,1]==spe1)
  num2 <- which(cladestotal3[,1]==spe2)
  
  phy1 <- cladestotal3[num1,9]
  phy1 <- unlist(strsplit(phy1, ";"))
  phy2 <- cladestotal3[num2,9]
  phy2 <- unlist(strsplit(phy2, ";"))
  
  if (isEmpty(phy1)==TRUE) {
    phy1 <- ""
  }
  arth <- any(c(phy1, phy2)=="Arthropoda")
  chor <- any(c(phy1, phy2)=="Chordata")
  if (arth==TRUE) {
    pim13data2u80[i,4] <- "Arthropoda"
  } else {
    pim13data2u80[i,4] <- "NO"
  }
  if (chor==TRUE) {
    pim13data2u80[i,5] <- "Chordata"
  } else {
    pim13data2u80[i,5] <- "NO"
  }
  if (chor==TRUE & arth==TRUE) {
    pim13data2u80[i,6] <- "Both"
  } else {
    pim13data2u80[i,6] <- "NO"
  }
  
  rm(spe1, spe2, num1, num2, phy1, phy2, arth, chor)
}

u80arth <- as.vector(pim13data2u80[,4])
u80chor <- as.vector(pim13data2u80[,5])
u80both <- as.vector(pim13data2u80[,6])

# Now for the actual calculations:

prop.table(table(o80arth))
prop.table(table(u80arth))
prop.table(table(o80chor))
prop.table(table(u80chor))
prop.table(table(o80both))
prop.table(table(u80both))

# Calculation done!

# --- Classes ---

#pim13data2
for(i in 1:length(pim13data2[,1])) {
  spe1 <- pim13data2[i,1]
  spe2 <- pim13data2[i,2]
  num1 <- which(cladestotal3[,1]==spe1)
  num2 <- which(cladestotal3[,1]==spe2)
  
  cls1 <- cladestotal3[num1,10]
  cls1 <- unlist(strsplit(cls1, ";"))
  cls2 <- cladestotal3[num2,10]
  cls2 <- unlist(strsplit(cls2, ";"))
  
  xtemp <- character(0)
  if (isEmpty(cls1)==TRUE) {
    cls1 <- ""
  }
  for (j in 1:length(cls1)) {
    xtemp[j] <- any(cls2==cls1[j])
  }
  if (any(xtemp==TRUE)==TRUE) {
    pim13data2[i,5] <- TRUE
  } else {
    pim13data2[i,5] <- FALSE
  }
  
  rm(spe1, spe2, num1, num2, cls1, cls2, xtemp)
}

# --- Orders ---

#pim13data2
for(i in 1:length(pim13data2[,1])) {
  spe1 <- pim13data2[i,1]
  spe2 <- pim13data2[i,2]
  num1 <- which(cladestotal3[,1]==spe1)
  num2 <- which(cladestotal3[,1]==spe2)
  
  ord1 <- cladestotal3[num1,11]
  ord1 <- unlist(strsplit(ord1, ";"))
  ord2 <- cladestotal3[num2,11]
  ord2 <- unlist(strsplit(ord2, ";"))
  
  xtemp <- character(0)
  if (isEmpty(ord1)==TRUE) {
    ord1 <- ""
  }
  for (j in 1:length(ord1)) {
    xtemp[j] <- any(ord2==ord1[j])
  }
  if (any(xtemp==TRUE)==TRUE) {
    pim13data2[i,6] <- TRUE
  } else {
    pim13data2[i,6] <- FALSE
  }
  rm(spe1, spe2, num1, num2, ord1, ord2, xtemp)
}

# --- Families ---

#pim13data2
for(i in 1:length(pim13data2[,1])) {
  spe1 <- pim13data2[i,1]
  spe2 <- pim13data2[i,2]
  num1 <- which(cladestotal3[,1]==spe1)
  num2 <- which(cladestotal3[,1]==spe2)
  
  fam1 <- cladestotal3[num1,12]
  fam1 <- unlist(strsplit(fam1, ";"))
  fam2 <- cladestotal3[num2,12]
  fam2 <- unlist(strsplit(fam2, ";"))
  
  xtemp <- character(0)
  if (isEmpty(fam1)==TRUE) {
    fam1 <- ""
  }
  for (j in 1:length(fam1)) {
    xtemp[j] <- any(fam2==fam1[j])
  }
  if (any(xtemp==TRUE)==TRUE) {
    pim13data2[i,7] <- TRUE
  } else {
    pim13data2[i,7] <- FALSE
  }
  
  rm(spe1, spe2, num1, num2, fam1, fam2, xtemp)
}

# Awesome, now to graph the results in a histogram:

# First, make a distribution of the sequence similarity with varying bin sizes:

x <- as.numeric(pim13data2[,3])
plot(hist(x, breaks=20)) # bin size 5
plot(hist(x, breaks=25)) # bin size 4
plot(hist(x, breaks=50)) # bin size 2
plot(hist(x, breaks=100)) # bin size 1

ggplot(pim13data2, aes(as.numeric(pim13data2$Similarity))) + 
  geom_histogram(col="black", fill="black", binwidth=0.05) +
  ylim(0,6000) +
  theme_classic()
ggplot(pim13data2, aes(as.numeric(pim13data2$Similarity))) + 
  geom_histogram(col="black", fill="black", binwidth=0.04) +
  theme_classic()
ggplot(pim13data2, aes(as.numeric(pim13data2$Similarity))) + 
  geom_histogram(col="black", fill="black", binwidth=0.02) +
  theme_classic()
ggplot(pim13data2, aes(as.numeric(pim13data2$Similarity))) + 
  geom_histogram(col="black", fill="black", binwidth=0.01) +
  theme_classic()

# --- Phyla ---

pim13data2phy <- data.frame(as.numeric(pim13data2[,3]))
colnames(pim13data2phy) <- c("Similarity")

phytrus <- which(pim13data2[,4]==TRUE)
pim13data2phytrus <- data.frame(as.numeric(pim13data2[phytrus,3]))
colnames(pim13data2phytrus) <- c("Similarity")

px <- as.vector(pim13data2phy[,1])
py <- as.vector(pim13data2phytrus[,1])

pxx <- hist(px, breaks=21)
pyy <- hist(py, breaks=21)

pzz0 <- pyy
pzz0$counts <- (pyy$counts/pxx$counts)
pzz0$counts[which(pzz0$counts=="NaN")] <- 0

plot(pzz0, col="red")

pzz1 <- data.frame(pzz0$breaks[2:21], pzz0$counts)
colnames(pzz1) <- c("Sequence Pair Similarity", "Frequency of Shared Phyla")

ggplot(pzz1, aes(y=`Frequency of Shared Phyla`, x=`Sequence Pair Similarity`)) + 
  geom_col(col="red", fill="red", size=1.5) +
  ylim(0,1) +
  theme_classic()

# --- Classes ---

pim13data2cls <- data.frame(as.numeric(pim13data2[,3]))
colnames(pim13data2cls) <- c("Similarity")

clstrus <- which(pim13data2[,5]==TRUE)
pim13data2clstrus <- data.frame(as.numeric(pim13data2[clstrus,3]))
colnames(pim13data2clstrus) <- c("Similarity")

cx <- as.vector(pim13data2cls[,1])
cy <- as.vector(pim13data2clstrus[,1])

cxx <- hist(cx, breaks=21)
cyy <- hist(cy, breaks=21)

czz0 <- cyy
czz0$counts <- (cyy$counts/cxx$counts)
czz0$counts[which(czz0$counts=="NaN")] <- 0

plot(czz0, col="purple")

czz1 <- data.frame(czz0$breaks[2:21], czz0$counts)
colnames(czz1) <- c("Sequence Pair Similarity", "Frequency of Shared Classes")

ggplot(czz1, aes(y=`Frequency of Shared Classes`, x=`Sequence Pair Similarity`)) + 
  geom_col(col="purple", fill="purple", size=1.5) +
  ylim(0,1) +
  theme_classic()

# --- Orders ---

pim13data2ord <- data.frame(as.numeric(pim13data2[,3]))
colnames(pim13data2ord) <- c("Similarity")

ordtrus <- which(pim13data2[,6]==TRUE)
pim13data2ordtrus <- data.frame(as.numeric(pim13data2[ordtrus,3]))
colnames(pim13data2ordtrus) <- c("Similarity")

ox <- as.vector(pim13data2ord[,1])
oy <- as.vector(pim13data2ordtrus[,1])

oxx <- hist(ox, breaks=21)
oyy <- hist(oy, breaks=21)

ozz0 <- oyy
ozz0$counts <- (oyy$counts/oxx$counts)
ozz0$counts[which(ozz0$counts=="NaN")] <- 0

plot(ozz0, col="blue")

ozz1 <- data.frame(ozz0$breaks[2:21], ozz0$counts)
colnames(ozz1) <- c("Sequence Pair Similarity", "Frequency of Shared Orders")

ggplot(ozz1, aes(y=`Frequency of Shared Orders`, x=`Sequence Pair Similarity`)) + 
  geom_col(col="blue", fill="blue", size=1.5) +
  ylim(0,1) +
  theme_classic()

# --- Families ---

# Little calculation for paper...

snoc <- which(pim13data2[,3]>0.65)
pim13data2snoc <- pim13data2[snoc,]
famtrus2 <- which(pim13data2snoc[,7]==TRUE)
length(famtrus2)/length(snoc)*100

fnoc <- which(pim13data2[,3]<0.65)
pim13data2fnoc <- pim13data2[fnoc,]
famtrus3 <- which(pim13data2fnoc[,7]==TRUE)
length(famtrus3)/length(fnoc)*100

# Calculation done!

pim13data2fam <- data.frame(as.numeric(pim13data2[,3]))
colnames(pim13data2fam) <- c("Similarity")

famtrus <- which(pim13data2[,7]==TRUE)
pim13data2famtrus <- data.frame(as.numeric(pim13data2[famtrus,3]))
colnames(pim13data2famtrus) <- c("Similarity")

fx <- as.vector(pim13data2fam[,1])
fy <- as.vector(pim13data2famtrus[,1])

fxx <- hist(fx, breaks=21)
fyy <- hist(fy, breaks=21)
fxxx <- hist(fx, breaks=2)
fyyy <- hist(fy, breaks=2)

fzz0 <- fyy
fzz0$counts <- (fyy$counts/fxx$counts)
fzz0$counts[which(fzz0$counts=="NaN")] <- 0

fzzz0 <- fyyy
fzzz0$counts <- (fyyy$counts/fxxx$counts)
fzzz0$counts[which(fzzz0$counts=="NaN")] <- 0

plot(fzz0, col="green")
plot(fzzz0, col="green")

fzz1 <- data.frame(fzz0$breaks[2:21], fzz0$counts)
colnames(fzz1) <- c("Sequence Pair Similarity", "Frequency of Shared Families")
fzzz1 <- data.frame(fzzz0$breaks[1:2], fzzz0$counts)
colnames(fzzz1) <- c("Sequence Pair Similarity", "Frequency of Shared Families")

ggplot(fzz1, aes(y=`Frequency of Shared Families`, x=`Sequence Pair Similarity`)) + 
  geom_col(col="green", fill="green", size=1.5) +
  ylim(0,1) +
  theme_classic()

# --- Making a Stacked Barplot with all four taxonyms ---

pzz2 <- data.frame(pzz1, "Phyla", stringsAsFactors = FALSE)
czz2 <- data.frame(czz1, "Classes", stringsAsFactors = FALSE)
ozz2 <- data.frame(ozz1, "Orders", stringsAsFactors = FALSE)
fzz2 <- data.frame(fzz1, "Families", stringsAsFactors = FALSE)

colnames(pzz2) <- c("Sequence Pair Similarity", "Frequency of Shared Taxonomy", "Shared Taxonomy")
pzz3 <- pzz2
pzz3[,2] <- pzz2[,2]-czz2[,2]

colnames(czz2) <- c("Sequence Pair Similarity", "Frequency of Shared Taxonomy", "Shared Taxonomy")
czz3 <- czz2
czz3[,2] <- czz2[,2]-ozz2[,2]

colnames(ozz2) <- c("Sequence Pair Similarity", "Frequency of Shared Taxonomy", "Shared Taxonomy")
ozz3 <- ozz2
ozz3[,2] <- ozz2[,2]-fzz2[,2]

colnames(fzz2) <- c("Sequence Pair Similarity", "Frequency of Shared Taxonomy", "Shared Taxonomy")

tzz0 <- rbind(pzz2, czz2, ozz2, fzz2)
tzz1 <- rbind(fzz2, pzz3, ozz3, czz3)
tzz1$`Shared Taxonomy` <- factor(tzz1$`Shared Taxonomy`,
                                 levels = c("Phyla", "Classes", "Orders", "Families")) # reorder factors so ggplot uses right order

cols <- character(80)
cols[1:20] <- "red"
cols[21:40] <- "purple"
cols[41:60] <- "blue"
cols[61:80] <- "green"

# --- Side by Side distribution ---

ggplot(tzz0, aes(fill=`Shared Taxonomy`, y=`Frequency of Shared Taxonomy`, x=`Sequence Pair Similarity`)) + 
  geom_col(position=position_dodge(), col=cols) +
  ylim(0,1) +
  scale_fill_manual(breaks=c("Phyla", "Classes", "Orders", "Families"),
                    values=c("Phyla"="red",
                             "Classes"="purple",
                             "Orders"="blue",
                             "Families"="green"),
                    name="") +
  theme_classic()

# --- Stacked distribution ---

ggplot(tzz1, aes(fill=`Shared Taxonomy`, y=`Frequency of Shared Taxonomy`, x=`Sequence Pair Similarity`, order=`Shared Taxonomy`)) + 
  geom_bar(position="stack", stat="identity", col=NA, width=0.05) +
  scale_fill_manual(breaks=c("Phyla", "Classes", "Orders", "Families"),
                    values=c("Phyla"="red",
                             "Classes"="purple",
                             "Orders"="blue",
                             "Families"="green"),
                    name="") +
  ylim(0,1) +
  theme_classic()

# ---------------------------------------------------------------------------------
#
# SKIP! THIS ANALYSIS IS FINISHED!
#
# Adressing Aaron's questions about the ~10% of species with greater than 80% sequence similarity that
# don't infect the same phyla:

sim80 <- which(pim13data2[,3]>0.8)
sim80data <- pim13data2[sim80,]
sim80phy <- which(sim80data[,4]==FALSE)
sim80phydata <- sim80data[sim80phy,]

sim80phydata <- data.frame(sim80phydata[,1:4], "Species 1 Host Phyla"=NA, "Species 2 Host Phyla"=NA, "Missing Phylum Data?"=NA, 
                           "Species 1 Clade"=NA, "Species 2 Clade"=NA, "Species 1 Accession"=NA, "Species 2 Accession"=NA)
colnames(sim80phydata) <- c("Species 1", "Species 2", "Similarity", "Shared Phylum?", "Species 1 Host Phyla", "Species 2 Host Phyla", "Missing Phylum Data?",
                            "Species 1 Clade", "Species 2 Clade", "Species 1 Accession", "Species 2 Accession")

for (i in 1:length(sim80phydata[,1])) {
  
  spe1 <- sim80phydata[i,1]
  spe2 <- sim80phydata[i,2]
  
  num1 <- which(cladestotal2[,1]==spe1)
  num2 <- which(cladestotal2[,1]==spe2)
  
  phy1 <- cladestotal2[num1,9]
  sim80phydata[i,5] <- phy1
  phy1 <- unlist(strsplit(phy1, ";"))
  
  phy2 <- cladestotal2[num2,9]
  sim80phydata[i,6] <- phy2
  phy2 <- unlist(strsplit(phy2, ";"))
  
  if (isEmpty(phy1==TRUE)) {
    sim80phydata[i,7] <- TRUE
  } else if (isEmpty(phy2)==TRUE) {
    sim80phydata[i,7] <- TRUE
  } else {
    sim80phydata[i,7] <- FALSE
  }
  
  cla1 <- cladestotal2[num1,7]
  sim80phydata[i,8] <- cla1
  
  cla2 <- cladestotal2[num2,7]
  sim80phydata[i,9] <- cla2
  
  acc1 <- cladestotal2[num1,8]
  sim80phydata[i,10] <- acc1
  
  acc2 <- cladestotal2[num2,8]
  sim80phydata[i,11] <- acc2
  
  rm(spe1, spe2, num1, num2, phy1, phy2, cla1, cla2, acc1, acc2)
}

write.csv(sim80phydata,"./Analysis/Fig6/80%sim_unsharedphyla_v2.csv", row.names = FALSE)

# Now, for less than 80% similarity:

sim00 <- which(pim13data2[,3]<0.8)
sim00data <- pim13data2[sim00,]
sim00phy <- which(sim00data[,4]==FALSE)
sim00phydata <- sim00data[sim00phy,]

sim00phydata <- data.frame(sim00phydata[,1:4], "Species 1 Host Phyla"=NA, "Species 2 Host Phyla"=NA, "Missing Phylum Data?"=NA, 
                           "Species 1 Clade"=NA, "Species 2 Clade"=NA, "Species 1 Accession"=NA, "Species 2 Accession"=NA)
colnames(sim00phydata) <- c("Species 1", "Species 2", "Similarity", "Shared Phylum?", "Species 1 Host Phyla", "Species 2 Host Phyla", "Missing Phylum Data?",
                            "Species 1 Clade", "Species 2 Clade", "Species 1 Accession", "Species 2 Accession")

for (i in 1:length(sim00phydata[,1])) {
  
  spe1 <- sim00phydata[i,1]
  spe2 <- sim00phydata[i,2]
  
  num1 <- which(cladestotal2[,1]==spe1)
  num2 <- which(cladestotal2[,1]==spe2)
  
  phy1 <- cladestotal2[num1,9]
  sim00phydata[i,5] <- phy1
  phy1 <- unlist(strsplit(phy1, ";"))
  
  phy2 <- cladestotal2[num2,9]
  sim00phydata[i,6] <- phy2
  phy2 <- unlist(strsplit(phy2, ";"))
  
  if (isEmpty(phy1==TRUE)) {
    sim00phydata[i,7] <- TRUE
  } else if (isEmpty(phy2)==TRUE) {
    sim00phydata[i,7] <- TRUE
  } else {
    sim00phydata[i,7] <- FALSE
  }
  
  cla1 <- cladestotal2[num1,7]
  sim00phydata[i,8] <- cla1
  
  cla2 <- cladestotal2[num2,7]
  sim00phydata[i,9] <- cla2
  
  acc1 <- cladestotal2[num1,8]
  sim00phydata[i,10] <- acc1
  
  acc2 <- cladestotal2[num2,8]
  sim00phydata[i,11] <- acc2
  
  rm(spe1, spe2, num1, num2, phy1, phy2, cla1, cla2, acc1, acc2)
}

# Now compare:

sim00ap <- which(sim00phydata[,5]=="Arthropoda" |  sim00phydata[,6]=="Arthropoda")
sim00cp <- which(sim00phydata[,5]=="Chordata" |  sim00phydata[,6]=="Chordata")

sim80ap <- which(sim80phydata[,5]=="Arthropoda" |  sim80phydata[,6]=="Arthropoda")
sim80cp <- which(sim80phydata[,5]=="Chordata" |  sim80phydata[,6]=="Chordata")

sim00cap1 <- which(sim00phydata[,5]=="Arthropoda" | sim00phydata[,5]=="Chordata") 
sim00cap2 <- which(sim00phydata[,6]=="Arthropoda" | sim00phydata[,6]=="Chordata") 
sim00phydatacap1 <- sim00phydata[sim00cap1,]
sim00phydatacap2 <- sim00phydata[sim00cap2,]
sim00cap <- sim00cap1[sim00cap1 %in% sim00cap2]
sim00phydatacap <- sim00phydata[sim00cap,]

sim80cap1 <- which(sim80phydata[,5]=="Arthropoda" | sim80phydata[,5]=="Chordata") 
sim80cap2 <- which(sim80phydata[,6]=="Arthropoda" | sim80phydata[,6]=="Chordata") 
sim80phydatacap1 <- sim80phydata[sim80cap1,]
sim80phydatacap2 <- sim80phydata[sim80cap2,]
sim80cap <- sim80cap1[sim80cap1 %in% sim80cap2]
sim80phydatacap <- sim80phydata[sim80cap,]

# Fraction of species pairs with hosts from different phyla that infect a Chordate host:

length(sim00cp)/length(sim00phydata[,1])*100 # Pairs with <80% similarity
length(sim80cp)/length(sim80phydata[,1])*100 # Pairs with >80% similarity

# Fraction of species pairs with hosts from different phyla that infect an Arthropod host:

length(sim00ap)/length(sim00phydata[,1])*100 # Pairs with <80% similarity
length(sim80ap)/length(sim80phydata[,1])*100 # Pairs with >80% similarity

# Fraction of Chordate / Arthropod pairs among species pairs that have hosts from different phyla:

length(sim00cap)/length(sim00phydata[,1])*100 # Pairs with <80% similarity
length(sim80cap)/length(sim80phydata[,1])*100 # Pairs with >80% similarity

# OKAY! RESUME!

# Variable cleanup

rm(px, py, pxx, pyy, pzz0, pzz1, phytrus, pim13dataphytrus, pim13dataphy)
rm(cx, cy, cxx, cyy, czz0, czz1, clstrus, pim13dataclstrus, pim13datacls)
rm(ox, oy, oxx, oyy, ozz0, ozz1, ordtrus, pim13dataordtrus, pim13dataord)
rm(fx, fy, fxx, fyy, fzz0, fzz1, famtrus, pim13datafamtrus, pim13datafam)
rm(pzz2, czz2, ozz2, fzz2, tzz1)

# ------------------------------------------------------------------------------
#
# Part 5: Comparing Host Environments within Species Pairs (Figure 6C)
#
# ------------------------------------------------------------------------------

# Same idea as previous section, except now we are comparing host environments. 

# First, remove species with missing environments data!!!

cladestotal3 <- cladestotal2
empt <- which(cladestotal3[,3]=="")
emptspe <- as.vector(cladestotal3[empt,1])
emptspe <- emptspe[-c(6,7,8,14)]

pim13data3 <- pim13data

for (i in 1:length(emptspe)) {
  emptrows1 <- which(pim13data3[,1]==emptspe[i])
  emptrows2 <- which(pim13data3[,2]==emptspe[i])
  emptrows <- c(emptrows1, emptrows2)
  pim13data3 <- pim13data3[-emptrows,]
  rm(emptrows, emptrows1, emptrows2)
}

# Now do analysis.

#pim13data3
for(i in 1:length(pim13data3[,1])) {
  spe1 <- pim13data3[i,1]
  spe2 <- pim13data3[i,2]
  num1 <- which(cladestotal3[,1]==spe1)
  num2 <- which(cladestotal3[,1]==spe2)
  
  env1 <- cladestotal3[num1,3]
  env1 <- unlist(strsplit(env1, ";"))
  env2 <- cladestotal3[num2,3]
  env2 <- unlist(strsplit(env2, ";"))
  
  xtemp <- character(0)
  if (isEmpty(env1)==TRUE) {
    env1 <- ""
  }
  for (j in 1:length(env1)) {
    xtemp[j] <- any(env2==env1[j])
  }
  if (any(xtemp==TRUE)==TRUE) {
    pim13data3[i,8] <- TRUE
  } else {
    pim13data3[i,8] <- FALSE
  }
  rm(spe1, spe2, num1, num2, env1, env2, xtemp)
}

# Awesome, now to graph the results in a histogram:

pim13data3env <- data.frame(as.numeric(pim13data3[,3]))
colnames(pim13data3env) <- c("Similarity")

envtrus <- which(pim13data3[,8]==TRUE)
pim13data3envtrus <- data.frame(as.numeric(pim13data3[envtrus,3]))
colnames(pim13data3envtrus) <- c("Similarity")

ex <- as.vector(pim13data3env[,1])
ey <- as.vector(pim13data3envtrus[,1])

exx <- hist(ex, breaks=21)
eyy <- hist(ey, breaks=21)

ezz0 <- eyy
ezz0$counts <- (eyy$counts/exx$counts)
ezz0$counts[which(ezz0$counts=="NaN")] <- 0

plot(ezz0, col="black")

ezz1 <- data.frame(ezz0$breaks[2:21], ezz0$counts)
colnames(ezz1) <- c("Sequence Pair Similarity", "Frequency of Shared Environments")

ggplot(ezz1, aes(y=`Frequency of Shared Environments`, x=`Sequence Pair Similarity`)) + 
  geom_col(col="black", fill="black", size=1.5) +
  ylim(0,1) +
  theme_classic()

# ------------------------------------------------------------------------------
#
# Part 6: Comparing Spore Volume within Species Pairs (Figure 6E)
#
# ------------------------------------------------------------------------------

# Calculate the % spore volume similarity by first averaging the spore volumes for each
# species in a pair, then dividing the smaller spore volume by the larger spore volume 
# in that pair, and multiplying by 100. 

# First, make all the empty sv rows NA:

cladestotal2[which(cladestotal2[,4]=="NA"),4] <- NA
cladestotal2[which(cladestotal2[,4]==""),4] <- NA

#pim13data
for(i in 1:length(pim13data[,1])) {
  spe1 <- pim13data[i,1]
  spe2 <- pim13data[i,2]
  num1 <- which(cladestotal2[,1]==spe1)
  num2 <- which(cladestotal2[,1]==spe2)
  
  sv1 <- cladestotal2[num1,4]
  sv1 <- unlist(strsplit(sv1, ";"))
  sv1 <- as.numeric(sv1)
  sv1a <- mean(sv1)
  sv2 <- cladestotal2[num2,4]
  sv2 <- unlist(strsplit(sv2, ";"))
  sv2 <- as.numeric(sv2)
  sv2a <- mean(sv2)
  
  svs <- c(sv1a, sv2a)
  svs <- svs[order(svs)]
  svdiff <- (svs[1]/svs[2]*100)
  
  pim13data[i,9] <- svdiff
  rm(spe1, spe2, num1, num2, sv1, sv2, sv1a, sv2a, svs, svdiff)
}

# Now need to calculate average % spore volume similarity for each bin:

pim13datasv <- data.frame(as.numeric(pim13data[,3]), as.numeric(pim13data[,9]))
colnames(pim13datasv) <- c("Similarity", "% Spore Volume Similarity")

svnas <- which(is.na(pim13datasv[,2])) # remove pairs without spore volume data (i.e. NAs)
pim13datasv <- pim13datasv[-svnas,]
rownames(pim13datasv) <- (1:length(rownames(pim13datasv)))

svx <- hist(pim13datasv[,1], breaks=21)
svxcounts <- svx$counts
svxbreaksstart <- svx$breaks[1:20]
svxbreaksend <- svx$breaks[2:21]

pim13datasv <- pim13datasv[order(pim13datasv[,1]),] # order by increasing seq similarity
rownames(pim13datasv) <- (1:length(rownames(pim13datasv)))

svy <- data.frame("Similarity Bin Start"=svxbreaksstart, "Similarity Bin End"=svxbreaksend,
                  "Average % Spore Volume Similarity"=NA, stringsAsFactors = FALSE)
colnames(svy) <- c("Similarity Bin Start", "Similarity Bin End", "Average % Spore Volume Similarity")

counter1 <- 0
counter2 <- numeric(0)
for (i in 1:length(svxcounts)) {
  if (svxcounts[i]==0) {
    svy[i,3] <- 0
  } else {
    counter2 <- (counter1 + svxcounts[i])
    svdiffs <- pim13datasv[counter1:counter2,2]
    avgsvdiff <- mean(svdiffs)
    svy[i,3] <- avgsvdiff
    counter1 <- (counter2)
    rm(svdiffs, avgsvdiff)
  }
}
rm(counter1, counter2)

# Awesome, now make a histogram.

svy <- svy[,2:3]
colnames(svy) <- c("Sequence Pair Similarity", "Average % Spore Volume Similarity")

ggplot(svy, aes(y=`Average % Spore Volume Similarity`, x=`Sequence Pair Similarity`)) + 
  geom_col(col="black", fill="black", size=1.5) +
  ylim(0,100) +
  theme_classic()

# Variable cleanup:

rm(svnas, svx, svxcounts, svxbreaksend, svxbreaksstart, svy)

# ------------------------------------------------------------------------------
#
# Part 7: Comparing Polar Tube Length within Species Pairs (Figure 6F)
#
# ------------------------------------------------------------------------------

# Sameas above, but for polar tube lengths. 

# First, make all the empty pt length rows NA:

cladestotal2[which(cladestotal2[,5]==""),5] <- NA

#pim13data
for(i in 1:length(pim13data[,1])) {
  spe1 <- pim13data[i,1]
  spe2 <- pim13data[i,2]
  num1 <- which(cladestotal2[,1]==spe1)
  num2 <- which(cladestotal2[,1]==spe2)
  
  pt1 <- cladestotal2[num1,5]
  pt1 <- unlist(strsplit(pt1, ";"))
  pt1 <- as.numeric(pt1)
  pt1a <- mean(pt1)
  pt2 <- cladestotal2[num2,5]
  pt2 <- unlist(strsplit(pt2, ";"))
  pt2 <- as.numeric(pt2)
  pt2a <- mean(pt2)
  
  pts <- c(pt1a, pt2a)
  pts <- pts[order(pts)]
  ptdiff <- (pts[1]/pts[2]*100)
  
  pim13data[i,10] <- ptdiff
  rm(spe1, spe2, num1, num2, pt1, pt2, pt1a, pt2a, pts, ptdiff)
}

# Now need to calculate average % spore volume similarity for each bin:

pim13datapt <- data.frame(as.numeric(pim13data[,3]), as.numeric(pim13data[,10]))
colnames(pim13datapt) <- c("Similarity", "% Polar Tube Length Similarity")

ptnas <- which(is.na(pim13datapt[,2])) # remove pairs without spore volume data (i.e. NAs)
pim13datapt <- pim13datapt[-ptnas,]
rownames(pim13datapt) <- (1:length(rownames(pim13datapt)))

ptx <- hist(pim13datapt[,1], breaks=21)
ptxcounts <- ptx$counts
ptxbreaksstart <- ptx$breaks[1:20]
ptxbreaksend <- ptx$breaks[2:21]

pim13datapt <- pim13datapt[order(pim13datapt[,1]),] # order by increasing seq similarity
rownames(pim13datapt) <- (1:length(rownames(pim13datapt)))

pty <- data.frame("Similarity Bin Start"=ptxbreaksstart, "Similarity Bin End"=ptxbreaksend,
                  "Average % Polar Tube Length Similarity"=NA, stringsAsFactors = FALSE)
colnames(pty) <- c("Similarity Bin Start", "Similarity Bin End", "Average % Polar Tube Length Similarity")

counter1 <- 0
counter2 <- numeric(0)
for (i in 1:length(ptxcounts)) {
  if (ptxcounts[i]==0) {
    pty[i,3] <- 0
  } else {
    counter2 <- (counter1 + ptxcounts[i])
    ptdiffs <- pim13datapt[counter1:counter2,2]
    avgptdiff <- mean(ptdiffs)
    pty[i,3] <- avgptdiff
    counter1 <- (counter2)
    rm(ptdiffs, avgptdiff)
  }
}
rm(counter1, counter2)

# Awesome, now make a histogram.

pty <- pty[,2:3]
colnames(pty) <- c("Sequence Pair Similarity", "Average % Polar Tube Length Similarity")

ggplot(pty, aes(y=`Average % Polar Tube Length Similarity`, x=`Sequence Pair Similarity`)) + 
  geom_col(col="black", fill="black", size=1.5) +
  ylim(0,100) +
  theme_classic()

# Variable cleanup:

rm(ptnas, ptx, ptxcounts, ptxbreaksend, ptxbreaksstart, pty)

# ------------------------------------------------------------------------------
#
# Part 8: Comparing Infected Tissues within Species Pairs (6D)
#
# ------------------------------------------------------------------------------

# Working directory:

setwd(anadir)

# First, read the Host Tissue Data and Host Tissue Dictionary .csv files, and convert them
# into useable R data frames:

TisDat1 <- read.csv("./Host_Tissue_Data.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
TisDat1 <- data.frame(TisDat1[!apply(is.na(TisDat1) | TisDat1 == "", 1, all),])
colnames(TisDat1) <- TisDat1[1,]
TisDat1 <- TisDat1[-1,]

TDspe <- as.vector(TisDat1$Species)
TDtis <- as.vector(TisDat1$Tissues)

TDspe <- gsub("unnamed ", "", TDspe)
TDspe <- gsub("Unnamed ", "", TDspe)
TDspe <- gsub("\\s*\\([^\\)]+\\)","", TDspe)
TDspe <- gsub("\\.", "", TDspe)       # "." normally refers to "any character", so, we escape it with "//" to
TDspe <- gsub(" ", "_", TDspe)        # search for an actual "."
TDspe <- gsub(";_", ";", TDspe)

TDtis <- gsub("unnamed ", "", TDtis)
TDtis <- gsub("Unnamed ", "", TDtis)
TDtis <- gsub("\\s*\\([^\\)]+\\)","", TDtis)
TDtis <- gsub("\\.", "", TDtis)       # "." normally refers to "any character", so, we escape it with "//" to
TDtis <- gsub(" ", "_", TDtis)        # search for an actual "."
TDtis <- gsub(";_", ";", TDtis)

TisDat2 <- data.frame(TDspe, TDtis, NA, stringsAsFactors = FALSE)
colnames(TisDat2) <- c("species", "tissues_raw", "tissues_sorted")

TisDict1 <- read.csv("./Host_Tissue_Dictionary.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
TisDict1 <- data.frame(TisDict1[!apply(is.na(TisDict1) | TisDict1 == "", 1, all),])
colnames(TisDict1) <- TisDict1[1,]
TisDict1 <- TisDict1[-1,]

TisDict2 <- data.frame("tissue_raw"=NA, "tissue_sorted"=NA)
colnames(TisDict2) <- c("tissue_raw", "tissue_sorted")
TDcols <- colnames(TisDict1)

counter1 <- 0
for (i in 1:length(colnames(TisDict1))) {
  rawtis <- TisDict1[,i]
  rawtis <- rawtis[rawtis != ""]
  sortis <- TDcols[i]
  for (j in 1:length(rawtis)) {
    counter1 <- (counter1+1)
    TisDict2[counter1,1] <- rawtis[j]
    TisDict2[counter1,2] <- sortis
  }
  rm(rawtis, sortis)
}
rm(counter1)

TDraw <- as.vector(TisDict2[,1])
TDsor <- as.vector(TisDict2[,2])

TDraw <- gsub("unnamed ", "", TDraw)
TDraw <- gsub("Unnamed ", "", TDraw)
TDraw <- gsub("\\s*\\([^\\)]+\\)","", TDraw)
TDraw <- gsub("\\.", "", TDraw)       # "." normally refers to "any character", so, we escape it with "//" to
TDraw <- gsub(" ", "_", TDraw)        # search for an actual "."
TDraw <- gsub(";_", ";", TDraw)

TDsor <- gsub("unnamed ", "", TDsor)
TDsor <- gsub("Unnamed ", "", TDsor)
TDsor <- gsub("\\s*\\([^\\)]+\\)","", TDsor)
TDsor <- gsub("\\.", "", TDsor)       # "." normally refers to "any character", so, we escape it with "//" to
TDsor <- gsub(" ", "_", TDsor)        # search for an actual "."
TDsor <- gsub(";_", ";", TDsor)

dup <- which(duplicated(TDraw)==TRUE)
TDraw <- TDraw[-dup]
TDsor <- TDsor[-dup]

TisDict2 <- data.frame("tissues_raw"=TDraw, "tissues_sorted"=TDsor, stringsAsFactors = FALSE)

# Awesome! We now have to go through every species in the TisDat2 data frame, and convert it's 
# infected tissue names into the streamlined tissue names from TisDict2. Lets do it:

for (i in 1:length(as.vector(TisDat2[,1]))) {
  rawtis <- TisDat2[i,2]
  if (rawtis=="") {
    TisDat2[i,3] <- ""
    rm(rawtis)
    } else {
    rawtis <- unlist(strsplit(rawtis, ";"))
    sortis <- character(0)
    for (j in 1:length(rawtis)) {
      num <- which(TDraw==rawtis[j])
      sortis[j] <- TDsor[num]
    }
    sortis <- unique(sortis)
    sortis <- paste(sortis, collapse=";")
    TisDat2[i,3] <- sortis
    rm(rawtis, sortis)
  }
}

# Cool. Next, we need to add all of this information to the cladestotal data frame.

cladestotal2 <- data.frame(cladestotal2, "tissues"=NA, stringsAsFactors=FALSE)

for (i in 1:length(files1)) {
  xtemp <- cladestotal2[i,1]
  num <- which(TisDat2[,1]==xtemp)
  ytemp <- TisDat2[num,3]
  cladestotal2[i,13] <- ytemp[1]
  rm(ytemp, xtemp, num)
}

cladestotal3 <- cladestotal2
empt <- which(cladestotal3[,13]=="")
napt <- which(is.na(cladestotal3[,13]))
empt <- c(empt, napt)
emptspe <- as.vector(cladestotal3[empt,1])
emptspe <- emptspe[-c(15,16,17,41)]

pim13data4 <- data.frame(pim13data, "Shared Tissue?"=NA, stringsAsFactors=FALSE)
colnames(pim13data4) <- c(colnames(pim13data), "Shared Tissue?")

# First, remove all the rows in the pim data frame that contain species with
# missing tissue data.

for (i in 1:length(emptspe)) {
  emptrows1 <- which(pim13data4[,1]==emptspe[i])
  emptrows2 <- which(pim13data4[,2]==emptspe[i])
  emptrows <- c(emptrows1, emptrows2)
  pim13data4 <- pim13data4[-emptrows,]
  rm(emptrows, emptrows1, emptrows2)
}

# Now do analysis.

#pim13data4
for(i in 1:length(pim13data4[,1])) {
  spe1 <- pim13data4[i,1]
  spe2 <- pim13data4[i,2]
  num1 <- which(cladestotal3[,1]==spe1)
  num2 <- which(cladestotal3[,1]==spe2)
  
  tis1 <- cladestotal3[num1,13]
  tis1 <- unlist(strsplit(tis1, ";"))
  tis2 <- cladestotal3[num2,13]
  tis2 <- unlist(strsplit(tis2, ";"))
  
  xtemp <- character(0)
  if (isEmpty(tis1)==TRUE) {
    tis1 <- ""
  }
  for (j in 1:length(tis1)) {
    xtemp[j] <- any(tis2==tis1[j])
  }
  if (any(xtemp==TRUE)==TRUE) {
    pim13data4[i,11] <- TRUE
  } else {
    pim13data4[i,11] <- FALSE
  }
  rm(spe1, spe2, num1, num2, tis1, tis2, xtemp)
}

# Now, make and plot the distribution.

pim13data4tis <- data.frame(as.numeric(pim13data4[,3]))
colnames(pim13data4tis) <- c("Similarity")

tistrus <- which(pim13data4[,11]==TRUE)
pim13data4tistrus <- data.frame(as.numeric(pim13data4[tistrus,3]))
colnames(pim13data4tistrus) <- c("Similarity")

tx <- as.vector(pim13data4tis[,1])
ty <- as.vector(pim13data4tistrus[,1])

txx <- hist(tx, breaks=21)
tyy <- hist(ty, breaks=21)

tzz0 <- tyy
tzz0$counts <- (tyy$counts/txx$counts)
tzz0$counts[which(tzz0$counts=="NaN")] <- 0

plot(tzz0, col="black")

tzz1 <- data.frame(tzz0$breaks[2:21], tzz0$counts)
colnames(tzz1) <- c("Sequence Pair Similarity", "Frequency of Shared Tissues")

ggplot(tzz1, aes(y=`Frequency of Shared Tissues`, x=`Sequence Pair Similarity`)) + 
  geom_col(col="black", fill="black", size=1.5) +
  ylim(0,1) +
  theme_classic()

# Done!

# Write cladestotal3 data frame to xlsx:

setwd(anadir)
write_xlsx(cladestotal3, "./cladestotal3.xlsx")


