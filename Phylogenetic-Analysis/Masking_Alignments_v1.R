# -----------------------------------------------------------------------------
#
# Alignment Masking - v1
#
# Brandon Murareanu - Created: 2020/07/17
#                     Last Edit: 2021/02/19 
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write some code that can mask a multiple sequence alignment.
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

library(BiocManager)
library(ape)
library(ggtree)
library(treeio)
library(tidytree)
library(ggnewscale)
library(ggplot2)
library(DECIPHER)

# Set working directory. Modify as need. For me, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Test_Trees/Rozella_Outgroup/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Test_Trees/Rozella_Outgroup"

# Big Tree 1 directory:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_1/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_1"

# Big Tree 2 directory:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_2/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_2"

# 230 Species Tree directory:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Test_Trees/230_Species/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Test_Trees/230_Species"

# For the Multi Trees Base Tree, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Multi_Trees/Base_Tree/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Multi_Trees/Base_Tree/"

# For 270 Species Tree, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_Tree/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_tree"

setwd(maindir)

# -----------------------------------------------------------------------------

# First, read the alignment into a DNA string from Biostrings.

MSAa <- readDNAStringSet("./MSA13.0_MAFFT_(FASTA).mafft", "fasta")

# Next, should be abe to mask the alignment with the following function from DECIPHER.

MSAb <- MaskAlignment(MSAa,
                        type = "sequences",
                        windowSize = 5,
                        threshold = 1,
                        maxFractionGaps = 0.2,
                        includeTerminalGaps = FALSE,
                        correction = FALSE,
                        showPlot = TRUE)

# Display only unmasked nucleotides for use in further analysis.

MSAb_not_masked <- as(MSAb, "DNAStringSet")
BrowseSeqs(MSAb_not_masked)

# Write the alignment to a file.

writeXStringSet(MSAb_not_masked, "./MSA13.0_MAFFT_(FASTA)_masked.mafft",
                format="fasta")

# -----------------------------------------------------------------------------



  
  