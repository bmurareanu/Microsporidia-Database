# -----------------------------------------------------------------------------
#
# Reading FASTA Files - v3
#
# Brandon Murareanu - Created: 2020/07/08
#                     Last Edit: 2021/02/19 
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write a script that can read multiple FASTA files and concatenate them
# into one large multi-layered FASTA file. This script differs from version 2
# in that the FASTA headers will be named according to accession number.
#
# -----------------------------------------------------------------------------
#
# First, load necessary packages.

if (! requireNamespace("ape", quietly=TRUE)) {
  install.packages("ape")
}

library(ape)

# Set the working directory. For the 1st Big Tree, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_1/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_1"

# Set the working directory. For the 2nd Big Tree, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_2/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_2"

# For the 230 species tree, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Test_Trees/230_Species/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Test_Trees/230_Species"

# For the Multi Trees Base Tree, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Multi_Trees/Base_Tree/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Multi_Trees/Base_Tree"

# For 270 Species Tree, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_Tree/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_tree"

setwd(sourcedir)

# -----------------------------------------------------------------------------

# Read all the file names (species names) into a vector:

files1 <- list.files()

# Now, read the FASTA files into a DNAbin object master list:

dnabin <- read.FASTA(files1[1])     # Initialize DNAbin master list.
accs <- character(length(files1))
for (i in 1:length(files1)) {
  temp <- read.FASTA(files1[i])
  dnabin[i] <- temp                 # read DNAbin.
  accs[i] <- names(temp)            # store accessions.
}

names(dnabin) <- accs               # name each object in the DNAbin master list
                                    # according to the accession number

# Read the DNAbin master list into a single multi-layered FASTA file:

setwd(maindir)
write.FASTA(dnabin, "multiFASTA")

# Create a data frame associating each accession number with each species name
# for future reference.

spenam <- data.frame(Species=files1,
                     Accessions=accs)

# ----------------------------------------------------------------------------
#
# Done!
