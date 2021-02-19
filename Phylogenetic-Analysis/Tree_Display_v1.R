# -----------------------------------------------------------------------------
#
# Tree Display - v1
#
# Brandon Murareanu - Created: 2020/10/22
#                     Last Edit: 2021/02/19
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write some code to edit the display of trees, like adding identifiers
# that can be converted into symbols signifying a certain trait about a species,
# say, environment.
#
# -----------------------------------------------------------------------------

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

# Set working directory. Modify as needed. For me, it's this:

sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_Tree/Source"
maindir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/270_Spe_tree"

setwd(maindir)

# -----------------------------------------------------------------------------------------------
#
# Part 1: Formatting Data for Tree Annotations
#
# -----------------------------------------------------------------------------------------------

# First, we are going to go through the "cladestotal" data frame, and assign numbers to each species
# corresponding to the environments that they are found in.

envs <- cladestotal$environments
spec <- as.vector(cladestotal$species)
envsID <- character(length(files1))

for (i in 1:length(files1)) {
  xtemp <- files1[i]
  num <- which(spec==xtemp)
  ytemp <- envs[num]
  ytemp <- unlist(strsplit(ytemp, ';'))
  ztemp <- numeric(4)
  
  brac <- any(ytemp=="brackish")
  if (brac==TRUE) {
    ztemp[1] <- 1
  } else {
    ztemp[1] <- 0
  }
  fres <- any(ytemp=="freshwater")
  if (fres==TRUE) {
    ztemp[2] <- 2
  } else {
    ztemp[2] <- 0
  }
  mari <- any(ytemp=="marine")
  if (mari==TRUE) {
    ztemp[3] <- 3
  } else {
    ztemp[3] <- 0
  }
  terr <- any(ytemp=="terrestrial")
  if (terr==TRUE) {
    ztemp[4] <- 4
  } else {
    ztemp[4] <- 0
  }
  
  qtemp <- paste(ztemp, collapse=" ")
  envsID[i] <- qtemp
}

taxanames <- paste(envsID, files1)

# Now, we are going to use the envsID to create a data frame to be the basis for a ggtree heatmap.

tre13.2 <- read.mrbayes("./Tree13.2.nex.con.tre")
tre13.2@phylo[["tip.label"]] <- files1

df <- data.frame("Species"=files1, "Brackish"="", "Freshwater"="", "Marine"="", "Terrestrial"="", stringsAsFactors=FALSE)
rownames(df) <- files1
df <- df[-1]

for (i in 1:length(envsID)) {
  xtemp <- envsID[i]
  ytemp <- unlist(strsplit(xtemp, ' '))
  df[i,1] <- ytemp[1]
  df[i,2] <- ytemp[2]
  df[i,3] <- ytemp[3]
  df[i,4] <- ytemp[4]
}

missingenvs <- rownames(filter(df, Brackish==0, Freshwater==0, Marine==0, Terrestrial==0))
nums <- integer(length(missingenvs))
for (i in 1:length(missingenvs)) {
  xtemp <- missingenvs[i]
  ytemp <- which(files1==xtemp)
  nums[i] <- ytemp
  rm(xtemp, ytemp)
}
df[nums,] <- 5

# Now, we are going to go through the "cladestotal" data frame again, assign numbers to each species
# corresponding to their multiple host environments tag, and create a data frame that will be the basis
# for a second ggtree heatmap:

mhtags <- cladestotal$`multiple host environments tag`
spec <- as.vector(cladestotal$species)
mhdf <- data.frame("Species"=files1, "mhtag1"="", "mhtag2"="", stringsAsFactors=FALSE)
rownames(mhdf) <- files1
mhdf <- mhdf[-1]

for (i in 1:length(files1)) {
  xtemp <- files1[i]
  num <- which(spec==xtemp)
  ytemp <- mhtags[num]
  ztemp <- numeric(2)
  
  if (ytemp=="Multiple_Hosts-Different_Environments") {
    ztemp[1] <- 6
    ztemp[2] <- 6
  } else if (ytemp=="Multiple_Hosts-Same_Environments") {
    ztemp[1] <- 6
    ztemp[2] <- 0
  } else {
    ztemp[1] <- 0
    ztemp[2] <- 0
  }
  mhdf[i,] <- ztemp
}

mhdf2 <- data.frame(label=rownames(mhdf),
                    mhtag1=mhdf[,1],
                    mhtag2=mhdf[,2])

# Awesome! Now to create the actual tree plot.

# -----------------------------------------------------------------------------------------------
#
# Part 2: Display Bayesian Tree with Annotations (Figure 5A)
#
# -----------------------------------------------------------------------------------------------

# First, initialize a separate R window to visualize plots:

win.graph()

# Load basic tree:

tree <- ggtree(tre13.2, branch.length="none", layout="circular", size=0.75)

# Group clades:

cladelist <- list(clade1=as.vector(clade1[,1]),
                  clade3=as.vector(clade3[,1]),
                  clade4=as.vector(clade4[,1]),
                  clade5=as.vector(clade5[,1]))
tree2 <- groupOTU(tree, cladelist) + aes(color=group) +
  scale_color_manual(values=c("black", "#F73CFA", "#37D132", "#35B9EE", "#F89714"))

# Add environment heatmap:

tree3 <- gheatmap(tree2, 
         df[, c("Brackish", "Freshwater", "Marine", "Terrestrial")], 
         width=0.2,
         offset=0.4, 
         hjust=0,
         color="black",
         colnames=FALSE)

# Add multiple hosts tags, colour scales, and display plot:

tree3 %<+% mhdf2 +
  geom_tippoint(shape=21, aes(x=x+5.8, fill=mhtag1), color="white", size=2) +
  geom_tippoint(shape=21, aes(x=x+6.45, fill=mhtag2), color="white", size=2) +
  scale_fill_manual(breaks=c("0", "1", "2", "3", "4", "5", "6"),
                    values=c("0"="white",
                             "1"="#76D7C4",
                             "2"="#AED6F1", 
                             "3"="#5DADE2",
                             "4"="#1E8449",
                             "5"="dark grey",
                             "6"="black"),
                    name="")

# Awesome! Now for the supplemental trees.

# -----------------------------------------------------------------------------------------------
#
# Part 3: Display Individual Clades with Name and Annotations (Figure 5 Supplemental 1)
#
# -----------------------------------------------------------------------------------------------

# Okay sweet. Now for the suplemental, we need to isolate each of the clades and display them
# by themselves, with labels, confidence values, and the annotations from the big tree.

tree <- ggtree(tre13.2, branch.length="none", size=0.75) 

# Group clades:

cladelist <- list(clade1=as.vector(clade1[,1]),
                  clade3=as.vector(clade3[,1]),
                  clade4=as.vector(clade4[,1]),
                  clade5=as.vector(clade5[,1]))
tree2 <- groupOTU(tree, cladelist) + aes(color=group) +
  scale_color_manual(values=c("black", "#F73CFA", "#37D132", "#35B9EE", "#F89714"))

# View node numbers, for reference in the next steps:

tree2 + geom_text(aes(label=node))

# Collapse appropriate clades based on node numbers:

# --- Clade 1 ---

# Subset data frames for annotations:

clade1names <- c(as.vector(cladelist$clade1), as.vector(noclade$species))
dfnames <- rownames(df)
c1dfnums <- numeric(length(clade1names))
for (i in 1:length(clade1names)) {
  xtemp <- character(0)
  xtemp <- clade1names[i]
  c1dfnums[i] <- which(xtemp==dfnames)
  rm(xtemp)
}
c1df <- df[c1dfnums,]
c1mhdf <- mhdf2[c1dfnums,]

# Collapse nodes:

tree3 <- ggtree::collapse(tree2, node=346)

# Add support values and species names:

tree4 <- tree3 + geom_tiplab(offset=3, size=3) + 
  geom_text(aes(label=prob_percent), hjust=-0.25, size=3, nudge_y=0.15) + 
  ggplot2::xlim(0, 30)

# Add environment heatmap:

tree5 <- gheatmap(tree4, 
         c1df[, c("Brackish", "Freshwater", "Marine", "Terrestrial")], 
         width=0.09,
         offset=-0.05, 
         hjust=0,
         color="black",
         colnames=FALSE)

# Add multiple hosts tags, colour scales, and display plot:

tree5 %<+% c1mhdf +
  geom_tippoint(shape=21, aes(x=x+2.4, fill=mhtag1), color="white", size=2) +
  geom_tippoint(shape=21, aes(x=x+2.7, fill=mhtag2), color="white", size=2) +
  scale_fill_manual(breaks=c("0", "1", "2", "3", "4", "5", "6"),
                    values=c("0"="white",
                             "1"="#76D7C4",
                             "2"="#AED6F1", 
                             "3"="#5DADE2",
                             "4"="#1E8449",
                             "5"="dark grey",
                             "6"="black"),
                    name="")

# --- Clade 3 ---

# Subset data frames for annotations:

clade3names <- as.vector(cladelist$clade3)
dfnames <- rownames(df)
c3dfnums <- numeric(length(clade3names))
for (i in 1:length(clade3names)) {
  xtemp <- character(0)
  xtemp <- clade3names[i]
  c3dfnums[i] <- which(xtemp==dfnames)
  rm(xtemp)
}
c3df <- df[c3dfnums,]
c3mhdf <- mhdf2[c3dfnums,]

# Collapse nodes:

tree3 <- ggtree::collapse(tree2, node=538)
tree3 <- ggtree::collapse(tree3, node=276)
tree3 <- ggtree::collapse(tree3, node=347)
tree3 <- ggtree::collapse(tree3, node=285)

# Add support values and species names:

tree4 <- tree3 + geom_tiplab(offset=2.2, size=3) + 
  geom_text(aes(label=prob_percent), hjust=-0.25, size=3, nudge_y=0.15) + 
  ggplot2::xlim(0, 30)

# Add environment heatmap:

tree5 <- gheatmap(tree4, 
                  c3df[, c("Brackish", "Freshwater", "Marine", "Terrestrial")], 
                  width=0.07,
                  offset=-0.05, 
                  hjust=0,
                  color="black",
                  colnames=FALSE)

# Add multiple hosts tags, colour scales, and display plot:

tree5 %<+% c3mhdf +
  geom_tippoint(shape=21, aes(x=x+1.825, fill=mhtag1), color="white", size=2) +
  geom_tippoint(shape=21, aes(x=x+2.05, fill=mhtag2), color="white", size=2) +
  scale_fill_manual(breaks=c("0", "1", "2", "3", "4", "5", "6"),
                    values=c("0"="white",
                             "1"="#76D7C4",
                             "2"="#AED6F1", 
                             "3"="#5DADE2",
                             "4"="#1E8449",
                             "5"="dark grey",
                             "6"="black"),
                    name="")

# --- Clade 4 ---

# Subset data frames for annotations:

clade4names <- as.vector(cladelist$clade4)
dfnames <- rownames(df)
c4dfnums <- numeric(length(clade4names))
for (i in 1:length(clade4names)) {
  xtemp <- character(0)
  xtemp <- clade4names[i]
  c4dfnums[i] <- which(xtemp==dfnames)
  rm(xtemp)
}
c4df <- df[c4dfnums,]
c4mhdf <- mhdf2[c4dfnums,]

# Collapse nodes:

tree3 <- ggtree::collapse(tree2, node=538)
tree3 <- ggtree::collapse(tree3, node=276)
tree3 <- ggtree::collapse(tree3, node=285)
tree3 <- ggtree::collapse(tree3, node=457)
tree3 <- ggtree::collapse(tree3, node=427)

# Add support values and species names:

tree4 <- tree3 + geom_tiplab(offset=2.2, size=3) + 
  geom_text(aes(label=prob_percent), hjust=-0.25, size=3, nudge_y=0.15) + 
  ggplot2::xlim(0, 30)

# Add environment heatmap:

tree5 <- gheatmap(tree4, 
                  c4df[, c("Brackish", "Freshwater", "Marine", "Terrestrial")], 
                  width=0.07,
                  offset=-0.05, 
                  hjust=0,
                  color="black",
                  colnames=FALSE)

# Add multiple hosts tags, colour scales, and display plot:

tree5 %<+% c4mhdf +
  geom_tippoint(shape=21, aes(x=x+1.825, fill=mhtag1), color="white", size=2) +
  geom_tippoint(shape=21, aes(x=x+2.05, fill=mhtag2), color="white", size=2) +
  scale_fill_manual(breaks=c("0", "1", "2", "3", "4", "5", "6"),
                    values=c("0"="white",
                             "1"="#76D7C4",
                             "2"="#AED6F1", 
                             "3"="#5DADE2",
                             "4"="#1E8449",
                             "5"="dark grey",
                             "6"="black"),
                    name="")

# --- Clade 5 ---

# Subset data frames for annotations:

clade5names <- as.vector(cladelist$clade5)
dfnames <- rownames(df)
c5dfnums <- numeric(length(clade5names))
for (i in 1:length(clade5names)) {
  xtemp <- character(0)
  xtemp <- clade5names[i]
  c5dfnums[i] <- which(xtemp==dfnames)
  rm(xtemp)
}
c5df <- df[c5dfnums,]
c5mhdf <- mhdf2[c5dfnums,]

# Collapse nodes:

tree3 <- ggtree::collapse(tree2, node=538)
tree3 <- ggtree::collapse(tree3, node=276)
tree3 <- ggtree::collapse(tree3, node=285)
tree3 <- ggtree::collapse(tree3, node=457)
tree3 <- ggtree::collapse(tree3, node=348)

# Add support values and species names:

tree4 <- tree3 + geom_tiplab(offset=2.2, size=3) + 
  geom_text(aes(label=prob_percent), hjust=-0.25, size=3, nudge_y=0.15) + 
  ggplot2::xlim(0, 30)

# Add environment heatmap:

tree5 <- gheatmap(tree4, 
                  c5df[, c("Brackish", "Freshwater", "Marine", "Terrestrial")], 
                  width=0.07,
                  offset=-0.05, 
                  hjust=0,
                  color="black",
                  colnames=FALSE)

# Add multiple hosts tags, colour scales, and display plot:

tree5 %<+% c5mhdf +
  geom_tippoint(shape=21, aes(x=x+1.825, fill=mhtag1), color="white", size=2) +
  geom_tippoint(shape=21, aes(x=x+2.05, fill=mhtag2), color="white", size=2) +
  scale_fill_manual(breaks=c("0", "1", "2", "3", "4", "5", "6"),
                    values=c("0"="white",
                             "1"="#76D7C4",
                             "2"="#AED6F1", 
                             "3"="#5DADE2",
                             "4"="#1E8449",
                             "5"="dark grey",
                             "6"="black"),
                    name="")

# Done!

# Awesome! Now we can save these graphs as .pdf or .svg files and edit in inkscape!

