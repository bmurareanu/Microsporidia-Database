# -----------------------------------------------------------------------------
#
# Reading Master List - v1
#
# Brandon Murareanu - Created: 2020/08/31 
#                     Last Edit: 2021/02/19
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write a script that can read the species master list, and compare it to
# the set of 18S seqeunces to determine which species are or are not in the database. 
#
# -----------------------------------------------------------------------------

# Set directories for the master list, and 18S sequence set:

listdir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/SSU_Sequences/Master_Lists"
sourcedir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Big_Tree_1/Source"

setwd(listdir)

# Read species from Master List .csv file:

MasterList1 <- read.csv("./Master_List_Species_2.csv", quote="\"", stringsAsFactors=FALSE, header=FALSE)
MasterList2 <- as.vector(MasterList1$V1)

# Substitute improper terms and characetrs to match the format of the 18S sequence set:

MasterList3 <- gsub("unnamed ", "", MasterList2)
MasterList3 <- gsub("Unnamed ", "", MasterList3)
MasterList3 <- gsub("\\s*\\([^\\)]+\\)","", MasterList3)
MasterList3 <- gsub("\\.", "", MasterList3)       # "." normally refers to "any character", so, we escape it with "//" to
MasterList3 <- gsub(" ", "_", MasterList3)        # search for an actual "."


# Now, query the Master List with the 18S Sequence Set to find which 18S sequences represent
# species names that are  not in the Master List:

setwd(sourcedir)
spe1 <- list.files()

xtemp <- logical(0)
for (i in 1:length(spe1)) {
  xtemp[i] <- (spe1[i] %in% MasterList3)
}

spefnum <- which(xtemp == FALSE)  
spef <- spe1[spefnum]

spetnum <- which(xtemp == TRUE)  
spet <- spe1[spetnum]

length(spef) # 785 species not in the Master List.
length(spet) # 242 species in the Master List.

# Now, save these spef and spet datasets as R objects to be loaded in any R console:

setwd(listdir)

SpeF1 <- as.data.frame(spef)
SpeT1 <- as.data.frame(spet)
MasterListc <- as.data.frame(MasterList3)

save(MasterListc, file="Master_List_cleaned.RData")
save(SpeF1, file="SpeF.RData")
save(SpeT1, file="SpeT.RData")

# Done!

# ----------------------------------------------------------------------------------------

# Just manually adding some species that should be in the list, but weren't, because of stuff
# like alternate spelling, different genus name but same species name, or just because they
# would be good to place in the phylogeny, etc.

# Alternate spelling or different genus name in the master list:

addit1 <- spef[c(1,2,780,19,38,111,125,134,136,142,146,150,591,594,706,714,715,716,717,718,720,721,741,743,753,754,755,756,760)]

# Not in master list, but would be good to place in the phylogeny:

addit2 <- spef[c(20,37)]

# Outgroup (Rozella allomycis):

addit3 <- spef[742]

# Combine all categories with the Spet list:

spet2 <- c(spet, addit1, addit2, addit3)

# ----------------------------------------------------------------------------------------

# Now, query FASTA files from source folder, and move the ones with species names corresponding
# to the spet2 list to a separate folder for alignment and phylogenetic analysis.

destdir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/Phylogenetic_Trees/Source2"
for (i in 1:length(spet2)) {
  sname <- spet2[i]
  dir1 <- paste0(sourcedir, "/", sname)
  dir2 <- paste0(destdir, "/", sname)
  file.copy(from=dir1, to=dir2)
}

# Check:

setwd(destdir)
list.files()

# Done!

# ----------------------------------------------------------------------------------------
