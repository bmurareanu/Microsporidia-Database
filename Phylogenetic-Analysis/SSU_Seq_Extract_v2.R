# -----------------------------------------------------------------------------
#
# SSU Sequence Extraction - v2
#
# Brandon Murareanu - Created: 2020/05/20 
#                     Last Edit: 2021/02/19 
#
# Reinke Lab - Microsporidia Database Project
#
# Goal: Write a script that can extract a single 18S sequence for each unique 
# species, as well as what host it was found in, from the results of an NCBI
# nucleotide database search. 
#
# -----------------------------------------------------------------------------
#
# Follow the instructions given in the shared project folder to get 18S 
# microsporidia sequences from the NCBI. Save the accession numbers into
# a file, then copy the accession numbers to an excel spreadsheet. You can 
# also add your own accession numbers to the spreadsheet as well. Then, save
# the spreadsheet as a csv file named "accessions.csv". This just makes it 
# easier to import as an R object for manipulation.

# Now, make sure the following packages are installed and loaded:

if (! requireNamespace("ape", quietly=TRUE)) {
  install.packages("ape")
}
if (! requireNamespace("rlist", quietly=TRUE)) {
  install.packages("rlist")
}

library(ape)
library(rlist)

# Set working directories:

SSUdir <- "C:/Users/*****/OneDrive/Documents/R_Projects/Microsporidia_Database/SSU_Sequences"
setwd(SSUdir)

# Read the accession numbers from "accessions.csv" into a table which we
# can manipulate:

acc0 <- read.table("./accessions.csv", quote="\"", stringsAsFactors=FALSE)
acc1 <- unique(acc0)

# Read the table into a list of characters:

acc2 <- as.list(acc1)$V1

# Now, in order to be able to extract host names from GenBank in addition to
# other information (i.e. sequences, species names), I made some modifications
# to the "read.GenBank" function from ape. The new and improved function is
# called "read.GenBank2", and is shown below:

read.GenBank2 <- function(access.nb, seq.names = access.nb, species.names = TRUE,
                          host.names = TRUE, as.character = FALSE)
{
  N <- length(access.nb)
  ## if more than 400 sequences, we break down the requests
  a <- 1L
  b <- if (N > 400) 400L else N
  fl <- tempfile()
  repeat {
    URL <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                  paste(access.nb[a:b], collapse = ","), "&rettype=fasta&retmode=text")
    X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
    cat(X, sep = "\n", file = fl, append = TRUE)
    if (b == N) break
    a <- b + 1L
    b <- b + 400L
    if (b > N) b <- N
  }
  res <- read.FASTA(fl)
  if (is.null(res)) return(NULL)
  attr(res, "description") <- names(res)
  if (length(access.nb) != length(res)) {
    names(res) <- gsub("\\..*$", "", names(res))
    failed <- paste(access.nb[! access.nb %in% names(res)], collapse = ", ")
    warning(paste0("cannot get the following sequences:\n", failed))
  } else names(res) <- access.nb
  
  if (as.character) res <- as.character(res)
  
  if (species.names) {
    a <- 1L
    b <- if (N > 400) 400L else N
    sp <- character(0)
    repeat {
      URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                   paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", sep = "")
      X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE, n = -1)
      sp <- c(sp, gsub(" +ORGANISM +", "", grep("ORGANISM", X, value = TRUE)))
      if (b == N) break
      a <- b + 1L
      b <- b + 400L
      if (b > N) b <- N
    }
    attr(res, "species") <- gsub(" ", "_", sp)
  }
  
  if (host.names) {
    a <- 1L
    b <- if (N > 400) 400L else N
    hst <- character(0)
    repeat {
      URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                   paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", sep = "")
      X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE, n = -1)
      Y <- which(X=="//")
      Y2 <- numeric(length(Y)+1)
      Y2[1] <- 1
      Y2[2:length(Y2)] <- Y
      hst <- character(length(Y))
      for (i in 1:length(Y)) {
        tempgb <- X[Y2[i]:Y2[i+1]]
        hst2 <- paste(grep("ORIGIN", tempgb, value = TRUE),
                      grep("/host=", tempgb, value = TRUE))
        if (hst2=="ORIGIN       ") {
          hst2 <- "No_host_recorded"
        } else {
          hst2 <- gsub(" ", "_", hst2)
          hst2 <- gsub("ORIGIN____________________________/host=\"", "", hst2)
          hst2 <- gsub("\"", "", hst2)
        }
        hst[i] <- hst2
      }
      if (b == N) break
      a <- b + 1L
      b <- b + 400L
      if (b > N) b <- N
    }
    attr(res, "host") <- hst
  }
  res
}

# We can use the "read.GenBank2" function to download the SSU sequences, 
# species names, and host names directly from GenBank:

int <- seq(1, length(acc2), 299) # but, need to split the requests into intervals
int <- c(int, length(acc2))      # of ~300 so that the NCBI can handle it.

acc2gen <- read.GenBank(acc2[1])   # initialize acc2gen DNAbin,
spetemp <- character(length(acc2)) # species names vector,
hostemp <- character(length(acc2)) # host names vector,
namtemp <- character(length(acc2)) # and accession numbers vector.

for (i in 1:(length(int)-1)) {                         # for every interval,
  gentemp <- read.GenBank2(acc2[int[i]:int[i+1]],      # get info from GenBank and
                           species.names=T,            # store in a temporary object.
                           host.names=T)
  spetemp[int[i]:int[i+1]] <- attr(gentemp, "species") # store species names,
  hostemp[int[i]:int[i+1]] <- attr(gentemp, "host")    # host names,
  namtemp[int[i]:int[i+1]] <- names(gentemp)           # and accession numbers,
  acc2gen[int[i]:int[i+1]] <- gentemp                  # then copy sequences 
                                                       # to acc2gen.
}

lac <- length(acc2gen)                      # because attributes are lost every
names(acc2gen) <- namtemp[1:lac]            # time you modify the acc2gen DNAbin,
attr(acc2gen, "species") <- spetemp[1:lac]  # need to copy over the attributes once
attr(acc2gen, "host") <- hostemp[1:lac]     # we have all the 18S sequences stored
                                            # in acc2gen. That's why we had to store 
                                            # species names, host names and accession 
                                            # numbers in vectors as we iterated 
                                            # through our intervals.

# This is good! We have all the 18S sequences stored in a single list. Now, we
# need some code that can take the longest sequence for a single species, save
# that as a FASTA file, and delete the rest of the sequences for that particular
# species. Then, repeat for all other unique species in the list. All of the 
# resulting FASTA files will be stored into a folder called "FASTA_files", which
# we might need to create:

maindir <- getwd()
FASTAdir <- "./FASTA_files_methods_test"

if (dir.exists(FASTAdir)) {  # checks if a folder called "FASTA_files" exists
  setwd(maindir)             # in the working directory.
} else {
  setwd(maindir)
  dir.create(FASTAdir)       # creates "FASTA_files" in the working directory
}                            # if it doesn't exist.


# Awesome! Now we run the "for loop" that extracts and stores our desired 18S
# sequences as FASTA files, based on the GenBank information we imported as 
# the "acc2gen" list:

snames <- attr(acc2gen, "species") # get the species names of every list entry.
unames <- unique(snames)           # get the unique species names in the list. 

for (i in 1:length(unames)) {
  sname <- attr(acc2gen, "species")[i]            # get species name for a single list entry.
  hname <- attr(acc2gen, "host")[i]               # get host name for a single list entry
  iname <- which(attr(acc2gen, "species")==sname) # get which entries are same species.
  temp <- acc2gen[iname]                          # store these entries into temporary list "temp".
  
  f <- numeric(length=length(temp))  # initialize vector.
  for (n in 1:length(temp)) {        # for all entries in "temp", 
    f[n] <- (length(temp[[n]]))      # get the sequence length.
  }
  
  lmax <- max(f)        # get the length of the longest sequence.
  n0 <- which(f==lmax)  # get which "temp" sequence it was.
  n1 <- n0[1]           # if multiple sequences of the same length, just pick
                        # the first one.
  n2 <- iname[n1]       # get which "acc2gen" sequence it corresponds to.
  genseq <- acc2gen[n2] # get the corresponding Genbank DNAbin object.
  
  speheader <- gsub(".", "",     # remove periods and slashes from the file names since
                    sname,       # they can cause problems when specifying the file type.
                    fixed=TRUE)
  speheader <- gsub("/", "&",     
                    speheader,       
                    fixed=TRUE)
  hosheader <- gsub(".", "",     
                    hname,       
                    fixed=TRUE)
  hosheader <- gsub("/", "&",     
                    hosheader,       
                    fixed=TRUE)
  fasheader <- paste("Species:", # pastes the species and host information together
                     speheader,  # to be used in the FASTA file header.
                     "Host:", 
                     hosheader)
  
  setwd(FASTAdir)                     # open "FASTA_files".
  write.FASTA(genseq,                 # save the sequence as a FASTA file with
              speheader,              # the species name as the file name,
              fasheader)              # the species and host names as the header.
  setwd(maindir)                      # close "FASTA_files".
  
  if (length(iname)==1) { # if there is only one sequence for a given species,
    return                # no need to delete other sequences since they don't
                          # exist, so just return.
  } else {
    iname <- iname[-n1]        # removes the element corresponding to the longest sequence.
    
    spetemp <- spetemp[-iname] # removes elements from species, host and accessions vector 
    hostemp <- hostemp[-iname] # that don't correspond to the longest sequence.
    namtemp <- namtemp[-iname] 
    
    acc2gen <- acc2gen[-iname] # removes all acc2gen DNAbin entries that do not correspond
                               # to the longest sequence. Removes attributes in the process.
    
    lac <- length(acc2gen)                      # this string of code just adds back
    names(acc2gen) <- namtemp[1:lac]            # the missing attributes.
    attr(acc2gen, "species") <- spetemp[1:lac]  
    attr(acc2gen, "host") <- hostemp[1:lac]
  }
}

# Done! The resulting 18S FASTA sequences should be in the "FASTA_files" folder,
# which can be opened by running the following code:

open.dir <- function(dir = getwd()){            
  if (.Platform['OS.type'] == "windows"){
    shell.exec(dir)
  } else {
    system(paste(Sys.getenv("R_BROWSER"), dir))
  }
}

setwd(FASTAdir)
open.dir()
setwd(maindir)

rm(list = ls()) # clear environment

# -----------------------------------------------------------------------------

# Done!
