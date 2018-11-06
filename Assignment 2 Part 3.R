#install.packages("tidyverse")
#install.packages("stringi")
#install.packages("rentrez")
#install.packages("seqinr")
#install.packages("Biostrings")
#install.packages("stringr")
#install.packages("strataG")
#install.packages("ape")
#source("https://bioconductor.org/biocLite.R")
#biocLite("DECIPHER")

library(stringi)
library(tidyverse)
library(rentrez)
library(seqinr)
library(Biostrings)
library(stringr)
library(strataG)
library(DECIPHER)
library(ape)

#2. Create a phylogenetic hypothesis for a taxonomic group of your choice.

#Dowload all Primate data from Bold
Primates.bold <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Primates&format=tsv")

#Filter for COI-5P marker and require neucleotides
Primates.bold.filtered <- Primates.bold %>%
  filter(str_detect(markercode, "COI-5P")) %>%
  filter(!is.na(nucleotides))

#Dataset can be smaller in this case since it is an exploratory study. Take 1 sample from each species
Primates.bold.sampled <- Primates.bold.filtered %>% 
  group_by(species_name) %>%
  sample_n(1)

#check the data for sequences much longer or shorter than the others, add lengths to the dataframe
sequences <- Primates.bold.sampled$nucleotides
sequenceLength <- lapply(sequences, nchar)
Primates.bold.sampled$sequence_length = NA
Primates.bold.sampled$sequence_length <- sequenceLength
rm(sequences, sequenceLength)
#large variation in sequence length in some cases, lots of sequences are short (600 aa's in comparison to others around 1500 aa's). May need to consider separationg

#create DNA Stringset for MSA
ss.Primates <- DNAStringSet(Primates.bold.sampled$nucleotides)
#names(ss.Primates) <- Primates.bold.sampled$genus_name #using genus instead of species here because species level information is too specific at this point

#make DNA Bin for mafft
dnaBin.Primates <- as.DNAbin(ss.Primates)

#run Mafft with 10 maximum iterations, create a fasta file of the sequences
mafft(dnaBin.Primates, maxiterate = 10, run.label = "Primates.MSA", opts = "--globalpair", delete.output = FALSE)

#read in MSA
Primates.MSA <- readDNAStringSet("Primates.MSA.aligned.fasta")

#check average length and number of gaps
mean(unlist(lapply(Primates.MSA, length)))
mean(unlist(lapply(Primates.MSA, str_count, ("-"))))
#mean length = 1635, mean number of gaps is 436.016. This should be acceptable considering the large variation in amino acid sequence legnth.

#Calculate distance matrix

#convert to a dna bin for clustering
dnaBIN.Primates.MSA <- as.DNAbin(Primates.MSA)

#calculate distance matrix
#Used TN93 method because it allows for different rates of transitions, transversions, and different rates of substitution. It also allows for the base pair frequencies to be highly variable (ref Jombart).
dMPrimates <- dist.dna(dnaBIN.Primates.MSA, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#Create a simple neighbour joining tree, only confirming that COI-5P can be used to accurately predict phylogeny at this point.
NJTreePrimates <- nj(dMPrimates)

#plot the tree
plot(NJTreePrimates, cex = .2)

#18 is an outlier, re-run without 18 and add Genus names for next plot
names(Primates.MSA) <- Primates.bold.sampled$genus_name
Primates.MSA <- Primates.MSA[-c(18), ]

#convert to a dna bin for clustering
dnaBIN.Primates.MSA.2 <- as.DNAbin(Primates.MSA)

#calculate distance matrix
#Used TN93 method because it allows for different rates of transitions, transversions, and different rates of substitution. It also allows for the base pair frequencies to be highly variable (ref Jombart).
dMPrimates.2 <- dist.dna(dnaBIN.Primates.MSA.2, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#Create a simple neighbour joining tree, only confirming that COI-5P can be used to accurately predict phylogeny at this point.
NJTreePrimates.2 <- nj(dMPrimates.2)

#plot the tree
plot(NJTreePrimates.2, cex = .2)
title("Primate Phylogenetic Hypothesis")

#OTUs were grouped according to expectations. So, the above method using COI-5P is clustering samples according to expectations from existing taxonomy.

#3: Seek Out Other Biological Data

#Pantheria is a database containing life history, ecology, and geography of a large number of species.

#Read in Pantheria from a tab deliminated text file
Pantheria <- read.table("PanTHERIA_1-0_WR05_Aug2008.txt", sep="\t", header=TRUE)

#Filter for primates and social group size greater than 0
Primates <- Pantheria %>%
  filter(str_detect(MSW05_Order, "Primates")) %>%
  filter(0 < X10.2_SocialGrpSize)

#4. Using the species identification, match up your sequence data with the other biological trait or geographical parameter you have chosen.

#Convert the Primates MSA into a dataframe
dfPrimates.MSA <- as.data.frame(Primates.MSA)

#Go back and take the outlier out of the origional dataframe so we can add the species to the MSA
Primates.bold.sampled <- Primates.bold.sampled[-c(18), ]

#add the species name line to MSA dataframe
dfPrimates.MSA$species_name <- Primates.bold.sampled$species_name

#add the nucleotids column for later use in phylogeny
dfPrimates.MSA$nucleotides <- Primates.bold.sampled$nucleotides

#merge the social group data into the dataframe
dfPrimates.MSA.Pantheria <- merge(x = dfPrimates.MSA, y = Primates, by.x = "species_name", by.y = "MSW05_Binomial")

#5. Prepare a visualization that maps your chosen biological trait or geographical parameter onto your phylogeny, such as by using maximum parsimony or maximum likelihood character mapping.

#Reanalyze with new subset

#Too much data to make a trr that can be analyzed, sample one of each Genus
dfPrimates.Genus <- dfPrimates.MSA.Pantheria %>% 
  group_by(MSW05_Genus) %>%
  sample_n(1)

#convert MSA to DNAStringset for clustering
ssPrimates.Social <- DNAStringSet(dfPrimates.Genus$x)
names(ssPrimates.Social) <- dfPrimates.Genus$species_name

dnaBIN.Primates.Social <- as.DNAbin(ssPrimates.Social)

#calculate distance matrix
dMPrimates.Social <- dist.dna(dnaBIN.Primates.Social, model = "raw")

#calculate initial neighbour joining tree
NJTreePrimates.Social <- nj(dMPrimates.Social)
NJTreePrimates.Social <- ladderize(NJTreePrimates.Social)

#create DNA set for Maximum parsimony reconstruction
Primate.dna <- as.phyDat(dfPrimates.MSA.Pantheria$x)

#create a tree by maximum parsimony
OptPars.Primates.Social <- optim.parsimony(NJTreePrimates.Social, Primate.dna)
#failed with error: Error in rep(1:nr, attr(x, "weight")) : invalid 'times' argument

#plot tree with Speices
plot(NJTreePrimates.Social, cex = .3)
title("Neighbour Joining Tree of Primates with Social Group Information")

#plot tree with Species and Social Group Size
plot(NJTreePrimates.Social, show.tip = FALSE)
title("Primate Phylogeny and Social Group Size")

#Visualize Social Group Data by Colour Scale: red -> blue = small -> large groups
myPalette <- colorRampPalette(c("red","yellow","green","blue"))
tiplabels(dfPrimates.Genus$X10.2_SocialGrpSize, bg=num2col(dfPrimates.Genus$X10.2_SocialGrpSize, col.pal=myPalette), cex=.3, width = .1, height = .1)

