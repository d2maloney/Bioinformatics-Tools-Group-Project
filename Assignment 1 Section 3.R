# This script was written by Daniel Maloney and was edited by Cailin Harris

# Question: Given a species of Echinodermata, does it form a community with a specific group of other Echinodermata species?

# Echinodermata is a phylum which includes a diverse group of species such as: starfish, sea urchins, and sea cucumbers (ref1).
# Some of these species have become increasingly imfamous due to their effect on the environment. 
# For example, the Crown of Thorns sea star is being studied due to its involvment in reef degradation (ref2).
# They are often found within complex communities.
# The purpose of this project is to discover if data from the BOLD database can be used to predict the communities these species form.
# An understanding of the dependencies between Echinodermata species will help in understanding how invasive Echinodermata species can be controlled.

# Reference 1: https://authors.library.caltech.edu/35244/
# Reference 2: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0047363

#install.packages("tidyverse")
#install.packages("stringi")
#install.packages("vegan")
#install.packages("iNEXT")
#install.packages("gridExtra")

# Loading the libraries
library(tidyverse)
library(vegan)
library(stringi)
library(iNEXT)
library(gridExtra)

# Download a database of Echinodermata
Echi <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Echinodermata&format=tsv")

#Basic Filters, bin uri are used to represent species, regions a being used to represent a communities
Echi.filtered <- Echi %>%
  filter(str_detect(bin_uri, "[:]")) %>%
  filter(!is.na(region))

#Find the species with the most records
Echi.top.bins <- Echi.filtered %>%
  group_by(bin_uri) %>%
  summarize(count = length(unique(processid))) %>%
  arrange(desc(count))

#At least one common species must be found in order to calculate a distance.
#Filter dataset for regions that contain the species with the most records: BOLD:ACF3333 to avoid missing values when calculating distance
#Compare the top 5 regions with the most records
Echi.with.control <- Echi.filtered %>%
  group_by(region) %>%
  mutate(count = sum(str_detect(bin_uri, "BOLD:ACF3333"), na.rm = TRUE)) %>%
  filter(count > 20)

#Generate a histogram showing the community of Echinodermata found in the same regions as BOLD:ACF3333
ggplot(Echi.with.control) + 
  geom_bar(mapping = aes(x = region, fill = bin_uri))

#Find the distance between these populations in terms of community diversity

#Create a community object from the dataset for use in Vegan
comm.Echi <- Echi.with.control %>%
  group_by(region, bin_uri) %>%
  count(bin_uri) %>%
  spread(bin_uri, n)

comm.Echi <- comm.Echi %>%
  remove_rownames %>%
  column_to_rownames(var="region")

#calculate distances between communities based on the number of collected samples of each species in each region 
Echi.dis <- vegdist(comm.Echi)

#perform clustering analysis to determine which communities are the most similar
clus <- hclust(Echi.dis, "single")
plot(clus)

####use iNext to create a rarification curve to estimate species diversity ----

#turn comm.Echi into a data frame
comm.Echi.dataFrame <- as.data.frame(comm.Echi)

#transpose the data so that the regions are now the columns and the accession numbers are the rows
comm.Echi.dataFrame2 <- t(comm.Echi.dataFrame)

#check the class to ensure it is in the correct format (matrix)
class(comm.Echi.dataFrame2)

#rename the regions so there are no spaces in the names
colnames(dataframe.test) <- c("Coats.Land", "Scotia.Arc", "Scotia.Island", "Terres.Australes", "Wilhelm.Land")

#separate them into new dataframes.
coats.land.data <- subset(dataframe.test, select = 1)
scotia.arc.data <- subset(dataframe.test, select = 2)
scotia.island.data <- subset(dataframe.test, select = 3)
terres.austales.data <- subset(dataframe.test, select = 4)
wilhelm.land.data <- subset(dataframe.test, select = 5)

# Now remove all of the missing data.
coats.land.data2 <- na.omit(coats.land.data)
scotia.arc.data2 <- na.omit(scotia.arc.data)
scotia.island.data2 <- na.omit(scotia.island.data)
terres.austales.data2 <- na.omit(terres.austales.data)
wilhelm.land.data2 <- na.omit(wilhelm.land.data)

#clean up the workspace
rm(coats.land.data)
rm(scotia.arc.data)
rm(scotia.island.data)
rm(terres.austales.data)
rm(wilhelm.land.data)

#create iNEXT objects for all of my regions.
coats.land.iNext <- iNEXT(coats.land.data2)
scotia.arc.iNext <- iNEXT(scotia.arc.data2)
scotia.island.iNext <- iNEXT(scotia.island.data2)
terres.australes.iNext <- iNEXT(terres.austales.data2)
wilhelm.land.iNext <- iNEXT(wilhelm.land.data2)

#create the rarefaction curves
p1 <- ggiNEXT(coats.land.iNext)
p2 <- ggiNEXT(scotia.arc.iNext)
p3 <- ggiNEXT(scotia.island.iNext)
p4 <- ggiNEXT(terres.australes.iNext)
p5 <- ggiNEXT(wilhelm.land.iNext)

#format plots onto one frame so comparison
grid.arrange(p1, p2, p3, p4, p5, nrow = 3)

####determining simpson's diversity for each region ----
ChaoSimpson(coats.land.data2, datatype = "abundance", B=200)
ChaoSimpson(scotia.arc.data2, datatype = "abundance", B=200)
ChaoSimpson(scotia.island.data2, datatype = "abundance", B=200)
ChaoSimpson(terres.austales.data2, datatype = "abundance", B=200)
ChaoSimpson(wilhelm.land.data2, datatype = "abundance", B=200)

####conclusions ----

# Using data associated with bin numbers from the BOLD database, a species of Echinodermata was selected and a group of species that it consistently forms communities with could be found.
# The species with a large number of records was selected, Promachocrinus kerguelensis (BOLD:ACF3333). The regions where it was collected most were determined.
# Over the five regions of Antarctica where P. kerguelensis was collected most, it was consistently found in communities containing only a small group of other Echinodermata (figure 1).
# A cluster dendrogram was then created to determine which regions had the most similar communities using single linkage clustering (figure 2).
# Sample size based rarification curves were created based on the number of records collected in each region (figure 3).
# A strong correltation was observed between the records collected from the species with these BOLD numbers in these communities.
# A concern is that several of these BOLD bin uri's were linked to the same species including: ACF3333, AAA0602, and ABZ8776.
# Further investigation may prove that these are truly different species.
# Terres Australes Francaises has the highest diversity (0.808), while Wilhelm II Land had the lowest diversity of 0.35.
# 
