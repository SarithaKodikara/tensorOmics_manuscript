library(readr)
library(magrittr)
library(CoDaSeq)
library(stringr)
library(dplyr)
library(mixOmics)
library(tidyverse)
library(RColorBrewer)
library(zCompositions)

# Function to remove low abundance features
low.removal <- function(
    data, # df of size n (sample) x p (features)
    percent=0.01 # cutoff chosen
) 
{
  keep.features = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.features]
  return(list(data.filter = data.filter, keep.features = keep.features))
}

#### Reading Data ####
# A supernatant is a liquid or medium which remains above a pellet after 
# centrifugation and is composed of lighter or smaller materials.
# More accurate to predict the pellet based on the supernatant rather 
# than the other way around
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4313957/
# https://pubmed.ncbi.nlm.nih.gov/28667629/

LCMS_supernatant <- read.csv('Data_org/supernant/LCMS_xcms.csv', row.names = 1) 
metaData_supernatant <- read.csv('Data_org/supernant/metadata.csv', row.names = 1)
LCMS_pellet <- read.csv('Data_org/pellet/LCMS_xcms.csv', row.names = 1) 
metaData_pellet <- read.csv('Data_org/pellet/metadata.csv', row.names = 1)

# Add row IDs to metadata
metaData_supernatant$id <- rownames(metaData_supernatant)
metaData_pellet$id <- rownames(metaData_pellet)

# Get metadata ids in both files using inner join on id
ids <- inner_join(metaData_supernatant, metaData_pellet, by = join_by("id"))%>%
  select(id)%>%unlist()
# metaData_filtered to only those ids
metaData <- metaData_supernatant[ids, ]


# Filter metadata for experimental groups
meta_filtered <- filter(metaData, 
                        group %in% c("Control", "Carbonate",
                                       "Chloride", "Phosphate"))
meta_filtered$number_of_days<-as.numeric(meta_filtered$number_of_days)

# Filter metadata for blanks
meta_filtered_bl <- filter(metaData, 
                           !(group %in% c("Control", "Carbonate",
                                            "Chloride", "Phosphate")))


# Calculate median values from blanks
LCMS_supernatant_bl <- apply(LCMS_supernatant[meta_filtered_bl$id, ], 2, median)
LCMS_pellet_bl <- apply(LCMS_pellet[meta_filtered_bl$id, ], 2, median)

# Filter LCMS data for experimental samples
LCMS_supernatant_fil_1 <- LCMS_supernatant[meta_filtered$id, ]
LCMS_pellet_fil_1 <- LCMS_pellet[meta_filtered$id, ]

# Impute zeros with blank medians
LCMS_supernatant_fil_imp <- LCMS_supernatant_fil_1 %>%
  mutate(across(
    names(LCMS_supernatant_bl),
    ~ ifelse(. == 0, LCMS_supernatant_bl[cur_column()], .)
  ))

LCMS_pellet_fil_imp <- LCMS_pellet_fil_1 %>%
  mutate(across(
    names(LCMS_pellet_bl),
    ~ ifelse(. == 0, LCMS_pellet_bl[cur_column()], .)
  ))

# Remove low abundance features
LCMS_supernatant_fil2_imp <- low.removal(LCMS_supernatant_fil_imp, percent = 0.01)$data.filter
LCMS_pellet_fil2_imp <- low.removal(LCMS_pellet_fil_imp, percent = 0.01)$data.filter

# Log transformation
LCMS_pellet_log <- log(LCMS_pellet_fil2_imp + 1)
LCMS_supernatant_log <- log(LCMS_supernatant_fil2_imp + 1)


save(
     LCMS_supernatant_log, LCMS_pellet_log,
     meta_filtered,
     file = "Data_filtered/anaerobicStudy.RData")
