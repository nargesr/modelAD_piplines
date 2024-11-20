## Request data from CLIMB
## Add necessary information

library(tidyverse)
source("https://raw.githubusercontent.com/nargesr/modelAD_piplines/refs/heads/main/CLIMB/source/climbGET.R")

createIndividualMetadata <- function(animal_IDs=NULL) {
  
  gts <- climbGET(ID, facet="genotypes", "animalID") %>%
    dplyr::select(climbID="animalID", genotype)
  
  
  ## Pull animal metadata from climb based on climbID
  inds <- climbGET(ID, facet="animals", "animalID") %>%
    dplyr::select(individualID="animalName", 
                  climbID="animalID", 
                  animalName, 
                  birthID="birthId", 
                  matingID="parentMatingId", 
                  sex, 
                  species,
                  generation,
                  dateBirth="dateBorn",
                  stockNumber="line")
  
  CLIMB_df = merge(inds, gts, by = "climbID")
  
  CLIMB_df['individualIdSource'] = 'UCI_TMF'
  CLIMB_df['materialOrigin'] = 'UCI_TMF'
  CLIMB_df['ageDeathUnits'] = 'months'
  CLIMB_df['rodentDiet'] = '7% standard'
  CLIMB_df['genotypeBackground'] = 'C57BL6J'
  
  CLIMB_df['microchipID'] = ''
  CLIMB_df['ageDeath'] = ''
  CLIMB_df['ageDeathUnits'] = ''
  CLIMB_df['brainWeight'] = ''
  CLIMB_df['rodentWeight'] = ''
  CLIMB_df['bedding'] = ''
  CLIMB_df['room'] = ''
  CLIMB_df['waterpH'] = ''
  CLIMB_df['treatmentDose'] = ''
  CLIMB_df['treatmentType'] = ''
  CLIMB_df['genotypeBackground'] = ''
  CLIMB_df['modelSystemName'] = ''
  CLIMB_df['officialName'] = ''
  
  
  CLIMB_df <- CLIMB_df[, c("individualID" ,
                           "climbID",
                           "microchipID",
                           "birthID",
                           "matingID",
                           "individualIdSource",
                           "materialOrigin",
                           "sex",
                           "species",
                           "generation",
                           "dateBirth",
                           "ageDeath",
                           "ageDeathUnits",
                           "brainWeight",
                           "rodentWeight",
                           "rodentDiet",
                           "bedding",
                           "room",
                           "waterpH",
                           "treatmentDose",
                           "treatmentType",
                           "stockNumber",
                           "genotype",
                           "genotypeBackground",
                           #"individualCommonGenotype"
                           "modelSystemName",
                           "officialName")]
  
  return(CLIMB_df)
}

