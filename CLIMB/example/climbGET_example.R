library(tidyverse)
source("../source/climbGET.R")

fqf_key <- read_csv("fastqfiles_key.csv")

## Pull animal metadata from climb based on climbID
inds <- climbGET(fqf_key$climbID, facet="animals", "animalID") %>%
  dplyr::select(animalName, "climbID"=animalId, sex, line, dateBorn, origin)

## Pull genotypes from climb based on climbID
gts <- climbGET(fqf_key$climbID, facet="genotypes", "animalID") %>%
  dplyr::select("climbID"=animalID, genotype)
