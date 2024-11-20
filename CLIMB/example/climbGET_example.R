library(tidyverse)
source("https://raw.githubusercontent.com/nargesr/modelAD_piplines/refs/heads/main/CLIMB/source/createIndividualMetadata.R")
#source("../source/createIndividualMetadata.R")

ID = c("13436",
            "13438",
            "13426",
            "13427",
            "13430",
            "13431",
            "13424",
            "13421",
            "13422",
            "13416")

CLIMB_df = createIndividualMetadata(animal_IDs=ID)

write.csv(CLIMB_df, "CLIMB.csv")

