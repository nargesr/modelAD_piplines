## A function to GET all information from a given facet for a number of records. 
## facet specifies the facet in climb from which information is queried.
## queryValues is a character vector of (one or more) query values from one field.
## queryField is a single character object specifying the field in facet that matches the query values.
## Output is a data frame with each row corresponds to one query value.
## Query values that yield no data in the response are returned as an empty row 
## with only the query value in the query field.
## This function does not currently support other querying parameters 
## (e.g., AnimalNameSearchOptions) or combining fields.
## To build a more complex search use the climbRequest function directly

source("https://raw.githubusercontent.com/nargesr/modelAD_piplines/refs/heads/main/CLIMB/source/climbRequest.R")

climbGET <- function(queryValues, facet, queryField) {
  qv <- c()
  itemsL <- list()
  for (qi in 1:length(queryValues)) {
    qvi <- queryValues[qi]
    cat("getting data from",facet,"facet for",queryField,"==",qvi,"\n")
    qL <- list(qvi) ; names(qL) <- queryField
    resp <- climbRequest("GET", endpointPath=paste0("api/",facet), queryList=qL)
    if(length(resp$data)==0) {
      resp$data <- NA_character_
      qv <- c(qv, qvi) } else {qv <- c(qv, rep(qvi,nrow(resp$data)))}
    itemsL[[qi]] <- resp$data
  }
  df <- as.data.frame(do.call(rbind, itemsL))
  df[[queryField]] <- qv
  df[] <- lapply(df, as.character)
  return(df)
}


