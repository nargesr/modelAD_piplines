## A generic API request following https://api.climb.bio/docs/index.html
## Authentication is done by retrieving a temporary token with getToken.R
## queries are provided as a list of key-value pairs,
## e.g, list(AnimalName="50101", AnimalNameSearchOptions="StartsWith")
## output is a climb_api object including the full response as returned by package 
## httr (https://httr.r-lib.org/reference/response.html) along with the following objects:
## parsed data is returned as dataframe with records in rows and fields in columns
## url, errors, and status code returned as strings

source("https://raw.github.com/TheJacksonLaboratory/ClimbR/master/getToken.R")
climbRequest <- function(method, endpointPath, queryList=NULL) {
  
  # check if there is a token in the environment
  istoken <- class(try(token, silent=TRUE))
  # if there is one already, check that it is valid 
  test <- 0
  if (istoken=='character') {
    test <- GET("https://api.climb.bio/api/Diagnostics", 
                add_headers(Authorization = token))$status_code}
  # if there isn't any token, or it isn't valid, get a new one
  if (istoken=='try-error' | test==401) getToken()

  # build url including endpoint path and queries
  url <- modify_url("https://api.climb.bio/", path=endpointPath, query=queryList)
  
  # send request
  resp <- VERB(method, url, add_headers(.headers = c(Authorization = token)))
  resp[[9]] <- NULL # remove response element that has credentials
  
  # parse response content 
  parsed <- jsonlite::fromJSON(content(resp, "text"), simplifyVector = FALSE)
  data <- parsed$data
  # make data structure consistent across unique vs paged lists responses
  if(is.null(data$totalItemCount)) {
    itemsL <- list(data)
  } else {itemsL <- data$items}
  # turn data object into dataframe, all strings, NA for missing values
  items <- lapply(itemsL, function(it){
    sapply(it, function(x){ifelse(is.null(x), NA_character_, as.character(x))})
  })
  items <- as.data.frame(do.call(rbind, items))

  structure(
    list(
      status_code = resp$status_code,
      url = resp$url,
      error = parsed$errors,
      data = items,
      response = resp
    ),
    class = "climb_api"
  )
}
