# ClimbR (Adopted from [this repo](https://github.com/TheJacksonLaboratory/ClimbR/tree/master)

R API client for [Climb.bio](https://api.climb.bio/docs/index.html).  
Still in development; not fully tested.  
Based on the [httr](https://CRAN.R-project.org/package=httr) package.  
climbRequest.R is a generic function, enabling any request with any query parameters.   
climbGET.R is a wrapper for climbRequest, to get all information available from a given facet based on one query field for a vector of corresponding query values.    
Authentication is done by retrieving a temporary token (valid for 60 minutes) with getToken.R  
See annotations in the scripts for details.  

