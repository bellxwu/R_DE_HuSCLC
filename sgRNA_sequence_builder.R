# Description: Automated script to nominate common exons between transcripts directly from ensembl API.
# Goal is to identify the common sequences within exons. 
# Date created: 2025.09.14
# Author: Bell Wu

# Criteria for selecting exons (ordered in priority):
# 1. Exons common to all transcripts
# 2. Early exons
# 3. Asymmetrical exons
# 4. Not exon 1

# Changelog:
# 2025.09.29 - functions for retrieving exons (get_early_freq_exons) converted into GetSameExons
# rewrite functions using GetSameExons

library(httr2)
library(jsonlite)
library(tidyverse)
library(rlang)

source("~/R_programming/R_analyses_scripts/GetSameExons.R") # new function for retrieving 

# 1.0 request data from REST API ------------------------------------------------
req = request("http://rest.ensembl.org")

# 1.1 build GET wrapper for ensembl ---------------------------------------------
# path = path to URL
get_ensembl = function(path, query = list()) {
  request("http://rest.ensembl.org") |> # set the main URL to request
    req_url_path(path) |> # use path to alter the path of the main URL
    req_headers("Accept" = "application/json") |> # accept only JSON files
    req_url_query(!!!query) |> # modify query components for GET URL
    req_perform() |>  # perform request
    resp_body_json(simplifyVector = TRUE)
}

# 1.2 create function to GET exons ---------------------------------------------
get_ensembl_exon = function(gene_id) {
  get_ensembl(paste0("overlap/id/", gene_id),
              query)
}

# 1.3 flag common exons of all transcripts --------------------------------------
tx_exons  <- split(exon_df$id, exon_df$Parent) # split by transcripts
common_id <- Reduce(intersect, tx_exons) # find the exon present in all IDs 

# create a a table of counts for each exons in transcript
t = table(unlist(tx_exons)) |> 
  data.frame() 
# 1.3.1 for single exon --------------------------------------------------------
# find the most frequent exon
exon_t = t |> 
  filter(Freq == max(t$Freq))
# find exon transcripts present in over 20 transcripts
# select one of the most common exons
exon_id = exon_t[1,1]
# function to look up exon transcript 
get_ensembl_exon_info = function(exon_id) {
  get_ensembl(path = paste0("lookup/id/", exon_id),
              query = list("expand" = 0))
}
# find the exon information
get_ensembl_exon_info(exon_id)

# 1.3.2 for multiple exons -----------------------------------------------------
exon_t = t |> 
  filter(Freq > 20)
# select the list of exon transcripts
exon_ids = exon_t$Var1
# function to look up multiple exons 
get_multiple_exon = function(path, exon_ids) {
  request("http://rest.ensembl.org") |> # set the main URL to request
    req_url_path(path) |> # use path to alter the path of the main URL
    req_headers("Accept" = "application/json") |> # accept only JSON files
    req_body_json(list("ids" = exon_ids)) |> # modify query components for GET URL
    req_perform() |>  # perform request
    resp_body_json(simplifyVector = TRUE)
}
# arguments for multiple exon
path = "lookup/id"
query = list("id" = exon)
# obtain list of multiple exons for gene
exons = get_multiple_exon(path = path, exon_ids = exon_ids) # test to get multiple exon data from IDs
exon_starts = lapply(exons, function(x) x[["start"]]) # function to select the same ids
# find earliest start (criteria 2)
early_exon = unlist(exon_starts)
early_exon = early_exon[order(early_exon)] # order by ascending order of exons

# 2.0 putting it all together for a function ------------------------------------
get_early_freq_exons = function(gene_id) {
  get_ensembl_exon = function(gene_id) {
    query = list("feature" = "exon")
    get_ensembl(paste0("overlap/id/", gene_id),
                query)
  }
  exon_df = get_ensembl_exon(gene_id = gene_id)
  tx_exons  <- split(exon_df$id, exon_df$Parent) # split by transcripts
  common_id <- Reduce(intersect, tx_exons) # find the exon present in all IDs 
  if (length(common_id) == 0) {
    # unlist and table how exon transcript appearance
    t = table(unlist(tx_exons)) |> 
      data.frame() 
    
    # select exon transcripts in top 10% freq
    exon_t = t |> 
      dplyr::filter(Freq >= quantile(t$Freq, 0.9)) # find above 90% quantile
    # create exon_ids
    exon_ids = exon_t$Var1
    colnames(exon_t)[colnames(exon_t) == "Var1"] = "exon_id" # rename column 
    # pull exon data from ensembl REST
    exons = request("http://rest.ensembl.org") |> # set the main URL to request
        req_url_path("lookup/id") |> # use path to alter the path of the main URL
        req_headers("Accept" = "application/json") |> # accept only JSON files
        req_body_json(list("ids" = exon_ids)) |> # modify query components for GET URL
        req_perform() |>  # perform request
        resp_body_json(simplifyVector = TRUE)
    
    # find exon starts to determine early exons
    exon_starts = lapply(exons, function(x) x[["start"]])
    early_exon = unlist(exon_starts)
    # add exon strand information to table
    # create a df for the exon id with strands
    strand_ids = data.frame(exon_id = exon_df$exon_id,
                            strands = exon_df$strand)
    strand_ids = strand_ids[!duplicated(strand_ids$exon_id), ] # remove duplicates
    # match the ids to the tabled results
    exon_t = dplyr::inner_join(exon_t, strand_ids, by = "exon_id") # join the dfs 
    # add exon starts to the table
    exon_t$Start = early_exon[match(exon_t$exon_id, names(early_exon))]
    # if there are two strand indicators, remove the transcripts of fewest strands
    count_strands = exon_t |> 
      dplyr::count(strands) |> 
      dplyr::filter(n == max(n))
    # remove fewest
    exon_t = exon_t |> 
      dplyr::filter(strands == count_strands$strands)
  
    # create if-else statement for if -1 or +1 strand
    if (unique(exon_t$strands == -1)) { # if reverse strand
    exon_t = exon_t[order(exon_t$Start, decreasing = TRUE), ] # order by ascending order of exons
    } else {
      exon_t = exon_t[order(exon_t$Start, decreasing = FALSE), ]
    }
    return(exon_t)
    
  } else { 
    return(common_id)
  }
}


# 3.0 finding exons for GENE1 --------------------------------------------------
path = "xrefs/symbol/homo_sapiens/GENE1"
GENE1 = get_ensembl(path)
gene_id = GENE1$id

GENE1 = GetSameExons("GENE1")

# 4.0 finding exons for GENE2 ---------------------------------------------------
path = "xrefs/symbol/homo_sapiens/GENE2"
GENE2 = get_ensembl(path)
gene_id = GENE2$id

GENE2 = GetSameExons("GENE2")
 
# 5.0 finding exons for GENE3 --------------------------------------------------
path = "xrefs/symbol/homo_sapiens/GENE3"
GENE3 = get_ensembl(path)
gene_id = GENE3$id

GENE3 = GetSameExons("GENE3")

# 6.0 finding exons for GENE4 --------------------------------------------------
path = "xrefs/symbol/homo_sapiens/PECAM1"
GENE4 = get_ensembl(path)
gene_id = GENE4$id

GENE4 = GetSameExons("GENE4")


