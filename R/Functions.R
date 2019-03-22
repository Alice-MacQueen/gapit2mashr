require(tidyverse)
require(purrr)
require(rlang)

#' Identify phenotype names from GAPIT results in a folder.
#'
#'
#'
#' @param path File path to the csv files that GAPIT has created.
#' @param model Model type used in GAPIT runs.
gapit_phenotypes_in_folder <- function(path = ".", model = "CMLM"){
  result_files <- list.files(path = path, pattern = "*GWAS.Results*")
  if(model == "CMLM"){
  phemiddle <- purrr::partial(str_sub, start = 12, end = -18)
  # Eventually this needs to accomodate models other than "CMLM", which means
  # modifying where this function ends based on the model used.
  gapit_phenotypes <- purrr::map(result_files, phemiddle) %>%
    unlist()
  } else stop("Sorry, models other than CMLM are not currently supported.")
}

#' Load GAPIT Results Tables
#'
#'
#'
#' @param path File path to the csv files that GAPIT has created.
#' @param phenotype Phenotype name used in GAPIT.
#' @param model Model type used in GAPIT runs.
#'
#' @note To create a vector of phenotype names, use the
#' \code{gapit_phenotypes_in_folder} function.
load_GAPIT_GWAS_all <- function(path, phenotype, model = "CMLM"){
  out <- list()
  out$Pred <- read_csv(paste0(path, "GAPIT.", model, ".", phenotype,
                              ".PRED.csv"))
  out$Results <- read_csv(paste0(path, "GAPIT.", model, ".", phenotype,
                                 ".GWAS.Results.csv"), col_names = TRUE,
                          col_types = "ciinninnnn")
  out$Effects <- read_csv(paste0(path, "GAPIT.", model, ".", phenotype,
                                 ".Df.tValue.StdErr.csv"), col_names = TRUE,
                          col_types = "ciiinnn")
  out$ROC <- read_csv(paste0(path, "GAPIT.", model, ".", phenotype, ".ROC.csv"))
  out$Log <- read_csv(paste0(path, "GAPIT.", model, ".", phenotype, ".Log.csv"))
  return(out)
}



