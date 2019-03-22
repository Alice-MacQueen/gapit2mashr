require(tidyverse)
require(purrr)
require(rlang)

#' Identify phenotype names from GAPIT results in a folder.
#'
#'
#'
#' @param path File path to the csv files that GAPIT has created, a character
#' string.
#' @param model Model type used in GAPIT runs, a character string.
gapit_phenotypes_in_folder <- function(path = ".", model = "CMLM"){
  result_files <- list.files(path = path, pattern = "*GWAS.Results*")
  if(model == "CMLM"){
  phemiddle <- purrr::partial(str_sub, start = 12, end = -18)
  # Eventually this needs to accomodate models other than "CMLM", which means
  # modifying where this function ends based on the model used.
  gapit_phenotypes <- purrr::map(result_files, phemiddle) %>%
    unlist()
  } else stop("Currently, models other than CMLM are not supported.")
}

#' Load GAPIT Results Tables
#'
#'
#'
#' @param path File path to the csv files that GAPIT has created, a character
#' string.
#' @param phenotype Phenotype name used in GAPIT, a character string.
#' @param model Model type used in GAPIT runs, a character string.
#'
#' @note To create a vector of phenotype names, use the
#' \code{gapit_phenotypes_in_folder} function.
load_GAPIT_GWAS_all <- function(path = ".", phenotype, model = "CMLM"){
  # Probably want to test this to make sure path works with file.path and with
  # and without a trailing "/".
  out <- list()
  out$Pred <- readr::read_csv(file.path(path, paste0("GAPIT.", model, ".", phenotype,
                              ".PRED.csv")))
  out$Results <- readr::read_csv(file.path(path, paste0("GAPIT.", model, ".", phenotype,
                                 ".GWAS.Results.csv")), col_names = TRUE,
                          col_types = "ciinninnnn")
  out$Effects <- readr::read_csv(file.path(path, paste0("GAPIT.", model, ".", phenotype,
                                 ".Df.tValue.StdErr.csv")), col_names = TRUE,
                          col_types = "ciiinnn")
  out$ROC <- readr::read_csv(file.path(path, paste0("GAPIT.", model, ".", phenotype,
                                    ".ROC.csv")))
  out$Log <- readr::read_csv(file.path(path, paste0("GAPIT.", model, ".", phenotype,
                                    ".Log.csv")))
  return(out)
}

gapit_top_effects <- function(df = df1$Effects, phenotype, numSNPs = numSNPs){
df2 <- df %>%
  mutate(abs_tvalue = abs(`t Value`)) %>%
  top_n(as.integer(numSNPs), abs_tvalue)
names(df2)[4] <- paste0(phenotype, "_DF")
names(df2)[5] <- paste0(phenotype, "_tvalue")
names(df2)[6] <- paste0(phenotype, "_stderror")
names(df2)[7] <- paste0(phenotype, "_effect")
return(df2)
}

#' Convert GAPIT output to mashr input dataframes.
#'
#'
#'
#' @param path File path to the csv files that GAPIT has created, a character
#' string. Defaults to the working directory.
#' @param phenotypes A character vector of phenotype names used in GAPIT. If NA,
#' will find these from the GAPIT Results files found in the path.
#' @param numSNPs The number of most significant SNPs selected from each GWAS.
#' Ideally this will give 1 million or fewer total cells in the resultant mash
#' dataframes. Defaults to 1000. For many users this will be far too few.
#' @param model Model type used in GAPIT runs, a character string. Defaults to
#' "CMLM".
#' @param S_hat One of \code{c("Hedge's G", "ones")}. If too many standard
#' errors have not been estimated in GAPIT, how should NA's be replaced? The
#' default is to estimate standard errors using Hedges' G, which...
#' (Hedges, Olkin 1985).
#' @param saveoutput Logical. Should the function's output also be saved to RDS
#' files?
#'
#' @note Hedges' g (Hedges & Olkin 1985 p. 86) is used here to calculate S_hat,
#' or the standard error in the effect size difference between the reference and
#' alternate allele, because it allows the calculation of both the effect size
#' of the alternate allele, and the confidence interval around the effect size.
#' This function uses the effect sizes provided by GAPIT to compute the
#' confidence interval calculation. The calculations are:
#' \code{d = 2r/sqrt(1-r^2)}
#' \code{d_unbiased = (1-(3/(4*(N-2)-1)))*d}
#' \code{sigma^2_d_i = (n_i^e + n_i^c)/n_i^e*n_i^c + d_i^2 / 2*(n_i^e + n_i^c)}
#' where r is the effect size, scaled between -1 and 1; n's are the sample sizes
#' of the two experimental groups; N is the total sample size.
gapit2mashr <- function(path = ".", phenotypes = NA, numSNPs = 1000,
                        model = "CMLM", S_hat = c("Hedge's G", "ones"),
                        saveoutput = FALSE){
  match.arg(stderr, c("Hedge's G", "ones")) # Fix this or make sure it works...
  if(is.na(phenotypes)){
    phe_col <- gapit_phenotypes_in_folder(path = path)
  } else {
    phe_col <- phenotypes
  }
  if(is.na(phe_col)){
    stop("Can't find any GAPIT Results files in this path.")
  }

  message(paste0("Starting part one: Making a data frame of all SNPs that are",
                 " in the top ", numSNPs, " SNPs for at least one phenotype."))

  df1 <- load_GAPIT_GWAS_all(path = path, phenotype = phe_col[1])
  df2 <- df1$Effects %>%  # fix this: make it a not-exported function.
    mutate(abs_tvalue = abs(`t Value`)) %>%
    top_n(as.integer(numSNPs), abs_tvalue)
  names(df2)[4] <- paste0(phe_col[1], "_DF")
  names(df2)[5] <- paste0(phe_col[1], "_tvalue")
  names(df2)[6] <- paste0(phe_col[1], "_stderror")
  names(df2)[7] <- paste0(phe_col[1], "_effect")
  big_effects_df <- df2 %>%
    dplyr::select(-abs_tvalue)

  for(i in seq_along(phe_col)[-1]){
    df1 <- load_GAPIT_GWAS_all(path = path, phenotype = phe_col[i])
    df2 <- df1$Effects %>%
      mutate(abs_tvalue = abs(`t Value`)) %>%
      top_n(as.integer(numSNPs), abs_tvalue)
    names(df2)[4] <- paste0(phe_col[i], "_DF")
    names(df2)[5] <- paste0(phe_col[i], "_tvalue")
    names(df2)[6] <- paste0(phe_col[i], "_stderror")
    names(df2)[7] <- paste0(phe_col[i], "_effect")
    big_effects_df <- df2 %>%
      dplyr::select(-abs_tvalue) %>%
      full_join(big_effects_df)
  }

  if(saveoutput == TRUE){
  saveRDS(big_effects_df, file = file.path(path, paste0("effects_", numSNPs,
                                        "SNPs_PartOneOutput.rds")))
  }

  message(paste0("Part One: data frame of SNPs to keep complete. Starting Part",
          " Two: Creating strong and random dataframes of B_hat and S_hat",
          " values for use in mashr."))

  df1 <- load_GAPIT_GWAS_all(path = "./", phenotype = phe_col[1])

  if(S_hat == "Hedge's G"){ # fix this: need support for "ones"
    if(sum(is.na(df1$Effects$`std Error`)) > length(df1$Effects$`std Error`)*.05){
      # If there are too many NA's for standard errors, derive new standard errors
      # using Hedge's G (which requires the MAF).
      standardization <- max(abs(df1$Results$effect), na.rm = TRUE)

      df3 <- df1$Results %>% # fix this: make it a not-exported function.
        mutate(Stand_effect = effect / standardization,
               Obs = maf * nobs,
               Obs2 = (1-maf) * nobs,
               d = ifelse(abs(Stand_effect) < 0.98,
                          (2 * Stand_effect) / sqrt(1 - Stand_effect^2),
                          4),
               d_unbiased = (1 - (3 / (4 * (nobs -2) -1))) * d,
               sigma2_d = ((Obs + Obs2) / (Obs * Obs2)) +
                 (d_unbiased^2 / (2*(Obs + Obs2))),
               stderr_d = sqrt(sigma2_d)) %>%
        mutate(Stand_effect = ifelse(is.na(Stand_effect) |
                                       is.infinite(Stand_effect),
                                     0,
                                     Stand_effect),
               stderr_d = ifelse(is.na(stderr_d) | is.infinite(stderr_d),
                                 10,
                                 stderr_d))  %>%
        dplyr::select(SNP, Stand_effect, stderr_d)
     names(df3)[2] <- paste0("Bhat_", phe_col[1])
      names(df3)[3] <- paste0("Shat_", phe_col[1])
      message(paste0("Hedge's G standard errors were used for", phe_col[1]))
  } else {
    # or if not many of the standard errors are NA's, just use them for Shats.
    standardization <- max(abs(df1$Effects$effect), na.rm = TRUE)

    df3 <- df1$Effects %>%  # fix this: make it a not-exported function.
      mutate(stderr_d = `std Error` / standardization,
             Stand_effect = effect / standardization) %>%
      dplyr::select(SNP, Stand_effect, stderr_d) %>%
      mutate(stderr_d = ifelse(is.na(stderr_d),
                               10,
                               stderr_d),
             Stand_effect = ifelse(is.na(Stand_effect),
                                   0,
                                   Stand_effect))
    names(df3)[2] <- paste0("Bhat_", phe_col[i])
    names(df3)[3] <- paste0("Shat_", phe_col[i])
    message(paste0("GAPIT standard errors were used for", phe_col[1]))
  }

  bhat_df <- big_effects_df %>%
    dplyr::select(SNP) %>%
    left_join(df3) %>%
    dplyr::select(SNP, starts_with("Bhat"))
  shat_df <- big_effects_df %>%
    dplyr::select(SNP) %>%
    left_join(df3) %>%
    dplyr::select(SNP, starts_with("Shat"))

  set.seed(1234)
  random_sample <- sample(1:nrow(df3), nrow(big_effects_df)) %>% sort()
  bhat_random <- df3[random_sample,] %>%
    dplyr::select(SNP, starts_with("Bhat"))
  shat_random <- df3[random_sample,] %>%
    dplyr::select(SNP, starts_with("Shat"))

  for(i in seq_along(phe_col)[-1]){
    df1 <- load_GAPIT_GWAS_all(path = "./", phenotype = phe_col[i])

    if(sum(is.na(df1$Effects$`std Error`)) > length(df1$Effects$`std Error`)*.05){
      # If there are too many NA's for standard errors, derive new standard
      # errors using Hedge's G (which requires the MAF).
      standardization <- max(abs(df1$Results$effect), na.rm = TRUE)

      df3 <- df1$Results %>%
        mutate(Stand_effect = effect / standardization,
               Obs = maf * nobs,
               Obs2 = (1-maf) * nobs,
               d = ifelse(abs(Stand_effect) < 0.98,
                          (2 * Stand_effect) / sqrt(1 - Stand_effect^2),
                          4),
               d_unbiased = (1 - (3 / (4 * (nobs -2) -1))) * d,
               sigma2_d = ((Obs + Obs2) / (Obs * Obs2)) +
                 (d_unbiased^2 / (2*(Obs + Obs2))),
               stderr_d = sqrt(sigma2_d)) %>%
        mutate(Stand_effect =
                 ifelse(is.na(Stand_effect) | is.infinite(Stand_effect),
                        0,
                        Stand_effect),
               stderr_d = ifelse(is.na(stderr_d) | is.infinite(stderr_d),
                                 10,
                                 stderr_d)) %>%
        dplyr::select(SNP, Stand_effect, stderr_d)
      names(df3)[2] <- paste0("Bhat_", phe_col[i])
      names(df3)[3] <- paste0("Shat_", phe_col[i])
      message(paste0("Hedge's G standard errors were used for ", phe_col[i]))
    } else {
      # or if not many of the standard errors are NA's, just use them for Shats.
      standardization <- max(abs(df1$Effects$effect), na.rm = TRUE)

      df3 <- df1$Effects %>%
        mutate(stderr_d = `std Error` / standardization,
               Stand_effect = effect / standardization) %>%
        dplyr::select(SNP, Stand_effect, stderr_d) %>%
        mutate(stderr_d = ifelse(is.na(stderr_d),
                                 10,
                                 stderr_d),
               Stand_effect = ifelse(is.na(Stand_effect),
                                     0,
                                     Stand_effect))
      names(df3)[2] <- paste0("Bhat_", phe_col[i])
      names(df3)[3] <- paste0("Shat_", phe_col[i])
      message(paste0("GAPIT standard errors were used for ", phe_col[i]))
    }

    bhat_df <- bhat_df %>%
      left_join(df3) %>%
      dplyr::select(SNP, starts_with("Bhat"))
    shat_df <- shat_df %>%
      left_join(df3) %>%
      dplyr::select(SNP, starts_with("Shat"))

    bhat_random <- bhat_random %>%
      left_join(df3[random_sample,]) %>%
      dplyr::select(SNP, starts_with("Bhat"))
    shat_random <- shat_random %>%
      left_join(df3[random_sample,]) %>%
      dplyr::select(SNP, starts_with("Shat"))

    B_hat_random <- column_to_rownames(bhat_random, "SNP")
    S_hat_random <- column_to_rownames(shat_random, "SNP")
    B_hat_strong <- column_to_rownames(bhat_df, "SNP")
    S_hat_strong <- column_to_rownames(shat_df, "SNP")
  }
}

  if(saveoutput == TRUE){
  saveRDS(B_hat_strong, file = file.path(path, paste0("B_hat_strong_df_", numSNPs,
                                                 "topSNPs.rds")))
  saveRDS(S_hat_strong, file = file.path(path, paste0("S_hat_strong_df_", numSNPs,
                                 "topSNPs.rds")))
  saveRDS(B_hat_random, file = file.path(path, paste0("B_hat_random_df_", numSNPs,
                                     "topSNPs.rds")))
  saveRDS(S_hat_random, file = file.path(path, paste0("S_hat_random_df_", numSNPs,
                                     "topSNPs.rds")))
  }
  return(list(SNP_df = big_effects_df, B_hat_strong = B_hat_strong,
              S_hat_strong = S_hat_strong, B_hat_random = B_hat_random,
              S_hat_random = S_hat_random))
}






