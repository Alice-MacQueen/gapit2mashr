require(testthat)

test_that("Phenotypes in Raw Data Folder", {
  ph1 <- gapit_phenotypes_in_folder(path = file.path("inst", "extdata"))
  ph2 <- gapit_phenotypes_in_folder(path = "./inst/extdata")
  expect_equivalent(ph1, ph2)
  expect_equivalent(ph1, c("Earliest_Year_CDBN", "LG"))
})

test_that("gapit top effects", {
  df1 <- load_GAPIT_GWAS_all(path = "./inst/extdata", phenotype = "LG")
  df_out <- gapit_top_effects_FDRpvalue(df = df1$Results, phenotype = "LG", numSNPs = 100)
  expect_equivalent(nrow(df_out), 100)
  expect_equivalent(names(df_out)[5], "LG_effect")
})

test_that("GAPIT Standard Errors", {
  df1 <- load_GAPIT_GWAS_all(path = file.path("inst", "extdata"),
                             phenotype = "LG")
  df_gapit <- s_hat_gapit(df = df1$Effects, phenotype = "LG")
  expect_equivalent(names(df_gapit), c("SNP", "Bhat_LG", "Shat_LG"))
})

test_that("Hedges' G Standard Errors", {
  df1 <- load_GAPIT_GWAS_all(path = file.path("inst", "extdata"),
                             phenotype = "LG")
  df_hedges <- s_hat_hedges_g(df = df1$Results, phenotype = "LG")
  expect_equivalent(names(df_hedges), c("SNP", "Bhat_LG", "Shat_LG"))
  # plot(x = df_gapit$Shat_LG, y = df_hedges$Shat_LG)
})

test_that("gapit2mashr works on example data", {
  df_out1 <- gapit2mashr(path = file.path("inst", "extdata"), numSNPs = 60,
                        S_hat = "Hedges' G")
  expect_equivalent(names(df_out1$B_hat_strong), names(df_out1$B_hat_random))
  expect_equivalent(nrow(df_out1$SNP_df), nrow(df_out1$B_hat_strong),
                    nrow(df_out1$B_hat_random), nrow(df_out1$S_hat_strong),
                    nrow(df_out1$S_hat_random), 116)
})

#df_out1$SNP_df %>%
# filter(Chromosome == 7 & Position > 34400000) %>%
#  select(SNP, LG_FDR_adj_pvalue, Earliest_Year_CDBN_FDR_adj_pvalue)
#df_out1$SNP_df %>%
#  filter(!is.na(LG_FDR_adj_pvalue) & !is.na(Earliest_Year_CDBN_FDR_adj_pvalue))


## Make small example dataframes.
#  df1 <- load_GAPIT_GWAS_all(path = file.path("inst", "extdata"),
#                             phenotype = "LG")
#
#  shared_region <- df1$Results %>%
#    filter(Chromosome == 7 & between(Position, 34000000, 38000000))
#  shared_eff <- df1$Effects %>%
#    filter(Chromosome == 7 & between(Position, 34000000, 38000000))
#  example_dataset <- sample(nrow(df1$Results), 400000) %>% sort()
#
#  small_lg_results <- shared_region %>%
#    full_join(df1$Results[example_dataset,])
#  write_csv(small_lg_results, "./GAPIT.CMLM.LG.GWAS.Results.csv")
#  small_lg_effects <- shared_eff %>%
#    full_join(df1$Effects[example_dataset,])
#  write_csv(small_lg_effects, "./GAPIT.CMLM.LG.Df.tValue.StdErr.csv")

#  df2 <- load_GAPIT_GWAS_all(path = file.path("inst", "extdata"),
#                             phenotype = "Earliest_Year_CDBN")
#
#  shared_region2 <- df2$Results %>%
#    filter(Chromosome == 7 & between(Position, 34000000, 38000000))
#  shared_eff2 <- df2$Effects %>%
#    filter(Chromosome == 7 & between(Position, 34000000, 38000000))
#
#  small_eyc_results <- shared_region2 %>%
#    full_join(df2$Results[example_dataset,])
#  write_csv(small_eyc_results, "./GAPIT.CMLM.Earliest_Year_CDBN.GWAS.Results.csv")
#  small_eyc_effects <- shared_eff2 %>%
#    full_join(df2$Effects[example_dataset,])
#  write_csv(small_eyc_effects, "./GAPIT.CMLM.Earliest_Year_CDBN.Df.tValue.StdErr.csv")
