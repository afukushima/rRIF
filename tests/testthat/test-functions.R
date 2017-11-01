context("Calculation of Regulatory impact factors (RIF)")


test_that("RIF1/2 calculation", {
  data(ToniData)
  data(ToniData.DEGs)
  data(ToniData.TFs)

  ## calculation of RIF1/2 by Reverter et al. Bioinformatics (2010)
  res <- rRIF(eSet = ToniData,
              formula = geno1~geno2,
              target.factor = "Genotype",
              DEGs = ToniData.DEGs,
              cor.method = "pearson",
              regulator.list = ToniData.TFs
              )

  ## check parameters
  para.check <- parameter.check(
         eSet = ToniData,
         formula = geno1~geno2,
         target.factor = "Genotype",
         DEGs = ToniData.DEGs,
         cor.method = "pearson",
         regulator.list = ToniData.TFs
         )
  expect_true(is.null(para.check))
  
  ## PIF
  formula <- geno1~geno2
  samp <- list(A=as.character(formula[[2]]), B=as.character(formula[[3]]))
  res.PIF <- calculatePIF(eSet = ToniData,
                          target.factor = "Genotype",
                          samp$A, samp$B)
  expect_identical(class(res.PIF), "numeric")
  
    
  ## test for calculated RIFs
  expect_identical(class(res), c("rRIF", "list"))
  
  ## RIF1
  expect_equivalent(round(res$RIF1[1], 6), 1.503376)
  expect_equivalent(round(res$RIF1[10], 6), 1.140441)
  expect_equivalent(round(res$RIF1[length(res$RIF1)], 7), 0.717939)
  ## RIF2
  expect_equivalent(round(res$RIF2[1], 6), -0.024975)
  expect_equivalent(round(res$RIF2[10], 7), -0.2112097)
  expect_equivalent(round(res$RIF2[length(res$RIF2)], 7), 0.3993724)

})
