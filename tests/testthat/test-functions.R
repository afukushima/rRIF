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
  expect_equivalent(round(res$RIF1[1], 2), -0.66)
  expect_equivalent(round(res$RIF1[10], 2), -1.51)
  expect_equivalent(round(res$RIF1[length(res$RIF1)], 2), -0.15)
  ## RIF2
  expect_equivalent(round(res$RIF2[1], 2), 0.2)
  expect_equivalent(round(res$RIF2[10], 2), -1.79)
  expect_equivalent(round(res$RIF2[length(res$RIF2)], 2), 0.19)

})
