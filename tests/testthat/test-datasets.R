context("Test for Datasets")

test_that("Datasets", {
    data(ToniData)
    data(ToniData.DEGs)
    data(ToniData.TFs)

    expect_equivalent(class(ToniData), "ExpressionSet")
    expect_identical(class(ToniData.DEGs), "character")
    expect_identical(class(ToniData.TFs), "character")
    ##
    expect_equal(dim(Biobase::exprs(ToniData)), c(1391, 12))
    expect_equal(length(ToniData.DEGs), 520)
    expect_equal(length(ToniData.TFs), 898)
})
