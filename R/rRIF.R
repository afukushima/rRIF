##' This function calculates Regulatory Impact Factors (RIF1 and RIF2)
##'
##' 
##' @title calculates RIF1 and RIF2 in transcriptomic data
##'
##' @param eSet an ExpressionSet
##' @param formula an object of class "formula" (e.g., genotype1~genotype2)
##' @param target.factor a factor of interest (e.g., "Genotype")
##' @param DEGs a list of DEGs (differentially expressed genes)
##' @param cor.method a list of correlation methods [i.e. c("pearson", "spearman", "kendall")]
##' @param regulator.list a list of transcriptional regulators (e.g., transcription factors)
##' @return "rRIF" class object
##' @export
##' @references Reverter A et al. Bioinformatics 26:896 (2010)
##' (\href{https://www.ncbi.nlm.nih.gov/pubmed/20144946}{PubMed})
##' @examples
##' data(ToniData)
##' data(ToniData.DEGs)
##' data(ToniData.TFs)
##'
##' res <- rRIF(eSet = ToniData,
##'             formula = geno1~geno2,
##'             target.factor = "Genotype",
##'             DEGs = ToniData.DEGs,
##'             cor.method="pearson",
##'             regulator.list = ToniData.TFs
##'             )
##' @author Kevin Rue-Albrecht \url{https://github.com/kevinrue/HudsonRIF}
##' modified by Atsushi Fukushima
rRIF <- function(eSet, formula, 
                target.factor, DEGs,
                cor.method = "spearman",
                regulator.list = ""
                ) {

    ## Checks the validity of user-defined variables
    ## cat("Checking input variables...", fill = TRUE)
    parameter.check(eSet, formula, target.factor, DEGs, cor.method, regulator.list)

    ## Set the sample info from the formula
    samp <- list(A = as.character(formula[[2]]), 
                 B = as.character(formula[[3]])
                 )

    ## Identifies the list of regulators to consider
    if(any(regulator.list == "")) {
        regulator.list = rownames(Biobase::exprs(eSet))[!rownames(
          Biobase::exprs(eSet)
          ) %in% DEGs]}

    ## Calculate the average expression of each gene in each condition
    EiAB <- calculateEiAB(eSet = eSet,
                          target.factor = target.factor,
                          samp$A, samp$B)

    ## Calculates the average abundance of each gene across conditions
    Ai <- (EiAB[, samp$A] + EiAB[, samp$B]) / 2

    ## Calculate the differential expression in log2foldchange
    ## for each gene in formula
    dEi <- EiAB[, samp$A] - EiAB[, samp$B]

    ## Calculate the PIF
    PIFi <- Ai * dEi
    
    ## Calculate the coexpression of all pairs of genes in each condition
    ## condition A
    rAij <- stats::cor(x = t(Biobase::exprs(
      eSet[, which(Biobase::pData(eSet)[, target.factor] == samp$A)])), 
      method = cor.method)
    rAij <- rAij[regulator.list, DEGs]
    ## condition B
    rBij <- stats::cor(x = t(Biobase::exprs(
      eSet[, which(Biobase::pData(eSet)[, target.factor] == samp$B)])), 
      method = cor.method)
    rBij <- rBij[regulator.list, DEGs]

    ## Calculate the difference in coexpression between the two conditions
    dCij <- rAij - rBij

    ## Calculate the RIF1/2
    RIF1 <- calculateRIF1(PIFi = PIFi, dCij = dCij, DEGs = DEGs)
    RIF2 <- calculateRIF2(EiAB = EiAB, rAij = rAij, rBij = rBij, 
                         DEGs = DEGs, samp = samp)
    res <- list(eSet = eSet, DEGs = DEGs, target.factor = target.factor, 
                  samp = samp, regulator.list = regulator.list, 
                  EiAB = EiAB, Ai = Ai, dEi = dEi, PIFi = PIFi, 
                  rAij = rAij, rBij = rBij, dCij = dCij, 
                  RIF1 = RIF1, RIF2 = RIF2)
    class(res) <- c("rRIF", "list")
    return(res)
}

## calculatePIF
##' This function calculates phenotype impact factor (PIF).
##'
##' @title calculates phenotype impact factor (PIF) from transcriptomic data
##'
##' @param eSet an ExpressionSet
##' @param target.factor a factor of interest
##' @param A a level of a factor (e.g., condition 1)
##' @param B a level of a factor (e.g., condition 2)
##' @return numeric
##' @export
##' @references Reverter A et al. Bioinformatics 26:896 (2010)
##' (\href{https://www.ncbi.nlm.nih.gov/pubmed/20144946}{PubMed})
##' @examples
##' data(ToniData)
##'
##' formula <- geno1~geno2
##' samp <- list(A=as.character(formula[[2]]), B=as.character(formula[[3]]))
##' target.factor <- "Genotype"
##' PIF <- calculatePIF(eSet=ToniData, target.factor=target.factor, samp$A, samp$B)
##' 
##' @author Kevin Rue-Albrecht \url{https://github.com/kevinrue/HudsonRIF}
##' cosmetical changes by Atsushi Fukushima
calculatePIF <- function(eSet = eSet,
                          target.factor = target.factor,
                          A, B) {
  ## Calculate the average expression of each gene in each condition
  EiAB <- calculateEiAB(eSet = eSet,
                        target.factor = target.factor,
                        A, B)
  
  ## Calculates the average abundance of each gene across conditions
  Ai <- (EiAB[, A] + EiAB[, B]) / 2
  
  ## Calculate the differential expression in log2foldchange 
  ## for each gene in formula
  dEi <- EiAB[, A] - EiAB[, B]
  
  ## Calculate the PIF
  PIFi <- Ai * dEi
  
  return(PIFi)
}

## This function calculates expression levels of each gene 
## in each data subset.
calculateEiAB <- function(eSet, target.factor, A, B)
{
  # Calculates the mean expression of features in samples from class A
  EiAB <- calculateEiA(eSet, target.factor, A)
  # Calculates the mean expression of features in samples from class B
  EiAB <- cbind(EiAB, calculateEiA(eSet, target.factor, B))
  # Relevant column names
  colnames(EiAB) <- c(A,B)
  # Returns the 2-column matrix
  return(EiAB)
}


## This function calculates expression levels of each gene 
## in one condition (A).
calculateEiA <- function(eSet, target.factor, A)
{
  # Calculates the mean expression of features in samples from a given class
  return(apply(X = Biobase::exprs(
    eSet[, which(Biobase::pData(eSet)[, target.factor] == A)]), 
    MARGIN = 1, FUN = "mean"))
}


##' This function calculates RIF1.
##'
##' @title calculates RIF1
##'
##' @param PIFi PIF
##' @param dCij a differential coexpression matrix
##' @param DEGs a list of DEGs (differentially expressed genes)
##' @return vector
##' @export
##' @examples
##' data(ToniData)
##'
##' #formula <- geno1~geno2
##' #samp <- list(A=as.character(formula[[2]]), B=as.character(formula[[3]]))
##' #target.factor <- "Genotype"
##' #EiA <- calculateEiA(eSet=ToniData, target.factor=target.factor, samp$A)
##'
##' @author Kevin Rue-Albrecht \url{https://github.com/kevinrue/HudsonRIF}
##' cosmetical changes by Atsushi Fukushima
calculateRIF1 <- function(PIFi, dCij, DEGs) {
  # Subset the PIF values to the DE genes
  DE.PIFi <- PIFi[DEGs]
  # Squares the values
  dCij.coefs <- dCij^2
  # The matrix product return the RIF values (absolute PIF differs from original formula)
  # Original Hudson formula using non-absolute PIF value
  res <- apply(X = ((dCij.coefs %*% DE.PIFi) / length(DEGs)), 
               MARGIN = 1, FUN = sum)
  ## convert to Z-score
  RIF1 <- as.vector(scale(res))
  names(RIF1) <- names(res)
  return(RIF1)
}

##' This function calculates RIF2.
##'
##' @title calculates RIF2
##'
##' @param EiAB an expression level
##' @param rAij a coexpression matrix in A
##' @param rBij a coexpression matrix in B
##' @param DEGs a list of DEGs (differentially expressed genes)
##' @param samp conditions
##' @return vector
##' @export
##' @examples
##' data(ToniData)
##'
##' #formula <- geno1~geno2
##' #samp <- list(A=as.character(formula[[2]]), B=as.character(formula[[3]]))
##' #target.factor <- "Genotype"
##' #EiA <- calculateEiA(eSet=ToniData, target.factor=target.factor, samp$A)
##'
##' @author Atsushi Fukushima
###################################################
calculateRIF2 <- function(EiAB, rAij, rBij, DEGs, samp) {
    EiAB.DE <- EiAB[DEGs,]
    res <- apply(X = ((rAij^2 %*% EiAB.DE[,samp$A]^2) - 
                      (rBij^2 %*% EiAB.DE[, samp$B]^2))/length(DEGs), 1, sum)
    ## convert to Z-score
    RIF2 <- as.vector(scale(res))
    names(RIF2) <- names(res)
    return(RIF2)
}

## This function validates all parameters
parameter.check <- function(eSet, formula, target.factor, DEGs, cor.method, regulator.list) {
  ## eSet is an expressionSet
  if(class(eSet) != "ExpressionSet") 
    stop("\"eSet\" (", class(eSet),") is not an ExpressionSet.", call. = FALSE)
  # formula should be a formula
  if(class(formula) != "formula")
    stop("\"formula\" (", class(formula),") is not a formula.", call. = FALSE)
  ## formula respect format A~B
  if(length(formula[[2]]) > 1 )
    stop("\"formula\" (", 
         paste(as.character(formula[c(2,3)]), collapse=" ~ ") ,
         ") contains more than one term on the left side.", call. = FALSE)
  if(length(formula[[3]]) > 1){
    stop("\"formula\" (", paste(as.character(
      formula[c(2,3)]), collapse=" ~ ") ,
      ") contains more than one term on the right side.", call. = FALSE)}
  ## target.factor is a valid column name in pData(eSet)
  if(!target.factor %in% names(Biobase::pData(eSet))){
    stop("\"target.factor\" (", target.factor ,
         ") is not a valid column name in \"pData(eSet)\"", call. = FALSE)}
  # A and B should be valid levels of target.factor
  samp <- list(A = as.character(formula[[2]]), 
               B = as.character(formula[[3]]))
  if(!samp$A %in% levels(as.factor(Biobase::pData(eSet)[,target.factor]))){
    stop("(", samp$A ,") is not a valid class level in \"pData(eSet)[,", 
         target.factor,"]\"", call. = FALSE)}
  if(!samp$B %in% levels(as.factor(Biobase::pData(eSet)[,target.factor]))){
    stop("(", samp$B ,") is not a valid class level in \"pData(eSet)[,", 
         target.factor,"]\"", call. = FALSE)}
  ## DEGs fully included in list of features
  if(sum(DEGs %in% rownames(Biobase::exprs(eSet))) != length(DEGs)){
    stop("\"DEGs\ contains ", sum(!DEGs %in% rownames(Biobase::exprs(eSet))) ,
         " feature(s) absent from \"rownames(exprs(eSet))\"", call. = FALSE)
  }
  if (!is.character(cor.method))
    stop("cor.method must specify a method of stats::cor (e.g., \"pearson\").")
  
}
