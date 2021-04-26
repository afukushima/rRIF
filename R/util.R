##' @export
print.rRIF <- function(x, ...) {
    print.data.frame(make.rRIFtable(x))
}
##' @export
summary.rRIF <- function(object, ...) {
    res <- make.rRIFtable(object)

    cat(paste("\n\n===== rRIF Summary =====\n"))
    cat(paste("Top 5 RIF1:\n"))
    utils::head(res[base::order(abs(res$RIF1), decreasing = TRUE),], 5)
}

make.rRIFtable <- function(x) {
    tgt <- names(x$RIF1)
    Ai <- round(x$Ai[tgt], digits = 2)
    dEi <- round(x$dEi[tgt], digits = 2)
    PIFi <- round(x$PIFi[tgt], digits = 2)
    RIF1 <- round(x$RIF1, digits = 2)
    RIF2 <- round(x$RIF2, digits = 2)
    RIFtable <- data.frame(cbind(Ai, dEi, PIFi, RIF1, RIF2))
    colnames(RIFtable) <- c("A", "DE", "PIF", "RIF1", "RIF2")
    return(RIFtable)
}

##' This function provides an interactive table of rRIF results
##'
##'
##' @title an interactive table of rRIF object
##'
##' @param rRIF "rRIF" class object
##' @param top specifies top X RIF
##' @return none
##' @export
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
##' iTable(res)
iTable <- function(rRIF, top = 30) {
    res <- make.rRIFtable(rRIF)
    if (top == "all" || top == "ALL")
        res <- res[base::order(abs(res$RIF1), decreasing = TRUE),]
    else
        res <- utils::head(res[base::order(abs(res$RIF1),
                                            decreasing = TRUE),], top)
    DT::datatable(res)
}

##' This function provides an interactive MA-plot of rRIF results
##' Differential expression of all genes plotted against the average abundance
##'
##' @title an interactive MA-plot of rRIF object
##'
##' @param rRIF "rRIF" class object
##' @return none
##' @export
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
##' iMAplot(res)
iMAplot <- function(rRIF) {
    # figure lables for 2D charts
    f <- list(
        family = "Courier New, monospace",
        size = 20,
        color = "black"
    )
    x <- list(
        title = "A",
        titlefont = f
    )
    y <- list(
        title = "M",
        titlefont = f
    )
    # Plot with symmetric Y-axis range
    # plot
    p <- plotly::plot_ly(type = "scatter", mode = "markers")
    p <- plotly::add_trace(p, x = rRIF$Ai, y = rRIF$dEi,
                            text = paste("GeneName: ", names(rRIF$Ai)),
                            color = rRIF$dEi, marker = list( size = rRIF$Ai ),
                            name = "")
    p <- plotly::layout(p, xaxis = x, yaxis = y, showlegend = TRUE)
    p
}

##' Plots for the coexpression of the specified gene with each of the DE.
##'
##' Modified Hudson::diffCoexPIF.plot()
##' (https://github.com/kevinrue/HudsonRIF)
##'
##' @title Co-expression relationships between the specified TF and the DE
##'
##' @param rRIF "rRIF" class object
##' @param reg.gene The regulatory gene/probeset to plot co-expression values
##' @param pch.max A scaling factor defining the size of the data point
##' @return none
##' @export
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
##' iMAplot(res)
diffCoexPIF.plot <- function(rRIF, reg.gene, pch.max = 0) {
    # Check the validity of user-defined variables
    paramCheck.diffCoexPIF.plot(rRIF, probe=reg.gene, pch.max=pch.max)
    # Filters and normalise the the maximum circle size if required
    if (pch.max == 0) {
        pch.max = max(rRIF$PIFi)
    }
    # Scales the data points to the maximal size defined
    PIF.DE = rRIF$PIFi[rRIF$DElist]/max(rRIF$PIFi) * pch.max
    # plots the coexpression values in each condition
    graphics::plot(x = rRIF$rBij[reg.gene, ], y = rRIF$rAij[reg.gene, ],
                xlim = c(-1, 1), ylim = c(-1, 1), cex = PIF.DE,
                xlab = rRIF$conds$B, ylab = rRIF$conds$A,
                pch = 20, main = reg.gene, col = "red")
    # adds the identity line (no coexpression change between conditions)
    graphics::abline(coef = c(0, 1), col = "grey", lwd = 2)
}


paramCheck.diffCoexPIF.plot <- function(rRIF, probe, pch.max) {
    # rRIF is a working output of the rRIF wrapper
    if(methods::is(rRIF)[1] != "rRIF"){
        stop("\"rRIF\" (", methods::is(rRIF),") is not an rRIF object.",
            call.=FALSE)
    }
    # Note: $eSet, $classCol, $regulator.list,
    # EiAB, Ai, dEi, dCij not necessary
    if(any(!c("target.factor", "DEGs", "PIFi",
            "rAij", "rBij", "RIF1", "RIF2") %in% names(rRIF))){
        stop("\"rRIF\" (", names(rRIF),") is not a valid rRIF object.
            One item of $target.factor, $DEGs, $PIFi,
            $rAij, $rBij or $RIF1/2 is missing.", call.=FALSE)
    }
    # Check that probe is a valid probeID
    if(!is.character(probe)){
        stop("\"probe\" (", probe,") is not a character object.", call.=FALSE)
    }
    if(! probe %in% rownames(rRIF$rAij)){
        stop("\"probe\" (", probe,") is not a valid probe ID in rAij rows.",
            call.=FALSE)
    }
    if(! probe %in% rownames(rRIF$rBij)){
        stop("\"probe\" (", probe,") is not a valid probe ID in rBij rows.",
            call.=FALSE)
    }
    # pch.max should be a positive integer
    if(!is.numeric(pch.max)){
        stop("\"pch.max\" (", pch.max,") is not a numeric object.",
            call.=FALSE)
    }
    if(pch.max < 0){
        stop("\"pch.max\" (", pch.max,") should be a positive integer.",
            call.=FALSE)
    }
}
