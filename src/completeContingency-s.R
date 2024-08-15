# Reads r and creates p and c
rrFromContingency <- function(kContin, fdr = "bonferroni", neg = FALSE) {
# r is the partial contingency table data frame
# Takes r and creates p
# p is a data frame containing riskFactor.kcd, disease.kcd, RR, pBalue, and pFdr.

# I have commented out the parts where it embellishes the output with information about the kcd codes.
# Embellishmen of the output should be performed in a separate wrapper function.
    time <- system.time(
        complete <- (function() {
            if(!(mode(fdr) == "character" && length(fdr) == 1)) stop("fdr must be string of length 1.")
            disease <- vapply(strsplit(row.names(kContin), "|", fixed = TRUE), function(x) { x[2] }, "")
            riskFactor <- vapply(strsplit(row.names(kContin), "|", fixed = TRUE), function(x) { x[1] }, "")
RR <- kContin$ED * (as.double(kContin$C) - kContin$E) / kContin$E / kContin$EcD
            pValue <- plnorm(RR ,sdlog = sqrt(1 / kContin$ED + 1 / kContin$EcD - 1 / kContin$E - 1 / (as.double(kContin$C) - kContin$E)), lower.tail = as.logical(neg))
            complete <- data.frame(riskFactor.kcd = riskFactor, disease.kcd = disease, RR, pValue, pFdr = p.adjust(pValue, fdr))
        })()
    )
    attr(complete, "time") <- time
    return(complete)
}
rr.igraph <- function(completeGraph) {
    library(igraph)
    pValues <- rrFromContingency(data.frame(edge_attr(completeGraph), row.names = as_ids(E(completeGraph)), check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE))
    edge_attr(completeGraph) <- c(edge_attr(completeGraph), pValues[,c("RR", "pValue", "pFdr")])
    completeGraph$time2 <- attr(pValues, "time")
    return(completeGraph)
}
