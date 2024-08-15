#f <- f()
#k <- k()
#j <- j()
###r <- m(f, k$kcd, j, sex = 1, age = 0)
###p <- cc(r, k)
###pr <- data.frame(p, r)
###prs <- subset(pr, RR > 4 & pFdr < 0.001 & ED > 8 & EcD > 8)
###s <- lapply(1:2, function(sex) { lapply(0:18, function(age) { r <- m(f, k$kcd, j, sex = sex, age = age); p <- cc(r, k); pr <- data.frame(p, r); subset(pr, RR > 4 & pFdr > 0.001 & ED > 8 & EcD > 8) }) })
## Retry
#st <- by(f, j[match(f$PERSON_ID, j$PERSON_ID), 3:4], m, k$kcd, j)
#save(st, file = "strata.nohup.RData")
## Calculation took about 2 hours up to this point.

# The following code stores the objects st, stu, and st1 in met/strata.RData
# It updates st1 in met/strata.RData if it already exists.
# It uses "st" stored in met/strata.RData if it is not already in the search space
# if (!exists("st")) load("met/strata.RData")
# pst <- lapply(st, cc, k)
# Next part applies a uniform cutoff for the stratified results and sticks p value and contingency tables together.
# Note that cc is a function that returns a data frame summarizing p values and RR.
# pst is a data frame containing p values and RR
# psts <- mapply(function(r, p) { 
#     m <- (p$RR > 4 & p$pFdr < 0.001 & r$ED > 7 & r$EcD > 7)
#     data.frame(p[m,], r[m,], stringsAsFactors = FALSE)
# }, st, pst, SIMPLIFY = FALSE)
# Prepare columns that label strata.
# le <- vapply(psts, function(x) { dim(x)[1] }, integer(1))
# s <- rep(rep(c(1, 2), 19), le)
# ag <- rep(rep(0:18, rep(2, 19)), le)
# # Assemble into one data frame.
# stu <- unlist.partdf(psts)
# stu$sex <- s
# stu$age <- ag
# stu is a data frame containing all the significant edges in each strata.

# This section selects the representative edges that survive stratification,
# and saves it in st1.
# st1 <- (function() {
#     # I select one strata per edge, 
#     # the one with the smallest pValue, 
#     # and if there are ties, the one with the largest RR.
#     key <- paste(stu$riskFactor.kcd, "(interacts with)", stu$disease.kcd)
#     splitdf <- split(stu[,c("pValue", "RR", "age", "sex")], key)
#     result <- unlist.partdf(lapply(splitdf, function(x) { 
#         minp <- min(x$pValue)
#         ties <- which(x$pValue == minp)
#         repix <- ties[which(x$RR[ties] == max(x$RR[ties]))]
#         if (length(repix) > 1) stop("Bad data")
#         return(x[repix,])
#     }))
#     sdkey <- names(splitdf)
#     riskFactor.kcd <- substr(sdkey, 1, 3)
#     disease.kcd <- substr(sdkey, 22, 24)
#     result <- data.frame(key = sdkey, riskFactor.kcd, disease.kcd, result, stringsAsFactors = FALSE)
#     return(result)
# })()
# 
# This section creates an igraph object, igstrat, representing the subnetwork that survives stratification
# igstrat <- (function() {
#     if(!exists("igs")) source("src/functions.R")
#     if(!exists("n")) source("src/nodes.R")
#     eid <- paste0(st1$riskFactor.kcd, "|", st1$disease.kcd)
#     row.names(st1) <- eid
#     igstrat <- subgraph.edges(igs, which(attr(E(igs), "vname") %in% eid))
#     eid <- intersect(attr(E(igstrat), "vnames"), eid)
#     V(igstrat)$com <- membership(cluster_walktrap(igstrat))
#     V(igstrat)$com0 <- n[names(V(igstrat)), "com"]
#     V(igstrat)$comSCol <- gray(0:28/29)[V(igstrat)$com]
#     V(igstrat)$comSCol[V(igstrat)$com == 5] <- substr(rainbow(4)[1], 1, 7)
#     V(igstrat)$comSCol[V(igstrat)$com == 15] <- substr(rainbow(4)[2], 1, 7)
#     V(igstrat)$comSCol[V(igstrat)$com == 4] <- substr(rainbow(4)[3], 1, 7)
#     E(igstrat)[eid]$strat.age <- st1[eid, "age"]
#     E(igstrat)[eid]$strat.sex <- st1[eid, "sex"]
#     E(igstrat)[eid]$strat.pValue <- st1[eid, "pValue"]
#     E(igstrat)[eid]$strat.RR <- st1[eid, "RR"]
#     igstrat
# })()
# This section creates a communities object, comstrat, representing the community structure of igstrat.

# Write to met/strata.RData
# save(igstrat, file = "met/strata.RData")
graphFromStrataRrs <- function(strataRrs, properKcds = NULL) {
    if (!is.null(properKcds)) {
        strataRrs <- strataRrs[strataRrs$riskFactor.kcd %in% properKcds & strataRrs$disease.kcd %in% properKcds,]
    }
    graphs <- lapply(split(strataRrs, list(strataRrs$sex, strataRrs$age)), graph_from_data_frame)
    names(graphs) <- levels(interaction(strataRrs$sex, strataRrs$age, sep = '-'))
    return(graphs)
}
cmh <- function(by) {
    pairs <- Reduce(`+`, lapply(by, function(tallyMatrix) {
        with(tallyMatrix, cbind(ED * (C - E) / C, EcD * E / C))
    }), matrix(rep(0, 520 * 520 * 2), nrow = 520 * 520))
    pairs[,1] / pairs[,2]
}
stratGraphsFromIncidences <- function(incidences, demog, syms = V(completeGraph)$name, ...) {
    incidences <- subset.incidences(incidences, syms = syms, ...)
    cat(dim(incidences))
    incidencesDemog <- incidenceDemogMerge(incidences, demog)
    cat(dim(incidencesDemog))
    by(incidences, with(incidencesDemog, list(SEX, AGE_GROUP)), graphFromIncidences)
}
