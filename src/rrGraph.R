#Todo: attach pid list to rrGraph object.
#Useful for updating rrGraph with new calculations from raw objects, e.g. onset.
#At the cost of object size.
is.rrGraph <- function(object) {
    if (!inherits(object, "igraph")) return(FALSE)
    library(igraph)
    if (!all("studyYear" %in% graph_attr_names(object))) return(FALSE)
    return(TRUE)
}
check.rrGraph <- function(object) {
    if (!is.rrGraph(object)) stop("Invalid rrGraph object.")
}
cutOff <- function(rrGraph, rr = 4, pFdr = 0.001) {
    # Applies a cutoff to the RR and pFdr columns in a graph to obtain a subgraph.
    library(igraph)
    rrGraph; check.rrGraph(rrGraph)
    if (!all(c("RR", "pFdr") %in% edge_attr_names(rrGraph))) stop("Columns RR and pFdr not found in rrGraph. Run rr.igraph first.")
    rrLimit <- rr
    pfdrLimit <- pFdr
    if ("com" %in% vertex_attr_names(rrGraph)) rrGraph <- delete_vertex_attr(rrGraph, "com")
    if ("comCol" %in% vertex_attr_names(rrGraph)) rrGraph <- delete_vertex_attr(rrGraph, "comCol")
    newGraph <- subgraph.edges(rrGraph, E(rrGraph)[RR > rrLimit & pFdr < pfdrLimit & !is.na(RR)])
    return(newGraph)
}
writeCytoscapeEdges <- function(rrGraph, file = "met/edges.csv", mac = FALSE) {
    library(igraph)
    source("src/function.R")
    edgeData <- data.frame0(edge_attr(rrGraph))
    vnames <- attr(E(rrGraph), "vnames")
    edgeData$source <- sub("\\|[[:alnum:]]+", "", vnames)
    edgeData$target <- sub("[[:alnum:]]+\\|", "", vnames)
    edgeData$interaction <- "risk for"
    edgeData$invRr <- 1 / edgeData$RR
    source("src/functions.R")
    exportCyto(edgeData, file, mac = mac)
}
writeCytoscapeNodes <- function(rrGraph, file = "met/nodes.csv", mac = FALSE) {
    library(igraph)
    source("src/functions.R")
    nodeData <- data.frame0(vertex_attr(rrGraph))
    numNodes <- gorder(rrGraph)
    equalized <- rank(c(nodeData$male, nodeData$female)) / 8 + 4
    nodeData$height <- equalized[1:numNodes]
    nodeData$width <- equalized[(numNodes + 1):(2 * numNodes)]
    source("src/functions.R")
    exportCyto(nodeData, file, mac = mac)
}
updateEng.rrGraph <- function(rrGraph, codeInfo = read.csv("src/kcd-names.csv", colClasses = "character")) {
    library(igraph)
    V(rrGraph)$eng <- codeInfo$eng[match(V(rrGraph)$name, codeInfo$kcd)]
    if (anyNA(V(rrGraph)$eng)) stop("Missing code description for codes: ", paste(V(rrGraph)$name[is.na(V(rrGraph)$eng)], collapse = " "))
    return(rrGraph)
}
updateVDemog <- function(rrGraph, incidences, demogData, years = 2002L:(rrGraph$studyYear - 1)) {
    # <demogData> should contain demographic info for the cohort used to generate rrGraph
    check.rrGraph(rrGraph)
    c(incidences, demogData)
    if (length(years) < 1) stop("Bad years specification for analyzing the demographic composition of diseases'(symbols') patients.")
    source("src/functions.R")
    symDemogSummary <- summarySymDemog(incidences, demogData, syms = V(rrGraph)$name, years = years)
    rowNames <- row.names(symDemogSummary)
    columnNames <- names(symDemogSummary)
    for(columnName in columnNames) vertex_attr(rrGraph, columnName) <- symDemogSummary[,columnName]
    return(rrGraph)
}
updateVOnset <- function(rrGraph, incidences, demogData) {
    incidences; demogData
    source("src/onset.R")
    V(rrGraph)$moa <- newInc(incidences, demogData, syms = V(rrGraph)$name)[V(rrGraph)$name]
    return(rrGraph)
}
updateCategory <- function(rrGraph, catData = read.csv("met/cat.csv")) {
    library(igraph)
    patterns <- c("[AB][0-9]{2}", "C[0-9]{2}|D[0-4][0-9]", "D[5-9][0-9]", "E[0-9]{2}", "F[0-9]{2}", "G[0-9]{2}", "H[0-5][0-9]", "H[6-9][0-9]", "I[0-9]{2}", "J[0-9]{2}", "K[0-9]{2}", "L[0-9]{2}", "M[0-9]{2}", "N[0-9]{2}", "O[0-9]{2}", "P[0-9]{2}", "Q[0-9]{2}", "R[0-9]{2}", "[ST][0-9]{2}", "[V-Y][0-9]{2}", "Z[0-9]{2}", "U[0-9]{2}")
    patternClassifier <- function(pattern, kcd) {
        searchPatterns <- paste0("^", patterns, "$")
        cat <- integer(length(kcd))
        return(grepl(pattern, kcd))
    }
    multiClassifier <- Vectorize(patternClassifier, vectorize.args = c("pattern"))
    intermediateResult <- multiClassifier(patterns, kcd = V(rrGraph)$name)
    browser()
    if(any(rowSums(intermediateResult > 0) > 1)) stop("Ambiguity in category designation.")
    V(rrGraph)$cat <- intermediateResult %*% 1:22
    V(rrGraph)$catCol <- substr(rainbow(22), 1, 7)[V(rrGraph)$cat]
    return(rrGraph)
}
updateClustering <- function(rrGraph, numMajor = 4) {
    # Adds or updates clustering information and color designation to a graph object.
    library(igraph)
    originalMembership <- membership(cluster_walktrap(rrGraph))
    sortedCommunityIds <- names(sort(table(originalMembership), decreasing = TRUE))
    conversionVector <- 1:length(sortedCommunityIds)
    names(conversionVector) <- sortedCommunityIds
    V(rrGraph)$com <- conversionVector[as.character(originalMembership)]
    V(rrGraph)$comCol <- c(substr(rainbow(numMajor), 1, 7), rep("#777777", length(sortedCommunityIds) - numMajor))[V(rrGraph)$com]
    return(rrGraph)
}

