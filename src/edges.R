#!/usr/bin/Rscript
#Package
if("Edges" %in% search()) detach(Edges)
Edges <- list(
writeStrat = function() {
    #There might not be any use for this function, but we'll see.

load("strata.RData")
attach(stu)
name <- paste(riskFactor.kcd, "(interacts with)", disease.kcd)
outdf <- data.frame(name, riskFactor.kcd, disease.kcd, RR, sex, age, stringsAsFactors = FALSE)
write.csv(outdf, file = "edges.cys.csv", row.names = FALSE, quote = FALSE)

},
filterStrat = function(stu) {
    #First filter stu for overlap with p.
    #Resolve ties for pValue and RR
    #Return as a useful data frame

    #Filter stu for overlap with igraph
    p <- subset(p, RR > 4 & pFdr < 0.001)
    namesep <- paste0(p$riskFactor.kcd, p$disease.kcd)
    nameses <- paste0(stu$riskFactor.kcd, stu$disease.kcd)
    stu <- subset(stu, nameses %in% namesep)
    #Resolve ties for pValue and RR
    nameses <- paste0(stu$riskFactor.kcd, stu$disease.kcd)
    stuSplit <- split(stu, nameses)
    stuSplit <- lapply(stuSplit, function(splitStu) {
        minp <- min(splitStu$pValue)
        splitStu <- subset(splitStu, pValue == minp)
        maxRR <- max(splitStu$RR)
        splitStu <- subset(splitStu, RR == maxRR)
        if (dim(splitStu)[1] > 1) cat("Tie remains.\n")
        return(splitStu)
    })
    stu <- unlist.partdf(stuSplit)
    #Return as a useful data frame
    return(stu)
},
writeEdges = function() {
    p <- subset(p, RR > 4 & pFdr < 0.001)
    row.names(p) <- paste0(p$riskFactor.kcd, p$disease.kcd)
    load("met/strata.RData")
    stu <- filterStrat(stu)
    p$age <- stu[row.names(p), "age"]
    write.csv(p, "met/edges.csv", row.names = FALSE)
},
loadEdges = function() {
    read.csv("met/edges.csv", stringsAsFactors = FALSE)
},
exportCytoEdges = function(dataFrame) exportCyto(dataFrame, "out/edges.tsv")
)
attach(Edges)
rm(Edges)
