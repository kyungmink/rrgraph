#Eventually I want this script to load functions into the visible scope
#that outputs result figures.
#For the single-column figures required by PNAS, figure width must be 1028 px.
#For the page-wide figures, figure width must be 2103 px.
#For both, pointsize must be 20 or cex must be at least 1.64042 if the pointsize is at the default 12.
#Multifigure panels might insert a cex of 0.83, so for this the pointsize must be 24.

# Results > Risk ratios ..
figRes11 <- function(newDev = TRUE) {
    intervalSize <- c(4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 15)
    ageMidPoints <- c(3, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5, 62.5, 67.5, 72.5, 77.5, 82.5, 92.5)
    load("met/j.RData")
    if(!exists("fetch")) load("met/fetch.RData")
    countsInAgeGroups <- table(j[j$PERSON_ID %in% fetch$PERSON_ID[fetch$year %in% c(2002, 2003)], c("SEX", "AGE_GROUP")])
    densityM <- countsInAgeGroups[1,] / intervalSize
    densityF <- countsInAgeGroups[2,] / intervalSize
    if(newDev) png("out/figRes11.png", width = 1536, height = 1536, pointsize = 64)
    plot(c(0, ageMidPoints, 110), c(0, densityM + densityF, 0), type = "n", xlab = "age", ylab = "number of patients", yaxs = "i")
    polygon(c(0, ageMidPoints, 100), c(0, densityM + densityF, 0), col = "lightgray")
    polygon(c(0, ageMidPoints, 100), c(0, densityF, 0), col = "darkgray")
    lines(c(2.5, ageMidPoints[2:18]), c(2382,3169,3435,3101,3662,3672,4096,4113,4123,3901,2855,2278,1889,1680,1253,767,432,233), lwd = 4)
    legend("topright", legend = c("male", "female", "2005 census"), col = c("lightgray", "darkgray", "black"), pch = c(15, 15, NaN), lwd = c(0, 0, 4), lty = c(0, 0, 1))
    if(newDev) dev.off()
}
figRes12 <- function(newDev = TRUE) {
    if(!exists("fetch")) load("met/fetch.RData")
    prevalences <- table(fetch$SICK_SYM[fetch$year %in% c(2002, 2003)])
    histo <- hist(prevalences, breaks = "Scott", plot = FALSE)
    if(newDev) png("out/figRes12.png", width = 1536, height = 1536, pointsize = 64)
    plot(histo$breaks[1:(length(histo$breaks) - 1)], histo$density, log = "xy", xlab = "patients", ylab = "p(patients)", xaxt = "n", yaxt = "n")
    axis(1, at = c(1e2, 1e3, 1e4, 1e5), labels = c(expression(10 ^ 2), expression(10 ^ 3), expression(10 ^ 4), expression(10 ^ 5)))
    axis(2, at = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c(expression(10 ^ -7), expression(10 ^ -6), expression(10 ^ -5), expression(10 ^ -4), expression(10 ^ -3)))
    lines(c(100, 500000), c(1/100, 1/500000))
    text(30000, 1e-4, expression(patients ^ -1))
    if(newDev) dev.off()
}
figRes13m <- function() {
    load("met/completeGraph.RData")
    library(igraph)
    RR <- E(completeGraph)$RR
    # Histogram of RR
    png("out/figRes13.png", 896, 896, pointsize = 48)
    hist(RR, xlab = "RR", main = NULL, breaks = "Scott")
    dev.off()
    # Histogram of log RR
    png("out/figRes14.png", 1536, 1536, pointsize = 64)
    hist(log(RR), xlab = "log RR", main = NULL, breaks = "Scott")
    Rs <- log(RR[!is.na(RR) & RR != 0])
    abline(v = mean(Rs))
    dev.off()
    # Normal Q-Q plot of log RR
    png("out/figRes15.png", 1536, 1536, pointsize = 64)
    qqnorm(Rs, main = "Normal Q-Q plot of log RR")
    dev.off()
}
figRes1m <- function() {
    png("out/figRes1m.png", 1028, 1028, pointsize = 24)
    par(mfrow = c(2, 2), mar = c(5, 5, 1, 1), oma = c(0, 1, 1, 0))
    figRes11(newDev = FALSE)
    mtext("A", adj = -3/8)
    figRes12(newDev = FALSE)
    mtext("B", adj = -3/8)
    load("met/completeGraph.RData")
    library(igraph)
    RR <- E(completeGraph)$RR
    # Histogram of log RR
    hist(log(RR), xlab = "log RR", main = NULL, breaks = "Scott")
    Rs <- log(RR[!is.na(RR) & RR != 0])
    abline(v = mean(Rs))
    mtext("C", adj = -3/8)
    # Normal Q-Q plot of log RR
    qqnorm(Rs, main = NULL)
    mtext("D", adj = -3/8)
    dev.off()
}
figRes2leg1 <- function(resampleRate = 1) {
    #For some reason, the pixel size of a user unit is not 1 even if all margins are 0 and xlim, ylim sizes are equal to the device size.
    #I have found that the pixel size of a user unit is smaller by approximately 9%.
    #So even if you supply a resampleRate that exacts the ratio of the pixel dimensions of a node to the virtual dimensions in Cytoscape, the node size in the legend are not consistent with the image.
    #In order to make this consistent, I correct the resampleRate.
    png("out/figRes2leg1.png", 210 * resampleRate, 160 * resampleRate)
    plotBottom = 16
    cexRate <- 1.3 * resampleRate
    par(mar = c(3, 4, 2, 3), oma = c(0, 0, 0, 0), cex = par("cex") * cexRate, cex.main = par("cex.main") * cexRate, cex.lab = par("cex.lab") * cexRate)
    plot.new()
    plot.window(c(0, 77.3 + plotBottom), c(-plotBottom, 77.3), asp = 1)
    indexSizes <- c(4.68, 10.4, 31.3, 77.3)
    null <- mapply(function(indexSize, color) {
        rect(77.3 - indexSize, 0, 77.3, indexSize, col = color, lty = "blank")
    }, rev(indexSizes), c("lightgray", "gray", "darkgray", "black"))
    axisLabels = c(0, expression(10 ^ 2), expression(10 ^ 3), expression(10 ^ 5))
    axis(1, at = 77.3 / 3 * 3:0, labels = axisLabels, lty = "blank", las = 2, hadj = 1 / 2, padj = 1 / 4 + 1 / 8, pos = -plotBottom)
    axis(4, at = 77.3 / 3 * 0:3, labels = axisLabels, lty = "blank", las = 2, hadj = 1 / 2, padj = 1 / 4 + 1 / 8, pos = 77.3 + plotBottom)
    null <- mapply(function(indexSize, i) {
        equiCoord <- 77.3 / 3 * i
        lines(c(77.3 - equiCoord, 77.3 - indexSize), c(-plotBottom, 0))
        lines(77.3 + c(0, plotBottom), c(indexSize, equiCoord))
    }, indexSizes, 0:3)
    mtext("female patients", side = 1, cex = par("cex.lab"), line = 2)
    mtext("Node size: prevalence", side = 3, cex = par("cex.main"), line = 1)
    mtext("male patients", side = 4, cex = par("cex.lab"), line = 2)
    dev.off()
}
figRes2leg2 <- function(resampleRate = 1, plot = TRUE) {
    toCanvas <- cbind(c(0, 0, 1), c(sqrt(3), 3, 1), c(2 * sqrt(3), 0, 1)) * 2
    #toCanvas <- cbind(c(0, 512, 512), c(512 / sqrt(3), 0, 512), c(512 * 2 / sqrt(3), 512, 512)) / 256
    toCanvas <- matrix(c(0, 2, 2, 2 / sqrt(3), 0, 2, 4 / sqrt(3), 2, 2), nrow = 3)
    toColorspace <- solve(toCanvas)
    pixelMatrix <- outer(0:511, 0:592, Vectorize(function(y, x) {
        canvasCoordinates <- matrix(c(x, y, 511), nrow = 3)
        colorspaceCoordinates <- floor(toColorspace %*% canvasCoordinates)
        if (any(colorspaceCoordinates < 0)) color <- "#FFFFFF"
        else color <- sprintf("#%02X%02X%02X", colorspaceCoordinates[1], colorspaceCoordinates[2], colorspaceCoordinates[3])
        if (nchar(color) != 7) browser()
        return(color)
    }))
    rasterized <- as.raster(pixelMatrix)
    if (plot) {
        #Factor for balancing raster image to text
        factor = 1/4 + 1/16
        png("out/figRes2leg2.png", 1508 * resampleRate * factor, 996 * resampleRate * factor)
        plot(rasterized, xlim = c(-512, 996), ylim = c(-128, 768), xpd = NA)
        #To make the font size consistent with figRes2leg1, we need a corrective factor for cex.
        par(oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), xpd = NA, cex = 1.3 * resampleRate)
        text(-384, 700, pos = 4, labels = "Node color: age composition", cex = 1.2)
        text(289, 500, labels = " 30 to 59", pos = 3)
        text(0, 0, labels = "29 or younger", pos = 1)
        text(577, 0, labels = " 60 or older", pos = 1)
        dev.off()
    }
    return(rasterized)
}

# Draw networks using Cytoscape.
# Styles and details are stored in met/rr.cys

figC02 <- function() {
    # Export data for catCol
    V(igs)$catCol <- n[V(igs)$name, "catCol"]
    data.frame(vertex_attr(igs), stringsAsFactors = FALSE) -> x1
    write.csv(x1, "r.csv", row.names = FALSE)
    savehistory("h.R")
}

## Legend for catCol
#source("src/cat.R")
figR04 <- function() {
    png("figR04.png", 874, 330, pointsize = 48)
    defaultMar <- par("mar")
    par(mar = c(0, 0, 0, 0), family = "mono")
    txt <- paste(format(ca$name), ca$kcd, format(ca$desc), ca$num)
    plot.new()
    legend("topleft", legend = txt, fill = substr(ca$col, 1, 7), box.lty = 0)
    dev.off()
}
#
##Print modularity score
#igs <- ig()
#ca <- updat(ca, igs)
#cat("Modularity score of the network divided into KCD categories =", attr(ca, "modularity"), ".\n")
##Legend for comCol with community sizes 
#(function() {
#    png("figR05.png", 874, 330)
#    defaultMar <- par("mar")
#    par(mar = c(0, 0, 0, 0), family = "mono")
#    txt <- paste(format(co$id), co$num)
#    plot.new()
#    legend("topleft", legend = txt, fill = substr(co$col, 1, 7), box.lty = 0)
#    dev.off()
#})()

resR5 <- function() {
    cat("Modularity score of the network divided into communities =", attr(ca, "modularity"), ".\n")
}
figRes23 <- function(rrGraph, incidences, demogData) {
    # Demographic characteristics of major clusters
    library(igraph)
    name <- V(rrGraph)$name
    com <- V(rrGraph)$com
    source("src/functions.R")
    s <- factorDemog(incidences, demogData, name, com)
    s <- s[s$factor < 5,]
    s$age <- age(s$AGE_GROUP)
    png("out/figRes23.png", 1028, 517, pointsize = 24)
    par(mfrow = c(1, 2), mar = c(5, 5, 1, 1), oma = c(0, 1, 1, 4))
    boxplot(s$age ~ s$f, xlab = "cluster", ylab = "age")
    mtext("A", adj = -3/8)
    ta <- matrix(table(s$SEX, s$f), nrow = 2)
    barplotData <- ta / rep(colSums(ta), rep(2, 4))
    barplot(barplotData, xlab = "cluster", names.arg = c(1, 2, 3, 4))
    mtext("B", adj = -3/8)
    legend(5, 1, c("male", "female"), pch = c(15, 15), col = gray.colors(2), xpd = NA)
    dev.off()
    return(barplotData)
}

#Community class composition
#source("src/cat.R")
#source("src/com.R")
figSup1 <- function(rrVertexInfo) {
    rrVertexInfo
#table cannot deal with the class roman
    x <- table(factor(as.integer(as.roman((rrVertexInfo$cat))), levels = 1:22L), rrVertexInfo$com)
x <- x[,1:4]
d <- dimnames(x)[[1]]
    plotData <- cbind(x / rep(colSums(x), rep(22, 4)), matrix(rep(NaN, 66L), nrow = 22))
    source("src/functions.R")
    png("out/figSup1.png", 768, 512)
    par(oma = c(0, 0, 0, 0), mar = c(4, 2, 1, 1/2))
    barplot(plotData, col = rainbow(22), xlab = "cluster")
    legend("topright", legend = paste(as.character(as.roman(1:22L)), kcdCategories, sep = ". "), pch = rep(15, 22), col = rainbow(22), cex = 9/8)
dev.off()
}

#Enrichment analysis
tabRes21 <- function(nodeInfo) {
source("src/nodes.R")
en <- enrichComCat(nodeInfo)
print(en)
}

#Power-law distribution of degree
figRes09 <- function(rrGraph, newDev = TRUE) {
    # rrGraph must be an incomplete Graph. Otherwise, otherwise, degree distribution is uniform.
    rrGraph
    library(igraph)
    if (graph.density(rrGraph) == 1) stop("figRes09:rrGraph must be an incomplete Graph.")
    source("src/functions.R")	
    if (newDev) png("out/figR09.png", 2048, 2048, pointsize = 64)
    plotDeg(rrGraph)
    if (newDev) dev.off()
}

#Distribution of strength
figR10m <- function(rrGraph, newDev = TRUE) {
    # Remember that rrGraph must be the complete graph for the purpose of your research.
    E(rrGraph)$weight <- E(rrGraph)$RR
    oStrength <- strength(rrGraph, mode = "out")
    iStrength <- strength(rrGraph, mode = "in")
    n <- rrGraph
    if (newDev) png("out/figR10.png", 1344, 1344, pointsize = 48)
    browser()
    hist(oStrength, main = "Histogram of out strength", xlab = "out strength")
    if (newDev) dev.off()
    if (newDev) png("out/figR11.png", 1344, 1344, pointsize = 48)
qqnorm(oStrength, main = "Normal Q-Q plot of out strength")
    if (newDev) dev.off()
    plot.new()
    if (newDev) png("out/figR12.png", 1344, 1344, pointsize = 48)
    hist(iStrength, main = "Histogram of in strength", xlab = "in strength")
    if (newDev) dev.off()
    if (newDev) png("out/figR13.png", 1344, 1344, pointsize = 48)
qqnorm(iStrength, main = "Normal Q-Q plot of in strength")
    if (newDev) dev.off()
}
figSup2 <- function(completeGraph, subGraph) {
    png("out/figSup2.png", 512 * 3, 512 * 2, pointsize = 28)
    par(mfrow = c(2, 3), mar = c(5, 5, 1, 1), oma = c(0, 1, 1, 0))
    source("src/functions.R")
    plotDeg(subGraph, log = "y", fit.lines = TRUE)
    mtext("A", adj = -1/4)
    E(completeGraph)$weight <- E(completeGraph)$RR
    oStrength <- strength(completeGraph, mode = "out")
    iStrength <- strength(completeGraph, mode = "in")
    hist(oStrength, xlab = "out-strength", main = "Histogram of out-strength")
    mtext("B", adj = -1/4)
    qqnorm(oStrength, main = "Normal Q-Q plot of out-strength")
    mtext("C", adj = -1/4)
    plotDeg(subGraph, cumulative = TRUE, fit.lines = TRUE, log = "y")
    mtext("D", adj = -1/4)
    hist(iStrength, xlab = "in-strength", main = "Histogram of in-strength")
    mtext("E", adj = -1/4)
    qqnorm(iStrength, main = "Normal Q-Q plot of in-strength")
    mtext("F", adj = -1/4)
    dev.off()
}

#Top strength tables
tabR01 <- function() {
    topStro <- n[head(sort(n$strengtho, decreasing = TRUE, index.return = TRUE)$ix), c("name", "eng", "strengtho")]
    write.csv(topStro, "out/tabR01.csv", row.names = FALSE)
    library(xtable)
    toLatex(xtable(topStro))
}
tabR02 <- function() {
    topStri <- n[head(sort(n$strengthi, decreasing = TRUE, index.return = TRUE)$ix), c("name", "eng", "strengthi")]
    write.csv(topStri, "out/tabR02.csv", row.names = FALSE)
    browser("topStri[2,2] <- \"Parkinson disease\"")
    library(xtable)
    toLatex(xtable(topStri))
}

#Scatter plots
#Correlation between degree and strength
figR14 <- function() {
    if (!exists("igs")) source("src/functions.R")
    if (!exists("n")) source("src/nodes.R")
    igs
    if (!("dego" %in% names(n))) n <- updateDeg(n, igs)
    png("out/figR14.png", 1344, 1344, pointsize = 48)
    plot(n$dego, n$strengtho, xlab = "out degree", ylab = "out strength", main = "Correlation between out degree and strength")
    dev.off()
    png("out/figR15.png", 1344, 1344, pointsize = 48)
    plot(n$degi, n$strengthi, xlab = "in degree", ylab = "in strength", main = "Correlation between in degree and strength")
    dev.off()
}

#Lack of correlation between out degree and in degree
figR16 <- function() {
    png("out/figR16.png", 1536, 1536, pointsize = 48)
    if (!exists("n")) source("src/nodes.R")
    if (!("dego" %in% names(n))) {
	n <- updateDeg(n, igs)
        write(n)
    }
    pearson <- cor(n$dego, n$degi, use = "na.or.complete")
    plot(n$dego, n$degi, main = "Lack of correlation between out degree and in degree.", xlab = "out degree", ylab = "in degree", log = "xy")
    legend("topleft", legend = c(paste("r =", round(pearson, digits = 2)), expression(r ^ 2 == 0.24)))
    dev.off()
    cat("Pearson correlation coefficient between out degree and in degree :", pearson)
}

#Correlation and regression between out strength and in strength
figRes31 <- function(completeGraph) {
    png("out/figRes31.png", width = 1024, height = 1024, pointsize = 32)
    source("src/functions.R")
    n <- data.frame0(vertex_attr(completeGraph))
    library(igraph)
    strengtho <- strength(completeGraph, mode = "out")
    strengthi <- strength(completeGraph, mode = "in")
    pearson <- cor(strengtho, strengthi)
    lmCom <- lm(strengthi ~ strengtho, n)
    plot(strengtho, strengthi, main = "Correlation between out strength and in strength.", xlab = "out strength", ylab = "in strength", pch = 21, bg = n$ageCol, bty = "L")
    legend("bottomright", legend = c(
        paste("r =", round(pearson, digits = 2)), 
    #    paste("r ^ 2 =", round(pearson ^ 2, digits = 2)),
        paste0("y = ", round(lmCom$coefficients[2], digits = 2), "x + ", round(lmCom$coefficients[1], digits = 2))
    ), bty = "n")
    lines(c(min(strengtho), max(strengtho)), predict(lmCom, data.frame0(strengtho = c(min(strengtho), max(strengtho)))))
    #lines(c(0, 1614.641), c(1614.641/0.65, 0), lty = 2)
    #lines(c(0, 770.7876), c(770.7876 / 0.65, 0), lty = 2)
    #lines(c(0,900), c(672.5773, 900 * 0.65 + 672.5773), lty = 2)
    #lines(c(0, 1100), c(50.35198, 1100 * 0.65 + 50.35198), lty = 2)
    #rasterImage(figRes2leg2(plot = FALSE), 400, 1100, 500, 1200)
    dev.off()
}
#Multipanel correlation figure
figSup3 <- function(completeGraph, subGraph) {
    library(igraph)
    degi <- degree(subGraph, mode = "in")
    dego <- degree(subGraph, mode = "out")
    strengthi <- strength(completeGraph, mode = "in")
    strengthiSub <- strengthi[match(V(subGraph)$name, V(completeGraph)$name)]
    strengtho <- strength(completeGraph, mode = "out")
    strengthoSub <- strengtho[match(V(subGraph)$name, V(completeGraph)$name)]
    png("out/figSup3.png", width = 1024, height = 1024, pointsize = 32)
    par(mfrow = c(2, 2), mar = c(5, 5, 1, 1), oma = c(1, 1, 1, 1))
    plot(degi, strengthoSub, xlab = "in-degree", ylab = "in-strength", log = "x")
    legend("topleft", legend = paste("r =", round(cor(degi, strengthiSub), 2)))
    plot(strengtho, strengthi, xlab = "out-strength", ylab = "in-strength")
    legend("topleft", legend = paste("r =", round(cor(strengtho, strengthi), 2)))
    plot(dego, degi, xlab = "out-degree", ylab = "in-degree", log = "xy")
    legend("topleft", legend = paste("r =", round(cor(dego, degi), 2)))
    plot(dego, strengthoSub, xlab = "out-degree", ylab = "out-strength", log = "x")
    legend("topleft", legend = paste("r =", round(cor(dego, strengthoSub), 2)))
    mtext("out-degree", side = 1, outer = TRUE)
    mtext("in-degree", side = 2, outer = TRUE)
    mtext("in-strength", side = 3, outer = TRUE)
    mtext("out-strength", side = 4, outer = TRUE)
    dev.off()
}

#Examples tables from out-in strength graph
tabRes33 <- function(completeGraph, coef = 0.72) {
    source("src/functions.R")
    library(igraph)
    n <- data.frame0(vertex_attr(completeGraph))
    n$strengtho <- strength(completeGraph, mode = "out")
    n$strengthi <- strength(completeGraph, mode = "in")
    #Upper end of the regression line
    outdf <- n[head(sort(n$strengtho + coef * n$strengthi, decreasing = TRUE, index.return = TRUE)$ix),]
    outdf <- outdf[,c("name", "eng", "strengtho", "strengthi")]
    #write.csv(outdf, "out/tabR03.csv", row.names = FALSE)
    outdf$set <- "upper end"
    outdfm <- outdf
    #Lower end of the regression line
    outdf <- n[head(sort(n$strengtho + coef * n$strengthi, index.return = TRUE)$ix),]
    outdf <- outdf[,c("name", "eng", "strengtho", "strengthi")]
    outdf$set <- "lower end"
    outdfm <- rbind(outdfm, outdf)
    #Greatest offset above regression line
    outdf <- n[head(sort(n$strengthi - n$strengtho * coef, decreasing = TRUE, index.return = TRUE)$ix),]
    outdf <- outdf[,c("name", "eng", "strengtho", "strengthi")]
    outdf$set <- "offset above"
    outdfm <- rbind(outdfm, outdf)
    #Greatest offset below regression line
    outdf <- n[head(sort(n$strengthi - n$strengtho * coef, index.return = TRUE)$ix),]
    outdf <- outdf[,c("name", "eng", "strengtho", "strengthi")]
    outdf$set <- "offset below"
    outdfm <- rbind(outdfm, outdf)
    outdfm$strengtho <- round(outdfm$strengtho)
    outdfm$strengthi <- round(outdfm$strengthi)
    library(xtable)
    toLatex(xtable(outdfm))
}

#Correlation between strength and age
#I can plot it several correlations in one graph but since they have the same axes but it gets cluttered.
#In order to plot in one graph the greatest values should be used to generate axes.
#The greatest values are in strengthi and moa.
#But for now I am drawing separate graphs for visualization.
figR18 <- function () {
    png("out/figR18.png", 1344, 1344, pointsize = 48)
    plot(n$strengtho, n$moa, xlab = "out strength", ylab = "mean age of onset", pch = 21, bg = n$comCol)
    dev.off()
}
figR19 <- function () {
    png("out/figR19.png", 1344, 1344, pointsize = 48)
    plot(n$strengthi, n$moa, xlab = "in strength", ylab = "mean age of onset", pch = 22, bg = n$comCol)
    dev.off()
}
figR20 <- function () {
    png("out/figR20.png", 1344, 1344, pointsize = 48)
    plot(n$strengtho, n$ageMean, pch = 23, xlab = "out strength", ylab = "mean age of affected patients", bg = n$comCol)
    dev.off()
}
figR21 <- function () {
    png("out/figR21.png", 1344, 1344, pointsize = 48)
    plot(n$strengthi, n$ageMean, pch = 24, xlab = "in strength", ylab = "mean age of affected patients", bg = n$comCol)
    dev.off()
}
#Multipanel figure for correlation between strength and age
figSup4 <- function(completeGraph) {
    library(igraph)
    strengtho <- strength(completeGraph, mode = "out")
    strengthi <- strength(completeGraph, mode = "in")
    png("out/figSup4.png", 1028, 1028, pointsize = 24)
    par(mfrow = c(2, 2), mar = c(0, 0, 1, 1), oma = c(5, 5, 1, 0), xpd = NA)
    plot(V(completeGraph)$ageMean, strengthi, ylab = "in-strength", xaxt = "n", xlab = "")
    plot(V(completeGraph)$moa, strengthi, xaxt = "n", yaxt = "n", ylab = "", xlab = "")
    plot(V(completeGraph)$ageMean, strengtho, xlab = "mean age of affected", ylab = "out-strength")
    plot(V(completeGraph)$moa, strengtho, xlab = "mean age of onset", yaxt = "n", ylab = "")
    dev.off()
}

#Just to clarify the difference between mean age of onset and mean age of affected patients.
figR22 <- function() {
    png("out/figR22.png", 512, 512)
    plot(n$ageMean, n$moa, xlab = "mean age of affected patients", ylab = "mean age of onset", pch = 21, bg = n$comCol)
    dev.off()
}

#They all segregate communities to some degree, but which segregates it best?
#Does any segregate communities better than out/in strength?
figR23 <- function() {
    png("out/figR23.png")
    plot(n$strengtho, n$strengthi, main = "Correlation between out strength and in strength.", xlab = "out strength", ylab = "in strength", pch = 21, bg = n$comcol)
    dev.off()
}
#Does out/in strenth segregate categories as well?
figR24 <- function() {
    png("out/figR24.png", width = 2048, height = 2046, pointsize = 48)
    if(!exists("n")) source("src/nodes.R")
    plot(n$strengtho, n$strengthi, main = "Does out strength and in strength segregate categories?", xlab = "out strength", ylab = "in strength", pch = 21, bg = n$catCol)
    dev.off()
}

#Graph layout using strength as coordinates
figRes321 <- function(completeGraph, rrGraph) {
    source("src/functions.R")
    V(rrGraph)$x <- strength(completeGraph, V(completeGraph)[V(rrGraph)$name], mode = "out")
    V(rrGraph)$y <- strength(completeGraph, V(completeGraph)[V(rrGraph)$name], mode = "in")
    V(rrGraph)$color <- V(rrGraph)$comCol
    png("out/figRes321.png", 514, 514, pointsize = 20)
    plot(rrGraph, xlab = "out strength", ylab = "in strength", vertex.size = 5, vertex.label = NA, edge.arrow.size = 0.5)
    axis(1, at = c(-1, 1), labels = c("", ""))
    axis(2, at = c(-1, 1), labels = c("", ""))
    dev.off()
}
figRes32 <- function(completeGraph, rrGraph) {
    source("src/functions.R")
    V(rrGraph)$x <- strength(completeGraph, V(completeGraph)[V(rrGraph)$name], mode = "out")
    V(rrGraph)$y <- strength(completeGraph, V(completeGraph)[V(rrGraph)$name], mode = "in")
    V(rrGraph)$color <- V(rrGraph)$comCol
    png("out/figRes32.png", 1536, 768, pointsize = 24)
    par(mfrow = c(1, 2), mar = c(5, 5, 1, 1), oma = c(0, 1, 1, 0))
    plot(rrGraph, xlab = "out strength", ylab = "in strength", vertex.size = 5, vertex.label = NA, edge.arrow.size = 0.5)
    axis(1, at = c(-1, 1), labels = c("", ""))
    axis(2, at = c(-1, 1), labels = c("", ""))
    mtext("A", adj = -1/4)
    V(rrGraph)$color <- V(rrGraph)$catCol
    plot(rrGraph, xlab = "out strength", ylab = "in strength", vertex.size = 5, vertex.label = NA, edge.arrow.size = 0.5)
    axis(1, at = c(-1, 1), labels = c("", ""))
    axis(2, at = c(-1, 1), labels = c("", ""))
    mtext("B", adj = -1/4)
    dev.off()
}

#Enrichment analysis of categories in the communities
#Duplicated in EnrichComCat
tabR7 <- function() {
    #Returns LaTeX table.

    if (!("rr" %in% search())) source("src/functions.R")
    if (!exists("n")) source("src/nodes.R")
    if (!exists("ca")) source("src/cat.R")
    library(xtable)

    mask <- !is.na(n$com)
    com <- n$com[mask]
    cat <- as.integer(as.roman(n$cat[mask]))
    coms <- sort.int(unique(com))
    cats <- sort.int(unique(cat))

    #Get fisher test data all at once.
    results <- outer(coms, cats, Vectorize(function(x, y) {
        fisher.test(com == x, cat == y, alternative = "greater")
    }, c("x", "y"), SIMPLIFY = FALSE))
    dimnames(results) <- list(community = coms, category = cats)
    
    #Format into table suitable for publishing
    enrichments <- lapply(1:4, function(x) {
        p.values <- vapply(results[x,], function(y) { 
            y$p.value 
        }, 0)
        significants <- which(p.values < 0.05)
        enriched <- cats[significants]
        list(text = paste0(ca[enriched, "name"], ".", ca[enriched, "desc"]), p.value = p.values[significants])
    })
    community = rep(1:4, vapply(enrichments, function(x) { length(x$text) }, 0))
    enrichment <- unlist(lapply(enrichments, function(x) { x$text }))
    p.value <- unlist(lapply(enrichments, function(x) { x$p.value }))
    tbl <- data.frame(community, enrichment, p.value)
    toLatex(xtable(tbl, digits = c(0, 1, 0, -1)))
}

#Section for stratification studies.
#Stratification results for import into Cytoscape
tabR8 <- function() {
    if (!exists("igstrat")) load("met/strata.RData")
    st1 <- as_data_frame(igstrat, "edges")
    names(st1)[match(c("from", "to", "strat.RR", "strat.pValue", "strat.age", "strat.sex"), names(st1))] <- c("riskFactor.kcd", "disease.kcd", "RR", "pValue", "age", "sex")
    ageGrpTxt <- c("infant", "1-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+")
    st1$age <- ageGrpTxt[st1$age + 1]
    sexTxt <- c("M", "F")
    st1$sex <- sexTxt[st1$sex]
    st1$pValue <- sprintf("%1.1e", st1$pValue)
    st1$pValue <- replace(st1$pValue, st1$pValue == "0.0e+00", "< 2.2e-16")
    library(xtable)
    capture.output(toLatex(xtable(st1[c("riskFactor.kcd", "disease.kcd", "RR", "pValue", "age", "sex")], digits = 1)), file = "out/tabR8.tex", type = "output")
}

#Calculation of the modularity score of the communities in the stratified sub-network.
resR1 <- function() {
    if (!exists("igstrat")) load("met/strata.RData")
    library(igraph)
    if(!exists("n")) source("src/nodes.R")
    modularity(igstrat, n[names(V(igstrat)), "com"])
}

# Detection of communities in the stratified network using the random walktrap algorithm detects a community structure with a modularity of 0.57.
resR2 <- function() {
    if (!exists("igstrat")) load("met/strata.RData")
    if (!exists("cluster_walktrap")) library(igraph)
    cluster_walktrap(igstrat)
}

# Comparison of newly discovered community structure and the community membership from the original communities gives an index of 1.22
resR3 <- function() {
    if (!exists("igstrat")) load("met/strata.RData")
    if (!exists("n")) source("src/nodes.R")
    compare(V(igstrat)$com0, V(igstrat)$com)
}

# Figure for comparison of communities in the stratification network
figR25 <- function() {
    if (!exists("igstrat")) load("met/strata.RData")
    library(igraph)
    com <- V(igstrat)$com
    com0 <- V(igstrat)$com0
    png(file = "out/figR25.png")
    ma <- Reduce(rbind2, lapply(1:3, function(x) {
        ta <- table(factor(com[com0 == x], levels = 1:29))
        as.vector(ta / sum(ta))
    }), NULL)
    barplot(ma, beside = TRUE, col = rainbow(3), names.arg = 1:29, legend.text = 1:3)
    dev.off()
    return(NULL)
}
# I need another method of visualization, one that highlights just the major communities,
# also showing the size of the communities.
figR26 <- function() {
    if (!exists("igstrat")) load("met/strata.RData")
    library(igraph)
    com <- V(igstrat)$com
    sortix <- sort(as.numeric(table(com)), decreasing = TRUE, index.return = TRUE)$ix
    com <- match(com, sortix)
    com0 <- V(igstrat)$com0
    ma <- Reduce(rbind2, lapply(1:20, function(x) {
        ta <- table(factor(com[com0 == x], levels = 1:29))
        as.vector(ta)
    }), NULL)
    png(file = "out/figR26.png", width = 2048, height = 2048, pointsize = 64)
    barplot(ma, col = c(rainbow(4), gray.colors(16)), names.arg = c(1:4, rep("", 25)), legend.text = c("chronic debilitation cluster", "women's disease cluster", "neoplasm-digestive system cluster", "infectious disease cluster"), xlab = "new cluster number", ylab = "new cluster size")
    dev.off()
    return(NULL)
}

# Comparison of transitivity (clustering coefficient) of the major communities before and after stratification.
resR4 <- function() {
    Reduce(rbind2, sapply(list(igs, igstrat), function(ig) {
        vapply(1:4, function(x) {
            transitivity(subgraph(ig, n[n$com == x, "name"]))
        }, 0)
    }), NULL)
}

# Comparison of transitivity (clustering coefficient) of the major communities before and after stratification.
tabR9 <- function() {
    if (!exists("igs")) source("src/functions.R")
    if (!exists("igstrat")) load("met/strata.RData")
    if (!exists("n")) source("src/nodes.R")
    vec <- sapply(list(igs, igstrat), function(ig) {
        vapply(1:4, function(x) {
            transitivity(induced_subgraph(ig, intersect(names(V(ig)), n[n$com == x, "name"])))
        }, 0)
    })
    restab <- matrix(vec, nrow = 4, dimnames = list(community = c(1, 2, 3, 4), stratification = c("before", "after")))
    library(xtable)
    toLatex(xtable(restab))
}
figStr <- function(graphs) {
    mapply(function(graph, name) {
        png(paste0('out/figStr', name, '.png'))
        plot(graph, layout = layout_with_fr, vertex.size = 8, vertex.label = NA,
            vertex.color = V(defaultGraph)$comCol[match(V(graph)$name, V(defaultGraph)$name)],
            edge.arrow.size = 1 / 2)
        dev.off()
    }, graphs, names(graphs))
}
