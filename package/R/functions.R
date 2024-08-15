library(igraph)
if (any(search() == "rr")) {
    detach(rr)
}
rr <- list(

poolSick = function(T20, T40, year = NULL, filter = NULL) {
    # This function is the main data cruncher of this project.
    # It extracts incidences from files in the raw dir.
    # T20 should be a list of T20 file paths in chronological order.
    # T40 should be a list of corresponding T40 file paths.
    # year (optional) should be an integer vector corresponding year values.
    # If year is given, it tries to give the year of onset in the return value.
    # It will take long, despite optimizations.
    # But if you intend to later calculate the mean age of onset, 
    # you will need to pass the year parameter to this function.
    # I can not think of any faster way to calculate the mean age of onset.
    # filter (optional) is a list of the form colName = vector.
    # You can filter using PERSON_ID or SICK_SYM.
    l <- length(T20)
    if (l != length(T40)) stop("T20, T40 differ in length")
    if (!is.null(year) && l != length(year)) stop("bad year length")
    kC <- mapply(function (T20, T40, year) {
        l <- lapply(c(k2 = T20, k4 = T40), read.csv, colClasses = "character")
        k2 <- l$k2
        k4 <- l$k4
        k2$PERSON_ID <- as.integer(k2$PERSON_ID)
        if (!is.null(filter) && "PERSON_ID" %in% names(filter)) {
            k2 <- subset(k2, PERSON_ID %in% filter$PERSON_ID)
        }
        # Using match is faster than using join in the plyr package,
        # because the join spends much time sorting which we do not need.
        m <- match(k4$KEY_SEQ, k2$KEY_SEQ)
        r <- data.frame(PERSON_ID = c(k2$PERSON_ID, k2$PERSON_ID[k2$SUB_SICK != ""], k2$PERSON_ID[m[!is.na(m)]]), SICK_SYM = c(k2$MAIN_SICK, k2$SUB_SICK[k2$SUB_SICK != ""], k4$SICK_SYM[!is.na(m)]), stringsAsFactors = FALSE)
        if (length(year) == length(T20)) {
            r$year = rep(year, dim(r)[[1]])
        }
        # Faster with filtering inside the loop than outside the loop?
        r$SICK_SYM <- substr(r$SICK_SYM, 1, 3)
        if (!is.null(filter) && "SICK_SYM" %in% names(filter)) {
            r <- subset(r, SICK_SYM %in% filter$SICK_SYM)
        }
        # Faster without the following reduction?
        r <- unique(r)
        cat(year, ", ", sep = "") #debug
        return(r)
    }, T20, T40, year, SIMPLIFY = FALSE)
    su <- function(l, name) { unlist(lapply(l, function(x) { x[[name]] })) }
    if (is.null(year)) {
        # Use the faster unique if incidence is not required.
        return(unique(data.frame(PERSON_ID = su(kC, "PERSON_ID"), SICK_SYM = su(kC, "SICK_SYM"), stringsAsFactors = FALSE)))
    }
    # 'join, by, aggregate, tapply, by' are all slow when working with multiple highly variable factors.
    # Internally, 'interaction' is a bottleneck when using highly variable factors,
    # and is used by 'split' which is used by 'tapply' which is used by 'by'.
    # Therefore it is more efficient to generate a single highly variable factor using 'paste0'.
    # Even then the main bottleneck is 'tapply'.
    # It is even more efficient to use the same internal mechanism that 'unique' uses.
    Reduce(function(l, r) {
        p <- list(l, r)
        class(p) = c(class(p), "partdf")
        u <- unlist(p)
        d <- duplicated(u[1:2])
        u[!d,]
    }, kC, data.frame(PERSON_ID = integer(), SICK_SYM = character(), year = integer()))
},
J = function(year = 2004L) {
    if (year < 2002L || year > 2013L) stop("Bad year specification.")
    f <- sprintf("raw/jk-%02d.csv", year - 2000L)
    table <- read.csv(f, colClasses = c("integer", "character", "integer", "integer"))
    table$PERSON_ID <- as.integer(table$PERSON_ID)
    return(table)
},
F = function() {
    f <- read.csv("met/fetch.csv", colClasses = c("character", "character", "integer"))
    f$PERSON_ID <- as.integer(f$PERSON_ID)
    f
},
`$.partdf` = function(l, name) {
    unlist(lapply(l, function(x) { x[[name]] }))
},
unlist.partdf = function(x, recursive, use.names) {
    n <- names(x[[1]])
    d <- data.frame(lapply(n, function(n) { `$.partdf`(x, n) }), stringsAsFactors = FALSE)
    names(d) <- n
    return(d)
},
Place = function(x, bins) {
    #If bins are lower bounds of bins sorted in ascending order,
    #this function returns the index of the bin that each x belongs to.
    #Remember to coerce factors into characters before use.
    r <- Map(function(x) {
        Position(function(bin) { bin > x }, bins, nomatch = length(bins) + 1) - 1
    }, x)
    return(unlist(r))
},
ig = function(x = p, RR = 4, pFdr = 0.001) {
    if(!is(x, "data.frame")) {
        stop("Bad pValues table. Perhaps you should p <- p().")
    }
    R <- RR
    pFd <- pFdr
    x <- subset(x, RR > R & pFdr < pFd)
    x$weight = x$RR
    graph_from_data_frame(x)
},
write = function(...) UseMethod("write"),
sizes = function(...) UseMethod("sizes"),
updat = function(...) UseMethod("updat"),
sizes.default = function(v) {
    to <- length(v)
    f <- integer()
    v <- as.vector(v)
    t <- table(v)
    if (is(v, "character")) {

    f <- as.integer(t)
    names(f) <- names(t)

    } else {
        f[as.integer(names(t))] <- as.integer(t)
    }
    if(sum(f, na.rm = TRUE) != to) {
        stop("Sum of sizes does not equal total.")
    }
    return(f)
},
sizes.communities = function(community) {
    get("sizes", "package:igraph")(community)
},
modularit.default = function(x) attr(x, "modularity"),
age = function(x) { 
        ifelse(x == 0, 0, ifelse(x == 1, 2.5, 5 * x - 2.5)) 
},
plotDeg = function(igraph, fit.lines = FALSE, cumulative = FALSE, log = "xy") {
    #fit.lines only supported for exponential distribution.
    #I.e. fit.lines ignored for anything but log = "y".
    degDistO <- degree_distribution(igraph, mode = "out", cumulative = cumulative)
    degDistI <- degree_distribution(igraph, mode = "in", cumulative = cumulative)
    if (cumulative) {
        type <- "l"
        main <- "Cumulative degree distribution"
    } else {
        type <- "p"
        main <- "Degree distribution"
    }
    #In-degree must be plotted first because max(in-degree) is greater.
    plot(degDistI[2:length(degDistI)], log = log, ylab = "p(k)", xlab = "degree k", type = type, pch = 17, lty = 2, main = main)
    if (cumulative) {
        lines(1:(length(degDistO) - 1), degDistO[2:length(degDistO)])
        if (fit.lines && log == "y") {
            p <- 1 / (gsize(igraph) / gorder(igraph) + 1)
            q <- 1 - p
            lines(c(0, 80), q ^ par("usr")[1:2], lty = 3)
            legend("bottomleft", legend = c("out", "in", "fit"), lty = c(1, 2, 3))
        } else {
            legend("bottomleft", legend = c("out", "in"), lty = c(1, 2))
        }
    } else {
        points(degDistO[2:length(degDistO)])
        if (fit.lines && log == "y") {
            p <- gorder(igraph) / gsize(igraph)
            q <- 1 - p
            lines(c(0, 80), q ^ par("usr")[1:2] * p, lty = 1)
            legend("bottomleft", legend = c("out", "in", "fit"), pch = c(1, 17, NaN), lty = c(0, 0, 1))
        } else {
            legend("bottomleft", legend = c("out", "in"), pch = c(1, 17))
        }
    }
}
)
T20Files <- c("raw/T20-02.csv", "raw/T20-03.csv", "raw/T20-04.csv", "raw/T20-05.csv", "raw/T20-06.csv", "raw/T20-07.csv", "raw/T20-08.csv", "raw/T20-09.csv", "raw/T20-10.csv", "raw/T20-11.csv", "raw/T20-12.csv", "raw/T20-13.csv")
T40Files <- c("raw/T40-02.csv", "raw/T40-03.csv", "raw/T40-04.csv", "raw/T40-05.csv", "raw/T40-06.csv", "raw/T40-07.csv", "raw/T40-08.csv", "raw/T40-09.csv", "raw/T40-10.csv", "raw/T40-11.csv", "raw/T40-12.csv", "raw/T40-13.csv")
kcdCategories <- c(
    "Infectious and parasitic",
    "Neoplasms",
    "Blood and the immune mechanism",
    "Endocrine",
    "Mental and behavioral",
    "Nervous system",
    "Eye and adnexa",
    "Ear and mastoid process",
    "Circulatory",
    "Respiratory",
    "Digestive",
    "Skin and subcutaneous tissue",
    "Musculoskeletal and connective tissue",
    "Genitourinary",
    "Pregnancy",
    "Perinatal",
    "Congnital malformations",
    "Symptoms",
    "Injury",
    "External causes",
    "Influencing factors",
    "Special codes")
manyTo1Merge <- function(many, one, drop = FALSE) {
    # Works only for data frames.
    # Compatibility for matrices was dropped because the data frame implementation of the matrix subsetting operation `[` is memory-intensive,
    # and kept running out of memory.
    # For data frames, manipulations as lists is more efficient.
    oneNames <- names(one)
    manyNames <- names(many)
    by <- intersect(manyNames, oneNames)[1]
    newOneNames <- setdiff(oneNames, by)
    many <- unclass(many)
    one <- unclass(one)
    matches <- match(many[[by]], one[[by]])
    onePart <- lapply(one[newOneNames], `[`, matches)
    combined <- c(many, onePart)
    if (drop == TRUE) combined <- combined[!is.na(matches),]
    names(combined) <- c(manyNames, newOneNames)
    class(combined) <- 'data.frame'
    return(combined)
}
data.frame0 <- function(...) data.frame(..., stringsAsFactors = FALSE, check.names = FALSE, fix.empty.names = FALSE)
exportCyto <- function(dataFrame, fileName, mac = TRUE) {
    if(mac) {
        write.table(dataFrame, fileName, sep = "\t", quote = FALSE, row.names = FALSE)
	return(NULL)
    }
    write.csv(dataFrame, row.names = FALSE)
}
subset.incidences <- function(incidences, syms = NULL, pids = NULL, years = NULL, row.names = FALSE) {
    if (!all(c('SICK_SYM', 'PERSON_ID', 'year') %in% names(incidences))) stop('Bad INCIDENCES parameter passed to subset.incidences.');
    ix <- Reduce(function(ix, searchPair) {
        haystack <- searchPair[[1]]
        needles <- searchPair[[2]]
        if(is.null(needles)) return(ix)
        mask <- haystack %in% needles
        foundIx <- .Internal(which(mask))
        if (is.null(ix)) return(foundIx)
        return(intersect(ix, foundIx))
    }, list(list(incidences$SICK_SYM, syms), list(incidences$PERSON_ID, pids), list(incidences$year, years)), NULL)
    returnData <- data.frame0(lapply(incidences, `[`, ix))
    if(row.names) row.names(returnData) <- row.names(incidences)[ix]
    return(returnData)
}
incidenceDemogMerge <- function(incidences, demogData) {
    # An often-used merge used for various purposes.
    # incidences; demogData
    return(manyTo1Merge(incidences, demogData))
}
summarySymDemog <- function(incidences, demogData, syms, pids = NULL, years = 2002:2003L) {
    #At first I wrote this function to get age-related info,
    #But the same merge is useful for getting various demographic info.
    incidences; demogData; syms
    incidences <- subset.incidences(incidences, syms = syms, pids = pids, years = years)
    fj <- incidenceDemogMerge(incidences, demogData)
    factorBySyms <- function(x) factor(x, levels = syms)
    young <- as.integer(table(factorBySyms(fj$SICK_SYM[fj$AGE_GROUP <= 6])))
    middle <- as.integer(table(factorBySyms(fj$SICK_SYM[fj$AGE_GROUP >= 7 & fj$AGE_GROUP <= 12])))
    old <- as.integer(table(factorBySyms(fj$SICK_SYM[fj$AGE_GROUP >= 13])))
    num <- as.vector(table(factorBySyms(fj$SICK_SYM)))
    ageCol <- mapply(function(young, middle, old, num) {
        col <- sprintf("#%02X%02X%02X", round(young / num * 255), round(middle / num * 255), round(old / num * 255))
        return(col)
    }, young, middle, old, num)
    #mean age of affected
    age <- age(fj$AGE_GROUP)
    ageMean <- tapply(age, factorBySyms(fj$SICK_SYM), mean)
    #counts and portion of sexes
    male <- as.vector(table(factorBySyms(fj$SICK_SYM[fj$SEX == 1])))
    female <- as.vector(table(factorBySyms(fj$SICK_SYM[fj$SEX == 2])))
    maleP = male / num
    return(data.frame(num, young, middle, old, ageCol, ageMean, male, female, maleP = male / (male + female), row.names = syms, stringsAsFactors = FALSE))
}
factorDemog <- function(incidences, demogData, syms, factor, years = 2002:2003L) {
    #Given kcd codes grouped by a factor factor,
    #Returns demographic profiles 
    #as a data frame with columns factor (the factor), PERSON_ID, SEX, AGE_GROUP
    #Formerly <demog>
    #It is more efficient to merge incidences with the factor first, as it can remove duplicate PERSON_IDs before merging with demogData.
    kf <- data.frame0(SICK_SYM = syms, factor = factor)
    m <- manyTo1Merge(subset.incidences(incidences, syms = syms, years = years), kf)
    m <- unique(m[c("PERSON_ID", "factor")])
    m <- manyTo1Merge(m, demogData, drop = TRUE)
    return(m)
}
attach(rr)
delayedAssign("p", read.csv("met/pValue.csv", stringsAsFactors = FALSE))
#Assigning promises for eligibility tables.
#Commenting this section out because the only use for all eligibility tables is so far only for calculating mean age of onset.
#(function() {
#    varNames <- sprintf("j%02d", 2L:13L)
#    mapply(function(varName, year) {
#        delayedAssign(varName, J(year), assign.env = .GlobalEnv)
#    }, varNames, 2002L:2013L, SIMPLIFY = FALSE)
#    return(0L)
#})()
delayedAssign("igs", ig())
delayedAssign("k", read.csv("met/selectedKcd.csv", colClasses = "character"))
delayedAssign("f", F())
