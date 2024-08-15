#Dependencies
#functions.R
if("classN" %in% search()) detach(classN)
classN <- list(
N = function(new = FALSE) {
    if (new) {
        #This will not work yet

datCom <- read.csv("pValue.csv")
library(igraph)
graSta <- graph_from_data_frame(subset(datCom, RR > 4 & pFdr < 0.001))
E(graSta)$weight <- E(graSta)$RR
clu <- cluster_walktrap(graSta)
#Cluster info should be loaded not calculated.
k <- k()
community <- membership(clu)[match(k$kcd, names(membership(clu)))]
#Onset info
load("onset.RData")
names(moa) <- k$kcd
#Onset info should not be modified.
save(inc, moa, file = "onset.RData")

    } else {
        n <- read.csv("met/nodes.csv", stringsAsFactors = FALSE)
    }
    row.names(n) <- n$name
    structure(n, class = c("N", "data.frame"))
},
updat.N = function(n, x) {
    #This function is used to update the information in an N object
    #with information from another class.
    if(is(x, "Ca")) {
    #Calculate categories.
        ca <- x
        n$cat <- as.roman(Place(as.character(n$name), as.character(ca$kcd)))
        n$catCol <- as.character(ca$col[n$cat])
        return(n)
    } else if (is(x, "Co")) {
        co <- attr(x, "communities")
        co <- x[as.character(membership(co)[match(n$name, co$names)]),]
        n$com <- co$id
        n$comCol <- co$col
        return(n)
    } else stop(paste("No method available for", class(x)))
},
write.N = function(n, file = "met/nodes.csv", ...) {
    #If possible, update nodes.csv with current information.
    write.csv(n, file, row.names = FALSE, ...)
},
enrichedCategories = function(n, kcd) {
    #kcd: a vector of kcd values
    #n: an N object
    #value: roman
    lk <- length(kcd)
    lo <- length(n$cat) - lk
    ta <- table(n$name %in% kcd, n$cat)
    ma <- mapply(function(v1, v2) {
        ma <- matrix(c(v1, lk - v1, v2, lo - v2), nrow = 2, dimnames = list(c("cat", "notCat"), c("com", "notCom")))
        fisher.test(ma, alternative = "greater")
    }, ta["TRUE",], ta["FALSE",], SIMPLIFY = FALSE)
    ma <- mapply(function(name, fisher) { ifelse(fisher$p.value < 0.05, name, NA) }, names(ma), ma)
    ma <- ma[!is.na(ma)]
    ma <- as.roman(ma)
    return(ma)
},
enrichComCat = function(n, coms = 1:4) {
    #kcd: a vector of kcd values
    #n: an N object
    #value: roman
    ta <- table(n$com, n$cat)
    cow <- colSums(ta)
    numNodes <- sum(cow)
    retur <-lapply(coms, function(com, tabl, catSize, totNum) {
        #For each community
        #Vector of category sizes inside the community.
        taTrue <- tabl[com,]
        numCom <- sum(taTrue)
        #Vector of category sizes outside the community.
        taFalse <- catSize - taTrue

        ma <- mapply(function(catCom, catNotCom, totCom, totNotCom) {
            notCatCom <- totCom - catCom
            notCatNotCom <- totNotCom - catNotCom
            ma <- matrix(c(catCom, notCatCom, catNotCom, notCatNotCom), nrow = 2, dimnames = list(c("cat", "notCat"), c("com", "notCom")))
            fisher.test(ma, alternative = "greater")
        }, taTrue, taFalse, MoreArgs = list(totCom = numCom, totNotCom = totNum - numCom), SIMPLIFY = FALSE)
        #Category names with NAs
        ma <- mapply(function(name, fisher) { ifelse(fisher$p.value < 0.05, name, NA) }, names(ma), ma)
        #Category names (character)
        ma <- ma[!is.na(ma)]
        #Category names (roman)
        ma <- as.roman(ma)
        return(ma)

    }, ta, cow, numNodes)
    return(retur)
},
modularity.default = function(x) attr(x, "modularity"),
demog = function(kcd, factor, demogData) {
    #Given kcd codes grouped by a factor factor,
    #Returns demographic profiles (sex and age in 2002) of patients who were 
    #eligible in 2002 and
    #afflicted in 2002 to 2003 with diseases in each subset kcd,
    #as a data frame with columns factor (the factor), PERSON_ID, SEX, AGE_GROUP
    inc <- subset(f, cohort == 1)
    kf <- data.frame(kcd, factor, stringsAsFactors = FALSE)
    m <- merge(inc, kf, by.x = "MAIN_SICK", by.y = "kcd", sort = FALSE)
    m <- unique(m[c("PERSON_ID", "factor")])
    m <- merge(m, j, sort = FALSE)
    return(m)
},
demogCom = function(N) {
    #Returns demographic profiles (sex and age in 2002) of patients who were
    #eligible in 2002 and 
    #afflicted in 2002 to 2003 with diseases in each disease community,
    #as a data frame with columns the factor (the community), PERSON_ID, SEX, AGE_GROUP
    demog(N$name, N$com)
},
updateStrength = function(N, igraph) {
    strengtho <- strength(igraph, mode = "out")
    strengthi <- strength(igraph, mode = "in")
    names <- names(strengtho)
    mapping <- match(N$name, names)
    N$strengtho <- strengtho[mapping]
    N$strengthi <- strengthi[mapping]
    return(N)
},
updateDeg = function(N, igraph) {
    dego <- degree(igraph, mode = "out")
    degi <- degree(igraph, mode = "in")
    names <- names(dego)
    mapping <- match(N$name, names)
    N$dego <- dego[mapping]
    N$degi <- degi[mapping]
    return(N)
},
updateEng = function(N) {
    mapping <- match(N$name, k$kcd)
    N$eng <- k$eng[mapping]
    return(N)
},
updateAge = function(N) {
    #Unlike before, I am using 2002-2003 for ageCol because 
    #I suspect that eligibility data for later years migght not represent the population at that time
    #Maybe I should be using 2002 only.
    #At first I wrote this function to get age-related info,
    #But the same merge is useful for getting various demographic info.
    f <- subset(f, cohort == 1 & MAIN_SICK %in% N$name)[,c("PERSON_ID", "MAIN_SICK")]
    fj <- merge(f, j, by = "PERSON_ID", sort = FALSE, stringsAsFactors = FALSE)
    young <- subset(fj, AGE_GROUP <= 6)
    middle <- subset(fj, AGE_GROUP >= 7 & AGE_GROUP <= 12)
    old <- subset(fj, AGE_GROUP >=13)
    young <- table(young$MAIN_SICK)
    middle <- table(middle$MAIN_SICK)
    old <- table(old$MAIN_SICK)
    ageNums <- list(young = young[N$name], middle = middle[N$name], old = old[N$name])
    ageNums <- lapply(ageNums, function(v) { ifelse(is.na(v), 0, v) })
    ageCol <- mapply(function(young, middle, old) {
        max <- max(young, middle, old)
        col <- sprintf("#%02X%02X%02X", round(young / max * 255), round(middle / max * 255), round(old / max * 255))
        return(col)
    }, ageNums$young, ageNums$middle, ageNums$old)
    N$ageCol <- ageCol
    #mean age of affected
    age <- age(fj$AGE_GROUP)
    meanAge <- tapply(age, fj$MAIN_SICK, mean)
    N$ageMean <- meanAge[N$name]
    #counts and portion of sexes
    fjm <- subset(fj, SEX == 1)
    fjf <- subset(fj, SEX == 2)
    N$male <- table(fjm$MAIN_SICK)[N$name]
    N$male[is.na(N$male)] <- 0
    N$female <- table(fjf$MAIN_SICK)[N$name]
    N$female[is.na(N$female)] <- 0
    N$maleP <- N$male / (N$male + N$female)
    return(N)
},
updateAgeCol = function(N) updateAge(N)
#Updates ageCol and mean age.
#Synonymous to updateAge
,
exportCytoNodes = function(N) exportCyto(N, "out/nodes.tsv")
)
attach(classN)
rm(classN)
n <- N()
