#! /usr/bin/Rscript
# A script for generating contingency tables per stratum.
# Usage: makeContingency.R --args [sex] [age]
# If neither sex nor age is given, the entire cohort is used.
# [sex] and [age] must be valid sex and age codes.
# Often the data reading portion is time consuming.
# For several strata, consider using makeContingency-h.R and makeContingency-t.R.
kFet <- read.csv("fetch.csv", colClasses = c("character", "character", "integer"))
kFet$PERSON_ID <- as.integer(kFet$PERSON_ID)
kCodes <- read.csv("selectedKcd.csv")[[1]]
# Do not filter kFet by kCodes as this changes the cohort.
# Do not:
# kFet <- subset(kFet, MAIN_SICK %in% kCodes)
c <- as.integer(commandArgs(trailingOnly = TRUE))
if (length(c) > 0) {
    j <- read.csv("jk-02.csv", colClasses = "numeric")
    j <- data.frame(lapply(j, function (x) { as.integer(x) }))
    kFet <- subset(kFet, PERSON_ID %in% j$PERSON_ID[j$SEX == c[1] & j$AGE_GROUP == c[2]])
}
kRem <- unique(kFet$PERSON_ID[kFet$cohort == 1])
repeat {
    if (!file.exists("pid")) {
        kPre <- integer()
    } else kPre <- as.integer(scan("pid"))
    kRem <- setdiff(kRem, kPre)
    cat(length(kRem), "left.\n")
    if (length(kRem) == 0) {
        break
    } else if (length(kRem) < 10000) {
        kPer <- kRem
    } else {
        kPer <- sample(kRem, 10000)
    }
    kWor <- subset(kFet, PERSON_ID %in% kPer)
    kIni <- subset(kWor, cohort == 1)
    kFol <- subset(kWor, cohort == 0)
    kLen <- length(kCodes)
    if (length(kPre) == 0) {
        kCohort <- rep(0, kLen)
        kE <- rep(0, kLen ^ 2)
        kEcD <- kED <- kE
    } else {
        kRr <- read.csv("rr.csv")
        kCohort <- kRr[["C"]][1:kLen]
        kE <- kRr[["E"]]
        kED <- kRr[["ED"]]
        kEcD <- kRr[["EcD"]]
    }
    dim(kE) <- dim(kED) <- dim(kEcD) <- c(kLen, kLen)
    kIniMas <- rep(FALSE, kLen)
    names(kIniMas) <- kCodes
    kFinMas <- kIniMas
    library(plyr)
    
    # Tally matrix algorithm. Because bootstrapping takes too long.
    for (i in 1:length(kPer)) {
        kFinMas[] <- FALSE
        kPerson <- kPer[i]
        if(i %in% 2^(1:20)) cat("Currently on", i, "th person\n")
        # Applying to cohort count
        kDel <- kIni$PERSON_ID == kPerson
        kIniMas <- kCodes %in% kIni$MAIN_SICK[kDel]
        kIni <- subset(kIni, !kDel)
        kCohort[!kIniMas] <- kCohort[!kIniMas] + 1
        # Applying to exposure count
        kE[!kIniMas, kIniMas] <- kE[!kIniMas, kIniMas] + 1
        # Applying to exposure disease count
        kDel <- kFol$PERSON_ID == kPerson
        kFinMas <- kCodes %in% kFol$MAIN_SICK[kDel]
        kFol <- subset(kFol, !kDel)
        kFinInc <- kFinMas & !kIniMas
        kED[kFinInc, kIniMas] <- kED[kFinInc, kIniMas] + 1
        # Applying to non-exposed disease count
        kEcD[kFinInc, !kIniMas] <- kEcD[kFinInc, !kIniMas] + 1
    }
    write.csv(data.frame(C = rep(kCohort, kLen), E = as.vector(kE), ED = as.vector(kED), EcD = as.vector(kEcD)), file = "rr.csv", row.names = FALSE)
    write(kPer, file = "pid", append = TRUE)
}
