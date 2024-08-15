#! /usr/bin/Rscript
mh <- function () {

kFet <- read.csv("fetch.csv", colClasses = c("character", "character", "integer"))
kFet$PERSON_ID <- as.integer(kFet$PERSON_ID)
kCodes <- read.csv("selectedKcd.csv", colClasses = "character")
# Do not filter kFet by kCodes as this changes the cohort.
# Do not:
# kFet <- subset(kFet, MAIN_SICK %in% kCodes)
j <- read.csv("jk-02.csv", colClasses = "numeric")
j <- data.frame(lapply(j, function (x) { as.integer(x) }))
r = list(f = kFet, k = kCodes, j = j)
browser()
r

}
