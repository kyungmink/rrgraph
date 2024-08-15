# The following algorithm is old, and based on looking up the age of patients at each year of incidence from raw jk files.
# This does not add to the accuracy of the patients ages, because they are still coded in 5-year blocks.
# The newInc function utilizes a single jk file taken at the time of cohort selection.
oldInc <- function() {
# inc stands for incidence.
# I am not counting the incidences in the first two years because they may not represent true incidences.
i <- subset(inc, year >= 2004)
system.time(m <- mapply(function(i, y) {
    j <- read.csv(paste0("raw/jk-", substr(y, 3, 4), ".csv"))
    j$PERSON_ID <- as.integer(j$PERSON_ID)
    r <- matrix(unlist(tapply(mapply(function(x) { ifelse(x == 0, 0, ifelse(x == 1, 2.5, 5 * x - 2.5)) }, j[match(i$PERSON_ID, j$PERSON_ID), "AGE_GROUP"]), factor(i$SICK_SYM, levels = k$kcd), function(x) {
        c(length(x), sum(x)) } )), nrow = 2)
    cat(y, ", ")
    return(r)
    }, split(i[1:2], i$year), 2004:2013))
nm <- matrix(rowSums(m, na.rm = TRUE), nrow = 2)
moa <- nm[2,] / nm[1,]
names(moa) <- k$kcd
save(inc, moa, file = "met/onset.RData")
}
newInc <- function(incidence, demogData, years = 2004:2013L, SICK_SYM = NULL, ...) {
    source("src/functions.R")
    #Makes use of driver functions in functions.R
    #Selects subset of data before calculation, because my machines cannot handle the whole thing.
    incidence <- incidence[incidence$SICK_SYM %in% SICK_SYM,]
    ip <- incidenceDemogMerge(incidence, demogData, years = years, ...)
    baseYear <- demogData$STND_Y[1]
    onsetAge <- age(ip$AGE_GROUP) + ip$year - baseYear
    as.vector(tapply(onsetAge, ip$SICK_SYM, mean, na.rm = TRUE))
}
