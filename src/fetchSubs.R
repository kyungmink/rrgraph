#!/usr/bin/Rscript
source("/home/kmko/rr/functions.R")
k2C <- c("T20-02.csv", "T20-03.csv")
k4C <- c("T40-02.csv", "T40-03.csv")
kC <- poolSick(k2C, k4C)
kC <- lapply(kC, function(x, mask) { x[mask] }, grepl("[A-Z][0-9][0-9]", kC$SICK_SYM))
kCo <- table(kC$SICK_SYM)
kCod <- names(kCo)[sort(as.vector(kCo), decreasing = TRUE, index.return = TRUE)$ix[1:520]]
# Do not attempt to filter kC for kCod. That changes the cohort.
source("/home/kmko/omim/functions.R")
write.csv(data.frame(kcd = kCod, eng = mapKcd(kCod)), file = "selectedKcd.csv", row.names = FALSE)
kN <- c("04", "05", "06", "07", "08", "09", "10", "11", "12", "13")
k2F <- paste(paste("T20-", kN, sep = ""), ".csv", sep = "")
k4F <- paste(paste("T40-", kN, sep = ""), ".csv", sep = "")
kF <- poolSick(k2F, k4F)
kF <- lapply(kF, function (x, mask) { x[mask] }, mask = (kF$SICK_SYM %in% kCod) & (kF$PERSON_ID %in% unique(kC$PERSON_ID)))
write.csv(data.frame(PERSON_ID = c(kC$PERSON_ID, kF$PERSON_ID), MAIN_SICK = c(kC$SICK_SYM, kF$SICK_SYM), cohort = c(rep(1, length(kC$PERSON_ID)), rep(0, length(kF[[1]]))), stringsAsFactors = FALSE), file = "fetch.csv", row.names = FALSE)
