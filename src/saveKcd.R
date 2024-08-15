#! /usr/bin/Rscript
kFet <- read.csv("/data/kmko/fetch.csv")
source("az.R")
source("../omim/functions.R")
kFet <- applyIndexToDataFrame(kFet, kFet$MAIN_SICK %in% kCodes)
library(plyr)
kCnt <- count(kFet$MAIN_SICK[kFet$cohort == 1])
kCodes <- kCnt[[1]][sort(kCnt[[2]], decreasing = TRUE, index.return = TRUE)$ix][1:520]
write.csv(data.frame(kcd = kCodes, eng = mapKcd(kCodes)), file = "/data/kmko/selectedKcd.csv", row.names = FALSE)
