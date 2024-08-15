#! /usr/bin/Rscript
# Reads rr.csv, selectedKcd.csv and writes pValue.csv.
# If --wd [working directory] is given, the script uses that working directory instead of the current one.
# If --fdr [method] is used, it uses [method] for FDR correction. Defaults to BH.
# If --neg TRUE is used, it calculates the lower tail for p values. This is used to find significant negative relationships.
kCom <- commandArgs(TRUE)
kO <- c(getwd(),"BH", FALSE)
names(kO) <- c("wd", "fdr", "neg")
for (word in names(kO)) {
    kO[word] <- ifelse(is.na(kMa <- match(paste(c("--", word), collapse = ""), kCom)), kO[word], kCom[kMa + 1])
}
kContin <- read.csv(file.path(kO["wd"], "rr.csv"))
kCodes <- read.csv(file.path(kO["wd"], "selectedKcd.csv"))
disease <- rep(kCodes$kcd, length(kCodes$kcd))
riskFactor <- rep(kCodes$kcd, rep(length(kCodes$kcd), length(kCodes$kcd)))
RR <- kContin$ED * (as.double(kContin$C) - kContin$E) / kContin$E / kContin$EcD
pValue <- plnorm(RR ,sdlog = sqrt(1 / kContin$ED + 1 / kContin$EcD - 1 / kContin$E - 1 / (as.double(kContin$C) - kContin$E)), lower.tail = as.logical(kO["neg"]))
complete <- data.frame(riskFactor.kcd = riskFactor, disease.kcd = disease, RR, pValue, pFdr = p.adjust(pValue, kO["fdr"]))
names(complete)[1] <- "kcd"
library(plyr)
complete <- join(complete, kCodes)
names(complete)[c(1, 2, 6)] <- c("riskFactor.kcd", "kcd", "riskFactor.eng")
complete <- join(complete, kCodes)
names(complete)[c(2, 7)] <- c("disease.kcd", "disease.eng")
write.csv(complete, file = file.path(kO["wd"], "pValue.csv"), row.names = FALSE)
